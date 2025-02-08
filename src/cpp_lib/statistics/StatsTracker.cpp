#include "StatsTracker.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

//

static void populateCdfQuantiles(map_type &sizeMap,
		const unsigned totalCount, std::vector<int> &quantiles) {
	const unsigned quantileNum(quantiles.size());
	const float pFactor(1 / static_cast<float>(totalCount));

	unsigned fillBase(0);
	unsigned cumulativeCount(0);

	for (auto map_it = (sizeMap.rbegin()); map_it != sizeMap.rend(); map_it++) {
		cumulativeCount += (map_it->second.count);
		assert(cumulativeCount <= totalCount);

		// update the hash map with cumulative prob value
		map_it->second.cprob = (cumulativeCount * pFactor);

		const unsigned fillNext = static_cast<unsigned>(rint(
				map_it->second.cprob * quantileNum));
		for (; fillBase < fillNext; fillBase++) {
			quantiles[fillBase] = map_it->first;
		}
	}
}

void SizeDistribution::calcStats() const {
#ifdef DEBUG_RPS
  log_os << "Calculating stats...\n"
         << "numOfSized=" << _sizeMap.size() << "\n";
#endif
	_isStatsComputed = true;
	if (_sizeMap.empty())
		return;

	populateCdfQuantiles(_sizeMap, _totalCount, _quantiles);
}

int SizeDistribution::quantile(const float prob) const {
	assert((prob >= 0.) && (prob <= 1.));

	static const int maxBin(_quantileNum - 1);
	if (!_isStatsComputed)
		calcStats();

	int bin(static_cast<int>(ceil(prob * _quantileNum) - 1));
	if (bin < 0)
		bin = 0;
	if (bin > maxBin)
		bin = maxBin;
	return _quantiles[bin];
}

float SizeDistribution::cdf(const int size) const {
	if (!_isStatsComputed)
		calcStats();

	// map uses greater<int> for comp, so lower bound is "first element not greater than" size, from a list
	// sorted high->low
	const map_type::const_iterator sizeIter(_sizeMap.lower_bound(size));
	if (sizeIter == _sizeMap.end())
		return 0;
	return sizeIter->second.cprob;
}

float SizeDistribution::pdf(const int size) const {
	if (!_isStatsComputed)
		calcStats();

	static const unsigned targetSampleSize(5);

	unsigned count(0);
	int minSize(size);
	int maxSize(size);

	bool isMinBound(false);
	bool isMaxBound(false);

	/// scheme: get the five closest (in bin space) samples and sum them divided by the range required to find
	/// them

	// map uses greater<int> for comp, so lower bound is "first element not greater than" size, from a list
	// sorted high->low
	map_type::const_iterator lowIter(_sizeMap.lower_bound(size));

	if (lowIter == _sizeMap.end()) {
		isMinBound = true;
	}

	map_type::const_iterator highIter(lowIter);

	if (highIter == _sizeMap.begin()) {
		isMaxBound = true;
	} else {
		--highIter;
	}

	for (unsigned sampleIndex(0); sampleIndex < targetSampleSize;
			++sampleIndex) {
		// determine whether fwd or rev pointer is closer to size:
		if (isMinBound && isMaxBound)
			break;

		bool isChooseLow(true);
		if (isMinBound) {
			isChooseLow = false;
		} else if (isMaxBound) {
			isChooseLow = true;
		} else {
			isChooseLow = (std::abs(lowIter->first - size)
					<= std::abs(highIter->first - size));
		}

		if (isChooseLow) {
			minSize = lowIter->first;
			count += lowIter->second.count;
			++lowIter;

			if (lowIter == _sizeMap.end())
				isMinBound = true;
		} else {
			maxSize = highIter->first;
			count += highIter->second.count;
			if (highIter == _sizeMap.begin()) {
				isMaxBound = true;
			} else {
				--highIter;
			}
		}
	}

	assert(maxSize >= minSize);

	return count
			/ (static_cast<float>(_totalCount)
					* static_cast<float>(1 + maxSize - minSize));
}

void SizeDistribution::filterObservationsOverQuantile(const float prob) {
	const int maxSize(quantile(prob));
	const map_type::iterator sizeBegin(_sizeMap.begin());
	map_type::iterator sizeEnd(_sizeMap.lower_bound(maxSize));

	for (map_type::iterator sizeIter(sizeBegin); sizeIter != sizeEnd;
			++sizeIter) {
		if (sizeIter->first <= maxSize) {
			sizeEnd = sizeIter;
			break;
		}
		_totalCount -= sizeIter->second.count;
	}
	_sizeMap.erase(sizeBegin, sizeEnd);

	_isStatsComputed = false;
}

bool SizeDistribution::isStatSetMatch(const SizeDistribution &pss2) {
	static const float cdfPrecision(0.001f);

	for (float prob(0.05f); prob < 1; prob += 0.1f) {
		// check if percentile values equal
		if (std::abs(quantile(prob) - pss2.quantile(prob)) >= 1) {
			return false;
		}

		// check the convergence of fragsize cdf
		const int fragSize(pss2.quantile(prob));
		if (std::abs(cdf(fragSize) - pss2.cdf(fragSize)) >= cdfPrecision) {
			return false;
		}
	}
	return true;
}

//**************************************************************************************/

void OrientTracker::finalize(const ReadCounter &readCounter) {
	if (_isFinalized)
		return;
	bool isMaxIndex(false);
	unsigned maxIndex(0);
	for (unsigned i(0); i < _orientCount.size(); ++i) {
		if ((!isMaxIndex) || (_orientCount[i] > _orientCount[maxIndex])) {
			isMaxIndex = true;
			maxIndex = i;
		}
	}

	assert(isMaxIndex);

	_finalOrient.setVal(maxIndex);

	{
		// make sure there's a dominant consensus orientation and that we have a minimum number of samples:
		static const float minMaxFrac(0.9f);

		if (!isOrientCountGood()) {
			std::cerr << "Too few high-confidence read pairs ("
					<< _totalOrientCount
					<< ") to determine pair orientation for " << _rgLabel
					<< "'\n" << "\tAt least " << getMinCount()
					<< " high-confidence read pairs are required to determine pair orientation.\n"
					<< readCounter << "\n";
		}

		const unsigned minMaxCount(
				static_cast<unsigned>(minMaxFrac * _totalOrientCount));
		if (_orientCount[maxIndex] < minMaxCount) {
			const unsigned maxPercent(
					(_orientCount[maxIndex] * 100) / _totalOrientCount);
			std::cerr << "Can't determine consensus pair orientation of "
					<< _rgLabel << ".\n"
					<< "' (" << maxPercent << "% of " << _totalOrientCount
					<< " total used read pairs)\n" << "\tThe fraction of '"
					<< "' among total high-confidence read pairs needs to be more than "
					<< minMaxFrac
					<< " to determine consensus pair orientation.\n"
					<< readCounter << "\n";
		}
	}

	_isFinalized = true;
}

void StatsTracker::addBufferedData() {
	for (const SimpleRead &srd : _buffer.getBufferedReads()) {
		const PAIR_ORIENT::index_t ori(srd._orient);
		addOrient(ori);

		// we define "high-confidence" read pairs as those reads passing all filters
		addHighConfidenceReadPairCount();

		if (ori != PAIR_ORIENT::Rp)
			continue;
		addInsertSize(srd._insertSize);
	}
}

bool StatsTracker::addObservation(PAIR_ORIENT::index_t ori, unsigned sz) {
	bool isNormal(true);

	_buffer.updateBuffer(ori, sz);

	if (_buffer.isBufferFull()) {
		// check abnormal fragment-size distribution in the buffer
		if (_buffer.isBufferNormal()) {
			addBufferedData();
			checkInsertSizeCount();
		} else {
			isNormal = false;
#ifdef DEBUG_RPS
        std::cerr << "The previous region (buffered) contains too many abnormal reads. "
                  << "abnormalCount=" << _buffer.getAbnormalCount()
                  << " observationCount=" << _buffer.getObservationCount() << "\n";
#endif
		}
		_buffer.clearBuffer();
	}

	return isNormal;
}

void StatsTracker::finalize() {
	if (_isFinalized)
		return;

	// add the remaining data in the buffer
	if (_buffer.isBufferNormal()) {
		addBufferedData();
	}
	_buffer.clearBuffer();

	// finalize pair orientation:
	_stats.relOrients = _orientInfo.getConsensusOrient(_stats.readCounter);

	if (_stats.relOrients.val() != PAIR_ORIENT::Rp) {
		std::cerr << "Unexpected consensus read orientation ("
				<< 	") for " << _rgLabel << "\n"
				<< "\tManta currently handles paired-end (FR) reads only.\n";
	}

	// finalize insert size distro:
	if (!isInsertSizeConverged()) {
		if (!isObservationCountGood()) {

			std::cerr << "Can't generate pair statistics for " << _rgLabel
					<< "\n"
					<< "\tTotal high-confidence read pairs (FR) used for insert size estimation: "
					<< insertSizeObservations() << "\n" << "\tAt least "
					<< getMinObservationCount()
					<< " high-confidence read pairs (FR) are required to estimate insert size.\n"
					<< _stats.readCounter << "\n";
		} else if (!isInsertSizeChecked()) {
			updateInsertSizeConvergenceTest();
		}

		if (!isInsertSizeConverged()) {
			std::cerr << "WARNING: read pair statistics did not converge for "
					<< _rgLabel << "\n"
					<< "\tTotal high-confidence read pairs (FR) used for insert size estimation: "
					<< insertSizeObservations() << "\n" << _stats.readCounter
					<< "\n";
		}
	}

	// final step before saving is to cut-off the extreme end of the fragment size distribution, this
	// is similar the some aligner's proper-pair bit definition of (3x the standard mean, etc.)
	static const float filterQuant(0.9995f);
	_stats.fragStats.filterObservationsOverQuantile(filterQuant);

	_isFinalized = true;
}

/// given an input integer, return an integer with all but the highest 4 decimal digits set to zero
///
/// this method is not written effeciently, and not intended for general integer truncation.
/// it is used as part of a simple compression scheme for the fragment sizes of the frag size
/// distribution
///
static unsigned getSimplifiedFragSize(unsigned fragmentSize) {
	unsigned fragSize(fragmentSize);
	unsigned steps(0); // reduce fragsize resolution for very large sizes:
	while (fragSize > 1000) {
		fragSize /= 10;
		steps++;
	}
	for (unsigned stepIndex(0); stepIndex < steps; ++stepIndex)
		fragSize *= 10;
	return fragSize;
}

/// This produces a useful result only when both reads align to the same
/// chromosome.
static PAIR_ORIENT::index_t getRelOrient(bam1_t *b) {
	pos_t pos1 = b->core.pos;
	bool is_fwd_strand1 = bam_is_fwd_strand(b);
	pos_t pos2 = b->core.mpos;
	bool is_fwd_strand2 = bam_is_mate_fwd_strand(b);

	if (!bam_is_filter(b)) {
		std::swap(pos1, pos2);
		std::swap(is_fwd_strand1, is_fwd_strand2);
	}
	return PAIR_ORIENT::get_index(pos1, is_fwd_strand1, pos2, is_fwd_strand2);
}

void StatsTracker::handleReadRecordBasic(bam1_t *b) {
	addReadCount(); //#记录该read的信息；结构体rgInfo记录了该slice的统计信息
	if (bam_is_paired(b)) {
		addPairedReadCount();
		if (bam_map_qual(b) == 0)
			addPairedLowMapqReadCount();
	} else
		addUnpairedReadCount();
}

RGT_RETURN StatsTracker::handleReadRecordCheck(bam1_t *b) {
	if (isInsertSizeConverged())
		return RGT_CONTINUE; //#过滤情况2
	const PAIR_ORIENT::index_t ori(getRelOrient(b));
	unsigned fragSize(0);
	if (ori == PAIR_ORIENT::Rp)
		fragSize = getSimplifiedFragSize(getFragSizeMinusSkip(b));
	if (!addObservation(ori, fragSize))
		return RGT_BREAK;
	if (!isInsertSizeChecked())
		return RGT_CONTINUE; //#过滤情况3，需要100K个read顺利pair了，才能通过
	// check convergence
	updateInsertSizeConvergenceTest();
	return RGT_NORMAL;
}
