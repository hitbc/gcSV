/*
 * StatsTracker.hpp
 *
 *  Created on: 2020-5-1
 *      Author: fenghe
 */

#pragma once

#include<iostream>
#include <cstring>
#include <cassert>
#include <functional>
#include <iosfwd>
#include <map>
#include <vector>

extern "C" {
#include "../../clib/bam_file.h"
}

/****************PART ONE: Size Distribution***************************/
/// \brief Accumulate size observations and provide cdf/quantile/smoothed-pdf for the distribution
///
struct SizeData {
	SizeData(unsigned initCount = 0, float initCprob = 0.) :
			count(initCount), cprob(initCprob) {
	}
	unsigned count;
	float cprob;
};
typedef std::map<int, SizeData, std::greater<int>> map_type;

struct SizeDistribution {
	SizeDistribution() :
			_isStatsComputed(false), _totalCount(0), _quantiles(_quantileNum, 0) {
	}

	/// \brief Implements the quantile function for this distribution
	int quantile(const float prob) const;	/// \return The size at which all sizes equal or less are observed with probability \p prob
	float cdf(const int x) const;/// \return Probability of observing value <= \p x
	float pdf(const int x) const;/// \return Probability of observing value \p x, (with a smoothing window)

	unsigned totalObservations() const {
		return _totalCount;
	}

	void addObservation(const int size) {
		_isStatsComputed = false;
		_totalCount++;
		_sizeMap[size].count++;
	}

	void filterObservationsOverQuantile(const float prob);/// filter high value outliers:
	bool isStatSetMatch(const SizeDistribution &pss2);/// compare distributions to determine stats convergence
	//get the count of read for a special ISIZE
	int getSizeCount(int isize) const{
		return _sizeMap[isize].count;
	}

private:
	void calcStats() const;
	static const int _quantileNum = 1000;
	mutable bool _isStatsComputed;
	unsigned _totalCount;
	mutable std::vector<int> _quantiles;
	mutable map_type _sizeMap;
};

//********************************PART TWO: Read PAIR DIRECTION**************************/
typedef int32_t pos_t;

namespace PAIR_ORIENT {

enum index_t { UNKNOWN, Fm, Fp, Rm, Rp, SIZE};

inline const char* label(const index_t i) {
	switch (i) {
	case Fm: return "Fm";
	case Fp: return "Fp";
	case Rm: return "Rm";
	case Rp: return "Rp";
	default:
		return "UNKNOWN";
	}
}

inline index_t get_index(const pos_t pos1, const bool is_fwd_strand1,
		const pos_t pos2, const bool is_fwd_strand2) {
	const bool is_read1_left(pos1 < pos2);

	if (is_fwd_strand1 != is_fwd_strand2) {
		// special-case very short fragments as innies:
		//
		// a few bases of overhang are allowed to account for random matches of
		// the reverse read to the primer
		//
		if (std::abs(pos1 - pos2) <= 2)
			return Rp;

		const bool left_strand(is_read1_left ? is_fwd_strand1 : is_fwd_strand2);
		return (left_strand ? Rp : Rm);
	} else {
		return ((is_read1_left == is_fwd_strand1) ? Fp : Fm);
	}
}

/// inefficient label to id lookup, returns SIZE for unknown string:
inline index_t get_index(const char *str) {
	for (int i(0); i < SIZE; ++i) {
		if (0 == strcmp(str, label(static_cast<index_t>(i))))
			return static_cast<index_t>(i);
	}
	return SIZE;
}
}  // END namespace PAIR_ORIENT

/// pair orientation status wrapper:
struct ReadPairOrient {
	ReadPairOrient() :
			_val(PAIR_ORIENT::UNKNOWN) {
	}

	PAIR_ORIENT::index_t val() const {
		return _val;
	}

	void setVal(const unsigned newVal) {
		assert(newVal < PAIR_ORIENT::SIZE);
		_val = static_cast<PAIR_ORIENT::index_t>(newVal);
	}

private:
	PAIR_ORIENT::index_t _val;
};

//********************************PART Three: ReadCounter**************************/
/// \brief Accumulate read statistics scanned for insert size estimation
struct ReadCounter {
	ReadCounter() :
			_totalReadCount(0), _totalPairedReadCount(0), _totalUnpairedReadCount(
					0), _totalPairedLowMapqReadCount(0), _totalHighConfidenceReadPairCount(
					0) {}

	unsigned totalReadCount() const { return _totalReadCount;}
	unsigned totalPairedReadCount() const {	return _totalPairedReadCount;}
	unsigned totalUnpairedReadCount() const {return _totalUnpairedReadCount;}
	unsigned totalPairedLowMapqReadCount() const {return _totalPairedLowMapqReadCount;}
	unsigned totalHighConfidenceReadPairCount() const {	return _totalHighConfidenceReadPairCount;}
	void addReadCount() {	_totalReadCount++;}
	void addPairedReadCount() {	_totalPairedReadCount++;}
	void addUnpairedReadCount() {	_totalUnpairedReadCount++;}
	void addPairedLowMapqReadCount() {	_totalPairedLowMapqReadCount++;}
	void addHighConfidenceReadPairCount() {	_totalHighConfidenceReadPairCount += 1;}

	friend std::ostream& operator<<(std::ostream &os, const ReadCounter &rs) {
		os << "\tTotal sampled reads: " + std::to_string(rs.totalReadCount()) + "\n"
				<< "\tTotal sampled paired reads: "
						+ std::to_string(rs.totalPairedReadCount()) + "\n"
				<< "\tTotal sampled paired reads passing MAPQ filter: "
						+ std::to_string(
								rs.totalPairedReadCount()
										- rs.totalPairedLowMapqReadCount()) + "\n"
				<< "\tTotal sampled high-confidence read pairs passing all filters: "
						+ std::to_string(rs.totalHighConfidenceReadPairCount())
						+ "\n";
		return os;
	}

private:
	///////////////////////////////////// data:
	unsigned _totalReadCount;
	unsigned _totalPairedReadCount;
	unsigned _totalUnpairedReadCount;
	unsigned _totalPairedLowMapqReadCount;
	unsigned _totalHighConfidenceReadPairCount;
};

//********************************PART Four: Unique Stats**************************/

/// Read pair insert stats can be computed for each sample or read group, this
/// class represents the statistics for one group:
///
struct UniqueStats {
public:
	SizeDistribution fragStats;
	ReadPairOrient relOrients;
	ReadCounter readCounter;
};

struct StatLabel {
	/// if isCopyPtrs then the strings are copied and alloced/de-alloced by
	/// the object, if false the client is responsible these pointers over
	/// the lifetime of the label:
	StatLabel(const char *bamLabelInit, const char *rgLabelInit,
			const bool isCopyPtrsInit = true) :
			isCopyPtrs(isCopyPtrsInit), bamLabel(
					(isCopyPtrs && (nullptr != bamLabelInit)) ?
							strdup(bamLabelInit) : bamLabelInit), rgLabel(
					(isCopyPtrs && (nullptr != rgLabelInit)) ?
							strdup(rgLabelInit) : rgLabelInit) {
		assert(nullptr != bamLabel);
		assert(nullptr != rgLabel);
	}

	StatLabel(const StatLabel &rhs) :
			isCopyPtrs(rhs.isCopyPtrs), bamLabel(
					isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel), rgLabel(
					isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel) {
	}

	StatLabel& operator=(const StatLabel &rhs) {
		if (this == &rhs)
			return *this;
		clear();
		isCopyPtrs = rhs.isCopyPtrs;
		bamLabel = (isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel);
		rgLabel = (isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel);
		return *this;
	}

public:
	~StatLabel() {
		clear();
	}

	/// sort allowing for nullptr string pointers in primary and secondary key:
	bool operator<(const StatLabel &rhs) const {
		const int scval(strcmp(bamLabel, rhs.bamLabel));
		if (scval < 0)
			return true;
		if (scval == 0) {
			return (strcmp(rgLabel, rhs.rgLabel) < 0);
		}

		return false;
	}

	friend std::ostream& operator<<(std::ostream &os,
			const StatLabel &rgl) {
		os << "read group '" << rgl.rgLabel << "' in bam file '"
				<< rgl.bamLabel;
		return os;
	}

private:
	void clear() {
		if (isCopyPtrs) {
			if (nullptr != bamLabel)
				free(const_cast<char*>(bamLabel));
			if (nullptr != rgLabel)
				free(const_cast<char*>(rgLabel));
		}
	}

	bool isCopyPtrs;

public:
	const char *bamLabel;
	const char *rgLabel;
};

/// track pair orientation so that a consensus can be found for a read group
struct OrientTracker {
	OrientTracker(const char *bamLabel, const char *rgLabel) :
			_isFinalized(false), _totalOrientCount(0), _rgLabel(bamLabel,
					rgLabel) {
		std::fill(_orientCount.begin(), _orientCount.end(), 0);
	}

	void addOrient(const PAIR_ORIENT::index_t ori) {
		static const unsigned maxOrientCount(100000);
		if (_totalOrientCount >= maxOrientCount)
			return;
		if (ori == PAIR_ORIENT::UNKNOWN)
			return;
		addOrientImpl(ori);
	}

	const ReadPairOrient& getConsensusOrient(const ReadCounter &readCounter) {
		finalize(readCounter);
		return _finalOrient;
	}

	unsigned getMinCount() {
		static const unsigned minCount(100);
		return minCount;
	}

	bool isOrientCountGood() {
		return (_totalOrientCount >= getMinCount());
	}

private:
	void addOrientImpl(const PAIR_ORIENT::index_t ori) {
		assert(!_isFinalized);
		assert(ori < PAIR_ORIENT::SIZE);

		_orientCount[ori]++;
		_totalOrientCount++;
	}
	void finalize(const ReadCounter &readCounter);

	bool _isFinalized;
	unsigned _totalOrientCount;
	const StatLabel _rgLabel;
	std::array<unsigned, PAIR_ORIENT::SIZE> _orientCount;
	ReadPairOrient _finalOrient;
};

struct SimpleRead {
	SimpleRead(PAIR_ORIENT::index_t ort, unsigned sz) :
			_orient(ort), _insertSize(sz) {
	}

	PAIR_ORIENT::index_t _orient;
	unsigned _insertSize;
};

struct ReadBuffer {
	ReadBuffer() :
			_abnormalRpCount(0), _observationRpCount(0) {
	}

	void updateBuffer(PAIR_ORIENT::index_t ort, unsigned sz) {
		_readInfo.emplace_back(ort, sz);
		if (ort == PAIR_ORIENT::Rp) {
			_observationRpCount++;
			if (sz >= 5000)
				_abnormalRpCount++;
		}
	}
	bool isBufferFull() const {
		return (_observationRpCount >= 1000);
	}
	bool isBufferNormal() const {
		if (_observationRpCount == 0)
			return false;
		return ((_abnormalRpCount / (float) _observationRpCount) < 0.01);
	}
	unsigned getAbnormalCount() {
		return _abnormalRpCount;
	}
	unsigned getObservationCount() {
		return _observationRpCount;
	}
	const std::vector<SimpleRead>& getBufferedReads() {
		return _readInfo;
	}
	void clearBuffer() {
		_abnormalRpCount = 0;
		_observationRpCount = 0;
		_readInfo.clear();
	}

private:
	unsigned _abnormalRpCount;
	unsigned _observationRpCount;
	std::vector<SimpleRead> _readInfo;
};

#define statsCheckCnt 100000
enum RGT_RETURN {
	RGT_CONTINUE, RGT_BREAK, RGT_NORMAL
};
struct StatsTracker {
	StatsTracker(const char *bamLabel = nullptr, const char *rgLabel =
			nullptr, const std::string &defaultStatsFilename = "") :
			_rgLabel(bamLabel, rgLabel), _orientInfo(bamLabel, rgLabel), _defaultStatsFilename(
					defaultStatsFilename) {
	}

	unsigned insertSizeObservations() const {	return _stats.fragStats.totalObservations();}
	unsigned getMinObservationCount() const {	static const unsigned minObservations(100);	return minObservations;	}
	bool isObservationCountGood() const {	return (insertSizeObservations() >= getMinObservationCount());	}
	void checkInsertSizeCount() {	if ((insertSizeObservations() % statsCheckCnt) == 0) _isChecked = true;	}
	bool isInsertSizeChecked() const {	return _isChecked; }
	void clearChecked() { _isChecked = false; }
	bool isInsertSizeConverged() const { return _isInsertSizeConverged; }
	bool isCheckedOrConverged() const {	return (_isChecked || isInsertSizeConverged()); }
	void updateInsertSizeConvergenceTest() { // check convergence
		if (_oldInsertSize.totalObservations() > 0) {
			_isInsertSizeConverged = _oldInsertSize.isStatSetMatch(
					_stats.fragStats);
		}
		_oldInsertSize = _stats.fragStats;
	}
	ReadBuffer& getBuffer() { return _buffer; }
	void addBufferedData();
	const OrientTracker& getOrientInfo() const {	return _orientInfo; }
	const UniqueStats& getStats() const {	assert(_isFinalized); return _stats; } /// getting a const ref of the stats forces finalization steps:
	void addReadCount() {	_stats.readCounter.addReadCount();}
	void addPairedReadCount() {	_stats.readCounter.addPairedReadCount();}
	void addUnpairedReadCount() {	_stats.readCounter.addUnpairedReadCount();}
	void addPairedLowMapqReadCount() {	_stats.readCounter.addPairedLowMapqReadCount();}
	void addHighConfidenceReadPairCount() {	_stats.readCounter.addHighConfidenceReadPairCount();}
	/// Add one observation to the buffer
	/// If the buffer is full, AND if the fragment size distribution in the buffer looks normal, add the buffered data;
	/// otherwise, discard the buffer and move to the next region
	bool addObservation(PAIR_ORIENT::index_t ori, unsigned sz);
	void finalize();
	void handleReadRecordBasic(bam1_t *b);
	RGT_RETURN handleReadRecordCheck(bam1_t *b);

private:
	void addOrient(const PAIR_ORIENT::index_t ori) {
		assert(!_isFinalized);
		_orientInfo.addOrient(ori);
	}
	void addInsertSize(const int size) {
		assert(!_isFinalized);
		_stats.fragStats.addObservation(size);
	}

	bool _isFinalized = false;
	const StatLabel _rgLabel;
	OrientTracker _orientInfo;

	bool _isChecked = false;
	bool _isInsertSizeConverged = false;
	SizeDistribution _oldInsertSize; // previous fragment distribution is stored to determine convergence

	ReadBuffer _buffer;
	UniqueStats _stats;
	const std::string _defaultStatsFilename;
};
