/*
 * StatOutputer.cpp
 *
 *  Created on: 2020年5月1日
 *      Author: fenghe
 */


//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "StatsManager.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

bool ReadPairDepthFilter::isFilterRead(bam1_t *b) {
		static const unsigned maxPosCount(1);

		if (b->core.tid != _lastTargetId) {
			_goodMates.clear();
			_lastTargetId = b->core.tid;
			_posCount = 0;
			_lastPos = b->core.pos;
		} else if (b->core.pos != _lastPos) {
			_posCount = 0;
			_lastPos = b->core.pos;
		}

		// Assert only two reads per fragment
		const unsigned readNum(bam_is_first(b) ? 1 : 2);
		assert(bam_is_second(b) == (readNum == 2));

		// Filter pairs with templateSize 0 (unknown)
		if (b->core.isize == 0)
			return true;

		// sample each read pair once by sampling stats from
		// downstream read only, or whichever read is encountered
		// second if the read and its mate start at the same position:
		const bool isDownstream(b->core.pos > b->core.mpos);
		const bool isSamePos(b->core.pos == b->core.mpos);

		if (isDownstream || isSamePos) {
			const int mateReadNo(bam_is_first(b) ? 2 : 1); //当read是first的时候，mate的No就是2；当read的No是2的时候，mate的No就是1
			mateMap_t::iterator i(_goodMates.find(getKey((char*) bam_qname(b), mateReadNo)));

			if (i == _goodMates.end()) {    //i == the last one: item not found
				if (isDownstream)
					return true;
			} else {
				_goodMates.erase(i);
				return false;
			}
		}

		// to prevent high-depth pileups from overly biasing the
		// read stats, we only take maxPosCount read pairs from each start
		// pos. by not inserting a key in goodMates, we also filter
		// the downstream mate:
		if (_posCount >= maxPosCount)
			return true;
		++_posCount;

		// crude mechanism to manage total set memory
		if (_goodMates.size() > maxMateSetSize)
			_goodMates.clear();

		// Ignore pairs where the upstream mate has a refskip, since we cannot
		// compute the correct insert size later when looking at the downstream mate
		// (Or we would have to save the total refskip length here)
		if (hasRefSkip(b))
			return true;

		_goodMates.insert(getKey(b));
		return true;
	}

bool CoreReadFilter::isFilterRead(bam1_t *b) {
	// filter common categories of undesirable reads:
	if (isReadFilteredCore(b))
		return true;
	if (isNonStrictSupplement(b))
		return true;
	if (!bam_is_mapped_pair(b))
		return true;
	if (bam_map_qual(b) == 0)
		return true;
	if (bam_isSASplit(b))
		return true;    // filter any split reads with an SA tag:
	if (readFilteredAlignment(b))
		return true; // remove alignments other than {X}M({Z}N{X2}M)?({Y}S)? (or reverse for reverse strand)
	if (pairFilter.isFilterRead(b))
		return true;  // filter out upstream reads and high depth regions:
	return false;
}


/// this struct exists for the sole purpose of xml output:
struct ReadGroupStatsExporter {

	std::string bamFile;
	std::string readGroup;
	UniqueStats groupStats;
};

void StatsManager::merge(const StatsManager &rhs) {
	const unsigned numGroups(rhs.size());
	for (unsigned i(0); i < numGroups; ++i) {
		const StatLabel &mkey(rhs.getKey(i));
		if (_group.test_key(mkey)) {
			std::cerr << "Can't merge stats set objects with repeated key: '"
					<< mkey << "'";
			exit(EXIT_FAILURE);
		}

		setStats(mkey, rhs.getStats(i));
	}
}

#ifdef READ_GROUPS
	#define USE_RG true
#else
	#define USE_RG false
#endif

void StatsManager::handleBamCramStats(const char *alignmentFilename, double *ave_depth)
{
	TrackerManager rgManager(USE_RG, alignmentFilename,	defaultStatsFilename);
	Bam_file bf = {0}; bam_file_open(alignmentFilename, referenceFilename, NULL, &bf);
	const bam_hdr_t &header(*(bf._hdr));

	std::vector<ChromInfo> chromList;
	for (int32_t i(0); i < header.n_targets; ++i)
		chromList.push_back(ChromInfo(i, header.target_len[i]));

	bool isStopEstimation(false);
	bool isActiveChrom(true);

	float total_depth = 0; int total_region_num = 0;

	while (isActiveChrom && (!isStopEstimation)) {
		isActiveChrom = false;

		for (auto &chrom : chromList) {
			if (isStopEstimation)
				break; // keep sampling until either the chromosome has been exhuasted or the current chunk has been sufficiently sampled
			bool isFinishedSlice(false);
			R_region region;
			region.chr_ID = chrom.chr_ID;

			while (!isFinishedSlice) {
				const int32_t startPos(chrom.highestPos + 1);
				if (startPos >= chrom.size)
					break;
				region.st_pos = startPos;
				region.ed_pos = chrom.size;
				resetRegion_ID(&bf, &region);
				uint64_t total_base_region = 0;

				while (bam_next(&bf)) {
					bam1_t *b = &(bf._brec);
					if (b->core.pos < startPos)
						continue;
					chrom.highestPos = b->core.pos;
					isActiveChrom = true;

					if(bam_is_secondary(b) || bam_is_supplementary(b) || bam_is_duplicate(b))
						continue;
					total_base_region += b->core.l_qseq;
					//total_base += b->core.l_qseq;

					StatsTracker &rgInfo(rgManager.getTracker(b));
					rgInfo.handleReadRecordBasic(b);
					if (rgManager.isFiltered(b))
						continue;
					RGT_RETURN r = rgInfo.handleReadRecordCheck(b);
					if (r == RGT_CONTINUE)
						continue;
					else if (r == RGT_BREAK) {
						chrom.highestPos += std::max(1, chrom.size / 100);
						break;
					} //else do nothing

					isFinishedSlice = rgManager.isFinishedSlice();
					if (!isFinishedSlice)
						continue;
					isStopEstimation = rgManager.isStopEstimation();
					// break from reading the current chromosome
					break;
				}
				if(isFinishedSlice == true){
					total_region_num += 1;
					float region_depth = total_base_region/(chrom.highestPos - startPos);
					total_depth += region_depth;
					if(false){
						fprintf(stderr, "%f @ [%d:%d-%d] isFinishedSlice: %d \n", region_depth, chrom.chr_ID, startPos, chrom.highestPos, isFinishedSlice);
					}
				}

				// move to next region if no read falling in the current region
				if (chrom.highestPos <= startPos) {
					chrom.highestPos += std::max(1, chrom.size / 100);
				}
			}
		}
	}

	*ave_depth = (double)total_depth/(total_region_num);

	for (auto &val : rgManager.getMap())
		setStats(val.first, val.second.getStats());
	bam_file_close(&bf);
}
