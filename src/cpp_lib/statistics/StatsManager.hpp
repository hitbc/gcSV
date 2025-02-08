/*
 * StatsManager.hpp
 *
 *  Created on: 2020年5月1日
 *      Author: fenghe
 */

#pragma once

#include <iostream>
#include <map>
#include <vector>

#include "StatsTracker.hpp"
#include <unordered_set>

//using function: StatsManager::simpleGetStats to simple get ISIZE

/// this object handles filtration of reads which are:
///
/// 1. not downstream or (if both reads start at same position) not the second in order in the bam file
/// 2. part of a high depth pileup region
///
struct ReadPairDepthFilter {
	bool isFilterRead(bam1_t *b);
private:
	std::string getKey(const char *initQname, const int initReadNo){
		std::string s(initQname); s += std::to_string(initReadNo);
		return s;
	}
	std::string getKey(bam1_t *b){
		std::string s((char *)bam_qname(b)); s += std::to_string(bam_read_no(b));
		return s;
	}
	static const unsigned maxMateSetSize = 100000;
	typedef std::unordered_set<std::string> mateMap_t;
	unsigned _posCount = 0;
	mateMap_t _goodMates{maxMateSetSize};
	int _lastTargetId = 0;
	int _lastPos = 0;
};

struct CoreReadFilter {
	bool isFilterRead(bam1_t *b);
private:
	ReadPairDepthFilter pairFilter;
};

/// \return The BAM record's auxillary RG tag value, or an empty string if no RG tag exists.
inline const char* getReadGroup(bam1_t *b) {
	char RG[5] = "RG";
	const char *rgStr(bam_get_string_tag(b, RG));
	return ((nullptr == rgStr) ? "" : rgStr);
}

/// manage the info structs for each RG
struct TrackerManager {
	typedef std::map<StatLabel, StatsTracker> RGMapType;

	TrackerManager(bool useRG, const std::string &statsBamFile,
			const std::string &defaultStatsFilename) :
			_useRG(useRG), _statsBamFile(statsBamFile), _defaultStatsFilename(
					defaultStatsFilename) {
		if(!_useRG) {//get default Tracker
			std::pair<RGMapType::iterator, bool> retval = _rgTracker.insert(
					std::make_pair(
							StatLabel(_statsBamFile.c_str(), ""),
							StatsTracker(_statsBamFile.c_str(), "",
									_defaultStatsFilename)));
			assert(retval.second);
			_defaultTracker = &(retval.first->second);
		}
	}

	StatsTracker& getTracker(bam1_t *b) {
		if (!_useRG)
			return *_defaultTracker;
		const char *readGroup(getReadGroup(b));
		StatLabel rgKey(_statsBamFile.c_str(), readGroup, false);
		RGMapType::iterator rgIter(_rgTracker.find(rgKey));
		if (rgIter == _rgTracker.end()) {
			std::pair<RGMapType::iterator, bool> retval;
			retval = _rgTracker.insert(
					std::make_pair(
							StatLabel(_statsBamFile.c_str(), readGroup),
							StatsTracker(_statsBamFile.c_str(), readGroup,
									_defaultStatsFilename)));
			assert(retval.second);
			rgIter = retval.first;
		}
		return rgIter->second;
	}

	/// check if all read groups have been sufficiently sampled for a slice
	/// for each read group, either 100k samples has been collected, or the insert size distrubution has
	/// converged
	bool isFinishedSlice() {
		for (auto &val : _rgTracker) {
			if (!val.second.isCheckedOrConverged())
				return false;
		}
		for (RGMapType::value_type &val : _rgTracker) {
			val.second.clearChecked();
		}
		return true;
	}

	/// test if all read groups have converged or hit other stopping conditions
	bool isStopEstimation() {
		static const unsigned maxRecordCount(5000000);
		for (RGMapType::value_type &val : _rgTracker) {
			if (!(val.second.isInsertSizeConverged()
					|| (val.second.insertSizeObservations() > maxRecordCount)))
				return false;
		}
		return true;
	}

	const RGMapType& getMap() {
		finalize();
		return _rgTracker;
	}

	bool isFiltered(bam1_t *b)
	{
		return filter.isFilterRead(b);
	}

private:
	TrackerManager(const TrackerManager & c)
	{
		xassert(0,"TrackerManager can`t be copy create!");
	}

	void finalize() {
		if (_isFinalized)
			return;
		for (RGMapType::value_type &val : _rgTracker) {
			val.second.finalize();
		}
		_isFinalized = true;
	}
	bool _useRG = false;
	bool _isFinalized = false;

	RGMapType _rgTracker;
	StatsTracker * _defaultTracker = nullptr;

	const std::string _statsBamFile;
	std::string _defaultStatsFilename;

	CoreReadFilter filter;
};


/// \brief Provides something like a map, but with sequential id numbers
/// assigned to each key starting from 0
///
/// The id numbers can be useful for faster lookup of the value, while
/// retaining the option of doing key lookup when required
///
template<typename K, typename V, typename COMPARE = std::less<K>>
struct id_map {
	/// \brief Update map with (key,value) and return id
	///
	unsigned insert(const K &key, const V &value) {
		const typename k2id_t::const_iterator i(_k2id.find(key));
		if (i == _k2id.end()) {
			const unsigned id(_id2kv.size());
			_k2id[key] = id;
			_id2kv.push_back(std::make_pair(key, value));
			return id;
		} else {
			_id2kv[i->second] = std::make_pair(key, value);
			return i->second;
		}
	}

	/// \brief Test if key exists in map
	bool test_key(const K &key) const {
		return (_k2id.find(key) != _k2id.end());
	}

	/// \brief Get id of inserted key
	unsigned get_optional_id(const K &key) const {
		const typename k2id_t::const_iterator i(_k2id.find(key));
		if (i == _k2id.end()) {
			return MAX_uint32_t;
		}
		return i->second;
	}

	/// \brief Get id of inserted key
	unsigned get_id(const K &key) const {
		const typename k2id_t::const_iterator i(_k2id.find(key));
		if (i == _k2id.end()) {
			std::cerr << "id_map.get_id(): invalid key" << std::endl;
		}
		return i->second;
	}

	/// \brief Get pre-existing key
	const K& get_key(const unsigned id) const {
		if (id >= _id2kv.size()) {
			std::cerr << "idmap.get_key(): invalid id" << std::endl;
		}
		return _id2kv[id].first;
	}

	/// \brief Get pre-existing key
	const V& get_value(const unsigned id) const {
		if (id >= _id2kv.size()) {
			std::cerr << "idmap.get_value(): invalid id" << std::endl;
		}
		return _id2kv[id].second;
	}

	bool empty() const {
		return _id2kv.empty();
	}

	unsigned size() const {
		return _id2kv.size();
	}

	void clear() {
		_k2id.clear();
		_id2kv.clear();
	}

private:
	typedef std::map<K, unsigned, COMPARE> k2id_t;

	k2id_t _k2id;
	std::vector<std::pair<K, V>> _id2kv;
};

/***************************PART TWO: Output************************************/

/// \brief manages multiple read_group_stats
///
#define default_insert_size_min 0
#define default_insert_size_max 1300
struct StatsManager {
	StatsManager(const char * referenceFilename_, std::string defaultStatsFilename_):
		referenceFilename(referenceFilename_), defaultStatsFilename(defaultStatsFilename_){}
	StatsManager(){}

//	static void simpleGetStats(
//			const std::string &referenceFilename,
//			const std::string &alignmentFilename,
//			FILE* output) {
//		// calculate fragment size statistics for all read groups in all bams
//		StatsManager rstats(referenceFilename.c_str(), "");
//		rstats.handleBamCramStats(alignmentFilename.c_str());
//		uint32_t minInsertLen = rstats.getInsertLen(alignmentFilename.c_str(), 0.01f);
//		uint32_t middleInsertLen = rstats.getInsertLen(alignmentFilename.c_str(), 0.5f);
//		uint32_t maxInsertLen = rstats.getInsertLen(alignmentFilename.c_str(), 0.99f);
//		fprintf(output, "%d %d %d\n", minInsertLen, middleInsertLen, maxInsertLen);
//	}

	typedef StatLabel KeyType;
	bool empty() const { return _group.empty(); }
	unsigned size() const { return _group.size(); }

	/// \brief get the index of a read group
	///
	/// the index can be used for fast lookup of the
	/// stats for that group
	///
	/// if the group does not exist, the returned value
	/// evaluates to false per boost::optional
	///
	/// Each read group is identified as a combination of a bam filename and
	/// an RG tag label. An empty label refers to the "default" read group
	/// for the file (all records that had no RG tag).
	unsigned getGroupIndex(const StatLabel &rgLabel) const { 	return _group.get_optional_id(rgLabel); 	}
	/// get stats associated with index
	const UniqueStats& getStats(const unsigned groupIndex) const { 	return _group.get_value(groupIndex); }
	const KeyType& getKey(const unsigned groupIndex) const {	return _group.get_key(groupIndex); }
	/// set stats for index
	void setStats(const StatLabel &rgLabel, const UniqueStats &rps) { _group.insert(rgLabel, rps); }

	unsigned getMinInsertLen(const std::string &bamFilename){
		unsigned idx = getGroupIndex(StatLabel(bamFilename.c_str(), ""));
		if(idx != MAX_uint32_t)		return getStats(idx).fragStats.quantile(0.01f);
		else	{
			std::cerr << "Warning: can`t properly get MIN insert_size for " << bamFilename <<
					", using default_insert_size_min" << default_insert_size_min << "\n";
			return default_insert_size_min;
		}
	}

	unsigned getMaxInsertLen(const std::string &bamFilename){
		unsigned idx = getGroupIndex(StatLabel(bamFilename.c_str(), ""));
		if(idx != MAX_uint32_t)	return getStats(idx).fragStats.quantile(0.99f);
		else{
			std::cerr << "Warning: can`t properly get MAX insert_size for " << bamFilename <<
					", using default_insert_size_max" << default_insert_size_max << "\n";
			return default_insert_size_max;
		}
	}

	unsigned getInsertLen(const char * bamFilename, float percent){
		unsigned idx = getGroupIndex(StatLabel(bamFilename, ""));
		if(idx != MAX_uint32_t)		return getStats(idx).fragStats.quantile(percent);
		else{
			if(percent < 0.5){
				std::cerr << "Warning: can`t properly get isize for " << bamFilename <<
						" @[" << percent <<
						"], using default isize_min [" << default_insert_size_min << "]\n";
				return default_insert_size_min;
			}
			else{
				std::cerr << "Warning: can`t properly get isize for " << bamFilename <<
						" @[" << percent <<
						"], using default isize_max[" << default_insert_size_max << "]\n";
				return default_insert_size_max;
			}
		}
	}

	void getBreakPoint_Distribution(const std::string &bamFilename, int read_len,
			std::vector<float> &DR_bp_distribution, std::vector<float> &SH_bp_distribution,
			std::vector<float> &UM_stPos_distribution, int &START_OFFSET_UM){
		unsigned idx = getGroupIndex(StatLabel(bamFilename.c_str(), ""));
		if(idx == MAX_uint32_t){
			std::cerr << "Warning: can`t properly get break point distribution for "
					<< bamFilename << "\n";
			return;
		}
		//for discordant /hard clip signals break point distribution
		int totalPairedReadCount = getStats(idx).readCounter.totalHighConfidenceReadPairCount();
		int max_len = getStats(idx).fragStats.quantile(0.99f);
		int max_probability_size = max_len -  2*read_len;
		DR_bp_distribution.clear();
		const SizeDistribution &fragStats = getStats(idx).fragStats;
		if(max_probability_size > 50){
			//dis.reserve(max_len);
			DR_bp_distribution.resize(max_probability_size, 0);
			for(int i = 1; i < max_probability_size; i++){
				int index_count = fragStats.getSizeCount(i + 2*read_len);
				float PI = (float)index_count/totalPairedReadCount;
				float PI_pos = PI/i;
				for(int j = 0; j < i; j++ ){
					DR_bp_distribution[j] += PI_pos;
				}
			}
			float sum_p = 0;
			for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++)
				sum_p += DR_bp_distribution[i];
			float levelup_rate = 1/sum_p;//sum will be 1
			for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++)
				DR_bp_distribution[i] *= levelup_rate;
		}
		else{
			max_probability_size = 50;
			DR_bp_distribution.resize(max_probability_size, 0.02);
		}
		//for Soft/hard clip signals break point distribution
		int max_probability_size_SH = 10;
		SH_bp_distribution.resize(max_probability_size_SH, 0.1);
		//for UM start position distribution:
		int min_len = getStats(idx).fragStats.quantile(0.03f);
			max_len = getStats(idx).fragStats.quantile(0.97f);
		int min_um_distribution = min_len -  read_len;
		int max_um_distribution = max_len -  read_len;
		int um_distribution_size = max_um_distribution - min_um_distribution;
		START_OFFSET_UM = min_um_distribution;
		UM_stPos_distribution.resize(um_distribution_size, 0);
		for(int i = 0; i < um_distribution_size; i++){
			int index_count = fragStats.getSizeCount(i + min_len);
			float PI = (float)index_count/totalPairedReadCount;
			UM_stPos_distribution[i] = PI;
		}
	}

	/// merge in the contents of another stats set object:
	void merge(const StatsManager &rhs);
	void save(const char *filename) const;
	void load(const char *filename){ //todo:::


	}
	bool isEmpty() {
		return _group.empty();
	}

	void handleBamCramStats(const char *alignmentFilename, double *ave_depth);

	struct ChromInfo {
		ChromInfo(int32_t chr_ID_, int32_t size_) :
				chr_ID(chr_ID_), size(size_) {
		}
		int32_t chr_ID;
		int32_t size;
		int32_t highestPos = size*0.2;
	};

private:
	void clear() {
		_group.clear();
	}

	id_map<KeyType, UniqueStats> _group;

	const char * referenceFilename;
	std::string defaultStatsFilename;
};

