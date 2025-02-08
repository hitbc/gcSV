/*
 * bam_stat.hpp
 *
 *  Created on: 2021年12月28日
 *      Author: fenghe
 */

#ifndef BAM_STAT_HPP_
#define BAM_STAT_HPP_

#include "../cpp_lib/statistics/StatsManager.hpp"
#include "../cpp_lib/JsonObject/CJsonObject.hpp"

struct BAM_STATUS{
	//global analysis variables
	//analysis results before getting all signals by sampling
	uint32_t minInsertLen;
	uint32_t middleInsertLen;
	uint32_t maxInsertLen;
	std::vector<float> isize_distribution;//from minInsertLen to maxInsertLen
	int analysis_read_length;
	double ave_read_len;
	double ave_read_depth;

	void json_load(const char * json_fn){
		//check file exist:


		neb::CJsonObject input_j;
		neb::load_from_file_to_json(json_fn, input_j);

		input_j.Get("minInsertLen", minInsertLen);
		input_j.Get("middleInsertLen", middleInsertLen);
		input_j.Get("maxInsertLen", maxInsertLen);

		isize_distribution.resize(input_j["isize_distribution"].GetArraySize());

		for (uint i = 0; i < isize_distribution.size(); ++i)
			input_j["isize_distribution"].Get(i, isize_distribution[i]);

		input_j.Get("minInsertLen", minInsertLen);
		input_j.Get("minInsertLen", minInsertLen);
		input_j.Get("minInsertLen", minInsertLen);

		input_j.Get("analysis_read_length", analysis_read_length);
		input_j.Get("ave_read_len", ave_read_len);
		input_j.Get("ave_read_depth", ave_read_depth);
	}

	//generate statistics and dump data to json
	void GEN_STAT_and_json_dump(const char * referenceFilename, const char * bamFile){
		generate_state(referenceFilename, bamFile, false);
		neb::CJsonObject output_j;

		output_j.Add("minInsertLen", minInsertLen);
		output_j.Add("middleInsertLen", middleInsertLen);
		output_j.Add("maxInsertLen", maxInsertLen);

		output_j.AddEmptySubArray("isize_distribution");
		for(uint i = 0; i < isize_distribution.size(); i++)
			output_j["isize_distribution"].Add(isize_distribution[i]);

		output_j.Add("analysis_read_length", analysis_read_length);
		output_j.Add("ave_read_len", ave_read_len);
		output_j.Add("ave_read_depth", ave_read_depth);

		std::cout << output_j.ToFormattedString();
	}

	void generate_read_length(const char * referenceFilename, const char * bamFile){
		//S1 :
		int max_ngz_read_len = 1000; //1000
		Bam_file c_b;
		memset(&c_b, 0, sizeof(Bam_file));
		bam_file_open(bamFile, referenceFilename, NULL, &c_b);

		bam_hdr_t* hdr = c_b._hdr;
		int total_read_number = 0;
		bam1_t b1 = {0};//BAM record for the first read in a pair
		int sam_rst1 = 0;
		uint64_t *read_length_analysis = (uint64_t *)xcalloc(max_ngz_read_len, sizeof(uint64_t));
		while (1){
			//load SAM 1 & 2
			do{
				sam_rst1 = sam_read1(c_b._hfp, hdr, &b1);
			} while( (sam_rst1 >= 0) && (bam_is_secondary(&b1) || bam_is_supplementary(&b1) || bam_is_duplicate(&b1)));
			if(sam_rst1 < 0)		break;
			total_read_number ++;
			if(total_read_number == 100000)
				break;
			//global analysis:
			if(b1.core.l_qseq < max_ngz_read_len)
				read_length_analysis[b1.core.l_qseq]++;
		}
		bam_file_close(&c_b);

		//output analysis results
		//part1: get normal read length
		analysis_read_length = -1;
		double total_read_len = 0;
		for(int i = 0; i < max_ngz_read_len; i++){
			total_read_len += i*read_length_analysis[i];
			if(read_length_analysis[i] > 0.6*total_read_number){
				analysis_read_length = i; break;
			}
		}
		ave_read_len = total_read_len/total_read_number;
		if(analysis_read_length == -1)
			analysis_read_length = ave_read_len;
		free(read_length_analysis);
	}

	void generate_isize_distribution(const char * referenceFilename, const char * bamFile){
		StatsManager rstats(referenceFilename, "");
		rstats.handleBamCramStats(bamFile, &ave_read_depth);
		minInsertLen = rstats.getInsertLen(bamFile, 0.01f);
		middleInsertLen = rstats.getInsertLen(bamFile, 0.5f);
		maxInsertLen = rstats.getInsertLen(bamFile, 0.99f);

		unsigned idx = rstats.getGroupIndex(StatLabel(bamFile, ""));
		int totalPairedReadCount = rstats.getStats(idx).readCounter.totalHighConfidenceReadPairCount() + 1;
		const SizeDistribution &fragStats = rstats.getStats(idx).fragStats;
		for(uint32_t i = minInsertLen; i < maxInsertLen; i++){
			int index_count = fragStats.getSizeCount(i);
			isize_distribution.emplace_back((float)index_count/totalPairedReadCount);
		}
	}

	void generate_state(const char * referenceFilename, const char * bamFile, bool will_output){
		fprintf(stderr, "Generating BAM status...\n");
		generate_read_length(referenceFilename, bamFile);
		//get stats: M2
		generate_isize_distribution(referenceFilename, bamFile);
		if(will_output)
			fprintf(stderr, "BAM/CRAM status: read length: [Normal: %d, AVE: %f] ISIZE: [MIN: %d MIDDLE:%d MAX: %d] ave_read_depth [%f]\n",
				analysis_read_length, ave_read_len,
				minInsertLen, middleInsertLen, maxInsertLen, ave_read_depth);
		if(will_output){
			int isize_idx = minInsertLen;
			fprintf(stderr, "ISIZE distribution\n");
			for(float & isize_rate: isize_distribution)
				fprintf(stderr, "ISIZE %d, rate %f%%\n", isize_idx++, isize_rate*100);
		}
	}

	void getBreakPoint_Distribution( std::vector<float> &DR_bp_distribution, std::vector<float> &SH_bp_distribution,
			std::vector<float> &UM_stPos_distribution, int &START_OFFSET_UM){
		//for discordant /hard clip signals break point distribution
		int min_probability_size = minInsertLen -  2*analysis_read_length; if(min_probability_size < 1) min_probability_size = 1;
		int max_probability_size = maxInsertLen -  2*analysis_read_length;
		DR_bp_distribution.clear();
		DR_bp_distribution.resize(max_probability_size, 0);

		//for DR signals
		for(int i = min_probability_size; i < max_probability_size; i++){
			float PI_per_position = ( isize_distribution[i + 2*analysis_read_length - minInsertLen])/i;
			for(int j = 0; j < i; j++ ) DR_bp_distribution[j] += PI_per_position;
		}
		float sum_p = 0; for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++) sum_p += DR_bp_distribution[i];
		float levelup_rate = 1/sum_p; for(unsigned int i = 0; i < DR_bp_distribution.size(); i ++) DR_bp_distribution[i] *= levelup_rate; //sum will be 1

		//for SR signals
		//for Soft/hard clip signals break point distribution
		int max_probability_size_SH = 10;
		SH_bp_distribution.resize(max_probability_size_SH, 0.0);
		float SH_bp_distribution_sum = 0;
		for(int i = 0; i < max_probability_size_SH; i++){
			SH_bp_distribution[i] = (max_probability_size_SH - i) * (max_probability_size_SH - i);
			SH_bp_distribution_sum += SH_bp_distribution[i];
		}
		for(int i = 0; i < max_probability_size_SH; i++) SH_bp_distribution[i] /= SH_bp_distribution_sum;
		//for UM start position distribution:
		START_OFFSET_UM = minInsertLen -  analysis_read_length;
		UM_stPos_distribution.resize(maxInsertLen - minInsertLen, 0);
		for(uint32_t i = 0; i < maxInsertLen - minInsertLen; i++)
			UM_stPos_distribution[i] = isize_distribution[i];

		//debug show distribution:
		if(0){
			fprintf(stderr, "DR_bp_distribution\n");
			for(float d: DR_bp_distribution)
				fprintf(stderr, "%f\n", d);

			fprintf(stderr, "SH_bp_distribution\n");
			for(float d: SH_bp_distribution)
				fprintf(stderr, "%f\n", d);

			fprintf(stderr, "UM_stPos_distribution\n");
			for(float d: UM_stPos_distribution)
				fprintf(stderr, "%f\n", d);
		}
	}
};



#endif /* BAM_STAT_HPP_ */
