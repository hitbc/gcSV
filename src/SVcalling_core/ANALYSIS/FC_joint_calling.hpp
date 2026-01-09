/*
 * FC_joint_calling.hpp
 *
 *  Created on: 2025-8-14
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_ANALYSIS_FC_JOINT_CALLING_HPP_
#define SVCALLING_CORE_ANALYSIS_FC_JOINT_CALLING_HPP_

#include <set>
#include <fstream>
#include <filesystem>

#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <tuple>
#include <unordered_set>



namespace fs = std::filesystem;

#include "../SVcalling_core/Contig_aligner.hpp"

#define BASE_IS_VNTR 0
#define BASE_IS_NOT_VNTR 1

#define BASE_EDGE_LEN 10

struct FC_joint_calling_handler
{

	static void masker_init(std::vector<uint8_t> &mask, int len)
	{
		mask.resize(len);
		for (int i = 0; i < len; i++)
		{
			mask[i] = BASE_IS_NOT_VNTR;
		}
	}
	static int region_overlap(int base_pos, int st1, int ed1, int st2, int ed2, int &st_r, int &ed_r, std::vector<uint8_t> &mask)
	{
		st_r = MAX(st1, st2);
		ed_r = MIN(ed1, ed2);
		for (int i = st_r; i < ed_r; i++)
			mask[i - base_pos] = BASE_IS_VNTR;
		if (st_r >= ed_r)
			return 0;
		else
			return ed_r - st_r;
	}

	static bool region_overlap_simple(int st1, int ed1, int st2, int ed2, int &st_r, int &ed_r)
	{
		st_r = MAX(st1, st2);
		ed_r = MIN(ed1, ed2);
		if (st_r >= ed_r)
			return false;
		else
			return true;
	}

	static int pos_2_region_distance(int pos, int region_st, int region_ed)
	{
		if (pos >= region_st && pos <= region_ed)
			return 0;
		int dis1 = ABS_U(pos, region_st);
		int dis2 = ABS_U(pos, region_ed);
		return MIN(dis1, dis2);
	}

	static int count_mask_size(std::vector<uint8_t> &mask)
	{
		int mask_size = 0;
		for (uint i = 0; i < mask.size(); i++)
		{
			if (mask[i] == BASE_IS_VNTR)
			{
				mask_size++;
			}
		}
		return mask_size;
	}

	static int count_zero_max_size(std::vector<uint8_t> &mask)
	{
		int max_length = 0;
		int current_length = 0;
		
		for (uint i = 0; i < mask.size(); i++)
		{
			if (mask[i] == BASE_IS_NOT_VNTR)
			{
				current_length++;
				max_length = std::max(max_length, current_length);
			}
			else
			{
				current_length = 0;
			}
		}
		return max_length;
	}

	// RM record
	struct RM_record
	{
		int sw_score;
		float perc_div;
		float perc_del;
		float perc_ins;
		// query_sequence name
		std::string sample_name;
		std::string fasta_region;
		int rm_q_st;
		int rm_q_ed;
		int position_in_query_distance_2_query_end; //{left}
		char match_mode;
		std::string matching_repeat;
		std::string repeat_class_family;
		int r_st;
		int r_ed;
		int position_in_repeat_distance_2_query_end; //{left}
		int RM_ID;
		std::string ori_string;
		void load_one_record(std::string &RM_i, std::vector<std::string> &tmp)
		{
			ori_string = RM_i;
			split_string_SKIP_blank(tmp, RM_i.c_str(), " ");
			if (tmp.size() < 15)
				return;
			sw_score = (int)atoi(tmp[0].c_str());
			perc_div = (float)atof(tmp[1].c_str());
			perc_del = (float)atof(tmp[2].c_str());
			perc_ins = (float)atof(tmp[3].c_str());
			sample_name = tmp[4];						// copy
			rm_q_st = (int)atoi(tmp[5].c_str());			// P_query:ST
			rm_q_ed = (int)atoi(tmp[6].c_str());			// P_query:ED
			position_in_query_distance_2_query_end = 0; //(int)atoi(item_value1[7].c_str());;//{left}
			match_mode = tmp[8][0];
			matching_repeat = tmp[9];
			repeat_class_family = tmp[10];
			r_st = (int)atoi(tmp[11].c_str());			 // position_in_repeat_bg
			r_ed = (int)atoi(tmp[12].c_str());			 // position_in_repeat_ed
			position_in_repeat_distance_2_query_end = 0; //???{left}
			RM_ID = (int)atoi(tmp[14].c_str());
		}
		std::string log_main()
		{
			std::ostringstream oss;
			oss << sample_name << " P_query: [" << rm_q_st << "-" << rm_q_ed << ", len " << (rm_q_ed - rm_q_st)
				<< "] mode " << match_mode << " REP:[" << matching_repeat << " " << repeat_class_family
				<< "] perc_div " << perc_div << "\n";
			return oss.str();
		}

		int query_overlap_len(int base_pos, int st, int ed, std::vector<uint8_t> &mask)
		{
			int st_r = 0;
			int ed_r = 0;
			return region_overlap(base_pos, st, ed, rm_q_st, rm_q_ed, st_r, ed_r, mask);
		}

		bool is_simple_repeat_or_satellite() {
			return repeat_class_family.find("Simple_repeat") != std::string::npos ||
				repeat_class_family.find("Satellite") != std::string::npos;
		}

		// Add this function to check if repeat_class_family contains SINE, LINE, or LTR
		bool isMobileElement() const
		{
			return (repeat_class_family.find("SINE") != std::string::npos ||
					repeat_class_family.find("LINE") != std::string::npos ||
					repeat_class_family.find("Retroposon") != std::string::npos ||
					repeat_class_family.find("DNA") != std::string::npos ||
					repeat_class_family.find("RC") != std::string::npos ||
					repeat_class_family.find("LTR") != std::string::npos);
		}

	};

	// TRF record
	struct TRF_record
	{
		std::string sample_name;
		int q_st;
		int q_ed;
		int trf_q_st_ext;
		int trf_q_ed_ext;
		float repeat_size;
		float nb_copies;
		float motif_size;
		float percMatches;
		float percIndels;
		float alignmentScore;
		float percA;
		float percC;
		float percG;
		float percT;
		float entropy;
		std::string motifSeq;
		std::string trSeq;
		std::string noTr5;
		std::string noTr3;
		std::string ori_string;
		int load_one_record(std::string &sample_name, std::string &TRF_i, std::vector<std::string> &tmp)
		{
			ori_string = TRF_i;
			this->sample_name = sample_name;
			split_string(tmp, TRF_i.c_str(), " ");
			if (tmp.size() < 17)
				return 0 ;
			q_st = (int)atoi(tmp[0].c_str()) - 1;
			q_ed = (int)atoi(tmp[1].c_str()) - 1;

			trf_q_st_ext = q_st - BASE_EDGE_LEN;
			trf_q_st_ext = MAX(0, trf_q_st_ext);
			trf_q_ed_ext = q_ed + BASE_EDGE_LEN;

			repeat_size = (float)atof(tmp[2].c_str());
			nb_copies = (float)atof(tmp[3].c_str());
			motif_size = (float)atof(tmp[4].c_str());
			percMatches = (float)atof(tmp[5].c_str()); // Percent of matches between adjacent copies overall.
			percIndels = (float)atof(tmp[6].c_str());
			alignmentScore = (float)atof(tmp[7].c_str());
			percA = (float)atof(tmp[8].c_str());
			percC = (float)atof(tmp[9].c_str());
			percG = (float)atof(tmp[10].c_str());
			percT = (float)atof(tmp[11].c_str());
			entropy = (float)atof(tmp[12].c_str());
			motifSeq = tmp[13];
			trSeq = tmp[14];
			noTr5 = tmp[15];
			noTr3 = tmp[16];
			return q_ed - q_st;
		}
		std::string log_main()
		{
			std::ostringstream oss;
			oss << sample_name << " R_region: [" << trf_q_st_ext << "-" << trf_q_ed_ext
				<< ", len " << (trf_q_ed_ext - trf_q_st_ext) << "] percMatches " << percMatches
				<< ", motifSeq " << motifSeq << " trSeq " << trSeq << "\n";
			return oss.str();
		}

		int query_overlap_len(int base_pos, int st, int ed, std::vector<uint8_t> &mask)
		{
			int st_r = 0;
			int ed_r = 0;
			return region_overlap(base_pos, st, ed, trf_q_st_ext, trf_q_ed_ext, st_r, ed_r, mask);
		}

		int pos_2_VNTR_distance(int pos)
		{
			return pos_2_region_distance(pos, q_st, q_ed);
		}
	};

	struct VNTR_Range
	{
		int id;
		int start;
		int end;
		bool is_vntr;

		bool with_in_region(int pos) { return (pos >= start && pos <= end); }
		bool is_region_overlap_the_SV(bool sv_is_insertion_viewed_by_cur_contig, int base_BP_core, int insert_bg, int insert_ed)
		{
			if (sv_is_insertion_viewed_by_cur_contig)
			{
				if (with_in_region(insert_bg))
					return true;
				if (with_in_region(insert_ed))
					return true;
			}
			else
			{
				if (with_in_region(base_BP_core))
					return true;
			}
			return false;
		}
	};

	struct CONTIG_COMPARE_RST_ITEM
	{
		// align results
		int a_mapped_len, b_mapped_len, match_base, gap_INS, gap_DEL;
		// summary result
		bool with_NO_VNTR_SV;
	};

	// contig list:
	struct CONTIG_INFO
	{
		std::string sample_name;
		int sample_contig_id;
		int FULL_COVER_R_N;
		int CLIP_AT_LEFT;
		int CLIP_AT_RIGHT;
		std::string contig;
		std::vector<uint8_t> contig_bin;
		std::vector<RM_record> rm_l;
		std::vector<TRF_record> trf_l;
		std::vector<uint8_t> trf_mask; // set to 1 when NOT mask by TRF annotation, set to 0 when masked by TRF
		std::vector<uint8_t> rm_mask;  // repeat masker mask, same as above;
		std::vector<VNTR_Range> vntr_region_l;
		std::vector<VNTR_Range> non_vntr_region_l;
		Contig_String_aligner mapping_2_ref; //
		CONTIG_COMPARE_RST_ITEM c_r;		 //
		bool right_aln_2_ref;

		void get_contig_anno_line(std::ostringstream &ss){
			
			ss << sample_name << "\t";
			ss << "VNTR-ANNO:" << "\t";
			
			for (size_t i = 0; i < vntr_region_l.size(); ++i) {
				if (i > 0) ss << ";";
				const VNTR_Range& range = vntr_region_l[i];
				ss << range.start << "," << range.end << "," << "VNTR";
			}
			ss << "\tNON-VNTR-ANNO:" << "\t";
			for (size_t i = 0; i < non_vntr_region_l.size(); ++i) {
				if (i > 0) ss << ";";
				const VNTR_Range& range = non_vntr_region_l[i];
				ss << range.start << "," << range.end << "," << "NON_VNTR";
			}

			//RM
			ss << "\tRM-ANNO:" << "\t";
			for (size_t i = 0; i < rm_l.size(); ++i) {
				if (i > 0) ss << ";";
				const RM_record& rm = rm_l[i];
				ss << rm.rm_q_st << "," << rm.rm_q_ed << "," << rm.repeat_class_family;
			}
			ss << "\tEND";
			ss << "\n";
		}

		void clear() {
			sample_name.clear();
			sample_contig_id = 0;
			FULL_COVER_R_N = 0;
			CLIP_AT_LEFT = 0;
			CLIP_AT_RIGHT = 0;
			contig.clear();
			contig_bin.clear();
			rm_l.clear();
			trf_l.clear();
			trf_mask.clear();
			rm_mask.clear();
			vntr_region_l.clear();
			non_vntr_region_l.clear();
			mapping_2_ref.init(); // Reinitialize the aligner
			memset(&c_r, 0, sizeof(CONTIG_COMPARE_RST_ITEM));
			right_aln_2_ref = false;
    	}
		std::string show_RM_annotation(std::string &header)
		{
			std::ostringstream oss;
			oss << "[" << header << "] " << sample_name << " RM anno: rm_l.size() " << rm_l.size() << "\n";
			for (uint i = 0; i < rm_l.size(); i++)
				oss << "[" << header << "] " << rm_l[i].log_main();
			return oss.str();
		}

		std::string show_TRF_annotation(std::string &header)
		{
			std::ostringstream oss;
			oss << "[" << header << "] " << sample_name << " TRF anno trf_l.size() " << trf_l.size() << "\n";
			for (uint i = 0; i < trf_l.size(); i++)
				oss << "[" << header << "] " << trf_l[i].log_main();
			return oss.str();
		}

		void store_the_RM_record(RM_record &r_s)
		{
			rm_l.emplace_back(r_s);
		}
		void store_the_TRF_record(TRF_record &t_s)
		{
			t_s.trf_q_ed_ext = MIN(t_s.trf_q_ed_ext, (int)contig.size());
			trf_l.emplace_back(t_s);
		}

		int search_RM(int base_pos, int st_pos, int ed_pos, std::vector<int> &overlap_list, int &max_idx)
		{
			masker_init(rm_mask, ed_pos - st_pos);
			max_idx = -1;
			overlap_list.clear();
			int max_match_len = 0;
			for (uint i = 0; i < rm_l.size(); i++)
			{
				int cur_overlap = rm_l[i].query_overlap_len(base_pos, st_pos, ed_pos, rm_mask);
				if (cur_overlap > 0)
					overlap_list.emplace_back(i);
				if (cur_overlap > max_match_len)
				{
					max_idx = i;
					max_match_len = cur_overlap;
				}
			}
			return count_mask_size(rm_mask);
		}

		int search_RM_MEI(int base_pos, int st_pos, int ed_pos, std::vector<int> &overlap_list, int &max_idx)
		{
			std::vector<uint8_t> mei_mask;
			masker_init(mei_mask, ed_pos - st_pos);
			max_idx = -1;
			overlap_list.clear();
			int max_match_len = 0;
			for (uint i = 0; i < rm_l.size(); i++)
			{
				if(rm_l[i].isMobileElement() == false){
					continue;
				}
				int cur_overlap = rm_l[i].query_overlap_len(base_pos, st_pos, ed_pos, mei_mask);
				if (cur_overlap > 0)
					overlap_list.emplace_back(i);
				if (cur_overlap > max_match_len)
				{
					max_idx = i;
					max_match_len = cur_overlap;
				}
			}
			return count_mask_size(mei_mask);
		}

		// search the region in contig whether the MEI is inside?
		bool search_REGION_WITH_MEI(int st_pos, int ed_pos, std::string &base_contig_sample_name)
		{
			std::ostringstream oss;
			oss << "[BEGIN]-------------------------------search_REGION_WITH_MEI----------------------------------------------------------\n";
			bool with_MEI = false;
			std::string repeat_class_family_cat;
			for (uint i = 0; i < rm_l.size(); i++)
			{
				int st_r = 0;
				int ed_r = 0;
				bool is_overlap = region_overlap_simple(rm_l[i].rm_q_st, rm_l[i].rm_q_ed, st_pos, ed_pos, st_r, ed_r);
				if (is_overlap)
				{
					// check if the annotation is MEI?
					oss << "RM anno in " << sample_name << ", ANNO_region_in_contig [" << st_r << "-" << ed_r << " len " << (ed_r - st_r) << "]\n";
					std::string MEI_STR = contig.substr(st_r, ed_r - st_r);
					oss << "THE ANNO string is " << MEI_STR << " " << rm_l[i].log_main() << "\n";
					repeat_class_family_cat += rm_l[i].repeat_class_family;
					repeat_class_family_cat += "&";
					oss << "\n";
					with_MEI = true;
				}
			}
			oss << "[search_REGION_WITH_MEI-SUMMARY];with_MEI=" << with_MEI 
			<< ";INS_region_in_contig=[" << st_pos << "-" << ed_pos << "];" << 
			"insert_contig=" << sample_name << ";base_contig=" << base_contig_sample_name << ";repeat_class_family_cat=" << repeat_class_family_cat << "\n";
			oss << "[END]-------------------------------search_REGION_WITH_MEI----------------------------------------------------------\n";
			fprintf(stdout, "%s", oss.str().c_str());
			return with_MEI;
		}

		int search_TRF(int base_pos, int st_pos, int ed_pos, std::vector<int> &overlap_list, int &max_idx)
		{
			masker_init(trf_mask, ed_pos - st_pos);
			max_idx = -1;
			int max_match_len = 0;
			for (uint i = 0; i < trf_l.size(); i++)
			{
				int cur_overlap = trf_l[i].query_overlap_len(base_pos, st_pos, ed_pos, trf_mask);
				if (cur_overlap > 0)
					overlap_list.emplace_back(i);
				if (cur_overlap > max_match_len)
				{
					max_idx = i;
					max_match_len = cur_overlap;
				}
			}
			return count_mask_size(trf_mask);
		}

		int get_insert_seq_VNTR_length(int base_pos, int st_pos, int ed_pos, int & non_VNTR_max_size)
		{
			std::vector<uint8_t> tmp_mask;
			masker_init(tmp_mask, ed_pos - st_pos);
			
			// 第一步：添加TRF数据到mask
			for (uint i = 0; i < trf_l.size(); i++)
				trf_l[i].query_overlap_len(base_pos, st_pos, ed_pos, tmp_mask);
			
			// 第二步：添加RM数据到mask（只添加简单重复或卫星重复）
			for (uint i = 0; i < rm_l.size(); i++){
				if(rm_l[i].is_simple_repeat_or_satellite())
					rm_l[i].query_overlap_len(base_pos, st_pos, ed_pos, tmp_mask);
			}
			if(true){
				// 第三步：处理MEI区域，将符合条件的简单重复MEI重新标记为非V碱基
				for (uint i = 0; i < rm_l.size(); i++) {
					if (rm_l[i].isMobileElement() && !rm_l[i].is_simple_repeat_or_satellite()) {
						// 检查该MEI区域是否被TRF覆盖80%以上且周期小于6
						if (false == isSimpleRepeatMEI(rm_l[i], base_pos)) {
							// 将该MEI区域重新标记为非V碱基
							unmarkMEIRegion(rm_l[i], base_pos, st_pos, ed_pos, tmp_mask);
						}
					}
				}
			}

			non_VNTR_max_size = count_zero_max_size(tmp_mask); 
			return count_mask_size(tmp_mask); 
		}

		// 辅助函数：判断是否为简单重复MEI
		bool isSimpleRepeatMEI(const RM_record& rm, int base_pos) const {
			int mei_start = rm.rm_q_st - base_pos;
			int mei_end = rm.rm_q_ed - base_pos;
			int mei_length = mei_end - mei_start;
			
			if (mei_length <= 0) return false;
			
			int trf_covered_length = 0;
			
			// 遍历所有TRF记录，计算覆盖该MEI区域的TRF总长度
			for (const auto& trf : trf_l) {
				int trf_start = trf.trf_q_st_ext - base_pos;
				int trf_end = trf.trf_q_ed_ext - base_pos;
				
				// 计算TRF与MEI的重叠区域
				int overlap_start = std::max(mei_start, trf_start);
				int overlap_end = std::min(mei_end, trf_end);
				int overlap_length = std::max(0, overlap_end - overlap_start);
				
				// 检查TRF周期是否小于6
				if (overlap_length > 0 && trf.motif_size < 6) {
					trf_covered_length += overlap_length;
				}
			}
			
			// 计算TRF覆盖比例
			float coverage_ratio = static_cast<float>(trf_covered_length) / mei_length;
			return coverage_ratio >= 0.8f; // 覆盖80%以上
		}
		
		// 辅助函数：将MEI区域重新标记为非V碱基
		void unmarkMEIRegion(const RM_record& rm, int base_pos, int st_pos, int ed_pos, std::vector<uint8_t>& mask) {
			int mei_start = rm.rm_q_st - base_pos;
			int mei_end = rm.rm_q_ed - base_pos;
			
			// 确保在目标区间内
			int unmark_start = std::max(mei_start, 0);
			int unmark_end = std::min(mei_end, ed_pos - st_pos);
			
			// 将该区域重新标记为0（非V碱基）
			for (int i = unmark_start; i < unmark_end; i++) {
				if (i < mask.size()) {
					mask[i] = BASE_IS_NOT_VNTR;
				}
			}
		}

		void findConsecutiveRanges()
		{
			vntr_region_l.clear();
			non_vntr_region_l.clear();
			if (trf_mask.empty())
				return;
			int current = trf_mask[0];
			int start = 0;
			int region_id = 0;
			for (int i = 1; i < (int)trf_mask.size(); ++i)
			{
				if (trf_mask[i] != current)
				{
					if (current == BASE_IS_VNTR)
						vntr_region_l.push_back({region_id++, start, i - 1, true});
					else
						non_vntr_region_l.push_back({region_id++, start, i - 1, false});
					current = trf_mask[i];
					start = i;
				}
			}
			// last region
			if (current == BASE_IS_VNTR)
				vntr_region_l.push_back({region_id++, start, (int)trf_mask.size() - 1, true});
			else
				non_vntr_region_l.push_back({region_id++, start, (int)trf_mask.size() - 1, false});
		}

		void get_VNTR_REGION()
		{
			masker_init(trf_mask, contig_bin.size());
			// add trf data to mask
			for (uint i = 0; i < trf_l.size(); i++)
				trf_l[i].query_overlap_len(0, 0, contig_bin.size(), trf_mask);
			// add rm data to mask
			for (uint i = 0; i < rm_l.size(); i++){
				if(rm_l[i].is_simple_repeat_or_satellite())
					rm_l[i].query_overlap_len(0, 0, contig_bin.size(), trf_mask);
			}	
			// get region list
			findConsecutiveRanges();
		}

		VNTR_Range *find_VNTR_or_non_VNTR_region_for_position(int position)
		{
			for (VNTR_Range &vr : vntr_region_l)
			{
				if (vr.with_in_region(position))
					return &vr;
			}
			for (VNTR_Range &vr : non_vntr_region_l)
			{
				if (vr.with_in_region(position))
					return &vr;
			}
			return NULL;
		}

		int dis_to_nearst_TRF(int pos)
		{
			int min_dis = MAX_int32t;
			for (uint i = 0; i < trf_l.size(); i++)
			{
				int cur_dis = trf_l[i].pos_2_VNTR_distance(pos);
				min_dis = MIN(cur_dis, min_dis);
			}
			return min_dis;
		}

		void store(std::string &sample_name, int FULL_COVER_R_N, int CLIP_AT_LEFT, int CLIP_AT_RIGHT, std::string &contig, int c_i)
		{
			this->sample_name = sample_name;
			this->sample_contig_id = c_i;
			this->FULL_COVER_R_N = FULL_COVER_R_N;
			this->CLIP_AT_LEFT = CLIP_AT_LEFT;
			this->CLIP_AT_RIGHT = CLIP_AT_RIGHT;
			std::swap(contig, this->contig);
			rm_l.clear();
		}

		bool is_High_depth(int min_full)
		{
			return (FULL_COVER_R_N >= min_full || (CLIP_AT_LEFT >= min_full && CLIP_AT_RIGHT >= min_full));
		}

		bool ksw2_aln_2_ref(std::vector<uint8_t> &ref_bin)
		{
			// size similarL
			mapping_2_ref.init();
			mapping_2_ref.reset_paramater(1, 6, 18, 3, 65, 1);
			int ref_len = ref_bin.size();
			int contig_len = contig_bin.size();
			int zdrop = MAX(ref_len, contig_len);
			mapping_2_ref.setZdrop(zdrop, zdrop);
			if(contig_len > 100000){
				right_aln_2_ref = false;
				std::cerr << "WARNING: skip contig, too long" << "ref_len: " << ref_len;
				std::cerr << " contig_len: " << contig_len;
				std::cerr << " zdrop: " << zdrop << std::endl;
				clear();
				return false;
			}
			mapping_2_ref.align_non_splice_normal_size(contig_bin, ref_bin);
			// ca.analysis_cigar(b_mapped_len, a_mapped_len, match_base);
			int contig_aln_len, ref_aln_len, match_base, MIN_GAP_LEN = 50, gap_INS, gap_DEL;
			mapping_2_ref.analysis_cigar_M2(contig_aln_len, ref_aln_len, match_base, MIN_GAP_LEN, gap_INS, gap_DEL);
			right_aln_2_ref = (contig_aln_len == contig_len && ref_aln_len == ref_len);
			return right_aln_2_ref;
		}

		int get_absolute_position_in_ref(int pos,int &mapping_type, int &event_len)
		{
			mapping_type = -1;
			event_len = -1;
			if (false == right_aln_2_ref)
				return -1;
			int ref_pos = mapping_2_ref.convert_query_position_2_target(pos, mapping_type, event_len);
			return ref_pos;
		}
	};

	struct CLUSTER
	{
		int cluster_ID;
		int main_contig_id = -1;
		std::vector<int> contig_id_list;
		void store(int id)
		{
			contig_id_list.emplace_back(id);
		}
		void new_c(int id, int cluster_ID)
		{
			this->cluster_ID = cluster_ID;
			contig_id_list.emplace_back(id);
			main_contig_id = id;
		}

		static inline int cmp_by_r_n(const CLUSTER &a, const CLUSTER &b)
		{
			// var basic
			return a.contig_id_list.size() > b.contig_id_list.size();
		}
		static inline int cmp_by_r_n_p(const CLUSTER *a, const CLUSTER *b)
		{
			// var basic
			return a->contig_id_list.size() > b->contig_id_list.size();
		}
		bool operator>(const CLUSTER &other) const
		{
			if (cluster_ID == 0 || other.cluster_ID == 0)
			{
				return cluster_ID < other.cluster_ID;
			}
			return contig_id_list.size() > other.contig_id_list.size();
		}
	};

	struct REGION_CONTIG
	{
		R_region r;
		std::vector<CONTIG_INFO> cis;

		int a_mapped_len, b_mapped_len, match_base, gap_INS, gap_DEL;
		bool with_NO_VNTR_SV;

		std::string QUERY_STR;
		std::string TARGET_STR;
	};

	void load_contig_list_one_region_fasta(const char *FC_Contig_L, std::vector<CONTIG_INFO> &contig_l)
	{
		std::vector<std::string> contig_dataset;
		load_string_list_from_file(FC_Contig_L, contig_dataset);
		int contig_n = contig_dataset.size() / 2;
		contig_l.clear();
		for (int c_i = 0; c_i < contig_n; c_i++)
		{
			int line_n = 2 * c_i;
			contig_l.emplace_back();
			std::string sample_name = contig_dataset[line_n].substr(1, contig_dataset[line_n].size() - 1);
			contig_l.back().store(sample_name, 99, 0, 0, contig_dataset[line_n + 1], c_i);
		}
	}

	void load_repeat_masker_result(const char *RM_fn, std::vector<RM_record> &RM_record_l)
	{ // repeat masker file names
		// load all strings
		std::vector<std::string> load_buff;
		std::vector<std::string> RM_record_in_f;
		load_string_list_from_file(RM_fn, RM_record_in_f);
		// skip the header line
		for (uint RM_r_idx = 3; RM_r_idx < RM_record_in_f.size(); RM_r_idx++)
		{
			RM_record_l.emplace_back();
			RM_record_l.back().load_one_record(RM_record_in_f[RM_r_idx], load_buff);
		}
	}

	void load_trf_result(const char *TRF_fn, std::vector<TRF_record> &TRF_record_l)
	{ // repeat masker file names
		// load all strings
		std::vector<std::string> load_buff;
		std::vector<std::string> TRF_record_in_f;
		load_string_list_from_file(TRF_fn, TRF_record_in_f);
		std::string sample_name;
		for (uint RM_r_idx = 0; RM_r_idx < TRF_record_in_f.size(); RM_r_idx++)
		{
			if (TRF_record_in_f[RM_r_idx].substr(0, 1) == "@")
			{
				sample_name = TRF_record_in_f[RM_r_idx].substr(1, TRF_record_in_f[RM_r_idx].size() - 1);
				continue;
			}
			TRF_record_l.emplace_back();
			int repeat_len = TRF_record_l.back().load_one_record(sample_name, TRF_record_in_f[RM_r_idx], load_buff);
			//if(repeat_len < 50)
			if(repeat_len < 30)
				TRF_record_l.pop_back();
		}
	}

	void store_RM_TRF_data_to_contig_l(std::vector<CONTIG_INFO> &contig_l, std::vector<RM_record> &RM_record_l, std::vector<TRF_record> &TRF_record_l)
	{
		// build index
		std::unordered_map<std::string, int> contig_name_2_ID;
		for (uint i = 0; i < contig_l.size(); i++)
			contig_name_2_ID[contig_l[i].sample_name] = i;
		// store the RM data
		for (uint i = 0; i < RM_record_l.size(); i++)
		{
			std::unordered_map<std::string, int>::iterator it = contig_name_2_ID.find(RM_record_l[i].sample_name);
			if(it == contig_name_2_ID.end())
				continue;
			//xassert(it != contig_name_2_ID.end(), "");
			contig_l[it->second].store_the_RM_record(RM_record_l[i]);
		}
		// store the TRF data
		for (uint i = 0; i < TRF_record_l.size(); i++)
		{
			std::unordered_map<std::string, int>::iterator it = contig_name_2_ID.find(TRF_record_l[i].sample_name);
			if(it == contig_name_2_ID.end())
				continue;
			//xassert(it != contig_name_2_ID.end(), "");
			contig_l[it->second].store_the_TRF_record(TRF_record_l[i]);
		}
	}

	int k_value_cal(const std::string &s, int MIN_K_size, int MAX_K_size, int K_step_LEN)
	{
		std::set<std::string> k_value_KC_Set;
		bool withloop = false;
		int kmerSizeF = 0;
		for (int cur_k_value = MIN_K_size; cur_k_value <= MAX_K_size; cur_k_value += K_step_LEN)
		{
			k_value_KC_Set.clear();
			const std::string &seq = s;
			const unsigned seqLen = seq.size();
			// track all words from the read, including repetitive k
			int kmer_number = seqLen - cur_k_value + 1;
			if (kmer_number <= 0)
				break;
			for (int j = 0; j < kmer_number; ++j)
			{
				const std::string word(seq.substr(j, cur_k_value));
				if (k_value_KC_Set.find(word) == k_value_KC_Set.end())
				{
					k_value_KC_Set.emplace(word);
				}
				else
				{
					withloop = true;
					kmerSizeF = cur_k_value;
					break;
				}
			}
			if (!withloop)
				break;
		}
		return kmerSizeF;
	}

	bool ksw2_aln(Contig_String_aligner &ca, CONTIG_INFO &a, CONTIG_INFO &b, CONTIG_COMPARE_RST_ITEM &c_r)
	{
		// size similarL
		int a_len = a.contig_bin.size();
		int b_len = b.contig_bin.size();
		int zdrop = MAX(b_len, a_len);
		ca.setZdrop(zdrop, zdrop);
		ca.align_non_splice_normal_size(b.contig_bin, a.contig_bin);
		// ca.analysis_cigar(b_mapped_len, a_mapped_len, match_base);
		ca.analysis_cigar_M2(c_r.b_mapped_len, c_r.a_mapped_len, c_r.match_base, 50, c_r.gap_INS, c_r.gap_DEL);
		bool right_aln = (c_r.a_mapped_len == a_len && c_r.b_mapped_len == b_len);
		return right_aln;
	}

#define SV_TYPE_NOT_KNOWN 0
#define SV_VNTR_INSERT_IN_VNTR 1
#define SV_NOVEL_INSERT_in_VNTR 2
#define SV_VNTR_INSERT_IN_NON_VNTR 3
#define SV_NOVEL_INSERT_IN_NON_VNTR 4

	struct SV_ITEM
	{
		// basic
		std::string insert_STR;
		std::string base_STR;
		int base_bg;
		int base_ed;

		int insert_bg;
		int insert_ed;
		bool is_INS;

		// analysis
		int ANATYPE;

		int insert_VNTR_LEN_SUM;
		int insert_non_VNTR_LEN;
		std::vector<uint8_t> trf_mask; // set to 1 when NOT mask by TRF annotation, set to 0 when masked by TRF
		int insert_overlap_MEI_len;//when over 0, SV overlap MEI

		//insertion_seq_with_non_VNTR_region
		bool insertion_seq_with_is_vntr;
		//position in base-contig
		int base_BP_core;
		bool base_is_VNTR_in_contig;
		int  base_VNTR_idx_in_contig;//when is VNTR, is the VNTR idx; when is NON-VNTR, is the NON-VNTR idx;
		
		//base- absolute position in ref
		int  base_absolute_BP_in_ref;
		int  base_absolute_BP_in_ref_event_len;//for nest SVs, the outter ins len 
		bool base_is_VNTR_in_ref;
		int  base_VNTR_idx_in_ref;//when is VNTR, is the VNTR idx; when is NON-VNTR, is the NON-VNTR idx;
		int distanceToNearestVNTRRegion_ref;
		bool base_is_also_INS_in_ref;
		int contig_3_1_INS_absolute_BP_in_ref;

		//insert - absolute position in ref
		int contig_bg_absolute_BP_in_ref;
		int contig_ed_absolute_BP_in_ref;
		
		std::string insert_trf_motifSeq;

		//insert annotation:
		int insert_trf_max_i;
		std::vector<int> insert_trf_overlap_l;
		int insert_rm_max_i;
		std::vector<int> insert_rm_overlap_l;

		CONTIG_INFO *base_contig;
		CONTIG_INFO *insert_contig;

		int	replacement_Distance;
		float replacement_Len_Diff_Ratio;

		std::string nest_string;

		/**
		 * @brief Checks if the SV type is not a VNTR insertion type
		 * @return true if the SV type is not SV_VNTR_INSERT_IN_VNTR or SV_VNTR_INSERT_IN_NON_VNTR, false otherwise
		 */
		bool isNonVntrInsertType() {
			return ANATYPE != SV_VNTR_INSERT_IN_VNTR;
			//return ANATYPE != SV_VNTR_INSERT_IN_VNTR && ANATYPE != SV_VNTR_INSERT_IN_NON_VNTR;
		}

		void store(bool is_INS, int a_i, int b_i, int c_length, CONTIG_INFO *contig_a, CONTIG_INFO *contig_b)
		{
			insert_STR.clear();
			base_STR.clear();
			this->is_INS = is_INS;
			// the a contig is first, the b is second, when CIGAR==INS, the base is contig a, insert is contig b;
			if(is_INS == true){
				insert_contig = contig_b;
				base_contig = contig_a;
				insert_bg = b_i;
				base_BP_core = a_i;
			}else{
				insert_contig = contig_a;
				base_contig = contig_b;
				insert_bg = a_i;
				base_BP_core = b_i;
			}
			insert_ed = insert_bg + c_length;

			int base_max_len = base_contig->contig_bin.size();
			base_bg = base_BP_core;
			base_bg = MAX(base_bg, 0);
			base_ed = base_BP_core + 1;
			base_ed = MIN(base_ed, base_max_len);
			xassert(base_bg <= base_ed, "");
			xassert(insert_bg <= insert_ed, "");

			xassert(base_ed <= base_contig->contig_bin.size(), "");
			xassert(insert_ed <= insert_contig->contig_bin.size(), "");

			uint8_t *insert_p = &(insert_contig->contig_bin[insert_bg]);
			for (int i = 0; i < c_length; i++)
				insert_STR += "ACGTNNN"[insert_p[i]];
			uint8_t *base_p = &(base_contig->contig_bin[0]);
			for (int i = base_bg; i < base_ed; i++)
				base_STR += "ACGTNNN"[base_p[i]];

			base_absolute_BP_in_ref = -1;
			base_is_also_INS_in_ref = false;
			contig_bg_absolute_BP_in_ref = -1;
			contig_ed_absolute_BP_in_ref = -1;

			insert_overlap_MEI_len = 0;
			contig_3_1_INS_absolute_BP_in_ref = -2;
		}

		void show_insert_STR_and_trf_mask(std::ostringstream &oss, std::string &var_name){
			oss << "[" << var_name << "] ";
			oss << "insert_STR: " << insert_STR << std::endl;
			oss << "[" << var_name << "] ";
			oss << "trf_mask  : ";
			for (size_t i = 0; i < trf_mask.size(); ++i) {
				oss << static_cast<int>(trf_mask[i]);
    		}
			oss << std::endl;
		}

		void show_basic(std::ostringstream &oss)
		{
			oss
				<< "base_contig: " << base_contig->sample_name << ", "
				<< "insert_contig: " << insert_contig->sample_name << ", "
				<< "ANATYPE: " << ANATYPE << ", "
				<< "is_INS: " << is_INS << ", "
				<< "base_R: [" << base_bg << "-" << base_ed << "], "
				<< "base_BP: " << base_BP_core << ", "
				<< "base_len: " << (base_ed - base_bg) << ", "
				<< "insert_R: [" << insert_bg << "-" << insert_ed << "], "
				<< "insert_len: " << (insert_ed - insert_bg) << ", "
				<< "VNTR_LEN: " << insert_VNTR_LEN_SUM << ", "
				<< "NO_VNTR_LEN: " << insert_non_VNTR_LEN << ", "
				<< "insert_STR: " << insert_STR << ", "
				<< "base_STR: " << base_STR << ", "
				<< "base_absolute_BP_in_ref: " << base_absolute_BP_in_ref << ", "
				<< "contig_bg_absolute_BP_in_ref: " << contig_bg_absolute_BP_in_ref << ", "
				<< "contig_ed_absolute_BP_in_ref: " << contig_ed_absolute_BP_in_ref
				<< "\n";
		}

		void show_ANNO_ONLY_OVERLAP_TRF(std::ostringstream &oss, std::string &var_name)
		{
			oss << "[" << var_name << "] ";
			oss << " insert_rm_overlap_l.size(): " << insert_rm_overlap_l.size();
			oss << " insert_trf_overlap_l.size(): " << insert_trf_overlap_l.size()
			    << std::endl;

			for (uint i = 0; i < insert_rm_overlap_l.size(); i++)
			{
				oss << "[" << var_name << "] [insert_rm]";
				oss << insert_contig->rm_l[insert_rm_overlap_l[i]].log_main();
			}
			for (uint i = 0; i < insert_trf_overlap_l.size(); i++)
			{
				oss << "[" << var_name << "] [insert_trf]";
				oss << insert_contig->trf_l[insert_trf_overlap_l[i]].log_main();
			}
			oss << "\n";
		}

		void show_ANNO_ALL(std::ostringstream &oss)
		{
			oss << "[INSERT]";
			if (-1 != insert_trf_max_i)
				oss << insert_contig->trf_l[insert_trf_max_i].log_main();
			else
				oss << "[NO VNTR], sample_name " << insert_contig->sample_name << "\n";
			oss << "\n";
		}

		void log_anno(std::ostringstream &oss, std::string &var_name)
		{
			oss << "[" << var_name << "] ";
			show_basic(oss);
			show_insert_STR_and_trf_mask(oss, var_name);
			show_ANNO_ONLY_OVERLAP_TRF(oss, var_name);
			//printf("%s", oss.str().c_str());
		}

		void log_simple(std::ostringstream &oss, std::string &var_name)
		{
			oss
				<< var_name << "\t"
				<< insert_contig->sample_name << "\t"
				<< "BG_POS_" << insert_bg << "\t"
				<< insert_STR << "\n";
		}
	};

	void store_SVs_from_cigar(Contig_String_aligner &ca, std::vector<SV_ITEM> &sv_l, CONTIG_INFO &a, CONTIG_INFO &b,
							  int MIN_GAP_LEN, int suggest_st_pos)
	{
		// get the insertion or deletion sequence:
		uint32_t *bam_cigar;
		int n_cigar = ca.get_cigar(&bam_cigar);
		int a_i = suggest_st_pos, b_i = 0;
		// for SUM insertions
		for (int cigar_ID = 0; cigar_ID < n_cigar; cigar_ID++)
		{
			int c_length = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch (c_type)
			{
			case 0:
			case 7:
				a_i += c_length;
				b_i += c_length;
				break; // M
			case 8:
				a_i += c_length;
				b_i += c_length;
				break; // X
			case 1:	//insertion
				if (cigar_ID != 0 && cigar_ID != n_cigar - 1 && c_length >= MIN_GAP_LEN)
				{
					sv_l.emplace_back(); // store the SVs
					sv_l.back().store(true, a_i, b_i, c_length, &a, &b);
				}
				b_i += c_length;
				break;
			case 2: // Deletions
				if (cigar_ID != 0 && cigar_ID != n_cigar - 1 && c_length >= MIN_GAP_LEN)
				{
					sv_l.emplace_back(); // store the SVs
					sv_l.back().store(false, a_i, b_i, c_length, &a, &b);
				}
				a_i += c_length;
				break;
			case 3:
			case 4:
				break; // S, print -
			default:
				fprintf(stdout, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}
		}

		// 合并函数：直接统计两个字符串之间共享的k-mer数量
		int countSharedKmers(const std::string &insert_sequence,
							 const std::string &ref_sequence,
							 int k = 10)
		{
			if (insert_sequence.length() < k || ref_sequence.length() < k)
			{
				return 0;
			}

			// 构建ref序列的k-mer索引
			std::unordered_set<std::string> ref_kmers;
			for (size_t i = 0; i <= ref_sequence.length() - k; i++)
			{
				std::string kmer = ref_sequence.substr(i, k);
				ref_kmers.insert(kmer);
			}

			// 统计insert序列中出现在ref中的k-mer数量
			int shared_count = 0;
			for (size_t i = 0; i <= insert_sequence.length() - k; i++)
			{
				std::string kmer = insert_sequence.substr(i, k);
				if (ref_kmers.find(kmer) != ref_kmers.end())
				{
					shared_count++;
				}
			}

			return shared_count;
		}

		/**
		 * Check if the given interval [insert_bg, insert_ed] is fully covered by any TRF record
		 * and meets the specified criteria (repeat_size > 50 and nb_copies <= 2).
		 *
		 * @param trf_l Vector of TRF records to search through
		 * @param insert_bg Start position of the interval to check
		 * @param insert_ed End position of the interval to check
		 * @return true if found a TRF record that fully covers the interval and meets criteria
		 *         false otherwise
		 */
		bool isCovered_BY_DUP(std::vector<TRF_record> &trf_l, int insert_bg, int insert_ed)
		{
			// Iterate through all TRF records
			for (auto &trf : trf_l)
			{
				// Check if current TRF record completely covers the given interval
				if (trf.trf_q_st_ext <= insert_bg && trf.trf_q_ed_ext >= insert_ed)
				{
					// Check if the record meets both criteria:
					// 1. repeat_size > 50 (tandem repeat length)
					// 2. nb_copies <= 2 (number of copies)
					if (trf.repeat_size >= 50 && trf.nb_copies <= 2.5)
					{
						printf("[V-INS-NON_V ---> NON_V-INS-NON_V as DUP]: insert_bg %d , insert_ed %d ; trf %s\n" ,  insert_bg, insert_ed, trf.log_main().c_str());
						return true; // Found a matching record
					}
				}
			}
			return false; // No matching record found
		}

#define sequencing_platform_ONT 0
#define sequencing_platform_HPRC 1

	void anno_type_analysis_for_all_SV(std::vector<SV_ITEM> &sv_l, CONTIG_COMPARE_RST_ITEM &c_r, 
		CONTIG_INFO& ref_contig,	
		bool &with_sv_is_NOT_vntr,
		int sequencing_platform)
	{
		for (uint i = 0; i < sv_l.size(); i++)
		{
			// search for the annotation
			SV_ITEM &sv = sv_l[i];
			// search the annotation
			// the repeat result of base:
			sv.insert_contig->search_TRF(sv.insert_bg, sv.insert_bg, sv.insert_ed, sv.insert_trf_overlap_l, sv.insert_trf_max_i);
			sv.trf_mask = sv.insert_contig->trf_mask;
			sv.insert_contig->search_RM(sv.insert_bg, sv.insert_bg, sv.insert_ed, sv.insert_rm_overlap_l, sv.insert_rm_max_i);
			// get VNTR length
			// sv.insert_VNTR_LEN_SUM = sv.insert_contig->search_TRF(sv.insert_bg, sv.insert_bg, sv.insert_ed, sv.insert_trf_overlap_l, sv.insert_trf_max_i);
			// sv.insert_non_VNTR_LEN = sv.insert_STR.size() - sv.insert_VNTR_LEN_SUM;

			// fowling not applied right now
			sv.insert_VNTR_LEN_SUM = sv.insert_contig->get_insert_seq_VNTR_length(sv.insert_bg, sv.insert_bg, sv.insert_ed, sv.insert_non_VNTR_LEN);
			sv.insert_non_VNTR_LEN = sv.insert_STR.size() - sv.insert_VNTR_LEN_SUM;

			{ // insert_overlap_MEI??
				std::vector<int> overlap_list;
				int max_idx;
				sv.insert_overlap_MEI_len = sv.insert_contig->search_RM_MEI(sv.insert_bg, sv.insert_bg, sv.insert_ed, overlap_list, max_idx);
			}

			bool insertion_seq_with_non_VNTR_region = false;
			if (sequencing_platform_ONT == sequencing_platform)
			{
				insertion_seq_with_non_VNTR_region = ((sv.insert_non_VNTR_LEN > 200) || (sv.insert_non_VNTR_LEN > sv.insert_STR.size() * 0.6));
			}
			else
			{
				insertion_seq_with_non_VNTR_region = ((sv.insert_non_VNTR_LEN > 200) || (sv.insert_non_VNTR_LEN > sv.insert_STR.size() * 0.2));
			}

			sv.insertion_seq_with_is_vntr = (insertion_seq_with_non_VNTR_region == false);

			// set type
			sv.ANATYPE = SV_TYPE_NOT_KNOWN;
			if (sv.base_is_VNTR_in_ref == true && insertion_seq_with_non_VNTR_region == false)
				sv.ANATYPE = SV_VNTR_INSERT_IN_VNTR;
			else if (sv.base_is_VNTR_in_ref == true && insertion_seq_with_non_VNTR_region == true)
				sv.ANATYPE = SV_NOVEL_INSERT_in_VNTR;
			else if (sv.base_is_VNTR_in_ref == false && insertion_seq_with_non_VNTR_region == false)
				sv.ANATYPE = SV_VNTR_INSERT_IN_NON_VNTR;
			else if (sv.base_is_VNTR_in_ref == false && insertion_seq_with_non_VNTR_region == true)
				sv.ANATYPE = SV_NOVEL_INSERT_IN_NON_VNTR;

			// when insertion is VNTR and is inserted in non-VNTR region
			if (sv.ANATYPE == SV_VNTR_INSERT_IN_NON_VNTR)
			{
				std::cout
					<< "SV_VNTR_INSERT_IN_NON_VNTR "
					<< sv.base_contig->sample_name.c_str() 
					<< sv.insert_contig->sample_name.c_str() 
					<< std::endl;

				bool is_nearby_VNTR = (sv.distanceToNearestVNTRRegion_ref > 0 && sv.distanceToNearestVNTRRegion_ref < 30);
				if (is_nearby_VNTR)
				{
					std::string vntr_sequence_insert = sv.insert_STR;
					int nearest_region_index = -1;
					findDistanceToNearestRegion(sv.base_absolute_BP_in_ref, ref_contig.vntr_region_l, nearest_region_index);
					// check if is similar VNTR -string ?
					xassert(nearest_region_index != -1, "");
					auto &cur_v = ref_contig.vntr_region_l[nearest_region_index];
					std::string vntr_sequence_ref = ref_contig.contig.substr(cur_v.start, cur_v.end - cur_v.start + 1);
					int shared_kmer_count = countSharedKmers(vntr_sequence_insert, vntr_sequence_ref, 10);
					if (shared_kmer_count >= 5)
					{
						sv.ANATYPE = SV_VNTR_INSERT_IN_VNTR;
						std::cout
							<< "Change color: " 
							<< ", nearest_region_index=" << nearest_region_index
							<< ", region_start=" << cur_v.start
							<< ", region_end=" << cur_v.end
							<< ", shared_kmer_count=" << shared_kmer_count
							<< ", ANATYPE=" << sv.ANATYPE
							<< ", insert_seq=" << vntr_sequence_insert
							<< ", ref_seq=" << vntr_sequence_ref
							<< std::endl;
					}
				}
			}
			
			// when insertion is still VNTR and is inserted in non-VNTR region
			//change DUP in NON_V region to pure NO_V INS
			if (sv.ANATYPE == SV_VNTR_INSERT_IN_NON_VNTR)
			{
				if (true == isCovered_BY_DUP(sv.insert_contig->trf_l, sv.insert_bg, sv.insert_ed))
					sv.ANATYPE = SV_NOVEL_INSERT_IN_NON_VNTR;
			}

			// when insertion is VNTR
			if (sv.isNonVntrInsertType())
				with_sv_is_NOT_vntr = true;
		}
	}

	void logAllVariants(std::vector<SV_ITEM>& sv_l, 
                   const std::string& CMP_header, 
                   int& var_ID, 
                   std::string& var_log, 
                   std::string& var_ANNO_log) 
	{
		for (uint i = 0; i < sv_l.size(); i++)
		{
			// search for the annotation
			// log for all vars-basic:
			var_ID++;
			std::string var_name = CMP_header + "-" + std::to_string(var_ID);
			{
				std::ostringstream oss;
				sv_l[i].log_simple(oss, var_name);
				var_log += oss.str();
			}
			// log for all vars-RM/TRF:
			{
				// search for the annotation
				std::ostringstream oss;
				sv_l[i].log_anno(oss, var_name);
				var_ANNO_log += oss.str();
			}
		}
	}

	struct absolute_BP_in_ref_ITEM
	{
		int absolute_position_in_ref;
		int relative_position_in_base_contig;

		VNTR_Range *ref_vntr_r_p;
		SV_ITEM *sv;

		std::string get_log()
		{
			std::ostringstream oss;
			oss << "[ref_non_vntr_absolute region id" << ref_vntr_r_p->id << " " << ref_vntr_r_p->start << "~" << ref_vntr_r_p->end << ", "
				<< (ref_vntr_r_p->is_vntr ? "VNTR" : "NON-VNTR") << "] "
				<< "relative_p_in_base " << relative_position_in_base_contig << " absolute_p_in_ref " << absolute_position_in_ref << "\n";
			return oss.str();
		}

		void set(int abs_pos, int rel_pos, VNTR_Range* vntr_range, SV_ITEM* sv_item)
		{
			absolute_position_in_ref = abs_pos;
			relative_position_in_base_contig = rel_pos;
			ref_vntr_r_p = vntr_range;
			sv = sv_item;
		}

	};

	void searchMEIRegionsForNovelInsertions(std::vector<SV_ITEM>& sv_l) {
		for (uint i = 0; i < sv_l.size(); i++)
		{
			if (sv_l[i].ANATYPE == SV_NOVEL_INSERT_in_VNTR)
			{
				sv_l[i].insert_contig->search_REGION_WITH_MEI(sv_l[i].insert_bg, sv_l[i].insert_ed, sv_l[i].base_contig->sample_name);
			}
		}
	}

	std::string generateMainContigLog(const std::string& CMP_header, 
                                 int64_t suggest_st_pos,
                                 Contig_String_aligner& ca, 
                                 CONTIG_INFO& target_contig, 
                                 CONTIG_INFO& query_contig) {  
		std::ostringstream oss;
		oss
			<< CMP_header << "\t"
			<< "suggest_st_pos" << "\t"
			<< std::to_string(suggest_st_pos) << "\t"
			<< ca.printCIGAR_core() << "\t"
			<< target_contig.sample_name << "\t"
			<< target_contig.contig << "\t"
			<< query_contig.sample_name << "\t"
			<< query_contig.contig << "\t"
			<< "\n";
		return oss.str();
	}

	bool findRegionForPosition(int absolute_position, 
                              const std::vector<VNTR_Range>& regions,
                              int& region_index)
    {
        for (uint j = 0; j < regions.size(); j++)
        {
            if (absolute_position >= regions[j].start && 
                absolute_position <= regions[j].end)
            {
                region_index = j;
                return true;
            }
        }
        return false;
    }

	int findDistanceToNearestRegion(int absolute_position, 
									const std::vector<VNTR_Range>& regions,
									int& nearest_region_index)
	{
		if (regions.empty()) {
			nearest_region_index = -1;
			return -1;
		}
		
		int min_distance = std::numeric_limits<int>::max();
		nearest_region_index = -1;
		
		for (int i = 0; i < regions.size(); i++) {
			const auto& region = regions[i];
			
			if (absolute_position >= region.start && absolute_position <= region.end) {
				nearest_region_index = i;
				return 0;
			}
			
			// 计算到区域边界的最近距离
			int distance_to_start = std::abs(absolute_position - region.start);
			int distance_to_end = std::abs(absolute_position - region.end);
			int current_min = std::min(distance_to_start, distance_to_end);
			
			if (current_min < min_distance) {
				min_distance = current_min;
				nearest_region_index = i;
			}
		}
		
		return min_distance;
	}

	std::vector<std::pair<int, char>> nest_parseCigar(std::string &cigar_string)
	{
		std::vector<std::pair<int, char>> operations;
		std::regex pattern("(\\d+)([MXID])");
		std::smatch matches;

		std::string remaining = cigar_string;
		while (std::regex_search(remaining, matches, pattern))
		{
			int length = std::stoi(matches[1]);
			char op = matches[2].str()[0];
			operations.push_back({length, op});
			remaining = matches.suffix();
		}
		return operations;
	}

	// Flip CIGAR operations: I becomes D, D becomes I, others remain unchanged
	std::vector<std::pair<int, char>>  nest_flipCigar(std::string &cigar_string)
	{
		auto operations = nest_parseCigar(cigar_string);
		std::vector<std::pair<int, char>> flipped_operations;

		for (const auto &op : operations)
		{
			int length = op.first;
			char operation = op.second;
			if (operation == 'I')
				flipped_operations.push_back({length, 'D'});
			else if (operation == 'D')
				flipped_operations.push_back({length, 'I'});
			else
				flipped_operations.push_back({length, operation});
		}

		return flipped_operations;
	}

	void count_match_base(int left_boundary, int right_boundary, std::vector<std::pair<int, char>> &genome_sequence, int &indel_count, int &mismatch_count, int &match_count_region)
	{
		bool in_insertion = false; 
		bool in_deletion = false;  

		for (int i = left_boundary; i <= right_boundary; i++)
		{
			char op = genome_sequence[i].second;

			if (op == 'I')
			{
				if (!in_insertion)
				{
					indel_count++; 
					in_insertion = true;
				}
				in_deletion = false;
			}
			else if (op == 'D')
			{
				if (!in_deletion)
				{
					indel_count++;
					in_deletion = true;
				}
				in_insertion = false; 
			}
			else
			{
				in_insertion = false;
				in_deletion = false;  
				if (op == 'X')
				{
					mismatch_count++;
				}
				else if (op == 'M')
				{
					match_count_region++;
				}
			}
		}
	}
	std::tuple<int, int, bool, double, int, int, int> nest_calculateMismatchRegion(std::string &cigar_string, int pos, bool pos_is_in_ref)
	{
		std::vector<std::pair<int, char>> operations;
		if (pos_is_in_ref)
			operations = nest_parseCigar(cigar_string);
		else
			operations =  nest_flipCigar(cigar_string);

		std::vector<std::pair<int, char>> genome_sequence;
		int current_pos = 0;

		for (const auto &op : operations)
		{
			int length = op.first;
			char operation = op.second;

			for (int i = 0; i < length; i++)
			{
				genome_sequence.push_back({current_pos, operation});
				if (operation == 'M' || operation == 'X' || operation == 'D')
					current_pos++;
			}
		}

		int pos_index = -1;
		for (size_t i = 0; i < genome_sequence.size(); i++)
		{
			if (genome_sequence[i].first == pos)
			{
				pos_index = i;
				break;
			}
		}
		if (pos_index == -1)
		{
			throw std::runtime_error("POS not found in the alignment");
		}

		int left_boundary = 0;
		int match_count = 0;
		for (int i = pos_index - 1; i >= 0; i--)
		{
			if (genome_sequence[i].second == 'M')
			{
				match_count++;
				if (match_count > 50)
				{
					left_boundary = i + match_count; 
					break;
				}
			}
			else
			{
				match_count = 0;
			}
		}

		int right_boundary = genome_sequence.size() - 1;
		match_count = 0;
		for (size_t i = pos_index + 1; i < genome_sequence.size(); i++)
		{
			if (genome_sequence[i].second == 'M')
			{
				match_count++;
				if (match_count > 50)
				{
					right_boundary = i - match_count; 
					break;
				}
			}
			else
			{
				match_count = 0;
			}
		}

		int indel_count = 0;
		int mismatch_count = 0;
		int match_count_region = 0;

		count_match_base(left_boundary, right_boundary, genome_sequence, indel_count, mismatch_count, match_count_region);

		double mismatch_rate = 0.0;
		int total_bases = indel_count + mismatch_count + match_count_region;
		if (total_bases > 0)
			mismatch_rate = static_cast<double>(indel_count + mismatch_count) / total_bases;

		int pos_left = genome_sequence[left_boundary].first;
		int pos_right = genome_sequence[right_boundary].first;
		bool many_mismatch = ((pos_right - pos_left) > 100 && mismatch_rate > 0.1);
		
		return std::make_tuple(
			pos_left,  // Start position of the region
			pos_right, // End position of the region
			many_mismatch,      // many_mismatch
			mismatch_rate,                         // Mismatch rate
			indel_count,                           // INDEL count
			mismatch_count,                        // Mismatch base count
			match_count_region                     // Match base count
		);
	}

	
	void storeAbsoluteBPForAllVariants(
		std::vector<SV_ITEM>& sv_l,
		CONTIG_INFO& target_contig,
		CONTIG_INFO& query_contig,
		CONTIG_INFO& ref_contig,
		std::vector<SV_ITEM> &target_contig_to_ref_sv_l,
		std::vector<SV_ITEM> &query_contig_to_ref_sv_l,
		std::string & cigar_str_3_to_2
	)
	{
		// S1:VNTR region list and ~VNTR region list for C1;
		for (uint i = 0; i < sv_l.size(); i++)
		{
			CONTIG_INFO& base_contig = (sv_l[i].is_INS == true)?target_contig:query_contig;
			CONTIG_INFO& insert_contig = (sv_l[i].is_INS == false)?target_contig:query_contig;
			int mapping_type = -1;
			int event_len = -1;
			if(sv_l[i].insert_contig->sample_contig_id == 0){//deletion to reference, the insert contig is reference 
				sv_l[i].base_absolute_BP_in_ref = sv_l[i].insert_bg;
			}
			else
				sv_l[i].base_absolute_BP_in_ref = base_contig.get_absolute_position_in_ref(sv_l[i].base_BP_core, mapping_type, event_len);	
			sv_l[i].base_is_also_INS_in_ref = (mapping_type == 1);// case 1: insertion
			sv_l[i].base_absolute_BP_in_ref_event_len = event_len;
			sv_l[i].contig_bg_absolute_BP_in_ref = insert_contig.get_absolute_position_in_ref(sv_l[i].insert_bg, mapping_type, event_len);
			sv_l[i].contig_ed_absolute_BP_in_ref = insert_contig.get_absolute_position_in_ref(sv_l[i].insert_ed, mapping_type, event_len);

			//for nest sv
			sv_l[i].contig_3_1_INS_absolute_BP_in_ref = -2;
			if(true == sv_l[i].base_is_also_INS_in_ref){
				//search contig2-2-ref SV list:
				int cur_sv_bg = sv_l[i].insert_bg;
				int cur_sv_ed = sv_l[i].insert_ed;
				sv_l[i].contig_3_1_INS_absolute_BP_in_ref = -1;

				std::vector<SV_ITEM> & contig_to_ref_sv_l = (sv_l[i].is_INS == true)?query_contig_to_ref_sv_l:target_contig_to_ref_sv_l;

				for (uint j = 0; j < contig_to_ref_sv_l.size(); j++){
					SV_ITEM &ref_sv = contig_to_ref_sv_l[j];
					if(ref_sv.is_INS == false)
						continue;
					if(std::max(ref_sv.insert_bg, cur_sv_bg) <= std::min(ref_sv.insert_ed, cur_sv_ed)){//overlap
						sv_l[i].contig_3_1_INS_absolute_BP_in_ref = ref_sv.base_bg;
					}
				}

				sv_l[i].nest_string = "NA";
				if(sv_l[i].contig_3_1_INS_absolute_BP_in_ref > 0){
					//count 
					//（1）判定（3到2）嵌套外层变异
					std::tuple<int, int, bool, double, int, int, int> result_3_2 = nest_calculateMismatchRegion(cigar_str_3_to_2, sv_l[i].base_bg, (sv_l[i].is_INS == true));
					bool many_mismatch_3_2 = std::get<2>(result_3_2);
					//（2）判定（2到ref）嵌套内存变异
					std::string cigar_2_ref = base_contig.mapping_2_ref.printCIGAR_core();
					std::tuple<int, int, bool, double, int, int, int> result_2_ref = nest_calculateMismatchRegion(cigar_2_ref, sv_l[i].base_absolute_BP_in_ref , true);
					bool many_mismatch_2_ref = std::get<2>(result_2_ref);
					//（3）判定（3到ref）3-1的变异
					std::string cigar_3_ref = insert_contig.mapping_2_ref.printCIGAR_core();
					std::tuple<int, int, bool, double, int, int, int> result_3_ref = nest_calculateMismatchRegion(cigar_3_ref, sv_l[i].contig_3_1_INS_absolute_BP_in_ref, true);
					bool many_mismatch_3_ref = std::get<2>(result_3_ref);

					//nest:string
					sv_l[i].nest_string = 
					std::string("is_many_mismatch;") + 
					((many_mismatch_3_2 || many_mismatch_2_ref || many_mismatch_2_ref) ? "true" : "false") + ";" +
					(many_mismatch_3_2 ? "true" : "false") + ";" +
					(many_mismatch_2_ref ? "true" : "false") + ";" +
					(many_mismatch_3_ref ? "true" : "false") + ";" +

					"result_3_2;" + 
					std::to_string(std::get<0>(result_3_2)) + ";" +
					std::to_string(std::get<1>(result_3_2)) + ";" +
					std::to_string(std::get<3>(result_3_2)) + ";" +
					"result_2_ref;" + 
					std::to_string(std::get<0>(result_2_ref)) + ";" +
					std::to_string(std::get<1>(result_2_ref)) + ";" +
					std::to_string(std::get<3>(result_2_ref)) + ";" +
					"result_3_ref;" + 
					std::to_string(std::get<0>(result_3_ref)) + ";" +
					std::to_string(std::get<1>(result_3_ref)) + ";" +
					std::to_string(std::get<3>(result_3_ref));
				}
			}

			///get VNTR idx
			///[VNTR_idx]: when is VNTR, is the VNTR idx; when is NON-VNTR, is the NON-VNTR idx;
            //for ref
			bool is_VNTR; int VNTR_idx; int distanceToNearestVNTRRegion;
			get_vntr_idx( is_VNTR, VNTR_idx, sv_l[i].base_absolute_BP_in_ref, sv_l[i], ref_contig, distanceToNearestVNTRRegion);
			sv_l[i].base_is_VNTR_in_ref = is_VNTR;
			sv_l[i].base_VNTR_idx_in_ref = VNTR_idx;
			sv_l[i].distanceToNearestVNTRRegion_ref = distanceToNearestVNTRRegion;

			//for base contig
			get_vntr_idx( is_VNTR, VNTR_idx, sv_l[i].base_BP_core, sv_l[i], *(sv_l[i].base_contig), distanceToNearestVNTRRegion);
			sv_l[i].base_is_VNTR_in_contig = is_VNTR;
			sv_l[i].base_VNTR_idx_in_contig = VNTR_idx;
		}
    }

    void get_vntr_idx(bool &is_VNTR, int &VNTR_idx, int bp_pos, SV_ITEM &sv_c, CONTIG_INFO &cur_contig, int &distanceToNearestVNTRRegion)
    {
		// Determine if the absolute position is in VNTR or non-VNTR region in reference
		xassert( !cur_contig.vntr_region_l.empty() || !cur_contig.non_vntr_region_l.empty() , "");
		is_VNTR = false;
		VNTR_idx = -1;
		distanceToNearestVNTRRegion = 0;
		// Search in VNTR regions first
		is_VNTR = findRegionForPosition(bp_pos, cur_contig.vntr_region_l, VNTR_idx);
		// If not found in VNTR regions, search in non-VNTR regions
		if (!is_VNTR){
			findRegionForPosition(bp_pos, cur_contig.non_vntr_region_l, VNTR_idx);
			int nearest_region_index = -1;
			distanceToNearestVNTRRegion = findDistanceToNearestRegion(bp_pos, cur_contig.vntr_region_l, nearest_region_index);
		}
    }

    void processAbsoluteBPAndVNTRRegions(std::vector<SV_ITEM>& sv_l,
										CONTIG_INFO& target_contig,
										CONTIG_INFO& query_contig,
										CONTIG_INFO& ref_contig,
										std::set<std::string>& non_vntr_sv_type_sum, bool is_ref_clu) 
	{
		// for each non VNTR region in B
		std::vector<absolute_BP_in_ref_ITEM> non_vntr_SV_in_ref_l; // store more detail info, include: base-contig ID, original non-VNTR region,
		// get_absolute_position_in_ref
		// convert all base-BP in non-VNTR region[contig position] into [ref position] and store into "absolute_BP_in_ref_l"
        get_non_vntr_SV_in_ref_l(non_vntr_SV_in_ref_l, sv_l, ref_contig);
        if (non_vntr_SV_in_ref_l.size() < 2)
			return;

		// Create a counter to track how many SVs fall into each VNTR/non-VNTR region in the reference
		std::map<VNTR_Range *, int> ref_NONE_VNTR_REGION_sv_counter;
        get_ref_NONE_VNTR_REGION_sv_counter(non_vntr_SV_in_ref_l, ref_NONE_VNTR_REGION_sv_counter);

        // Iterate through all counted VNTR regions
		std::map<VNTR_Range *, int>::iterator it = ref_NONE_VNTR_REGION_sv_counter.begin();
		for (; it != ref_NONE_VNTR_REGION_sv_counter.end(); it++)
		{
			// Only process regions with at least 2 SVs (minimum threshold for significance)
			bool sv_in_NON_VNTR_or_VNTR_region_is_enough = (it->second >= 2);
			if (false == sv_in_NON_VNTR_or_VNTR_region_is_enough)
				continue;
				
			// Process all SVs that belong to this significant region
			for (uint i = 0; i < non_vntr_SV_in_ref_l.size(); i++)
			{
				if (non_vntr_SV_in_ref_l[i].ref_vntr_r_p == it->first)
				{
					// Create a unique key for this region using start, end coordinates and type
					std::string REF_VNTR_R_K;
                    get_REF_VNTR_R_K(REF_VNTR_R_K, non_vntr_SV_in_ref_l[i].ref_vntr_r_p);
					//store data by KEY
					non_vntr_sv_type_sum.emplace(REF_VNTR_R_K);
                    // Output debug logs for this SV
					if(false) fprintf(stderr, "%s", non_vntr_SV_in_ref_l[i].get_log().c_str());
				}
			}
		}
    }

   void get_REF_VNTR_R_K(std::string &REF_VNTR_R_K, VNTR_Range * ref_vntr_r_p)
    {
        REF_VNTR_R_K =
            std::to_string(ref_vntr_r_p->start) + "_" +
            std::to_string(ref_vntr_r_p->end) + "_";
        REF_VNTR_R_K += ((ref_vntr_r_p->is_vntr) ? ("VNTR") : ("NON-VNTR"));
		REF_VNTR_R_K += (ref_vntr_r_p->id == -1) ? "FULL" : std::to_string(ref_vntr_r_p->id);
    }

    void get_ref_NONE_VNTR_REGION_sv_counter(std::vector<FC_joint_calling_handler::absolute_BP_in_ref_ITEM> &non_vntr_SV_in_ref_l, std::map<FC_joint_calling_handler::VNTR_Range *, int> &ref_NONE_VNTR_REGION_sv_counter)
    {
        {
            // Iterate through all non-VNTR SVs in the reference list
            for (uint i = 0; i < non_vntr_SV_in_ref_l.size(); i++)
            {
                // Search for the VNTR position for each reference position
                if (non_vntr_SV_in_ref_l[i].ref_vntr_r_p != NULL)
                {
                    // Count SVs per VNTR region - increment existing count or initialize to 1
                    if (ref_NONE_VNTR_REGION_sv_counter.find(non_vntr_SV_in_ref_l[i].ref_vntr_r_p) != ref_NONE_VNTR_REGION_sv_counter.end())
                        ref_NONE_VNTR_REGION_sv_counter[non_vntr_SV_in_ref_l[i].ref_vntr_r_p]++;
                    else
                        ref_NONE_VNTR_REGION_sv_counter[non_vntr_SV_in_ref_l[i].ref_vntr_r_p] = 1;
                }
            }
        }
    }

    void get_non_vntr_SV_in_ref_l(std::vector<FC_joint_calling_handler::absolute_BP_in_ref_ITEM> &non_vntr_SV_in_ref_l, std::vector<FC_joint_calling_handler::SV_ITEM> &sv_l, FC_joint_calling_handler::CONTIG_INFO &ref_contig)
    {
        {
            non_vntr_SV_in_ref_l.clear();
            for (uint i = 0; i < sv_l.size(); i++)
            {
                // sv is in non-VNTR regions
                if (sv_l[i].base_is_VNTR_in_ref == false)
                {
                    non_vntr_SV_in_ref_l.emplace_back();
                    VNTR_Range *ref_v = NULL;
                    if (sv_l[i].base_VNTR_idx_in_ref != -1)
                        ref_v = &(ref_contig.non_vntr_region_l[sv_l[i].base_VNTR_idx_in_ref]);
                    non_vntr_SV_in_ref_l.back().set(sv_l[i].base_absolute_BP_in_ref, sv_l[i].base_BP_core, ref_v, &(sv_l[i]));
                }
            }
        }
    }

    struct ALL_var_store{
		int contig_id;
		int target_contig_id;
		int target_clu_id;
		std::vector<SV_ITEM> sv_l;
	};

	void storeVariantsToFinalMap(std::map<std::string, ALL_var_store> &var_final,
							int contig_query_id,
							int clu_contig_target_id,
							int target_clu_id,
							std::vector<SV_ITEM> &sv_l)
	{
		std::string s = "CON_" + std::to_string(contig_query_id) + "_" + std::to_string(target_clu_id);
		ALL_var_store vtmp;
		vtmp.contig_id = contig_query_id;
		vtmp.target_contig_id = clu_contig_target_id;
		vtmp.target_clu_id = target_clu_id;
		var_final[s] = vtmp;
		std::swap(var_final[s].sv_l, sv_l);
	}

	/**
	 * Calculates the minimum distance between two structural variants (SV_ITEM)
	 * The calculation method depends on whether the SVs have the same insertion type
	 * 
	 * @param sv1 First structural variant
	 * @param sv2 Second structural variant
	 * @return Minimum distance between the two SVs
	 */
	int calculate_SV_Distance(const SV_ITEM &sv1, const SV_ITEM &sv2)
	{
		// If SVs have different insertion types
		if (sv1.is_INS != sv2.is_INS)
		{
			// Calculate distances between base position of one SV and insert positions of the other
			int dist1 = std::abs(sv2.insert_bg - sv1.base_bg); // sv2.insert_bg to sv1.base_bg
			int dist2 = std::abs(sv2.insert_ed - sv1.base_bg); // sv2.insert_ed to sv1.base_bg
			int dist3 = std::abs(sv1.insert_bg - sv2.base_bg); // sv1.insert_bg to sv2.base_bg
			int dist4 = std::abs(sv1.insert_ed - sv2.base_bg); // sv1.insert_ed to sv2.base_bg
			
			// Return the minimum of all calculated distances
			return std::min({dist1, dist2, dist3, dist4});
		}
		else // If SVs have the same insertion type
		{
			// Calculate distances between base positions and all insert position combinations
			int dist0 = std::abs(sv2.base_bg - sv1.base_bg);     // Base to base distance
			int dist1 = std::abs(sv2.insert_bg - sv1.insert_bg); // Insert begin to insert begin
			int dist2 = std::abs(sv2.insert_bg - sv1.insert_ed); // Insert begin to insert end
			int dist3 = std::abs(sv2.insert_ed - sv1.insert_bg); // Insert end to insert begin
			int dist4 = std::abs(sv2.insert_ed - sv1.insert_ed); // Insert end to insert end
			
			// Return the minimum of all calculated distances
			return std::min({dist0, dist1, dist2, dist3, dist4});
		}
	}

	/**
	 * Calculate the normalized length difference ratio between two SVs
	 * This computes the relative length difference normalized by the smaller SV length
	 * @param sv1 First SV item
	 * @param sv2 Second SV item  
	 * @return Normalized length difference ratio (absolute difference divided by min length)
	 */
	float calculate_SV_Length_Difference_Ratio(const SV_ITEM &sv1, const SV_ITEM &sv2)
	{
		// Calculate lengths based on insert string size
		int len1 = sv1.insert_STR.size();
		int len2 = sv2.insert_STR.size();
		
		// Return normalized difference ratio
		return static_cast<float>(std::abs(len1 - len2)) / std::min(len1, len2);
	}

	bool compare_by_annotation(int sequencing_platform, Contig_String_aligner &ca, CONTIG_INFO &target_contig, CONTIG_INFO &query_contig,
							   CONTIG_INFO &ref_contig, CONTIG_COMPARE_RST_ITEM &c_r,
							   std::set<std::string> &non_vntr_sv_type_sum,
							   std::string &CMP_header, int &var_ID, bool is_ref_clu, std::string &var_log, std::string &var_ANNO_log, std::string &main_contig_log, 
							std::map<std::string, ALL_var_store> & var_final, int contig_query_id,int  clu_contig_target_id, int  target_clu_id, int query_clu_id_tmp,
							std::vector<std::vector<SV_ITEM>> &sv_to_ref_l
						)
	{
		//S1: mapping and SV calling of the two contigs
		bool with_sv_is_NOT_vntr = false;
		bool right_aln = ksw2_aln(ca, target_contig, query_contig, c_r);
		if (false == right_aln)
			return -1;
		int MIN_GAP_LEN = 50;

		std::string cigar_str = ca.printCIGAR_core();
		
		std::vector<SV_ITEM> sv_l;
		//not adjust cigar:
		int suggest_st_pos = ca.adjustCIGAR();

		store_SVs_from_cigar(ca, sv_l, target_contig, query_contig, MIN_GAP_LEN, suggest_st_pos);
		//S2: set type for SVs
		// store the absolute_BP_in_ref for all vars
		storeAbsoluteBPForAllVariants(sv_l, target_contig, query_contig, ref_contig, sv_to_ref_l[clu_contig_target_id], sv_to_ref_l[contig_query_id], cigar_str);
		anno_type_analysis_for_all_SV(sv_l, c_r, ref_contig, with_sv_is_NOT_vntr, sequencing_platform);

        replacement_SV_analysis(sv_l);

        // logging all variants and contigs
		main_contig_log += generateMainContigLog(CMP_header, suggest_st_pos, ca, target_contig, query_contig);
		//var logs
		logAllVariants(sv_l, CMP_header, var_ID, var_log, var_ANNO_log);

		// search the region in contig whether the MEI is inside?;
		// get the SVs in current non-vntr region
		if (false)
			 searchMEIRegionsForNovelInsertions(sv_l);

		// but in following code:
		processAbsoluteBPAndVNTRRegions(sv_l, target_contig, query_contig, ref_contig, non_vntr_sv_type_sum, is_ref_clu);
		//store all vars for later using
		storeVariantsToFinalMap(var_final, contig_query_id, clu_contig_target_id, target_clu_id, sv_l);
		return with_sv_is_NOT_vntr;
    }

	void replacement_SV_analysis(std::vector<FC_joint_calling_handler::SV_ITEM> &sv_l)
	{
		// Analyze each SV to find closest replacement candidate
		for (size_t i = 0; i < sv_l.size(); ++i)
		{
			bool current_is_INS = sv_l[i].is_INS;
			int min_distance = INT_MAX;
			int closest_index = -1; // Initialize to -1 indicating no candidate found

			// Check previous SV (if exists)
			if (i > 0)
			{
				bool prev_is_INS = sv_l[i - 1].is_INS;
				if (current_is_INS != prev_is_INS)
				{
					int dis = calculate_SV_Distance(sv_l[i - 1], sv_l[i]);
					if (dis < min_distance)
					{
						min_distance = dis;
						closest_index = i - 1;
					}
				}
			}

			// Check next SV (if exists)
			if (i < sv_l.size() - 1)
			{
				bool next_is_INS = sv_l[i + 1].is_INS;
				if (current_is_INS != next_is_INS)
				{
					int dis = calculate_SV_Distance(sv_l[i + 1], sv_l[i]);
					if (dis < min_distance)
					{
						min_distance = dis;
						closest_index = i + 1;
					}
				}
			}

			// Initialize replacement metrics to -1 (no replacement candidate)
			sv_l[i].replacement_Distance = -1;
			sv_l[i].replacement_Len_Diff_Ratio = -1;

			// If replacement candidate found, calculate metrics
			if (closest_index != -1)
			{
				sv_l[i].replacement_Distance = min_distance;
				sv_l[i].replacement_Len_Diff_Ratio = calculate_SV_Length_Difference_Ratio(sv_l[i], sv_l[closest_index]);
			}
		}
	}

	void show_all_contig(std::vector<CONTIG_INFO> &contig_l, std::string &s)
	{
		for (uint contig_i = 0; contig_i < contig_l.size(); contig_i++)
		{
			fprintf(stdout, "contig[%d:%s][%d:%d:%d] string %s \n",
					contig_i + 1, contig_l[contig_i].sample_name.c_str(), contig_l[contig_i].FULL_COVER_R_N,
					contig_l[contig_i].CLIP_AT_LEFT, contig_l[contig_i].CLIP_AT_RIGHT, contig_l[contig_i].contig.c_str());
		}

		std::ostringstream oss;
		for (uint contig_i = 0; contig_i < contig_l.size(); contig_i++)
		{
			oss << "contig[" << contig_i + 1 << ":" << contig_l[contig_i].sample_name
				<< "][" << contig_l[contig_i].FULL_COVER_R_N << ":"
				<< contig_l[contig_i].CLIP_AT_LEFT << ":"
				<< contig_l[contig_i].CLIP_AT_RIGHT << "] string "
				<< contig_l[contig_i].contig << " \n";
		}
		s = oss.str();
	}

	void store_bin_contig_all(std::vector<std::vector<uint8_t>> &contig_bin_l, std::vector<CONTIG_INFO> &contig_l)
	{
		// store the bin contigs
		contig_bin_l.resize(contig_l.size());
	}

	template <typename T, typename Compare>
	class PointerSorter
	{
	private:
		std::vector<T> data_pool;	
		std::vector<T *> ptr_array; 
		Compare comp;				
	public:
		void addItem(const T &item)
		{
			data_pool.push_back(item);
			ptr_array.clear();
			// re-build index
			for (uint i = 0; i < data_pool.size(); i++)
				ptr_array.push_back(&data_pool[i]);
		}
		int get_size() { return ptr_array.size(); }
		const std::vector<T *> &getSortedPointers() const { return ptr_array; }
		const std::vector<T>   &getData() const { return data_pool; }
		void sortPointers()
		{
			std::sort(ptr_array.begin(), ptr_array.end(), [this](const T *a, const T *b)
					  { return comp(*a, *b); });
		}
	};

	void store_CSV_main_log(std::string &work_dir, std::string &file_name, PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter, std::vector<CONTIG_INFO> &contig_l)
	{
		// Construct full file path
		fs::path dir_path(work_dir);
		fs::path file_path = dir_path / file_name;
		// Create directory if it doesn't exist
		if (!fs::exists(dir_path))
		{
			fs::create_directories(dir_path);
		}
		// Open log file in append mode
		std::ofstream log_file(file_path, std::ios::out | std::ios::trunc);
		if (!log_file.is_open())
		{
			std::cerr << "Error: Could not open log file " << file_path << std::endl;
			return;
		}
		// Write cluster information to file
		int cluster_size = cluster_sorter.get_size();
		log_file << "Cluster list information:\n";
		log_file << "Total number of clusters: " << cluster_size << "\n";
		// Log size of each cluster
		log_file << "\nclusters AC:\n";
		for (int i = 0; i < cluster_size; ++i)
			log_file << "Cluster " << i << " has " << cluster_sorter.getSortedPointers()[i]->contig_id_list.size() << " elements\n";
		log_file << "\nclusters detail:\n";
		for (int i = 0; i < cluster_size; ++i)
		{
			log_file << "Cluster " << i << " has " << cluster_sorter.getSortedPointers()[i]->contig_id_list.size() << " elements\n";
			CLUSTER *c = cluster_sorter.getSortedPointers()[i];
			log_file << "Contig_N:" << std::setw(4) << c->contig_id_list.size() << "\t";
			for (int contig_id : c->contig_id_list)
			{
				if (contig_id > 0)
				{
					log_file
						<< "contig_ID " << std::setw(4) << contig_id << " "
						<< "contig_len " << std::setw(6) << contig_l[contig_id].contig_bin.size() << " "
						<< "SR " << std::setw(3) << contig_l[contig_id].FULL_COVER_R_N << " "
						<< "CL " << std::setw(3) << contig_l[contig_id].CLIP_AT_LEFT << " "
						<< "CR " << std::setw(3) << contig_l[contig_id].CLIP_AT_RIGHT << "\n";
				}
				else
				{
					log_file
						<< "contig_ID " << std::setw(4) << contig_id << " "
						<< "contig_len " << std::setw(6) << contig_l[contig_id].contig_bin.size() << " "
						<< "REF REF REF\n";
				}
			}
			log_file << "\n";
		}

		// Close the file
		log_file.close();
		std::cout << "Contig information logged successfully" << work_dir << file_name << std::endl;
	}

	void store_CLUSTER_main_log(std::string &work_dir, std::string &file_name, PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter, std::vector<CONTIG_INFO> &contig_l)
	{
		// Create directory if it doesn't exist
			fs::create_directories(work_dir);
		// Open log file
		std::ofstream log_file(fs::path(work_dir) / file_name, std::ios::out | std::ios::trunc);
		if (!log_file.is_open())
		{
			std::cerr << "Error: Failed to open log file" << std::endl;
			return;
		}
		// Write header
		log_file << "Contig_Name\tCluster_ID\tAlignment_Info\n";
		log_file << "----------------------------------------\n";
		// Write each contig's information
		for (int i = 0; i < cluster_sorter.get_size(); ++i)
		{
			log_file << "Cluster " << i << " has " << cluster_sorter.getSortedPointers()[i]->contig_id_list.size() << " elements\n";
			CLUSTER *c = cluster_sorter.getSortedPointers()[i];
			for (int i : c->contig_id_list)
			{
				std::string cigar;
				log_file
					<< "CLU_MAIN_" << i << "\t"
					<< contig_l[i].sample_name << "\t"
					<< contig_l[i].mapping_2_ref.printCIGAR_core() << "\t"
					<< contig_l[i].contig << "\t"
					<< "\n";
			}
			log_file << "\n";
		}

		// Close the file
		log_file.close();
		std::cout << "Contig information logged successfully" << work_dir << file_name << std::endl;
	}

	void store_main_log_simple(std::string &work_dir, std::string &file_name, std::string &data)
	{
		// Create directory if it doesn't exist
		if (!fs::exists(work_dir))
			fs::create_directories(work_dir);
		// Open log file
		std::ofstream log_file(fs::path(work_dir) / file_name, std::ios::out | std::ios::trunc);
		if (!log_file.is_open())
		{
			std::cerr << "Error: Failed to open log file" << std::endl;
			return;
		}
		// Write data
		log_file << data;
		// Close the file
		log_file.close();
		std::cout << "Contig information logged successfully" << work_dir << file_name << std::endl;
	}

	void generateClusterAnnotationLogs(PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter,
									  std::vector<CONTIG_INFO> &contig_l,
									  std::string &cluster_ANNO_log_RM,
									  std::string &cluster_ANNO_log_TRF)
	{
		int cluster_size = cluster_sorter.get_size();
		for (int i = 0; i < cluster_size; ++i)
		{
			CLUSTER *c = cluster_sorter.getSortedPointers()[i];
			std::string header = "CLU_ANNO_" + std::to_string(c->cluster_ID) + "_";
			for (int contig_id : c->contig_id_list)
			{
				CONTIG_INFO &c_q = contig_l[contig_id];
				std::string header1 = header + std::to_string(contig_id) + "_RM";
				cluster_ANNO_log_RM += c_q.show_RM_annotation(header1);
			}
			for (int contig_id : c->contig_id_list)
			{
				CONTIG_INFO &c_q = contig_l[contig_id];
				std::string header1 = header + std::to_string(contig_id) + "_TRF";
				cluster_ANNO_log_TRF += c_q.show_TRF_annotation(header1);
			}
		}
	}

	void analyzeSubRegionClusters(std::set<std::string> &non_vntr_sv_type_sum,
								 PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter,
								 std::vector<CONTIG_INFO> &contig_l,
								 CONTIG_INFO &ref_contig,
								 std::map<std::string, ALL_var_store> &var_final,
								 std::string &sub_region_detail)
	{
		//main region:
		{
			VNTR_Range full{-1, 0, ref_contig.contig.size(), 0};
			non_vntr_sub_region_analysis_core(cluster_sorter, contig_l, var_final, 
				full, 
				ref_contig, sub_region_detail, non_vntr_sv_type_sum);
		}
		for (uint i = 0; i < ref_contig.vntr_region_l.size(); i++)
		{
			non_vntr_sub_region_analysis_core(cluster_sorter, contig_l, var_final, 
				ref_contig.vntr_region_l[i], 
				ref_contig, sub_region_detail, non_vntr_sv_type_sum);
		}
		for (uint i = 0; i < ref_contig.non_vntr_region_l.size(); i++)
		{
			non_vntr_sub_region_analysis_core(cluster_sorter, contig_l, var_final, 
				ref_contig.non_vntr_region_l[i], 
				ref_contig, sub_region_detail, non_vntr_sv_type_sum);
		}
		fprintf(stderr, "%s\n", sub_region_detail.c_str());
    }

    void non_vntr_sub_region_analysis_core(PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter,
		 std::vector<CONTIG_INFO> &contig_l, 
		 std::map<std::string, ALL_var_store> &var_final,
		VNTR_Range &vr, 
	 	CONTIG_INFO &ref_contig, 
		std::string &sub_region_detail,
		std::set<std::string> &non_vntr_sv_type_sum)
    {
		std::string sv_type_count_s;
		//get the SV type count
        //sv_type_count_core(var_final, vr, sv_type_count_s);

        std::string REF_VNTR_R_K;
		get_REF_VNTR_R_K(REF_VNTR_R_K, &vr);
		bool is_CSV_sub_region = (non_vntr_sv_type_sum.find(REF_VNTR_R_K) != non_vntr_sv_type_sum.end());
        int totol_cluster_number = 0;
        int home_ref_contig_number = 0;
        std::map<int, std::vector<int>> sub_cluster_all_contig;
		std::map<int, std::vector<int>> sub_cluster_all_ori_clu;
		std::vector<int> type_VINV_GREEN_sv_base_absolute_BP_in_ref;
		int clu_n_with_ONLY_Insert_VNTR_svs = 0;
        //get the cluster number
		cluster_number_count_core(cluster_sorter, contig_l, var_final, vr, sub_cluster_all_contig, sub_cluster_all_ori_clu, type_VINV_GREEN_sv_base_absolute_BP_in_ref, clu_n_with_ONLY_Insert_VNTR_svs, sv_type_count_s);
		std::string non_vntr_sv_type_sum_size = " non_vntr_sv_type_sum_size: " + std::to_string(non_vntr_sv_type_sum.size());
		int vr_len = vr.end - vr.start;
		//store sub region results
        store_sub_region_result(non_vntr_sv_type_sum_size, sv_type_count_s, ref_contig, REF_VNTR_R_K.c_str(), sub_cluster_all_contig, sub_cluster_all_ori_clu, sub_region_detail, is_CSV_sub_region,
				 type_VINV_GREEN_sv_base_absolute_BP_in_ref, clu_n_with_ONLY_Insert_VNTR_svs, vr_len
			);
    }

    void cluster_number_count_core(PointerSorter<CLUSTER, std::greater<CLUSTER>> &cluster_sorter, 
		std::vector<CONTIG_INFO> &contig_l, std::map<std::string, ALL_var_store> &var_final,
		 VNTR_Range &vr, 
		 std::map<int, std::vector<int>> &sub_cluster_all_contig, 
		 std::map<int, std::vector<int>> &sub_cluster_all_ori_clu, 
		 std::vector<int> & type_VINV_GREEN_sv_base_absolute_BP_in_ref, int &clu_n_with_ONLY_Insert_VNTR_svs, std::string &sv_type_count_s)
    {
		clu_n_with_ONLY_Insert_VNTR_svs = 0; 
		int sv_type_count[5] = {0, 0, 0, 0, 0};

		//store yellow
		for (auto it = var_final.begin(); it != var_final.end(); ++it)
        {
            const std::string &key = it->first;
            ALL_var_store &value = it->second;
            for (size_t sv_id = 0; sv_id < value.sv_l.size(); ++sv_id)
            {
                SV_ITEM &sv_item = value.sv_l[sv_id];
                if (vr.with_in_region(sv_item.base_absolute_BP_in_ref) && false == sv_item.isNonVntrInsertType() )
                {
                    sv_type_count[sv_item.ANATYPE]++;
                }
            }
        }

        // for all culsters, if the main contig of a main-cluster belong to a cluster in the sub-region, all other contig belong to it, too.
        for (int cluster_ID = 0; cluster_ID < cluster_sorter.get_size(); cluster_ID++)
        {
            bool is_ref_clu = (cluster_ID == 0);
            int main_contig_id = cluster_sorter.getData()[cluster_ID].main_contig_id;
            xassert(cluster_ID == cluster_sorter.getData()[cluster_ID].cluster_ID, "");

            if (is_ref_clu == true)
                xassert(main_contig_id == 0, "REF is used for cluster 0");
            CONTIG_INFO &c_t = contig_l[main_contig_id];
			if(c_t.right_aln_2_ref == false)
				continue;
            bool is_a_new_sub_cluster = true;
            int combine_cluster_id = -1;


			//check and analysis for type-3 SVs
			//only used for NON-vntr regions
			std::vector<int> TMP_type_VINV_GREEN_sv_base_absolute_BP_in_ref;	
			bool with_blue_insertion = false;
			if(vr.is_vntr == false){
				// check the sub-cluster it belong:
				for (int target_c_idx = 0; target_c_idx < cluster_ID; target_c_idx++)
				{
					std::string s = "CON_" + std::to_string(main_contig_id) + "_" + std::to_string(target_c_idx);
					std::vector<SV_ITEM> &sv_l = var_final[s].sv_l;
					bool with_non_vntr_var = false;
					for (uint i = 0; i < sv_l.size(); i++)
					{
						//skip when the SV is a pure deletion to ref
						if (vr.with_in_region(sv_l[i].base_absolute_BP_in_ref) && sv_l[i].isNonVntrInsertType()){
							with_non_vntr_var = true;
							sv_type_count[sv_l[i].ANATYPE]++;
							if(sv_l[i].ANATYPE == SV_VNTR_INSERT_IN_NON_VNTR){
								TMP_type_VINV_GREEN_sv_base_absolute_BP_in_ref.emplace_back(sv_l[i].base_absolute_BP_in_ref);
							}
							if(sv_l[i].ANATYPE == SV_NOVEL_INSERT_IN_NON_VNTR){
								with_blue_insertion = true;
							}
						}
					}
					// when is same sub-cluster
					if (with_non_vntr_var == false)
					{
						is_a_new_sub_cluster = false;
						combine_cluster_id = target_c_idx;
						break;
					}
				}
			}
            
			//store the SVs to log
			if (is_a_new_sub_cluster)
            {
                // store all contig to the new cluster
                sub_cluster_all_contig[cluster_ID] = cluster_sorter.getData()[cluster_ID].contig_id_list;
				sub_cluster_all_ori_clu[cluster_ID].push_back(cluster_ID);
				if( vr.is_vntr == false){
					//check and analysis for type-3 SVs
					if(false == with_blue_insertion)
						clu_n_with_ONLY_Insert_VNTR_svs ++;
					type_VINV_GREEN_sv_base_absolute_BP_in_ref.insert(
						type_VINV_GREEN_sv_base_absolute_BP_in_ref.end(),
						TMP_type_VINV_GREEN_sv_base_absolute_BP_in_ref.begin(),
						TMP_type_VINV_GREEN_sv_base_absolute_BP_in_ref.end());
				}

            }
            else
            {
                // store all contig to the old cluster
                xassert(combine_cluster_id != -1, "");
                const std::vector<int> &list_to_combine = cluster_sorter.getData()[cluster_ID].contig_id_list;
                sub_cluster_all_contig[combine_cluster_id].insert(sub_cluster_all_contig[combine_cluster_id].end(), list_to_combine.begin(), list_to_combine.end());
				sub_cluster_all_ori_clu[combine_cluster_id].push_back(cluster_ID);
            }

        }

		sv_type_count_s = "sv_type_count: [" +
					std::to_string(sv_type_count[0]) + ", " +
					std::to_string(sv_type_count[1]) + ", " +
					std::to_string(sv_type_count[2]) + ", " +
					std::to_string(sv_type_count[3]) + ", " +
					std::to_string(sv_type_count[4]) + "]";
    }

    void sv_type_count_core(std::map<std::string, FC_joint_calling_handler::ALL_var_store> &var_final, FC_joint_calling_handler::VNTR_Range &vr, std::string &sv_type_count_s)
    {
        // cout all vars
        int sv_type_count[5] = {0, 0, 0, 0, 0};
        for (auto it = var_final.begin(); it != var_final.end(); ++it)
        {
            const std::string &key = it->first;
            ALL_var_store &value = it->second;
            for (size_t sv_id = 0; sv_id < value.sv_l.size(); ++sv_id)
            {
                SV_ITEM &sv_item = value.sv_l[sv_id];
                if (vr.with_in_region(sv_item.base_absolute_BP_in_ref))
                {
                    sv_type_count[sv_item.ANATYPE]++;
                }
            }
        }
        sv_type_count_s = "sv_type_count: [" +
                          std::to_string(sv_type_count[0]) + ", " +
                          std::to_string(sv_type_count[1]) + ", " +
                          std::to_string(sv_type_count[2]) + ", " +
                          std::to_string(sv_type_count[3]) + ", " +
                          std::to_string(sv_type_count[4]) + "]";
    }

    void store_sub_region_result(std::string &non_vntr_sv_type_sum_size,
		std::string &sv_type_count_s, 
		FC_joint_calling_handler::CONTIG_INFO &ref_contig,
		const char *input, 
		std::map<int, std::vector<int>> &sub_cluster_all_contig, 
		std::map<int, std::vector<int>> &sub_cluster_all_ori_clu,
		std::string &sub_region_detail, bool is_CSV_sub_region, std::vector<int> & type_VINV_GREEN_sv_base_absolute_BP_in_ref, int &clu_n_with_ONLY_Insert_VNTR_svs,
		int vr_len)
    {
		int used_hap_number = 0;
		for (const auto& pair : sub_cluster_all_contig)
			used_hap_number += pair.second.size();
        std::ostringstream debug_oss;
        debug_oss << "AT_subregion " << ref_contig.sample_name << "_" << input
				  << " : cluster_size " << sub_cluster_all_contig.size()
                  << ", hom-ref_size " << sub_cluster_all_contig[0].size() 
				  << " is_CSV_sub_region: " << is_CSV_sub_region << " "
				  << sv_type_count_s << non_vntr_sv_type_sum_size << " "
				  << "used_hap_number: " << used_hap_number << " "
				  << "vr_len: " << vr_len << " "
				  << "cluster_size_with_ONLY_Insert_VNTR_svs: " << clu_n_with_ONLY_Insert_VNTR_svs << " ";
		{
			debug_oss << "type_VINV_GREEN_sv_base_absolute_BP_in_ref: [";
			for (size_t i = 0; i < type_VINV_GREEN_sv_base_absolute_BP_in_ref.size(); ++i) {
				debug_oss << type_VINV_GREEN_sv_base_absolute_BP_in_ref[i];
				if (i < type_VINV_GREEN_sv_base_absolute_BP_in_ref.size() - 1) {
					debug_oss << ", ";
				}
			}
			debug_oss << "]";
		}
		debug_oss <<"\n";

        debug_oss << "Debug: sub_cluster_detail contents:" << std::endl;
        for (const auto &[key, vec] : sub_cluster_all_contig)
        {
            debug_oss << "  " << key << ": [";
            for (size_t i = 0; i < vec.size(); ++i)
            {
                debug_oss << vec[i];
                if (i < vec.size() - 1)
                    debug_oss << ", ";
            }
            debug_oss << "]" << std::endl;
        }

		debug_oss << "Debug: sub_cluster_all_ori_clu contents:" << std::endl;
		for (const auto &[key, vec] : sub_cluster_all_ori_clu) {
			debug_oss << "  " << key << ": [";
			for (size_t i = 0; i < vec.size(); ++i) {
				debug_oss << vec[i];
				if (i < vec.size() - 1)
					debug_oss << ", ";
			}
			debug_oss << "]" << std::endl;
		}
        sub_region_detail += debug_oss.str();
    }

	// 嵌套查找函数
	void findNestedStructures(std::map<std::string, ALL_var_store> var_final, std::string& result_record, CONTIG_INFO &ref_contig) {
		        // cout all vars
        
        for (auto it = var_final.begin(); it != var_final.end(); ++it)
        {
            const std::string &key = it->first;
            ALL_var_store &value = it->second;
            for (size_t sv_id = 0; sv_id < value.sv_l.size(); ++sv_id)
            {
                SV_ITEM &sv_item = value.sv_l[sv_id];
				//skip not nest vars
				if(sv_item.base_is_also_INS_in_ref == false){
					continue;
				}
				result_record += "Region: " + ref_contig.sample_name + 
								", key: " + key +
								", sv_id: " + std::to_string(sv_id) +
								", base_Sample: " + sv_item.base_contig->sample_name + 
								", insert_Sample: " + sv_item.insert_contig->sample_name + 
								", VNTRIdxInRef: " + std::to_string(sv_item.base_VNTR_idx_in_ref) + 
								", VNTRIdxInContig: " + std::to_string(sv_item.base_VNTR_idx_in_contig) + 
								", IsInsertionSeqVNTR: " + (sv_item.insertion_seq_with_is_vntr ? "true" : "false") + 
								", IsVNTRInRef: " + (sv_item.base_is_VNTR_in_ref ? "true" : "false") + 
								", IsVNTRInContig: " + (sv_item.base_is_VNTR_in_contig ? "true" : "false") +
								", VNTR_LEN: " + std::to_string(sv_item.insert_VNTR_LEN_SUM) + 
								", Non_VNTR_LEN: " + std::to_string(sv_item.insert_non_VNTR_LEN) + 
								", Total_STR_LEN: " + std::to_string(sv_item.insert_STR.size()) +
								", Non_VNTR_rate: " + std::to_string(static_cast<float>(sv_item.insert_non_VNTR_LEN)/(sv_item.insert_STR.size())) + 
								", Base_absolute_BP_in_ref: " + std::to_string(sv_item.base_absolute_BP_in_ref) +
								", Base_is_INS_in_ref: " + (sv_item.base_is_also_INS_in_ref ? "Yes" : "No") +
								", base_absolute_BP_in_ref_event_len: " + std::to_string(sv_item.base_absolute_BP_in_ref_event_len) +
								", Contig_bg_BP_in_ref: " + std::to_string(sv_item.contig_bg_absolute_BP_in_ref) +
								", Contig_ed_BP_in_ref: " + std::to_string(sv_item.contig_ed_absolute_BP_in_ref) +
								", insert_overlap_MEI_len: " + std::to_string(sv_item.insert_overlap_MEI_len) +
								", replacement_Distance: " + std::to_string(sv_item.replacement_Distance) +
								", replacement_Len_Diff_Ratio: " + std::to_string(sv_item.replacement_Len_Diff_Ratio) +
								", insert_bg: " + std::to_string(sv_item.insert_bg) +
								", insert_ed: " + std::to_string(sv_item.insert_ed) +
								", base_bg: " + std::to_string(sv_item.base_bg) +
								", contig_3_1_INS_absolute_BP_in_ref: " + std::to_string(sv_item.contig_3_1_INS_absolute_BP_in_ref) +
								", nest_string: " + sv_item.nest_string +
								
								+ "\n";
				;
            }
        }
	}

    void contig_clustering(Contig_String_aligner &ca, std::vector<CONTIG_INFO> &contig_l, int sequencing_platform, std::string &work_dir, int max_used_contig, 
		std::vector<std::vector<SV_ITEM>> &sv_to_ref_l)
	{
		// show basic data
		int contig_n = contig_l.size();
		// The clusters
		// std::vector<CLUSTER> cls;
		PointerSorter<CLUSTER, std::greater<CLUSTER>> cluster_sorter;
		// classification:
		int total_used_hap = 0;
		int var_ID = 0;
		//reference contig:
		CONTIG_INFO &ref_contig = contig_l[0];
		//get VNTR region for all contigs
		for (int contig_query_id = 0; contig_query_id < contig_n && contig_query_id < max_used_contig; contig_query_id++)
			 contig_l[contig_query_id].get_VNTR_REGION();

		//store log when new cluster generated
		std::string var_log_full_NC = "var_log\n";
		std::string var_ANNO_log_full_NC = "var_ANNO_log\n";
		std::string contig_log_full_NC = "main_contig_log\n";
		//store log all the times
		std::string var_log_full = "var_log\n";
		std::string var_ANNO_log_full = "var_ANNO_log\n";
		std::string contig_log_full = "main_contig_log\n";

		//store all the variants
		std::map<std::string, ALL_var_store> var_final;
		std::set<std::string> non_vntr_sv_type_sum;

		for (int contig_query_id = 0; contig_query_id < contig_n && contig_query_id < max_used_contig; contig_query_id++)
		{
			CONTIG_INFO &c_q = contig_l[contig_query_id];
			if (contig_query_id > 0 && !c_q.is_High_depth(0))
				continue;
			cluster_sorter.sortPointers();
			total_used_hap++;
			// store the contig to clusters
			bool is_new_cluster = true;
			std::string var_log;
			std::string var_ANNO_log;
			std::string contig_log;
			for (int c_idx = 0; c_idx < cluster_sorter.get_size(); c_idx++)
			{
				bool is_ref_clu = (c_idx == 0);
				int clu_contig_target_id = cluster_sorter.getSortedPointers()[c_idx]->main_contig_id;
				if(is_ref_clu == true)
					xassert(clu_contig_target_id == 0, "REF is used for cluster 0");

				CONTIG_INFO &c_t = contig_l[clu_contig_target_id];
				// KSW2 comparing
				CONTIG_COMPARE_RST_ITEM c_r;
				memset(&c_r, 0, sizeof(CONTIG_COMPARE_RST_ITEM));
				
				int target_clu_id = cluster_sorter.getSortedPointers()[c_idx]->cluster_ID;
				int query_clu_id_tmp = cluster_sorter.get_size();
				std::string cmp_header = "CLU_CMP_" + std::to_string(query_clu_id_tmp) + "_" + std::to_string(target_clu_id);
				bool with_sv_is_NOT_vntr = compare_by_annotation(sequencing_platform, ca, c_t, c_q, ref_contig, c_r, non_vntr_sv_type_sum,
																 cmp_header, var_ID, is_ref_clu,
																 var_log, var_ANNO_log, contig_log, var_final,
																 contig_query_id, clu_contig_target_id, target_clu_id, query_clu_id_tmp, sv_to_ref_l
																);
				bool is_same_hap = (with_sv_is_NOT_vntr == false);
				// main contig log
				if (true == is_same_hap)
				{
					is_new_cluster = false;
					cluster_sorter.getSortedPointers()[c_idx]->store(contig_query_id);
					break;
				}
			}
			//store log all the times
			{
				var_log_full += var_log;
				var_ANNO_log_full += var_ANNO_log;
				contig_log_full += contig_log;
			}
			//store log when new cluster generated
			if (is_new_cluster){
				var_log_full_NC += var_log;
				var_ANNO_log_full_NC += var_ANNO_log;
				contig_log_full_NC += contig_log;
				store_new_cluster(contig_query_id, cluster_sorter);
			}
                
		}
		// show classification results
		cluster_sorter.sortPointers();

		//found CSV in non-vntr region, re-calculate class number for all ref-based-non-vntr region: 
		std::string sub_region_detail = "sub_region_detail\n";
		analyzeSubRegionClusters(non_vntr_sv_type_sum, cluster_sorter, contig_l, ref_contig, var_final, sub_region_detail);
		std::string nestedStructures_detail = "nestedStructures_detail\n";
		findNestedStructures(var_final, nestedStructures_detail, ref_contig);

		//嵌套查找

		//LOGs and output
		std::string cluster_ANNO_log_RM = "cluster_ANNO_log_RM\n";
		std::string cluster_ANNO_log_TRF = "cluster_ANNO_log_TRF\n";
		generateClusterAnnotationLogs(cluster_sorter, contig_l, cluster_ANNO_log_RM, cluster_ANNO_log_TRF);

		// TIRE 1
		// LOG:CSV.main.log
		std::string file_name = "CSV.main.log";
		store_CSV_main_log(work_dir, file_name, cluster_sorter, contig_l);
		//LOG: sub_region_detail
		file_name = "sub_region_detail.log";
		store_main_log_simple(work_dir, file_name, sub_region_detail);

		// TIRE 2
		// LOG:cluster.main.log
		file_name = "cluster.main.log";
		store_CLUSTER_main_log(work_dir, file_name, cluster_sorter, contig_l);
		file_name = "cluster.ANNO.RM.log";
		store_main_log_simple(work_dir, file_name, cluster_ANNO_log_RM);
		file_name = "cluster.ANNO.TRF.log";
		store_main_log_simple(work_dir, file_name, cluster_ANNO_log_TRF);
		//contig log:
		std::string contig_VNTR_AND_ANNO_log;
		{
			std::ostringstream ss;
			for (int i = 0; i < contig_n && i < max_used_contig; i++)
 				contig_l[i].get_contig_anno_line(ss);
			contig_VNTR_AND_ANNO_log = ss.str();
		}
		file_name = "contig_VNTR_AND_ANNO.log";
		store_main_log_simple(work_dir, file_name, contig_VNTR_AND_ANNO_log);

		// TIRE 3
		// LOG:var log
		file_name = "var.main_NC.log";
		store_main_log_simple(work_dir, file_name, var_log_full_NC);
		// LOG:ANNO log
		file_name = "var.anno_TRF_RM_NC.log";
		store_main_log_simple(work_dir, file_name, var_ANNO_log_full_NC);
		// LOG: contig log
		file_name = "contig.main_NC.log";
		store_main_log_simple(work_dir, file_name, contig_log_full_NC);

		// LOG:var log
		file_name = "var.main.log";
		store_main_log_simple(work_dir, file_name, var_log_full);
		// LOG:ANNO log
		file_name = "var.anno_TRF_RM.log";
		store_main_log_simple(work_dir, file_name, var_ANNO_log_full);
		// LOG: contig log
		file_name = "contig.main.log";
		store_main_log_simple(work_dir, file_name, contig_log_full);

		// LOG: nestedStructures_detail log
		file_name = "nestedStructures.log";
		store_main_log_simple(work_dir, file_name, nestedStructures_detail);
		

    }

    void store_new_cluster( int contig_query_id, 
		FC_joint_calling_handler::PointerSorter<FC_joint_calling_handler::CLUSTER, std::greater<FC_joint_calling_handler::CLUSTER>> &cluster_sorter)
    {
        fprintf(stdout, "HAP is new!\n");
        CLUSTER n;
        n.new_c(contig_query_id, cluster_sorter.get_size());
        cluster_sorter.addItem(n);
    }

	int align_main(int argc, char *argv[])
	{
		if (argc < 8) {
			std::cerr << "Usage: " << argv[0] << " <filename> <match> <mismatch> <gap_open> <gap_extend> <gap_open2> <gap_extend2>" << std::endl;
			std::cerr << "       " << " <filename>: The first line is read as the target sequence " << ";" << "The second line is read as the query sequence" << std::endl;
			return 1;
		}

		const char* filename = argv[1];
		int match_D = atoi(argv[2]);
		int mismatch_D = atoi(argv[3]);
		int gap_open_D = atoi(argv[4]);
		int gap_ex_D = atoi(argv[5]);
		int gap_open2_D = atoi(argv[6]);
		int gap_ex2_D = atoi(argv[7]);

		// Read sequences from file
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error: Could not open file " << filename << std::endl;
			return 1;
		}

		std::string target_seq, query_seq;
		std::getline(file, target_seq);
		std::getline(file, query_seq);
		file.close();

		// Convert sequences to binary format
		std::vector<uint8_t> target_bin, query_bin;
		store_bin_contig(target_seq, target_bin);
		store_bin_contig(query_seq, query_bin);

		// Initialize aligner with parameters
		Contig_String_aligner ca;
		ca.init();
		ca.reset_paramater(match_D, mismatch_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D);

		// Perform alignment
		int zdrop = MAX(target_bin.size(), query_bin.size());
		ca.setZdrop(zdrop, zdrop);
		ca.align_non_splice_normal_size(query_bin, target_bin);

		// Print alignment results
		std::cout << ca.printCIGAR_core() << std::endl;
		return 0;
	}

    // AS input, the list of contigs, stored in FASTA format
	int ana_main(int argc, char *argv[])
	{
		std::string work_dir = argv[1];
		const char *FC_Contig_fasta_L = argv[2];
		const char *RM_L = argv[3];
		const char *TRF_L = argv[4];

		int sequencing_platform = 0;
		if (strcmp(argv[5], "HPRC") == 0) {
			sequencing_platform = sequencing_platform_HPRC;
		} else if (strcmp(argv[5], "ONT") == 0) {
			sequencing_platform = sequencing_platform_ONT;
		} else {
			xassert(0, "dataset must be one of HPRC or ONT");
		}

		Contig_String_aligner ca;
		ca.init();
		ca.reset_paramater(1, 4, 18, 3, 65, 1);

		std::vector<CONTIG_INFO> contig_l;
		std::vector<RM_record> RM_record_l;
		std::vector<TRF_record> TRF_record_l;
		std::vector<std::vector<SV_ITEM>> sv_to_ref_l;

		int max_used_contig = 3000000;

		load_contig_list_one_region_fasta(FC_Contig_fasta_L, contig_l);
		load_repeat_masker_result(RM_L, RM_record_l);
		load_trf_result(TRF_L, TRF_record_l);
		store_RM_TRF_data_to_contig_l(contig_l, RM_record_l, TRF_record_l);
		// align all contigs to reference, skip the reference
		xassert(contig_l[0].sample_name.find("ref") != std::string::npos, "The first contig must BE reference");
		for (uint i = 0; i < contig_l.size(); i++)
			store_bin_contig(contig_l[i].contig, contig_l[i].contig_bin);
		int mapping_right = 0;
		for (uint i = 0; i < contig_l.size() && i < max_used_contig; i++){
			if(contig_l[i].ksw2_aln_2_ref(contig_l[0].contig_bin)){
				mapping_right++;
				//keep var list to ref for each contig
				sv_to_ref_l.emplace_back();
				store_SVs_from_cigar(contig_l[i].mapping_2_ref, sv_to_ref_l.back(), contig_l[0], contig_l[i], 50, 0);	
			}
		}
		std::cout << "LOG: mapping_right contig number: " << mapping_right << " ,total contig number: " << contig_l.size() << std::endl;

		// analysis:
		contig_clustering(ca, contig_l, sequencing_platform, work_dir, max_used_contig, sv_to_ref_l);
		return 0; // skip others
	}
};

#endif /* SVCALLING_CORE_ANALYSIS_FC_JOINT_CALLING_HPP_ */
