/*
 * T2T_analysis_handler.hpp
 *
 *  Created on: 2025-8-27
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_ANALYSIS_ASM_ANALYSIS_HANDLER_HPP_
#define SVCALLING_CORE_ANALYSIS_ASM_ANALYSIS_HANDLER_HPP_

#include <iostream>
#include <iomanip>

struct ASM_ANALYSIS_HANDLER
{

	struct BAM_RST_ITEM
	{
		std::string LRS_bamFileName;
		std::string selected_read_name_forward;
		std::string selected_read_name_reverse;

		int pos_AT_read_forward = 0;
		int pos_AT_read_reverse = 0;
		int total_read_len_forward = 0;
		int total_read_len_reverse = 0;
		int ref_pos_forward = 0;
		int ref_pos_reverse = 0;

		bool will_be_used = false;

		int get_the_mapping_end_position_at_ref(bam1_t *br, int &region_AT_read_ed, int &total_read_len)
		{
			region_AT_read_ed = -1;
			int seq_i = 0;
			int ref_i = 0;
			total_read_len = 0;
			int cigar_end = br->core.n_cigar;
			uint32_t *bam_cigar = bam_get_cigar(br);
			for (int i = 0; i < cigar_end; i++)
			{
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));

				if (type == align_t::CIGAR_SOFT_CLIP || type == align_t::CIGAR_HARD_CLIP)
				{
					if (i == 0)
						seq_i += length;
					total_read_len += length;
				}
				else if (type == align_t::CIGAR_INSERT)
				{
					seq_i += length;
					total_read_len += length;
				}
				else if (type == align_t::CIGAR_DELETE)
				{
					ref_i += length;
				}
				else if (type == align_t::CIGAR_MATCH || type == align_t::CIGAR_SEQ_MATCH || type == align_t::CIGAR_SEQ_MISMATCH)
				{
					ref_i += length;
					seq_i += length;
					total_read_len += length;
				}
				else
				{
					// DO nothing
				}
			}
			region_AT_read_ed = seq_i;
			return br->core.pos + ref_i;
		}

		int get_the_mapping_begin_position_at_ref(bam1_t *br, int &region_AT_read_st, int &total_read_len)
		{
			region_AT_read_st = -1;
			int seq_i = 0;
			int ref_i = 0;
			total_read_len = 0;
			int cigar_end = br->core.n_cigar;
			uint32_t *bam_cigar = bam_get_cigar(br);
			for (int i = 0; i < cigar_end; i++)
			{
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
				if (type == align_t::CIGAR_SOFT_CLIP || type == align_t::CIGAR_HARD_CLIP)
				{
					if (i == 0)
					{
						seq_i += length;
						region_AT_read_st = seq_i;
					}
					total_read_len += length;
				}
				else if (type == align_t::CIGAR_INSERT)
				{
					seq_i += length;
					total_read_len += length;
				}
				else if (type == align_t::CIGAR_DELETE)
				{
					ref_i += length;
				}
				else if (type == align_t::CIGAR_MATCH || type == align_t::CIGAR_SEQ_MATCH || type == align_t::CIGAR_SEQ_MISMATCH)
				{
					ref_i += length;
					seq_i += length;
					total_read_len += length;
				}
				else
				{
					// DO nothing
				}
			}
			return br->core.pos;
		}

		int get_the_nearst_position_in_BAM(Bam_file *c_b, R_region &region, bool search_forward, std::string &selected_read_name, int &pos_AT_read, int &f_total_read_len)
		{
			selected_read_name.clear();
			resetRegion_ID(c_b, &region); // reset region
			int read_local_c = 0;
			int nearst_position = MAX_int32t;
			if (false == search_forward)
				nearst_position = -1;
			while (bam_next(c_b))
			{
				read_local_c++;
				// new sig for this read
				bam1_t *br = &(c_b->_brec);
				// filter://basic filter
				if (bam_is_secondary(br))
					continue;
				if (br->core.qual < 1)
					continue;
				// supplementary-min length check
				if (bam_is_supplementary(br) && br->core.l_qseq < 1200)
					continue;
				if (br->core.l_qseq < 200000)
					continue;
				// from the core position to a bigger position, to find the begin of first read;
				if (search_forward)
				{
					int region_AT_read_st;
					;
					int total_read_len;
					int the_bgn = get_the_mapping_begin_position_at_ref(br, region_AT_read_st, total_read_len);
					if (nearst_position > the_bgn)
					{
						nearst_position = the_bgn;
						selected_read_name.clear();
						selected_read_name.append(bam_get_qname(br));
						pos_AT_read = region_AT_read_st;
						f_total_read_len = total_read_len;
					}
				}
				else
				{
					int region_AT_read_ed;
					int total_read_len;
					int the_end = get_the_mapping_end_position_at_ref(br, region_AT_read_ed, total_read_len);
					if (nearst_position < the_end)
					{
						nearst_position = the_end;
						selected_read_name.clear();
						selected_read_name.append(bam_get_qname(br));
						pos_AT_read = region_AT_read_ed;
						f_total_read_len = total_read_len;
					}
				}
			}
			return nearst_position;
		}
		void closeBAM()
		{
		}

		void open_BAM(const char *LRS_BAM_F)
		{
			LRS_bamFileName.append(LRS_BAM_F);
		}
		void mark_bam(int chrID, int CSV_core, int &MIN_LEFT, int &MAX_RIGHT)
		{
			int ref_pos_len = ref_pos_forward - ref_pos_reverse;
			int read_pos_len = pos_AT_read_forward - pos_AT_read_reverse;
			int dis_A = ABS_U(ref_pos_reverse, CSV_core);
			int dis_B = ABS_U(ref_pos_forward, CSV_core);

			if (!selected_read_name_reverse.empty() && selected_read_name_reverse.compare(selected_read_name_forward) == 0 && dis_A < 1500000 && dis_B < 1500000)
			{
				xassert(total_read_len_forward == total_read_len_reverse, "");
				will_be_used = true;

				fprintf(stdout,
						"[SUMMARY] [SUCCESS] CORE_POS %d:%d "
						"ref_pos [%d-%d, len %d], "
						"pos_AT_read  [%d-%d, len %d]\n",
						chrID, CSV_core,
						ref_pos_reverse, ref_pos_forward, ref_pos_len,
						pos_AT_read_reverse, pos_AT_read_forward, read_pos_len);
			}
			else
			{
				fprintf(stdout,
						"[SUMMARY] [FAIL] CORE_POS %d:%d \n",
						chrID, CSV_core);
			}

			if (will_be_used)
			{
				MIN_LEFT = MIN(MIN_LEFT, ref_pos_reverse);
				MAX_RIGHT = MAX(MAX_RIGHT, ref_pos_forward);
			}
		}
		void search_BAM(int chrID, int CSV_core, int search_region[11], char *fa_fn_in)
		{
			BAM_handler LRS_read;

			if (!LRS_bamFileName.empty())
				LRS_read.init(LRS_bamFileName.c_str(), NULL, NULL, fa_fn_in);

			Bam_file *c_b = &(LRS_read.file);
			// search forward?
			{
				for (int level = 0; level < 11; level++)
				{
					// check the END;
					// set the region;
					R_region region;
					region.chr_ID = chrID;
					region.st_pos = CSV_core + 5000;
					region.ed_pos = CSV_core + search_region[level];
					ref_pos_forward = get_the_nearst_position_in_BAM(c_b, region, true, selected_read_name_forward, pos_AT_read_forward, total_read_len_forward);
					if (MAX_int32t != ref_pos_forward)
					{
						fprintf(stdout, "Search FORWARD: Right END = %d, read name %s  pos_AT_read %d , total_read_len %d \n", ref_pos_forward, selected_read_name_forward.c_str(), pos_AT_read_forward, total_read_len_forward);
						break;
					}
				}
			}

			// search reverse
			{
				for (int level = 0; level < 11; level++)
				{
					// check the END;
					// set the region;
					R_region region;
					region.chr_ID = chrID;
					region.st_pos = CSV_core - search_region[level];
					region.ed_pos = CSV_core - 5000;
					region.st_pos = MAX(0, region.st_pos);
					ref_pos_reverse = get_the_nearst_position_in_BAM(c_b, region, false, selected_read_name_reverse, pos_AT_read_reverse, total_read_len_reverse);
					if (-1 != ref_pos_reverse)
					{
						fprintf(stdout, "Search REVERSE Right END = %d, read name %s  pos_AT_read %d , total_read_len %d \n", ref_pos_reverse, selected_read_name_reverse.c_str(), pos_AT_read_reverse, total_read_len_reverse);
						break;
					}
				}
			}
			// close bam
			if (!LRS_bamFileName.empty())
				LRS_read.destroy();
		}
	};

	int run(int argc, char *argv[])
	{
		// input
		char *fa_fn_in = argv[1];
		char *LRS_bamFileName_l = argv[2];
		int chrID = (int)atoi(argv[3]);
		int CSV_core = (int)atoi(argv[4]);

		faidx_t *c_ref_idx = reference_index_load(fa_fn_in);
		BAM_handler LRS_read;
		// load reference index
		FILE *try_open = xopen(fa_fn_in, "r");
		fclose(try_open);

		// bam_hdr_t * cur_header = NULL;
		if (LRS_bamFileName_l != NULL)
			LRS_read.init(LRS_bamFileName_l, NULL, NULL, fa_fn_in);
		int search_region[11] = {5000, 10000, 20000, 40000, 80000, 200000, 400000, 800000, 2000000, 4000000, 10000000};
		// init:
		std::vector<std::string> all_bam_list;
		load_string_list_from_file(LRS_bamFileName_l, all_bam_list);
		std::vector<BAM_RST_ITEM> rst_list;
		// load all bam files????
		int MIN_LEFT = MAX_int32t;
		int MAX_RIGHT = -1;
		for (uint sid = 0; sid < all_bam_list.size(); sid++)
		{ //??????
			//			if(sid == 2)		break;
			rst_list.emplace_back();
			BAM_RST_ITEM &cur_rst = rst_list.back();
			cur_rst.open_BAM(all_bam_list[sid].c_str());
			cur_rst.search_BAM(chrID, CSV_core, search_region, fa_fn_in);
			cur_rst.mark_bam(chrID, CSV_core, MIN_LEFT, MAX_RIGHT);
			cur_rst.closeBAM();
		};
		fprintf(stderr, "%s\t%d\t%d\tORI_%d_%d_DIS_%d_%d\n ", faidx_iseq(c_ref_idx, chrID), MIN_LEFT, MAX_RIGHT, chrID, CSV_core, CSV_core - MIN_LEFT, MAX_RIGHT - CSV_core);
		return 0;
	}

	struct T2T_SA_INFO
	{
		int sa_idx;
		int SA_chrID;
		int SA_POS;
		std::vector<path_segment> sa_cigar_l;
		int SA_MAPQ;
		int SA_number_of_mismatch;
		int sa_in_read_position_bg = 0;
		int sa_direction;
		int mapping_len_ref = 0;
		int mapping_len_read = 0;

		void char_cigar_to_path(const char *string_cigar, std::vector<path_segment> &cigar_l)
		{
			cigar_l.clear();
			int length = 0;
			while ((*string_cigar) != 0)
			{
				char c = *string_cigar;
				string_cigar++;
				if (c <= '9' && c >= '0')
					length = (length * 10) + c - '0';
				else
				{
					int type = cigar_code_to_segment_type(c);
					xassert(type != CIGAR_NONE, "Wrong Cigar!\n");
					cigar_l.emplace_back();
					cigar_l.back().type = type;
					cigar_l.back().length = length;
					length = 0;
				}
			}
		}

		void get_the_mapping_length_at_ref()
		{
			mapping_len_ref = 0;
			mapping_len_read = 0;
			for (uint i = 0; i < sa_cigar_l.size(); i++)
			{
				int length = sa_cigar_l[i].length;
				int type = sa_cigar_l[i].type;
				if (type == align_t::CIGAR_SOFT_CLIP || type == align_t::CIGAR_HARD_CLIP)
				{
					// DO NOTHING
				}
				else if (type == align_t::CIGAR_INSERT)
				{
					mapping_len_read += length;
				}
				else if (type == align_t::CIGAR_DELETE)
				{
					mapping_len_ref += length;
				}
				else if (type == align_t::CIGAR_MATCH || type == align_t::CIGAR_SEQ_MATCH || type == align_t::CIGAR_SEQ_MISMATCH)
				{
					mapping_len_read += length;
					mapping_len_ref += length;
				}
				else
				{
					// DO nothing
				}
			}
		}

		void store(int sa_idx, std::vector<std::string> &item_tmp, faidx_t *c_ref_idx, int total_read_len, bool main_is_fwd)
		{
			this->sa_idx = sa_idx;
			// chr5,47154842,+,87624S77370M374I1749543S,1,12794
			SA_chrID = faidx_get_chrID(c_ref_idx, item_tmp[0].c_str(), NULL, 0);
			SA_POS = (int)atoi(item_tmp[1].c_str());
			bool is_fwd = ((item_tmp[2][0] == '+') ? true : false);
			sa_direction = ((main_is_fwd == is_fwd) ? 1 : 0);
			char_cigar_to_path(item_tmp[3].c_str(), sa_cigar_l);
			SA_MAPQ = (int)atoi(item_tmp[4].c_str());
			SA_number_of_mismatch = (int)atoi(item_tmp[5].c_str());
			// analysis cigar:
			sa_in_read_position_bg = 0;
			if (!sa_cigar_l.empty() && (sa_cigar_l[0].type == CIGAR_SOFT_CLIP || sa_cigar_l[0].type == CIGAR_HARD_CLIP))
			{
				sa_in_read_position_bg = sa_cigar_l[0].length;
			}
			get_the_mapping_length_at_ref();
		}

		void store_main(int SA_chrID, int SA_POS, uint32_t *bam_cigar, int n_cigar, int total_read_len)
		{
			// chr5,47154842,+,87624S77370M374I1749543S,1,12794
			this->SA_chrID = SA_chrID;
			this->SA_POS = SA_POS;
			sa_direction = 1;
			for (int i = 0; i < n_cigar; i++)
			{
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
				sa_cigar_l.emplace_back();
				sa_cigar_l.back().length = length;
				sa_cigar_l.back().type = type;
			}
			get_the_mapping_length_at_ref();
			SA_MAPQ = 99;
			// analysis cigar:
			sa_in_read_position_bg = 0;
			if (!sa_cigar_l.empty() && (sa_cigar_l[0].type == CIGAR_SOFT_CLIP || sa_cigar_l[0].type == CIGAR_HARD_CLIP))
			{
				sa_in_read_position_bg = sa_cigar_l[0].length;
			}
		}

		std::string toString(int Read_gap_size, int Ref_gap_size) const
		{
			//std::string cigar_str;
			//for (uint i = 0; i < sa_cigar_l.size(); i++)
			//	cigar_str += std::to_string(sa_cigar_l[i].length) + pathSTR[sa_cigar_l[i].type];

			return
				// main
				"ref_chrID:" + std::to_string(SA_chrID) + ", " +
				"read_st:" + std::to_string(sa_in_read_position_bg) + ", " +
				"read_ed:" + std::to_string(sa_in_read_position_bg + mapping_len_read) + ", " +
				"ref_st:" + std::to_string(SA_POS) + ", " +
				"ref_ed:" + std::to_string(SA_POS + mapping_len_ref) + ", " +
				"sa_direction:" + std::to_string(sa_direction) + ", " +
				"MAPQ:" + std::to_string(SA_MAPQ) + ", " +
				/// additional
				"read_mapping_len:" + std::to_string(mapping_len_read) + ", " +
				"read_mapping_diff:" + std::to_string(SA_POS - sa_in_read_position_bg) + ", " +
				"ref_mapping_len:" + std::to_string(mapping_len_ref) + ", " +
				"sa_idx:" + std::to_string(sa_idx) + ", " +
				"Read_gap_size:" + std::to_string(Read_gap_size) + ", " +
				"Ref_gap_size:" + std::to_string(Ref_gap_size) + ", " +
				// detail
				"number_of_mismatch:" + std::to_string(SA_number_of_mismatch) + ", " +
				//"cigar_str:" + cigar_str  + ", " +
				"END" 
				;
		}

		static inline int cmp_by_r_pos_mapDIR(const T2T_SA_INFO &a, const T2T_SA_INFO &b)
		{
			// var basic
			if (a.sa_direction != b.sa_direction)
				return a.sa_direction < b.sa_direction;
			return a.sa_in_read_position_bg < b.sa_in_read_position_bg;
		}

		static inline int cmp_by_r_pos(const T2T_SA_INFO &a, const T2T_SA_INFO &b)
		{
			// var basic
			return a.sa_in_read_position_bg < b.sa_in_read_position_bg;
		}
	};

	int get_total_read_len(uint32_t *bam_cigar, int n_cigar)
	{
		int total_read_len = 0;
		// get the reference sequence:
		for (int i = 0; i < n_cigar; i++)
		{
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			if (type == align_t::CIGAR_SOFT_CLIP || type == align_t::CIGAR_MATCH || type == align_t::CIGAR_SEQ_MATCH || type == align_t::CIGAR_SEQ_MISMATCH || type == align_t::CIGAR_INSERT || type == align_t::CIGAR_HARD_CLIP)
			{
				total_read_len += length;
			}
			else
			{
				// DO nothing
			}
		}
		return total_read_len;
	}

	struct Mapping
	{
		double a_start, a_end;
		double b_start, b_end;
		int ref_chr_id;
		Mapping(double a_start, double a_end, double b_start, double b_end, int ref_chr_id)
		{
			this->a_start = a_start;
			this->a_end = a_end;
			this->b_start = b_start;
			this->b_end = b_end;
			this->ref_chr_id = ref_chr_id;
		}
	};

	struct Result
	{
		struct ref_region
		{
			int st;
			int ed;
			int chrID;
			ref_region(int st, int ed, int chrID)
			{
				this->st = st;
				this->ed = ed;
				this->chrID = chrID;
			}
		};
		double a_start, a_end;
		std::vector<ref_region> b_mappings; // <start, end>
		double ratio;
		bool has_mapping;
	};

	void calculate_mappings_depth(int total_read_len, std::vector<T2T_SA_INFO> &sa_l)
	{

		double L = total_read_len;
		std::vector<Mapping> mappings;
		// store all mapping infos
		for (uint i = 0; i < sa_l.size(); i++)
		{
			if (sa_l[i].mapping_len_read < 30000)
				continue;
			if (sa_l[i].SA_MAPQ < 15)
				continue;
			if (sa_l[i].SA_chrID > 24)
				continue;
			mappings.emplace_back(
				sa_l[i].sa_in_read_position_bg, sa_l[i].sa_in_read_position_bg + sa_l[i].mapping_len_read,
				sa_l[i].SA_POS, sa_l[i].SA_POS + sa_l[i].mapping_len_ref,
				sa_l[i].SA_chrID);
		}

		std::vector<Result> results;
		std::vector<std::pair<double, int>> points; // <position, type: 0=start, 1=end>

		// 收集所有关键点并排序
		for (const auto &map : mappings)
		{
			points.emplace_back(map.a_start, 0);
			points.emplace_back(map.a_end, 1);
		}

		// 添加起点和终点
		points.emplace_back(0.0, -1);
		points.emplace_back(L, -1);

		// 排序关键点
		std::sort(points.begin(), points.end(),
				  [](const auto &a, const auto &b)
				  {
					  return a.first < b.first;
				  });

		// 处理每个区间
		for (size_t i = 0; i < points.size() - 1; ++i)
		{
			double start = points[i].first;
			double end = points[i + 1].first;

			if (start >= end)
				continue;

			Result res;
			res.a_start = start;
			res.a_end = end;
			res.has_mapping = false;

			// 检查该区间在哪些映射中
			for (const auto &map : mappings)
			{
				if (start >= map.a_start && end <= map.a_end)
				{
					// 计算B线段上的映射区间
					double ratio = (map.b_end - map.b_start) / (map.a_end - map.a_start);
					double b_seg_start = map.b_start + ratio * (start - map.a_start);
					double b_seg_end = map.b_start + ratio * (end - map.a_start);

					res.b_mappings.emplace_back(b_seg_start, b_seg_end, map.ref_chr_id);
					res.ratio = ratio;
					res.has_mapping = true;
				}
			}

			results.push_back(res);
		}
		for (const auto &res : results)
		{
			std::cout << "READ_REGION: [" << std::fixed << std::setprecision(6) << res.a_start / 1000000 << ", " << res.a_end / 1000000 << " len " << (res.a_end - res.a_start) / 1000000 << "] " << "-> ";
			if (!res.has_mapping)
			{
				std::cout << "GAP";
			}
			else
			{
				std::cout << "RefRegion: total: " << res.b_mappings.size() << " ";
				for (size_t i = 0; i < res.b_mappings.size(); ++i)
				{
					if (i > 0)
						std::cout << " | ";
					std::cout << std::fixed << std::setprecision(6) << "[" << res.b_mappings[i].chrID << ":" << (float)res.b_mappings[i].st / 1000000 << ", "
							  << (float)res.b_mappings[i].ed / 1000000 << "] (Ratio:"
							  << std::fixed << std::setprecision(2) << res.ratio << ")";
				}
			}
			std::cout << std::endl;
		}
	}

	void get_sa_list(bam1_t *&br, int &total_read_len, std::vector<T2T_SA_INFO> &sa_l, faidx_t *c_ref_idx, bool main_is_fwd)
	{
		{
			// get the total read length
			int n_cigar = br->core.n_cigar;
			uint32_t *bam_cigar = bam_get_cigar(br);
			total_read_len = get_total_read_len(bam_cigar, n_cigar);
			sa_l.emplace_back();
			sa_l.back().store_main(br->core.tid, br->core.pos, bam_cigar, n_cigar, total_read_len);
			char *SA_tag_char = bam_get_string_tag(br, "SA");
			if (SA_tag_char != NULL)
			{
				std::vector<std::string> item_tmp;
				split_string(item_tmp, SA_tag_char, ";");
				std::vector<std::string> item_tmp_1;
				for (uint sa_idx = 0; sa_idx < item_tmp.size(); sa_idx++)
				{
					split_string(item_tmp_1, item_tmp[sa_idx].c_str(), ",");
					sa_l.emplace_back();
					sa_l.back().store(sa_idx, item_tmp_1, c_ref_idx, total_read_len, main_is_fwd);
				}
			}
		}
	}

	int t2t_split_alignment_analysis(int argc, char *argv[])
	{
		// input
		char *fa_fn_in = argv[1];
		char *LRS_bamFileName = argv[2];

		faidx_t *c_ref_idx = reference_index_load(fa_fn_in);
		BAM_handler LRS_read;
		LRS_read.init(LRS_bamFileName, NULL, NULL, fa_fn_in);
		Bam_file *c_b = &(LRS_read.file);
		int read_ID = -1;
		while (bam_next(c_b))
		{
			++read_ID;
			// new sig for this read
			bam1_t *br = &(c_b->_brec);
			// filter://basic filter
			if (bam_is_secondary(br))
				continue;
			if (bam_is_supplementary(br))
				continue;
			// from the core position to a bigger position, to find the begin of first read;
			// SA_ANALYSIS
			bool main_is_fwd = bam_is_fwd_strand(br);

			std::vector<T2T_SA_INFO> sa_l;
			int total_read_len = 0;
			get_sa_list(br, total_read_len, sa_l, c_ref_idx, main_is_fwd);
			
			if (total_read_len < 500000)
				continue;

			std::sort(sa_l.begin(), sa_l.end(), T2T_SA_INFO::cmp_by_r_pos);
			if (false == dump_sa_info_to_file(br, total_read_len, read_ID, sa_l, c_ref_idx))
				return -1;

			// not used right now
			if (false)
				calculate_mappings_depth(total_read_len, sa_l);
		}
		return 0;
	}

	int dump_sa_info_to_file(bam1_t *&br, int total_read_len, int read_ID, std::vector<ASM_ANALYSIS_HANDLER::T2T_SA_INFO> &sa_l, faidx_t *c_ref_idx)
	{
		std::unordered_map<int, int> chr_mapping_lengths;
		int max_chrID = -1;
		int max_length = 0;
		get_main_chrID(sa_l, chr_mapping_lengths, max_length, max_chrID);
		if (max_chrID == -1)
			return true;

		std::string filename = std::string(bam_get_qname(br)) + ".txt";
		FILE *output_file = fopen(filename.c_str(), "w");
		if (output_file == nullptr)
		{
			fprintf(stderr, "Error: Cannot open file %s for writing\n", filename.c_str());
			return false;
		}

		//INFO and logs
		fprintf(output_file, "Read mapping full: NAME %s, ", bam_get_qname(br));
		fprintf(output_file, "READ_BEGIN, read_length %d (%.3fM)ID, %d\n", total_read_len, (float)total_read_len / 1000000, read_ID);
		fprintf(output_file, "[STATS] Chromosome mapping lengths:\t");
		for (const auto &entry : chr_mapping_lengths)
		{
			const char *chrName = faidx_iseq(c_ref_idx, entry.first);
			fprintf(output_file, "ChrID: %s, Total Mapping Length: %d;\t",
					chrName, entry.second);
		}
		fprintf(output_file, "\n");
		const char *main_chrName = faidx_iseq(c_ref_idx, max_chrID);
		int main_chrLength = faidx_seq_len(c_ref_idx, main_chrName);

		//CORE for figure:
		fprintf(output_file, "[F_INFO] %s,%d,%s,%d \n", main_chrName, main_chrLength, bam_get_qname(br), total_read_len );

		std::string s = toString_sa_l(sa_l, true);
		fprintf(output_file, "%s ", s.c_str());

		fclose(output_file);
		return true;
	}

	void get_main_chrID(std::vector<ASM_ANALYSIS_HANDLER::T2T_SA_INFO> &sa_l, std::unordered_map<int, int> &chr_mapping_lengths, int &max_length, int &max_chrID)
	{
		for (uint i = 0; i < sa_l.size(); i++)
		{
			if (pass_core_filter(sa_l[i]) == false)
				continue;
			chr_mapping_lengths[sa_l[i].SA_chrID] += sa_l[i].mapping_len_read;
		}

		for (const auto &pair : chr_mapping_lengths)
		{
			if (pair.second > max_length)
			{
				max_length = pair.second;
				max_chrID = pair.first;
			}
		}
	}

	std::string toString_sa_l(std::vector<ASM_ANALYSIS_HANDLER::T2T_SA_INFO> &sa_l, bool is_core)
	{
		std::string result;

		int pre_end_read = 0;
		int pre_ref_end = 0;
		for (uint i = 0; i < sa_l.size(); i++)
		{
			if (is_core && pass_core_filter(sa_l[i]) == false)
				continue;
			if (is_core)
				result += "CORE: " + std::to_string(i) + ", ";
			else
				result += "ALL: " + std::to_string(i) + ", ";
			int Read_gap_size = sa_l[i].sa_in_read_position_bg - pre_end_read;
			int Ref_gap_size = sa_l[i].SA_POS - pre_ref_end;
			result += sa_l[i].toString(Read_gap_size, Ref_gap_size) + "\n";
			pre_end_read = sa_l[i].sa_in_read_position_bg + sa_l[i].mapping_len_read;
			pre_ref_end = sa_l[i].SA_POS + sa_l[i].mapping_len_ref;
		}
		return result;
	}

	bool pass_core_filter(T2T_SA_INFO &sa)
	{
		if (sa.mapping_len_read < 10000)
			return false;
		if (sa.SA_MAPQ < 15)
			return false;
		if (sa.SA_chrID > 24)
			return false;
		return true;
	}
};

#endif /* SVCALLING_CORE_ANALYSIS_ASM_ANALYSIS_HANDLER_HPP_ */
