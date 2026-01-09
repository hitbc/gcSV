/*
 * VNTR_analysis.hpp
 *
 *  Created on: 2025-4-30
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_ANALYSIS_VNTR_ANALYSIS_HPP_
#define SVCALLING_CORE_ANALYSIS_VNTR_ANALYSIS_HPP_

#include <vector>
#include <algorithm>
#include <map>
#include "../NGS_merging.hpp"
#include "cpp_lib/cpp_utils.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/vcf_lib.h"
}

struct vntr_analysis_handler{
	struct SV_in_repeat{
		int SV_id;
		int repeatID;
		int SV_len;
		static inline int cmp_by_repeatID(const SV_in_repeat &a, const SV_in_repeat &b){
			if(a.repeatID != b.repeatID)
				return a.repeatID < b.repeatID;
			else if(a.SV_len != b.SV_len)
				return a.SV_len < b.SV_len;
			else
				return a.SV_id < b.SV_id;
		}
	};

	struct BED_ITEM{
		int chrID;
		int pos;
		int end;

		std::string info1;
		std::string info2;

		int repeat_len;
		std::string repeat_str;
		//for VNTR count:
		int MIN_GCA_sv_len;
		int max_sv_len;

		int repeat_range;

		int bed_original_index = 0;
		int bed_repeat_range_index = 0;

		void load_data(std::vector<std::string> &item_value, faidx_t *c_ref_idx, std::vector<std::string> &item_value_buff){
			chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			pos = atoi(item_value[1].c_str());
			end = atoi(item_value[2].c_str());
			info1 = item_value[3];
			info2 = item_value[4];
			split_string(item_value_buff, info2.c_str(), ")");
			repeat_len = item_value_buff[0].size() - 1;
			repeat_str.clear();
			repeat_str.append(item_value_buff[0].c_str() + 1);
		}

		void init_SV_min_max(){
			MIN_GCA_sv_len = MAX_int32t;
			max_sv_len = -MAX_int32t;
			repeat_range = -1;
		}

		void store_SV_min_max(int SV_len){
			MIN_GCA_sv_len = MIN(MIN_GCA_sv_len, SV_len);
			max_sv_len = MAX(max_sv_len, SV_len);
		}

		void store_repeat_range(){
			if(MIN_GCA_sv_len == MAX_int32t)
				repeat_range = -1;
			else
				repeat_range = (max_sv_len - MIN_GCA_sv_len) / repeat_len;
		}
		void show(){
			fprintf(stderr, "BED1: %d\t%d\t%d\t%s\t%s\t%d\t", chrID, pos, end, info1.c_str(), info2.c_str(), repeat_len);
			fprintf(stderr, "BED2: %d\t%d\t%d\t%d\t%d\n", MIN_GCA_sv_len, max_sv_len, repeat_range, bed_original_index, bed_repeat_range_index);
		}

		static inline int cmp_by_repeat_range(const BED_ITEM &a, const BED_ITEM &b){
			if(a.repeat_range != b.repeat_range)
				return a.repeat_range < b.repeat_range;
			else if(a.chrID != b.chrID)
				return a.chrID < b.chrID;
			else
				return a.pos < b.pos;
		}

		static inline int cmp_by_bed_original_index(const BED_ITEM &a, const BED_ITEM &b){
			if(a.bed_original_index != b.bed_original_index)
				return a.bed_original_index < b.bed_original_index;
			else xassert(0, "");
		}
	};
	#define BED_IDX_BIN_SIZE 1000

	struct BED_IDX{
		struct BLOCK_IDX{
			std::vector<int> overlap_bed_id;
		};

		struct CHR_IDX{
			int chrID;
			int chr_len;
			std::vector<BLOCK_IDX> pos_idx;
			void initBlock(int block_size){
				int total_block_length = chr_len/block_size + 1;
				pos_idx.resize(total_block_length);
			}
		};

		std::vector<CHR_IDX> chr_idx;
		void init(faidx_t *c_ref_idx ){
			//load reference index
			int N_seq = faidx_nseq(c_ref_idx);
			N_seq = MIN(24, N_seq);
			for(int i = 0; i < N_seq; i++){
				chr_idx.emplace_back();
				chr_idx.back().chrID = i;
				chr_idx.back().chr_len = faidx_seq_len(c_ref_idx, faidx_iseq(c_ref_idx, i));
				chr_idx.back().initBlock(BED_IDX_BIN_SIZE);
			}
		}
		void store_bed_item(int bed_id, BED_ITEM & b){
			int bg = b.pos/BED_IDX_BIN_SIZE;
			int ed = b.end/BED_IDX_BIN_SIZE;
			for(int i = bg; i <= ed; i++){
				chr_idx[b.chrID].pos_idx[i].overlap_bed_id.emplace_back(bed_id);
			}
		}

		int search_bed_full_cover_SV(int sv_chr, int sv_pos, int sv_end, std::vector<BED_ITEM> & bed_list){
			int bg = sv_pos/BED_IDX_BIN_SIZE;
			for(int bed_i: chr_idx[sv_chr].pos_idx[bg].overlap_bed_id){
				if(sv_pos >= bed_list[bed_i].pos && sv_end <= bed_list[bed_i].end)
					return bed_i;
			}
			return -1;
		}
	};


	int bed_in_repeat_masker_region(int argc, char *argv[]){
		bool print_log = true;
		//parameters
		char * ref_fn = argv[1];
		char * analysis_bed_fn_in = argv[2];//separate by ','
		char * repeat_masker_region_bed_fn = argv[3];

		faidx_t *c_ref_idx = NULL;
		BED_IDX bed_idx;
		std::vector<BED_ITEM> bed_l;
		//load reference index
		{
			c_ref_idx = reference_index_load(ref_fn);
			bed_idx.init(c_ref_idx);
		}
		//load bed file
		{
			std::vector<std::string> bed_record_l;
			load_string_list_from_file_MAX_line(repeat_masker_region_bed_fn, bed_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &bed_r :bed_record_l){
				split_string(item_value, bed_r.c_str(), "\t");
				bed_l.emplace_back();
				bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
				if(bed_l.back().chrID >= 24){
					bed_l.pop_back();
				}else
					bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
			}
		}
		//load input bed file
		{
			std::vector<std::string> bed_record_l;
			load_string_list_from_file_MAX_line(analysis_bed_fn_in, bed_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &bed_r :bed_record_l){
				split_string(item_value, bed_r.c_str(), "\t");
				if(item_value.size() < 3){
					fprintf(stderr, "BED_in_range: #\n");
				}
				int chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
				int pos = atoi(item_value[1].c_str());
				int end = atoi(item_value[2].c_str());
				if(chrID < 0){
					fprintf(stderr, "BED_in_range: %s\n", item_value[0].c_str());
				}else{
					int bed_search_idx = bed_idx.search_bed_full_cover_SV(chrID, pos, end, bed_l);
					//output:
					if(bed_search_idx == -1)
						fprintf(stderr, "BED_in_range: %d\t%d\t%d\tNA\n", chrID, pos, end);
					else
						fprintf(stderr, "BED_in_range: %d\t%d\t%d\t%s\t%s\n", chrID, pos, end, bed_l[bed_search_idx].info1.c_str(), bed_l[bed_search_idx].info2.c_str());
				}
			}
		}
		return 0;
	}
	int SV_in_simple_repeat(int argc, char *argv[]){
		bool print_log = true;
		//parameters
		char * ref_fn = argv[1];
		char * vcf_fn_in = argv[2];//separate by ','
		char * vntr_region_bed_fn = argv[3];
		char * sample_region_fn = argv[4];

		std::vector<Merging_SV_info> SV_info;
		std::vector<BED_ITEM> bed_l;
		std::vector<SV_in_repeat> sv_repeat;

		std::vector<int> sample_region;

		faidx_t *c_ref_idx = NULL;
		BED_IDX bed_idx;

		//load sample region file:
		{
			std::vector<std::string> sample_record_l;
			load_string_list_from_file_MAX_line(sample_region_fn, sample_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &s_r :sample_record_l){
				split_string(item_value, s_r.c_str(), "\t");
				int main_region_id = atoi(item_value[6].c_str()) - 1;
				sample_region.emplace_back(main_region_id);
			}
		}

		//load reference index
		{
			c_ref_idx = reference_index_load(ref_fn);
			bed_idx.init(c_ref_idx);
		}

		//load vcf binary
		{
			std::vector<std::string> vcf_record_l;
			load_string_list_from_file_MAX_line(vcf_fn_in, vcf_record_l, -1);//todo::
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &vcf_r :vcf_record_l){
				split_string(item_value, vcf_r.c_str(), "\t");
				SV_info.emplace_back();
				SV_info.back().load_data(item_value, item_value_buff);
			}
		}

		//load bed file
		{
			std::vector<std::string> bed_record_l;
			load_string_list_from_file_MAX_line(vntr_region_bed_fn, bed_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &bed_r :bed_record_l){
				split_string(item_value, bed_r.c_str(), "\t");
				bed_l.emplace_back();
				bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
				if(bed_l.back().chrID >= 24){
					bed_l.pop_back();
				}else
					bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
			}
		}

		//search vcf in bed
		{
			for(int SV_id = 0; SV_id < SV_info.size(); SV_id++){
				Merging_SV_info &cur_sv_info = SV_info[SV_id];
				int sv_chr = cur_sv_info.chrID;
				int sv_pos = cur_sv_info.pos;
				int sv_end = (cur_sv_info.SV_len > 0)?sv_pos:(sv_pos - cur_sv_info.SV_len);
				int bed_search_idx = bed_idx.search_bed_full_cover_SV(sv_chr, sv_pos, sv_end, bed_l);
				if(bed_search_idx != -1){
					//analysis if the SV is a true VNTR;
					std::string repeat_string = bed_l[bed_search_idx].repeat_str;
					int max_dis = bed_l[bed_search_idx].repeat_len * 5;
					std::string check_str;
					if(cur_sv_info.ALT_str.size() > cur_sv_info.REF_str.size())
						check_str = cur_sv_info.ALT_str;
					else
						check_str = cur_sv_info.REF_str;

					bool pass_filter = false;
					{
						int total_gap_length = 0;
						int found_pos = 0;
						int repeat_number = 0;
						while(found_pos != -1){
							int found = check_str.find(repeat_string, found_pos + 1);
							//fprintf(stderr, "%ld\t", found);
							int gap_len = 0;
							if(found != -1){
								repeat_number++;
								gap_len = found-found_pos;
							}
							else
								gap_len = check_str.size()-found_pos;
							if(gap_len > max_dis){
								//fprintf(stderr, "GAP:%d, %d, %d\n", found_pos, found, gap_len);
								total_gap_length += gap_len;
							}
							found_pos = found;
						}
						int final_sv_len = check_str.size() - 1 - total_gap_length;
						fprintf(stderr, "Final repeat_number: %d total_gap_length %d final_sv_len %d %s %s %s\n", repeat_number, total_gap_length, final_sv_len, cur_sv_info.REF_str.c_str(), cur_sv_info.ALT_str.c_str(), repeat_string.c_str());
						if(repeat_number != 0 && final_sv_len >= 30 && final_sv_len*2 > check_str.size()){
							pass_filter = true;
						}else
							pass_filter = false;
					}

					if(pass_filter){
						sv_repeat.emplace_back();
						sv_repeat.back().SV_id = SV_id;
						sv_repeat.back().repeatID = bed_search_idx;
						sv_repeat.back().SV_len = cur_sv_info.SV_len;
					}
	//				if(12422321 == cur_sv_info.pos){
	//					fprintf(stderr, " ");
	//				}
				}else{
					//cur_sv_info.show();
					//fprintf(stderr, "NO data\n");
				}
			}

			for(BED_ITEM & b : bed_l)
				b.init_SV_min_max();

			//analysis:
			std::sort(sv_repeat.begin(), sv_repeat.end(), SV_in_repeat::cmp_by_repeatID);
			int old_repeatID = -1;
			for(SV_in_repeat & sr : sv_repeat){
				if(old_repeatID != sr.repeatID){
					bed_l[sr.repeatID].show();
				}
				old_repeatID = sr.repeatID;
				SV_info[sr.SV_id].show();
				//set the max and min SV size:
				bed_l[sr.repeatID].store_SV_min_max(SV_info[sr.SV_id].SV_len);
			}

			//set the repeat_range and the original_index for each BED
			int bed_original_index = 0;
			for(BED_ITEM & b : bed_l){
				b.store_repeat_range();
				b.bed_original_index = bed_original_index;
				bed_original_index ++;
			}

			//sort by repeat range
			std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_repeat_range);
			//set the repeat range index
			int bed_repeat_range_index = 0;
			for(BED_ITEM & b : bed_l){
				if(b.repeat_range != -1){
					b.bed_repeat_range_index = bed_repeat_range_index;
					bed_repeat_range_index ++;
				}
			}

			//sort by original index
			std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_bed_original_index);

			std::vector<int> sample_region_counter;
			sample_region_counter.resize(5);
			//set the repeat range index for all SVs
			old_repeatID = -1;
			for(SV_in_repeat & sr : sv_repeat){
				if(old_repeatID != sr.repeatID){
					fprintf(stderr, "bed_l show final\t");
					bed_l[sr.repeatID].show();
				}
				old_repeatID = sr.repeatID;
				int repeat_unit_n = SV_info[sr.SV_id].SV_len/bed_l[sr.repeatID].repeat_len;
				fprintf(stderr, "SV_in_range: %d\t%d\t%d\t", bed_l[sr.repeatID].bed_original_index, bed_l[sr.repeatID].bed_repeat_range_index, repeat_unit_n);
				SV_info[sr.SV_id].count_sample_region(sample_region, sample_region_counter);
				SV_info[sr.SV_id].show();
			}
		}
		return 0;
	}

	int vntr_analysis(int argc, char *argv[]){
		bool print_log = true;
		int SV_MERGE_IDX_B_SIZE = 30;
		//parameters
		char * ref_fn = argv[1];
		char * single_sample_vcf_list = argv[2];//separate by ','
		char * vntr_region_bed_fn = argv[3];
		char * sample_region_fn = argv[4];

		faidx_t *c_ref_idx = NULL;
		std::vector<Merging_SV_info> final_SV_l;
		std::vector<BED_ITEM> bed_l;
		std::vector<SV_in_repeat> sv_repeat;

		std::vector<std::string> sample_file_name;//vcf file name for each samples
		std::vector<int> sample_region;

		BED_IDX bed_idx;//simple repeat BED file
		//char *joint_vcf_fn = argv[1];
		std::vector<SV_store_index> SV_idx;
		NGS_SV_MERGING_HANDLER ngs_sv_merging_h;
		ngs_sv_merging_h.init();
			//index
			//load sample list
		//S1: load reference index
		{
			c_ref_idx = reference_index_load(ref_fn);
			bed_idx.init(c_ref_idx);
			ngs_sv_merging_h.init_SV_index(SV_idx, c_ref_idx);
		}

		//load SVs from all single sample vcf data:
		{
			//S2: load sample file names
			load_string_list_from_file_MAX_line(single_sample_vcf_list, sample_file_name, 1000);
			//load single sample vcf
			int total_SV_num = 0;
			char *c_sv_type = (char *)malloc(1000);
			uint32_t * GT_output_buff = (uint32_t *) xcalloc(10,sizeof(uint32_t));
			for(uint sample_ID = 0; sample_ID < sample_file_name.size(); sample_ID++){
				//if(sample_ID > 10) break;
				fprintf(stderr, "Current handle %d: %s\n",sample_ID, sample_file_name[sample_ID].c_str());
				BCF_FILE input_vcf;

				VCF_open_read(&input_vcf, sample_file_name[sample_ID].c_str());
				bcf_hdr_t *vcf_header = input_vcf.header;
				int new_SV_num = 0;
				//load each SVs
				do
				{
					bcf1_t *c_r = &( input_vcf.r);
					bcf_unpack(c_r, BCF_UN_INFO);
					int sv_chr_ID = c_r->rid;
					//vcf_read.header->
					//faidx_iseq(c_ref_idx, SV_CHR_ID);
					int chrID = faidx_get_chrID(c_ref_idx, bcf_hdr_id2name(vcf_header, sv_chr_ID), NULL, 0);
					if(chrID >= 24)
						continue;
					vcf_get_sv_type(vcf_header, c_r, c_sv_type);
					if(strcmp(c_sv_type, "BND") == 0)
						continue;
					if(*c_r->d.flt != 0 && *c_r->d.flt != 13)
						continue;
					int SV_POS = c_r->pos;
					int SV_length = 0;
					vcf_get_sv_LENGTH(input_vcf.header, c_r, &SV_length);
					std::string REF_str; std::string ALT_str;
					REF_str.append(c_r->d.allele[0]);
					ALT_str.append(c_r->d.allele[1]);
					int mdat;
					//dat = 0/0	0/2	2/2	2/4	4/4	4/6	5/6	5/7	6/6	7/7
					//GT  = ./.	./0	0/0	0/1	1/1	1/2	1/2	1|2	2/2	2|2
					bcf_get_genotypes(vcf_header,c_r,(void**)&GT_output_buff, &mdat);
					//store the SVs to SV list
					bool is_new = true;
					bool high_complex_region = (SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.size() > 100);
					if(high_complex_region){
						fprintf(stderr, " ");
					}
					for(Merging_SV_info & sv: SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs){
						if(sv.sameWith(SV_POS, ALT_str, REF_str,high_complex_region, SV_length)){
							sv.supportSAMPLE_ID.emplace_back();
							sv.supportSAMPLE_ID.back().store(sample_ID, GT_output_buff[0], GT_output_buff[1], c_r->qual);
							is_new = false;
							break;
						}
					}
					if(is_new){
						new_SV_num ++;
						SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.emplace_back();
						Merging_SV_info & sv = SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.back();
						sv.store_basic(SV_POS, ALT_str, REF_str, SV_length);
						sv.supportSAMPLE_ID.emplace_back();
						sv.supportSAMPLE_ID.back().store(sample_ID, GT_output_buff[0], GT_output_buff[1], c_r->qual);
					}
				}while(VCF_next(&input_vcf));//read one
				total_SV_num += new_SV_num;
				fprintf(stderr, "\n new_SV_num %d total_SV_num %d\n", new_SV_num, total_SV_num);
				bcf_close(input_vcf.file);
			}

			//output the final result
			{
				//std::vector<SV_INFO_vcf_add_alt_string> SV_info;
				for(uint chrID = 0; chrID < SV_idx.size(); chrID++){
					for(uint b_id = 0; b_id < SV_idx[chrID].b.size(); b_id++){
						//for each block
						std::vector<Merging_SV_info> SV_l = SV_idx[chrID].b[b_id].SVs;
						// sort:
						std::sort(SV_l.begin(), SV_l.end(), Merging_SV_info::cmp_by_pos);
						final_SV_l.emplace_back();
						for(uint SV_idx_in_B = 0; SV_idx_in_B < SV_l.size(); SV_idx_in_B++){
							final_SV_l.emplace_back();
							std::swap(final_SV_l.back(), SV_l[SV_idx_in_B]);
						}
					}
				}
			}
		}

		//load sample region file:
		{
			std::vector<std::string> sample_record_l;
			load_string_list_from_file_MAX_line(sample_region_fn, sample_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &s_r :sample_record_l){
				split_string(item_value, s_r.c_str(), "\t");
				int main_region_id = atoi(item_value[6].c_str()) - 1;
				sample_region.emplace_back(main_region_id);
			}
		}

		//load vntr bed file
		{
			std::vector<std::string> bed_record_l;
			load_string_list_from_file_MAX_line(vntr_region_bed_fn, bed_record_l, -1);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &bed_r :bed_record_l){
				split_string(item_value, bed_r.c_str(), "\t");
				bed_l.emplace_back();
				bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
				if(bed_l.back().chrID >= 24){
					bed_l.pop_back();
				}else
					bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
			}
		}

		//search vcf in bed
		{
			for(int SV_id = 0; SV_id < final_SV_l.size(); SV_id++){
				Merging_SV_info &cur_sv_info = final_SV_l[SV_id];
				int sv_chr = cur_sv_info.chrID;
				int sv_pos = cur_sv_info.pos;
				int sv_end = (cur_sv_info.SV_len > 0)?sv_pos:(sv_pos - cur_sv_info.SV_len);
				int bed_search_idx = bed_idx.search_bed_full_cover_SV(sv_chr, sv_pos, sv_end, bed_l);
				if(bed_search_idx != -1){
					//analysis if the SV is a true VNTR;
					std::string repeat_string = bed_l[bed_search_idx].repeat_str;
					int max_dis = bed_l[bed_search_idx].repeat_len * 5;
					std::string check_str;
					if(cur_sv_info.ALT_str.size() > cur_sv_info.REF_str.size())
						check_str = cur_sv_info.ALT_str;
					else
						check_str = cur_sv_info.REF_str;

					bool pass_filter = false;
					{
						int total_gap_length = 0;
						int found_pos = 0;
						int repeat_number = 0;
						while(found_pos != -1){
							int found = check_str.find(repeat_string, found_pos + 1);
							//fprintf(stderr, "%ld\t", found);
							int gap_len = 0;
							if(found != -1){
								repeat_number++;
								gap_len = found-found_pos;
							}
							else
								gap_len = check_str.size()-found_pos;
							if(gap_len > max_dis){
								//fprintf(stderr, "GAP:%d, %d, %d\n", found_pos, found, gap_len);
								total_gap_length += gap_len;
							}
							found_pos = found;
						}
						int final_sv_len = check_str.size() - 1 - total_gap_length;
						fprintf(stderr, "Final repeat_number: %d total_gap_length %d final_sv_len %d %s %s %s\n", repeat_number, total_gap_length, final_sv_len, cur_sv_info.REF_str.c_str(), cur_sv_info.ALT_str.c_str(), repeat_string.c_str());
						if(repeat_number != 0 && final_sv_len >= 30 && final_sv_len*2 > check_str.size()){
							pass_filter = true;
						}else
							pass_filter = false;
					}

					if(pass_filter){
						sv_repeat.emplace_back();
						sv_repeat.back().SV_id = SV_id;
						sv_repeat.back().repeatID = bed_search_idx;
						sv_repeat.back().SV_len = cur_sv_info.SV_len;
					}
				}else{/*DO NOTHING*/}
			}

			for(BED_ITEM & b : bed_l)
				b.init_SV_min_max();

			//analysis:
			std::sort(sv_repeat.begin(), sv_repeat.end(), SV_in_repeat::cmp_by_repeatID);
			int old_repeatID = -1;
			for(SV_in_repeat & sr : sv_repeat){
				if(old_repeatID != sr.repeatID){
					bed_l[sr.repeatID].show();
				}
				old_repeatID = sr.repeatID;
				final_SV_l[sr.SV_id].show();
				//set the max and min SV size:
				bed_l[sr.repeatID].store_SV_min_max(final_SV_l[sr.SV_id].SV_len);
			}

			//set the repeat_range and the original_index for each BED
			int bed_original_index = 0;
			for(BED_ITEM & b : bed_l){
				b.store_repeat_range();
				b.bed_original_index = bed_original_index;
				bed_original_index ++;
			}

			//sort by repeat range
			std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_repeat_range);
			//set the repeat range index
			int bed_repeat_range_index = 0;
			for(BED_ITEM & b : bed_l){
				if(b.repeat_range != -1){
					b.bed_repeat_range_index = bed_repeat_range_index;
					bed_repeat_range_index ++;
				}
			}

			//sort by original index
			std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_bed_original_index);

			std::vector<int> sample_region_counter;
			sample_region_counter.resize(5);
			//set the repeat range index for all SVs
			old_repeatID = -1;
			for(SV_in_repeat & sr : sv_repeat){
				if(old_repeatID != sr.repeatID){
					fprintf(stderr, "bed_l show final\t");
					bed_l[sr.repeatID].show();
				}
				old_repeatID = sr.repeatID;
				int repeat_unit_n = final_SV_l[sr.SV_id].SV_len/bed_l[sr.repeatID].repeat_len;
				fprintf(stderr, "SV_in_range: %d\t%d\t%d\t", bed_l[sr.repeatID].bed_original_index, bed_l[sr.repeatID].bed_repeat_range_index, repeat_unit_n);
				final_SV_l[sr.SV_id].count_sample_region(sample_region, sample_region_counter);
				final_SV_l[sr.SV_id].show();
			}
		}

		fprintf(stderr, "\nEND\n");
		return 0;
	}
};



#endif /* SVCALLING_CORE_ANALYSIS_VNTR_ANALYSIS_HPP_ */
