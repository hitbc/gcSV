/*
 * NGS_merging.hpp
 *
 *  Created on: 2025年1月8日
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_NGS_MERGING_HPP_
#define SVCALLING_CORE_NGS_MERGING_HPP_

#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include "cpp_lib/cpp_utils.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/vcf_lib.h"
}

struct Merging_SV_info{
	struct SAMPLE_info{
		int ID;
		uint8_t GT1;
		uint8_t GT2;
		int16_t QUAL;
		void store(	int ID,	int GT1, int GT2, int16_t QUAL){
			this->ID = ID;
			this->GT1 = GT1;
			this->GT2 = GT2;
			this->QUAL = QUAL;
		}
	};

	int chrID;
	int pos;
	int SV_len;
	std::string ALT_str;
	std::string REF_str;
	std::vector<SAMPLE_info> supportSAMPLE_ID;//todo:: add info, like GT and QUAL

	void store_basic(int pos, std::string &ALT_str, std::string &REF_str, int SV_len){
		this->pos = pos;
		this->SV_len = SV_len;
		this->ALT_str = ALT_str;
		this->REF_str = REF_str;
	}

	void load_data(std::vector<std::string> &item_value, std::vector<std::string> &buff){
		split_string(buff, item_value[0].c_str(), ":");
		chrID = atoi(buff[0].c_str());
		pos = atoi(buff[1].c_str());
		SV_len = atoi(item_value[1].c_str());
		REF_str = item_value[2];
		ALT_str = item_value[3];
		for(uint i = 4; i < item_value.size(); i++){
			split_string(buff, item_value[i].c_str(), ":");
			int base = 0;
			if(i == 4)
				base = 1;
			supportSAMPLE_ID.emplace_back();
			supportSAMPLE_ID.back().store(atoi(buff[base].c_str()), atoi(buff[base + 1].c_str()), atoi(buff[base + 2].c_str()), atoi(buff[base + 3].c_str()));
		}
	}

	static inline int cmp_by_pos(const Merging_SV_info &a, const Merging_SV_info &b){
		if(a.pos != b.pos)
			return a.pos < b.pos;
		return a.SV_len < b.SV_len;
	}

	void show(){
		fprintf(stderr, "SV: %d\t%d\t%d\t%s\t%s\n", chrID, pos, SV_len, REF_str.c_str(), ALT_str.c_str());
	}

	void count_sample_region(std::vector<int> &sample_region, std::vector<int> &sample_region_counter){
		for(uint i = 0;  i< sample_region_counter.size(); i++){
			sample_region_counter[i] = 0;
		}

		for(SAMPLE_info &s : supportSAMPLE_ID){
			sample_region[s.ID];
			if(s.GT1 == 4)
				sample_region_counter[sample_region[s.ID]]++;
			if(s.GT2 == 4)
				sample_region_counter[sample_region[s.ID]]++;
		}
		for(uint i = 0;  i< sample_region_counter.size(); i++){
			fprintf(stderr, "SPOP%d=%d;", i, sample_region_counter[i]);
		}

	}


	bool sameWith(int cmp_pos, 	std::string &cmp_ALT_str, std::string &cmp_REF_str, bool high_complex_region, int cmp_SV_length){
		if(high_complex_region){
			float size_diff = ((float)cmp_ALT_str.size())/ALT_str.size();
			return ((cmp_pos == pos ) && (size_diff > 0.9) && (size_diff < 1.1));
		}
		else if(ALT_str.size() < 400 && cmp_ALT_str.size() < 400)
			return ((cmp_pos == pos ) && (cmp_SV_length == SV_len));
		else{
			float size_diff = ((float)cmp_SV_length)/SV_len;
			return ((cmp_pos == pos ) && (size_diff > 0.9) && (size_diff < 1.1));
		}
	}

	bool sameWithM2(int cmp_pos, 	std::string &cmp_ALT_str, std::string &cmp_REF_str, bool high_complex_region, int SV_length){
		if(high_complex_region){
			float size_diff = ((float)cmp_ALT_str.size())/ALT_str.size();
			return ((cmp_pos == pos ) && (size_diff > 0.9) && (size_diff < 1.1));
		}
		else if(ALT_str.size() < 100 && cmp_ALT_str.size() < 100)
			return ((cmp_pos == pos ) && (cmp_ALT_str == ALT_str) && (cmp_REF_str == REF_str));
		else if(ALT_str.size() < 400 && cmp_ALT_str.size() < 400)
			return ((cmp_pos == pos ) && (cmp_ALT_str.size() == ALT_str.size()) && (cmp_REF_str.size() == REF_str.size()));
		else{
			float size_diff = ((float)cmp_ALT_str.size())/ALT_str.size();
			return ((cmp_pos == pos ) && (size_diff > 0.9) && (size_diff < 1.1));
		}
	}
};

#define SV_MERGE_IDX_B_SIZE 30
struct SV_merge_block{
	int pos_begin;
	std::vector<Merging_SV_info> SVs;
};

struct SV_store_index{
	uint32_t ref_length;
	const char * chrName;
	int chrID;
	int total_block_length = 0;
	std::vector<SV_merge_block> b;

	SV_store_index(int chrID, uint32_t ref_length, const char * chrName){
		this->ref_length = ref_length;
		this->chrName = chrName;
		this->chrID = chrID;
	}

	SV_store_index(){
		this->ref_length = 0;
		this->chrName = NULL;
		this->chrID = 0;
	}

	void initBlock(int block_size){
		int total_block_length = ref_length/block_size + 1;
		b.resize(total_block_length);
		int position_bg_cur = 0;
		for(SV_merge_block & r : b){
			r.pos_begin = position_bg_cur;
			position_bg_cur += block_size;
		}
	}
};

struct NGS_SV_MERGING_HANDLER{
	void static init_SV_index(std::vector<SV_store_index> &SV_idx, faidx_t *c_ref_idx){
		//load reference index
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int i = 0; i < N_seq; i++){
			char fn[1024];
			sprintf(fn, "%d_dump_file.bin", i);
			SV_idx.emplace_back();
			SV_idx.back().chrID = i;
			SV_idx.back().ref_length = faidx_seq_len(c_ref_idx, faidx_iseq(c_ref_idx, i));
			SV_idx.back().initBlock(SV_MERGE_IDX_B_SIZE);
		}
	}

	int SV_merging(int argc, char *argv[]){
		//char *joint_vcf_fn = argv[1];
		faidx_t *c_ref_idx = reference_index_load("/media/fenghe/Data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa");
		const char *joint_vcf_fn = "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/filltag_S5.vcf";
		std::vector<SV_store_index> SV_idx;
		std::vector<std::string> sample_info;
		//index
		init_SV_index(SV_idx, c_ref_idx);
		//load sample list
		{
			std::vector<std::string> joint_vcf_data;
			load_string_list_from_file_MAX_line(joint_vcf_fn, joint_vcf_data, 1000);
			size_t line_num = joint_vcf_data.size();//skip the header line
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_2;
			std::map <int, int> SVLEN_counter;
			for(uint64_t i = 0; i < line_num; i++){
				//show header
				const char * line_data = joint_vcf_data[i].c_str();
				if(line_data[0] == '#' && line_data[1] == '#'){}
				else if(line_data[0] == '#'){
					split_string(item_value, joint_vcf_data[i].c_str(), "\t");//skip the header line
					for(uint j = 9; j < item_value.size(); j ++){
						sample_info.emplace_back(item_value[j]);
					}
					break;
				}
			}
		}

		//recover from single sample vcf data:
		{
			//load single sample vcf
			int total_SV_num = 0;
			char *c_sv_type = (char *)malloc(1000);

			uint32_t * GT_output_buff = (uint32_t *) xcalloc(10,sizeof(uint32_t));

			for(uint sample_ID = 0; sample_ID < sample_info.size(); sample_ID++){
				//if(sample_ID > 10) break;
				char vcf_fn[1024];
				sprintf(vcf_fn, "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/sorted_vcf/%s_sv_spa.sort.vcf.gz", sample_info[sample_ID].c_str());
				fprintf(stderr, "Current handle %d: %s\n",sample_ID, vcf_fn);
				BCF_FILE input_vcf;

				VCF_open_read(&input_vcf, vcf_fn);
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
						for(uint SV_idx_in_B = 0; SV_idx_in_B < SV_l.size(); SV_idx_in_B++){
							Merging_SV_info &sv =SV_l[SV_idx_in_B];
							fprintf(stdout, "%d:%d\t%d\t%s\t%s\tSA:", chrID, sv.pos, sv.SV_len ,sv.REF_str.c_str(), sv.ALT_str.c_str());//chrID, pos, name
							for(Merging_SV_info::SAMPLE_info &sid :sv.supportSAMPLE_ID){
								fprintf(stdout, "%d:%d:%d:%d\t", sid.ID,sid.GT1,sid.GT2,sid.QUAL);
							}
							fprintf(stdout, "\n");
						}
					}
				}
			}
		}
		fprintf(stderr, "\nEND\n");
		return 0;
	}
};



#endif /* SVCALLING_CORE_NGS_MERGING_HPP_ */
