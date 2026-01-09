/*
 * analysis.hpp
 *
 *  Created on: 2025-4-28
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_ANALYSIS_HPP_
#define SVCALLING_CORE_ANALYSIS_HPP_

#include <vector>
#include <algorithm>
#include <map>
#include <algorithm>

#include "cpp_lib/cpp_utils.hpp"
#include "../ReadHandler.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/vcf_lib.h"
}


struct Alignment_Entropy_Analysis{

	faidx_t * c_ref_idx;
	std::vector<R_region> region_l;

	struct BP_signal{
		int read_id;
		int position;
	};
	std::vector<BP_signal> BP_l;

	struct Read_INFO{
		int read_st;
		int read_ed;
	};
	std::vector<Read_INFO> RD_l;

	struct Entropy_window{
		std::set<int> BP_read_id;
		int total_read_n;
		int st_pos;
		int ed_pos;
		int sv_hap_n_INS;
		int sv_hap_n_DEL;
	};
	std::vector<Entropy_window> Entropy_window_l;

	void load_reion_list(const char * region_bed_in){
		region_l.clear();
		//load region list[BED file]
		std::vector<std::string> load_line;
		std::vector<std::string> item_value;
		load_string_list_from_file(region_bed_in, load_line);
		size_t line_num = load_line.size();//skip the header line
		//skip the head line
		for(uint64_t i = 1; i < line_num; i++){
			split_string(item_value, load_line[i].c_str(), "\t");//skip the header line
			region_l.emplace_back();
			int chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			if(chrID <= 23){
				region_l.back().chr_ID =chrID;
				region_l.back().st_pos = atoi(item_value[1].c_str());
				region_l.back().ed_pos = atoi(item_value[2].c_str());
			}
		}
	}

	static void cigar2Path(const uint32_t* bam_cigar, const unsigned n_cigar, path_segment* path){
		for (unsigned int i = 0; i < n_cigar; ++i){
			path[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			path[i].type   = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		}
	}

	bool load_read_seq(int st_pos, int end_pos, std::vector<char>& seq_S, bam1_t* _bp){
		uint8_t *bam_seq = bam_get_seq(_bp);
		end_pos = MIN(end_pos, _bp->core.l_qseq);
		seq_S.resize(end_pos - st_pos);
		char *seq = &(seq_S[0]);
		int n_count = 0;
		while(st_pos < end_pos)
		{
			uint8_t c_c = bam_seqi(bam_seq, st_pos);
			switch(c_c)
			{
			case 1: *seq++ ='A'; break;
			case 2: *seq++ ='C'; break;
			case 4: *seq++ ='G'; break;
			case 8: *seq++ ='T'; break;
			default: *seq++ =0; n_count ++; break;
			}
			st_pos++;
		}
		if(n_count > 10)
			return false;
		return true;
	}

	BAM_handler LRS_read;
	std::vector<char> seq_S;
	bool print_log;

	void collect_bam_signals_Entropy(int MIN_sv_len){
		//handler each region:
		//char *ref = fai_fetch(c_ref_idx, reg, &load_len);
		Bam_file *c_b = &(LRS_read.file);
		seq_S.clear();
		//init the BP & RD list
		BP_l.clear();
		RD_l.clear();
		//load all reads and analysis:
		int read_id = 0;
		while (bam_next(c_b)) {
			//new sig for this read
			bam1_t *br = &(c_b->_brec);
			//if (print_log)fprintf(stderr, "Read Begin NAME %s \n", bam_get_qname(br));
			//filter:
			//basic filter
			if (bam_is_secondary(br))
				continue;
			if (br->core.qual < 4)
				continue;
			//supplementary-min length check
			if (bam_is_supplementary(br) && br->core.l_qseq < 1200)
				continue;
			read_id++;
			//get seq,  get qual:
			load_read_seq(0, br->core.l_qseq, seq_S, br);
			//uint8_t * qual = bam_get_qual(br);
			//get cigar

			std::vector<path_segment> cigar_P;
			cigar_P.resize(br->core.n_cigar);
			cigar2Path( bam_get_cigar(br), br->core.n_cigar, &(cigar_P[0]));
			int seq_i = 0;
			int ref_i = 0;
			//get the reference sequence:
			//char *tseq = ref;//??????
			//char *qseq = &(seq_S[0]);
			for (uint i = 0; i < br->core.n_cigar; i++) {
				path_segment *p = &(cigar_P[i]);
				if (p->type == align_t::CIGAR_SOFT_CLIP)			//soft clip:
				{
					seq_i += p->length;
				} else if (p->type == align_t::CIGAR_INSERT) {
					if (p->length >= (unsigned) MIN_sv_len) {
						BP_l.emplace_back();
						BP_l.back().read_id = read_id;
						BP_l.back().position = ref_i + br->core.pos;
					}
					seq_i += p->length;
				} else if (p->type == align_t::CIGAR_DELETE) {
					if (p->length >= (unsigned) MIN_sv_len) {
						BP_l.emplace_back();
						BP_l.back().read_id = read_id;
						BP_l.back().position = ref_i + br->core.pos;
						BP_l.emplace_back();
						BP_l.back().read_id = read_id;
						BP_l.back().position = ref_i + br->core.pos + p->length;
					}
					ref_i += p->length;
				} else if (p->type == align_t::CIGAR_MATCH
						|| p->type == align_t::CIGAR_SEQ_MATCH
						|| p->type == align_t::CIGAR_SEQ_MISMATCH) {
					for(uint i = 0; i < p->length; i++, ref_i++, seq_i++){
						//Do nothing
					}
				} else if (p->type == align_t::CIGAR_HARD_CLIP) {
					//DO nothing
				} else {
					//DO nothing
				}
			}

			RD_l.emplace_back();
			RD_l.back().read_st = br->core.pos;
			RD_l.back().read_ed = br->core.pos + ref_i;
		}
	}

	struct VCF_BP_INFO{
		int chrID;
		int position;
		int hap_n;
		int int_type;
		static inline int cmp_by_pos(const VCF_BP_INFO &a, const VCF_BP_INFO &b){
			//var basic
			if(a.chrID != b.chrID) 	return a.chrID < b.chrID;
			return a.position < b.position;
		}
	};
	struct VCF_BP_INFO_IDX{
		int chrID;
		int st_position;
		int ed_position;
		int VCF_BP_INFO_idx_st;
		int VCF_BP_INFO_idx_ed;
		void add(int idx){
			if(VCF_BP_INFO_idx_st == -1){
				VCF_BP_INFO_idx_st = idx;
				VCF_BP_INFO_idx_ed = idx + 1;
			}else{
				VCF_BP_INFO_idx_st = MIN(VCF_BP_INFO_idx_st, idx);
				VCF_BP_INFO_idx_ed = MAX(VCF_BP_INFO_idx_ed, idx + 1);
			}
		}
	};

	std::vector<VCF_BP_INFO> vcf_bp_info_l;
	int vcf_bp_info_idx_window_len;
	std::vector<std::vector<VCF_BP_INFO_IDX>> vcf_bp_info_idx_l;

	int load_all_vcfs(const char * GS_VCF_FileName, int min_SV_LEN, int window_step_len){
		//load vcf files and index
		BCF_FILE vcf_read;//vcf for read
		VCF_open_read(&vcf_read, GS_VCF_FileName);//open for read
		//hts_idx_t * vcf_idx = bcf_index_load2(GS_VCF_FileName, NULL);
		bcf_hdr_t *bcf_header = vcf_read.header;
		char *c_sv_type = (char *)malloc(1000);
		//int SV_CHR_ID; int SV_POS; int SV_length; int SV_END;
		char *GT[1]; char GT_1[3]; GT[0] = GT_1;

        //int tid = bcf_hdr_name2id(bcf_header, chr_name);
        //if ( tid ==-1 ) return -1;    // the sequence not present in this file
        //hts_itr_t *reader_itr = bcf_itr_queryi(bcf_idx, tid, r.st_pos, r.ed_pos+1);

		//init the index
		int N_seq = faidx_nseq(c_ref_idx);
		vcf_bp_info_idx_window_len = 10000;
		N_seq = MIN(26, N_seq);
		//m-alloc for all data
		for(int i = 0; i < N_seq; i++){
			const char * chrName = faidx_iseq(c_ref_idx, i);
			int chrLength = faidx_seq_len(c_ref_idx, chrName);
			vcf_bp_info_idx_l.emplace_back();
			int total_window_size = (chrLength/vcf_bp_info_idx_window_len) + 1;
			std::vector<VCF_BP_INFO_IDX> & cur_idx = vcf_bp_info_idx_l.back();
			cur_idx.resize(total_window_size);

			for(int j = 0; j < total_window_size; j++){
				cur_idx[j].chrID = i;
				cur_idx[j].st_position = j*vcf_bp_info_idx_window_len;
				cur_idx[j].ed_position = (j+1)*vcf_bp_info_idx_window_len;
				cur_idx[j].VCF_BP_INFO_idx_st = cur_idx[j].VCF_BP_INFO_idx_ed = -1;
			}
		}

		do{
			bcf1_t *c_r = &( vcf_read.r);
			if(c_r->d.flt != NULL)
				*c_r->d.flt = 0;
			//unpack the vcf data to get the Filter
			bcf_unpack(c_r, BCF_UN_INFO);
			//vcf filters
			if(c_r->d.flt != NULL && *c_r->d.flt != 0)//filter: PASS
				continue;

			const char * bcf_chr_name = bcf_hdr_id2name(bcf_header, c_r->rid);
			int fa_chr_ID = faidx_get_chrID(c_ref_idx, bcf_chr_name, NULL, -1);
			int SV_POS = c_r->pos; int SV_length = -1; int SV_END = -1;
			//vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
			vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
			vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

			//if(SV_POS < 893700) continue;
			//skip * SVs:
			if(c_r->d.allele[1][0] == '*')		continue;

			if(SV_length < min_SV_LEN)
				continue;
			int int_type = 0;
			if(strcmp(c_sv_type, "INS") == 0)	{SV_END = SV_POS; int_type = 0;}
			else								{SV_END = SV_POS + SV_length; int_type = 1;}

			vcf_get_sv_GT(bcf_header, c_r, GT);
			//store data:
			int hap_number = 0;
			if(GT[0][0] == 4)			hap_number ++;
			if(GT[0][1] == 5)			hap_number ++;
			//store the data:

			vcf_bp_info_l.emplace_back();
			vcf_bp_info_l.back().chrID = fa_chr_ID;
			vcf_bp_info_l.back().position = SV_POS;
			vcf_bp_info_l.back().hap_n = hap_number;
			vcf_bp_info_l.back().int_type = int_type;
			if(SV_END % window_step_len != SV_POS % window_step_len){
				vcf_bp_info_l.emplace_back();
				vcf_bp_info_l.back().chrID = fa_chr_ID;
				vcf_bp_info_l.back().position = SV_END;
				vcf_bp_info_l.back().hap_n = hap_number;
				vcf_bp_info_l.back().int_type = int_type;
			}
		}while(VCF_next(&vcf_read));//read one
		//close VCFs
		bcf_close(vcf_read.file);
		//sort all the records
		std::sort(vcf_bp_info_l.begin(), vcf_bp_info_l.end(), VCF_BP_INFO::cmp_by_pos);
		//build index for all the vcf records
		for(uint i = 0; i < vcf_bp_info_l.size(); i++){
			int chrID = vcf_bp_info_l[i].chrID;
			int idx_wb_id = vcf_bp_info_l[i].position / vcf_bp_info_idx_window_len;
			vcf_bp_info_idx_l[chrID][idx_wb_id].add(i);
		}

		for(int i = 0; i < N_seq; i++){
			//const char * chrName = faidx_iseq(c_ref_idx, i);
			//int chrLength = faidx_seq_len(c_ref_idx, chrName);
			//vcf_bp_info_idx_l.emplace_back();
			//int total_window_size = (chrLength/vcf_bp_info_idx_window_len) + 1;
			std::vector<VCF_BP_INFO_IDX> & cur_idx = vcf_bp_info_idx_l[i];
			for(uint j = 1; j < cur_idx.size(); j++){
				if(cur_idx[j].VCF_BP_INFO_idx_st == -1){
					cur_idx[j].VCF_BP_INFO_idx_st = cur_idx[j - 1].VCF_BP_INFO_idx_ed;
					cur_idx[j].VCF_BP_INFO_idx_ed = cur_idx[j].VCF_BP_INFO_idx_st;
				}
			}
		}

		return 0;
	}

	void search_vcfs_in_region(R_region & r, std::vector<VCF_BP_INFO> &vcf_bp_info_in_region){
		//search the index:
		std::vector<VCF_BP_INFO_IDX> & cur_idx = vcf_bp_info_idx_l[r.chr_ID];
		int the_bg_idx = cur_idx[r.st_pos / vcf_bp_info_idx_window_len].VCF_BP_INFO_idx_st;
		int the_ed_idx = cur_idx[(r.ed_pos /vcf_bp_info_idx_window_len) + 1].VCF_BP_INFO_idx_ed;
		//search all vcf record
		vcf_bp_info_in_region.clear();
		for(int i = the_bg_idx; i < the_ed_idx; i++){
			if(vcf_bp_info_l[i].position >= r.st_pos && vcf_bp_info_l[i].position <= r.ed_pos) {
				vcf_bp_info_in_region.emplace_back();
				vcf_bp_info_in_region.back().chrID = vcf_bp_info_l[i].chrID;
				vcf_bp_info_in_region.back().position = vcf_bp_info_l[i].position;
				vcf_bp_info_in_region.back().hap_n = vcf_bp_info_l[i].hap_n;
				vcf_bp_info_in_region.back().int_type = vcf_bp_info_l[i].int_type;
			}
		}
	}

	void score_M1(double &total_entropy, double &score1, int &total_HET){
		total_entropy = 0;
		total_HET = 0;
		int total_HET_INS = 0;
		int total_HET_DEL = 0;
		//if(r.st_pos == 205208976){fprintf(stderr, " ");}
		for(uint E_idx = 0; E_idx < Entropy_window_l.size(); E_idx++){
			double P = ((float)Entropy_window_l[E_idx].BP_read_id.size()) / ((float)Entropy_window_l[E_idx].total_read_n);
			double Q = 1 - P;
			double cur_entropy = (P * log2(P + 0.000001)) + (Q * log2(Q + 0.000001));
			if(cur_entropy != 0){
				total_entropy += cur_entropy;
			}
			if(Entropy_window_l[E_idx].sv_hap_n_INS == 1){
				total_HET_INS += 1;
			}
			if(Entropy_window_l[E_idx].sv_hap_n_DEL == 1){
				total_HET_DEL += 1;
			}
			if(Entropy_window_l[E_idx].sv_hap_n_DEL == 1 || Entropy_window_l[E_idx].sv_hap_n_INS == 1){
				total_HET +=1;
			}
		}
		score1 = -total_entropy - total_HET;
	}

	void score_KLD(double &total_KLD){
		total_KLD = 0;
		int total_HET_INS = 0;
		int total_HET_DEL = 0;
		//if(r.st_pos == 205208976){fprintf(stderr, " ");}
		for(uint E_idx = 0; E_idx < Entropy_window_l.size(); E_idx++){
			double P = ((float)Entropy_window_l[E_idx].BP_read_id.size()) / ((float)Entropy_window_l[E_idx].total_read_n);
			double Q = 1 - P;
			double Pa = 0;
			int answer_SV_N = MAX(Entropy_window_l[E_idx].sv_hap_n_INS, Entropy_window_l[E_idx].sv_hap_n_DEL);
			if(answer_SV_N == 0){		Pa = 0; }
			else if(answer_SV_N == 1){	Pa = 0.5; }
			else						Pa = 1;
			double Qa = 1 - Pa;
			double KLD = P*log((P+0.01)/(Pa + 0.01)) + Q * log((Q + 0.01)/(Qa + 0.01));
			total_KLD += KLD;
		}
	}

	void analysis_final(R_region & r, std::vector<VCF_BP_INFO> & vcf_bp_info_in_region){
		//calculate theEntropy
		int step_len = 10;
		Entropy_window_l.clear();
		for(int st_pos = r.st_pos; st_pos < r.ed_pos; st_pos += step_len){
			Entropy_window_l.emplace_back();
			int read_number = 0;
			for(uint i = 0; i < RD_l.size(); i++){
				if(RD_l[i].read_st < st_pos && RD_l[i].read_ed > st_pos + step_len)
					read_number++;
			}
			//store read number
			Entropy_window_l.back().total_read_n = read_number;
			Entropy_window_l.back().st_pos = st_pos;
			Entropy_window_l.back().ed_pos = st_pos + step_len;
			Entropy_window_l.back().sv_hap_n_INS = 0;
			Entropy_window_l.back().sv_hap_n_DEL = 0;
		}
		//store the signals to Entropy_window_l
		int base_pos = r.st_pos;
		for(uint i = 0; i < BP_l.size(); i++){
			int w_id = (BP_l[i].position - base_pos) / step_len;
			if(w_id >= 0 && w_id < (int)Entropy_window_l.size())
				Entropy_window_l[w_id].BP_read_id.emplace(BP_l[i].read_id);
		}
		//store the SVs to the  Entropy_window_l
		for(VCF_BP_INFO & v:vcf_bp_info_in_region){
			int w_id = (v.position - base_pos) / step_len;
			if(w_id >= 0 && w_id < (int)Entropy_window_l.size()){
				if(v.int_type == 0)		Entropy_window_l[w_id].sv_hap_n_INS += v.hap_n ;
				else					Entropy_window_l[w_id].sv_hap_n_DEL += v.hap_n ;
			}
		}

		//entropy calculate
		if(false){
			double total_entropy = 0;
			double score1 = 0; int total_HET = 0;
			score_M1(total_entropy, score1,total_HET);
			//output
			fprintf(stdout, "%s:%d-%d\t"
				"%f\t%d\t"
				"%ld\t"
				"%f\n", faidx_iseq(c_ref_idx, r.chr_ID),
				r.st_pos, r.ed_pos,
				total_entropy, total_HET , Entropy_window_l.size(),
				score1);
		}

		//KLD:
		if(true){
			double total_KLD = 0;
			score_KLD(total_KLD);
			fprintf(stdout, "%s:%d-%d\t"
				"%f\t%ld\n",
				faidx_iseq(c_ref_idx, r.chr_ID), r.st_pos, r.ed_pos,
				total_KLD, Entropy_window_l.size());
		}

	}

	int ana_run_Entropy(int argc, char *argv[]){
		char * region_bed_in = argv[1];//separate by ','
		char * fa_fn_in = argv[2];
		char * GS_VCF_FileName = argv[3];
		char * LRS_bamFileName = argv[4];
		int MIN_sv_len = atoi(argv[5]);
		bool print_log = false;
		//load reference index
		FILE* try_open = xopen(fa_fn_in, "r");
		fclose(try_open);
		c_ref_idx = reference_index_load(fa_fn_in);
		//load region list
		load_reion_list(region_bed_in);
		//load bam, for each region:
		//bam_hdr_t * cur_header = NULL;
		if(LRS_bamFileName != NULL)
			LRS_read.init(LRS_bamFileName, NULL, NULL, fa_fn_in);
		//cur_header = LRS_read.file._hdr;
		//load all the vcfs
		load_all_vcfs(GS_VCF_FileName, MIN_sv_len, 10);
		//handler each region:
		Bam_file *c_b = &(LRS_read.file);
		std::vector<char> seq_S;
		std::vector<VCF_BP_INFO> vcf_bp_info_in_region;
		for(uint r_i = 0; r_i < region_l.size(); r_i++){
			//load reference in the region:
			char reg[1024];
			const char * chr_name = faidx_iseq(c_ref_idx, region_l[r_i].chr_ID);
			sprintf(reg, "%s:%d-%d", chr_name, region_l[r_i].st_pos, region_l[r_i].ed_pos);

			//S1: load all reads and analysis:
			resetRegion_ID(c_b, &region_l[r_i]);	//reset region
			collect_bam_signals_Entropy(MIN_sv_len);
			//S2: analysis the vcfs in the region
			search_vcfs_in_region(region_l[r_i], vcf_bp_info_in_region);
			//S3: analysis the final results
			analysis_final(region_l[r_i], vcf_bp_info_in_region);
		}
		//close bam
		if(LRS_bamFileName != NULL)
			LRS_read.destroy();
		return 0;
	}

};

struct SMALL_VAR_DENSITY{
	faidx_t * c_ref_idx;
	std::vector<R_region> region_l;
	BAM_handler LRS_read;
	std::vector<char> seq_S;
	bool print_log;

	void load_reion_list(const char * region_bed_in){
		region_l.clear();
		//load region list[BED file]
		std::vector<std::string> load_line;
		std::vector<std::string> item_value;
		load_string_list_from_file(region_bed_in, load_line);
		size_t line_num = load_line.size();//skip the header line
		//skip the head line
		for(uint64_t i = 1; i < line_num; i++){
			split_string(item_value, load_line[i].c_str(), "\t");//skip the header line
			region_l.emplace_back();
			int chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			if(chrID <= 23){
				region_l.back().chr_ID =chrID;
				region_l.back().st_pos = atoi(item_value[1].c_str());
				region_l.back().ed_pos = atoi(item_value[2].c_str());
			}
		}
	}

	bool load_read_seq(int st_pos, int end_pos, std::vector<char>& seq_S, bam1_t* _bp){
		uint8_t *bam_seq = bam_get_seq(_bp);
		end_pos = MIN(end_pos, _bp->core.l_qseq);
		seq_S.resize(end_pos - st_pos);
		char *seq = &(seq_S[0]);
		int n_count = 0;
		while(st_pos < end_pos)
		{
			uint8_t c_c = bam_seqi(bam_seq, st_pos);
			switch(c_c)
			{
			case 1: *seq++ ='A'; break;
			case 2: *seq++ ='C'; break;
			case 4: *seq++ ='G'; break;
			case 8: *seq++ ='T'; break;
			default: *seq++ =0; n_count ++; break;
			}
			st_pos++;
		}
		if(n_count > 10)
			return false;
		return true;
	}

	static void cigar2Path(const uint32_t* bam_cigar, const unsigned n_cigar, path_segment* path){
		for (unsigned int i = 0; i < n_cigar; ++i){
			path[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			path[i].type   = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		}
	}

	void collect_bam_signals_SMALL_VAR_DEN(
		int chrID, int region_begin, int region_end, //input
		int &N_read, int &N_SNP, int &N_match_base, int &N_INS, int &N_DEL, int &len_INS, int &len_DEL
	){
		//init
		N_read = 0;
		N_SNP = 0;
		N_match_base = 0;
		N_INS = 0;
		N_DEL = 0;
		len_INS = 0;
		len_DEL = 0;

		//set the region:
		R_region cur_R;
		cur_R.chr_ID = chrID;  cur_R.st_pos = region_begin; cur_R.ed_pos = region_end;
		//set the bam region
		Bam_file *c_b = &(LRS_read.file);
		resetRegion_ID(c_b, &cur_R);	//reset region
		//loading the reference sequence:
		int load_ref_len = 0;
		const char * chrName = faidx_iseq(c_ref_idx, chrID);
		char reg[1024];
		sprintf(reg, "%s:%d-%d", chrName, cur_R.st_pos, cur_R.ed_pos);
		char *ref = fai_fetch(c_ref_idx, reg, &load_ref_len);
		//handler each region:
		seq_S.clear();
		//init the BP & RD list
		//load all reads and analysis:
		int read_id = 0;
		while (bam_next(c_b)) {
			//new sig for this read
			bam1_t *br = &(c_b->_brec);
			//if (print_log)fprintf(stderr, "Read Begin NAME %s \n", bam_get_qname(br));
			//filter:
			//basic filter
			if (bam_is_secondary(br))
				continue;
			if (br->core.qual < 4)
				continue;
			//supplementary-min length check
			if (bam_is_supplementary(br) && br->core.l_qseq < 1200)
				continue;
			read_id++;
			N_read++;
			//get seq,  get qual:
			load_read_seq(0, br->core.l_qseq, seq_S, br);
			//uint8_t * qual = bam_get_qual(br);
			//get cigar

			int region_begin_local = region_begin - br->core.pos;
			int region_end_local   = region_end - br->core.pos;

			std::vector<path_segment> cigar_P;
			cigar_P.resize(br->core.n_cigar);
			cigar2Path( bam_get_cigar(br), br->core.n_cigar, &(cigar_P[0]));
			int seq_i = 0;
			int ref_i = 0;
			//get the reference sequence:
			char *tseq = ref;//??????
			char *qseq = &(seq_S[0]);
			for (uint i = 0; i < br->core.n_cigar; i++) {
				path_segment *p = &(cigar_P[i]);
				if (p->type == align_t::CIGAR_SOFT_CLIP)			//soft clip:
				{
					seq_i += p->length;
				} else if (p->type == align_t::CIGAR_INSERT) {
					if(ref_i > region_begin_local && ref_i < region_end_local){//within
						N_INS ++;
						len_INS += p->length;
					}
					seq_i += p->length;
				} else if (p->type == align_t::CIGAR_DELETE) {
					int DEL_END = ref_i + p->length;
					if(ref_i < region_end_local && region_begin_local < DEL_END){//overlap
						int overlap_length = MIN(region_end_local - ref_i, DEL_END - region_begin_local);
						xassert(overlap_length >= 0, "");
						N_DEL ++;
						len_DEL += overlap_length;
					}
					ref_i += p->length;
				} else if (p->type == align_t::CIGAR_MATCH
						|| p->type == align_t::CIGAR_SEQ_MATCH
						|| p->type == align_t::CIGAR_SEQ_MISMATCH) {
					for(uint i = 0; i < p->length; i++, ref_i++, seq_i++){
						if(ref_i > region_begin_local && ref_i < region_end_local){//within
							int ref_idx_in_ref = br->core.pos - cur_R.st_pos + ref_i + 1;
							xassert(ref_idx_in_ref >= 0 && ref_idx_in_ref < load_ref_len, "");
							if(tseq[ref_idx_in_ref] == qseq[seq_i]) //when match
								N_match_base ++;
							else
								N_SNP ++;
						}
					}
				} else if (p->type == align_t::CIGAR_HARD_CLIP) {
					//DO nothing
				} else {
					//DO nothing
				}
			}
		}
	}

	int ana_run_SMALL_VAR_DENSITY(int argc, char *argv[]){
		char * region_bed_in = argv[1];//separate by ','
		char * fa_fn_in = argv[2];
		char * LRS_bamFileName = argv[3];
		//load reference index
		FILE* try_open = xopen(fa_fn_in, "r");
		fclose(try_open);
		c_ref_idx = reference_index_load(fa_fn_in);
		//load region list
		load_reion_list(region_bed_in);
		//load bam, for each region:
		//bam_hdr_t * cur_header = NULL;
		if(LRS_bamFileName != NULL)
			LRS_read.init(LRS_bamFileName, NULL, NULL, fa_fn_in);
		//cur_header = LRS_read.file._hdr;
		//load all the vcfs
		//load_all_vcfs(GS_VCF_FileName, MIN_sv_len, 10);
		//handler each region:
		std::vector<char> seq_S;
		for(uint r_i = 0; r_i < region_l.size(); r_i++){
			R_region &r = region_l[r_i];
			//analysis P1:
			{
				int N_read = 0; int N_SNP = 0;	int N_match_base = 0; int N_INS = 0; int N_DEL = 0; 	int len_INS = 0; int len_DEL = 0;
				collect_bam_signals_SMALL_VAR_DEN(r.chr_ID, r.st_pos, r.ed_pos, N_read, N_SNP, N_match_base, N_INS, N_DEL, len_INS, len_DEL);
				fprintf(stdout, "R_FULL [%s:%d-%d]\t", faidx_iseq(c_ref_idx, r.chr_ID), r.st_pos, r.ed_pos);//regions
				int region_len = r.ed_pos - r.st_pos;
				float gap_excluded_identity = (float)N_match_base/(N_match_base + N_SNP + 0.001);
				float gap_compressed_identity = (float)N_match_base/(N_match_base + N_SNP + N_INS + N_DEL + 0.001);
				float AVE_MIS_BASE = ((float)N_SNP + N_INS + N_DEL)/N_read;
				float AVE_M_BASE = ((float)N_match_base)/N_read;;
				fprintf(stdout, "I1=%f;T2=%f\t", gap_excluded_identity, gap_compressed_identity);//identity
				fprintf(stdout, "AVE_MIS_BASE=%f;AVE_M_BASE=%f\t", AVE_MIS_BASE, AVE_M_BASE);//
				fprintf(stdout, "R_LEN=%d;N_read=%d;N_SNP=%d;N_match_base=%d;N_INS=%d;N_DEL=%d;len_INS=%d;len_DEL=%d\n",
						region_len, N_read, N_SNP, N_match_base, N_INS, N_DEL, len_INS,len_DEL);
			}
			{
				int N_read = 0; int N_SNP = 0;	int N_match_base = 0; int N_INS = 0; int N_DEL = 0; 	int len_INS = 0; int len_DEL = 0;
				collect_bam_signals_SMALL_VAR_DEN(r.chr_ID, r.st_pos + 500, r.ed_pos - 500,N_read, N_SNP, N_match_base, N_INS, N_DEL, len_INS, len_DEL);
				fprintf(stdout, "R_PART [%s:%d-%d]\t", faidx_iseq(c_ref_idx, r.chr_ID), r.st_pos, r.ed_pos);//regions
				int region_len = r.ed_pos - r.st_pos -1000;
				float gap_excluded_identity = (float)N_match_base/(N_match_base + N_SNP + 0.001);
				float gap_compressed_identity = (float)N_match_base/(N_match_base + N_SNP + N_INS + N_DEL + 0.001);
				float AVE_MIS_BASE = ((float)N_SNP + N_INS + N_DEL)/N_read;
				float AVE_M_BASE = ((float)N_match_base)/N_read;;
				fprintf(stdout, "I1=%f;T2=%f\t", gap_excluded_identity, gap_compressed_identity);//identity
				fprintf(stdout, "AVE_MIS_BASE=%f;AVE_M_BASE=%f\t", AVE_MIS_BASE, AVE_M_BASE);//
				fprintf(stdout, "R_LEN=%d;N_read=%d;N_SNP=%d;N_match_base=%d;N_INS=%d;N_DEL=%d;len_INS=%d;len_DEL=%d\n",
				region_len, N_read, N_SNP, N_match_base, N_INS, N_DEL, len_INS,len_DEL);
			}
		}
		//close bam
		if(LRS_bamFileName != NULL)
			LRS_read.destroy();
		return 0;
	}
};


struct COMBINE_VCF{

	struct VCF_COM_Record{
		std::string sample;
		bcf1_t 		r;
		VCF_COM_Record(){
			memset(&r, 0, sizeof(bcf1_t));
		}
		static inline int cmp_by_pos(const VCF_COM_Record &a, const VCF_COM_Record &b){
			//var basic
			if(a.r.rid != b.r.rid) 	return a.r.rid < b.r.rid;
			if(a.r.pos != b.r.pos) 	return a.r.pos < b.r.pos;
			return a.sample.compare(b.sample);
		}
	};

	int combine_sort_vcf(int argc, char *argv[]){
		char * vcf_fn_in = argv[1];//separate by ','
		char * vcf_fn_out = argv[2];

		//get bam file list
		std::vector<std::string> vcf_files_names;
		split_string(vcf_files_names, vcf_fn_in, ",");
		std::vector<BCF_FILE> bam_list;
		for(std::string &bam_fn:  vcf_files_names){
			bam_list.emplace_back();
			memset(&bam_list.back(), 0, sizeof(BCF_FILE));
			VCF_open_read(&bam_list.back(),bam_fn.c_str());
			//bcf_hdr_append(bam_list.back().header, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
		}

		std::vector<VCF_COM_Record> vcf_list;
		uint32_t load_size = 0;

		for(BCF_FILE & bcf_f: bam_list){
			std::string sample_name(bcf_f.header->samples[0]);
			do{
				if(vcf_list.size() < load_size + 1)
					vcf_list.emplace_back();
				 vcf_list[load_size].sample = sample_name;
			}while(VCF_next_dump(&bcf_f, &(vcf_list[load_size].r)) && vcf_list[load_size].r.rid == 0 && ++load_size);
		}

		std::sort(vcf_list.begin(), vcf_list.end(), VCF_COM_Record::cmp_by_pos);

		//write::
		//open write file
		BCF_FILE vcf_out;
		VCF_open_write(&vcf_out, vcf_fn_out, false);

		bcf_hdr_t *write_header = bam_list[0].header;
		//bcf_hdr_append(write_header, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
		//fprintf(stderr, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
	    vcf_hdr_write(vcf_out.file, write_header);

	    //int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);

		//append sample info for each line:
		for(VCF_COM_Record & vcf_r: vcf_list){
		    //int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);
			if ( !(vcf_r.r.unpacked & BCF_UN_ALL) ) bcf_unpack(&vcf_r.r, BCF_UN_ALL);
			bcf_add_id(write_header, &vcf_r.r, vcf_r.sample.c_str());
			//bcf_update_info(write_header, &vcf_r.r, "SAMPLE",vcf_r.sample.c_str() , vcf_r.sample.size(), BCF_HT_STR);
		    vcf_write(vcf_out.file, write_header, &vcf_r.r);
		}
		//close::
		for(BCF_FILE & bcf_f: bam_list){
			bcf_close(bcf_f.file);
		}
		return 0;
	}

	typedef struct
	{
		char ID[64];
		char type[32];
		int rid;
		int sample;
		int st;
		int ed;
	}VCF_ANALYSIS;
	kvec_T(VCF_ANALYSIS, VCF_A_L);
	//compare two vcf record, return overlap rate
	float vcf_compare_one(VCF_ANALYSIS *v1, VCF_ANALYSIS *v2)
	{
		int MIN_ED = MIN(v1->ed, v2->ed);
		int MAX_ST = MAX(v1->st, v2->st);
		int overlap = MIN_ED - MAX_ST;
		if(overlap < 0)
			return 0;
		int length1 = v1->ed - v1->st;
		int length2 = v2->ed - v2->st;
		int MAX_LENGTH = MAX(length1, length2);
		return ((float)overlap)/MAX_LENGTH;
	}
	void print_one_vcf_analysis(FILE *stream, VCF_ANALYSIS *c_va) {
		fprintf(stream, "%d\t"
				"%s\t"
				"%d\t"
				"%d\t"
				"%d\t"
				"%s\n", c_va->sample, c_va->type, c_va->rid, c_va->st, c_va->ed,
				c_va->ID);
	}
	//to compare data from one sample, one a sv type and one chromosome
	void vcf_compare_sample_type_chromosome(VCF_ANALYSIS *vl1, int n1, VCF_ANALYSIS *vl2, int n2,
			FILE* flog, int * overlap_90, int *overlap_50)
	{
		int overlap_90_percent = 0;
		int overlap_50_percent = 0;

		//for v1
		for(int i = 0; i < n1; i++)
		{
			VCF_ANALYSIS *v1 = vl1 + i;
			float max_overlap = 0;
			VCF_ANALYSIS c_MAX;
			for(int j = 0; j < n2; j++)
			{
				float overlap = vcf_compare_one(v1, vl2 + j);
				if(overlap > max_overlap)
				{
					max_overlap = overlap;
					memcpy(&c_MAX, vl2 + j, sizeof(VCF_ANALYSIS));
				}
			}
			print_one_vcf_analysis(flog, v1);
			if(max_overlap > 0.01)
			{
				fprintf(flog, "\t\t%f\t", max_overlap);
				print_one_vcf_analysis(flog, &c_MAX);
				if(max_overlap > 0.9)
					overlap_90_percent++;
				if(max_overlap > 0.5)
					overlap_50_percent++;
			}
		}
		*overlap_90 = overlap_90_percent;
		*overlap_50 = overlap_50_percent;
	}

	int get_one_vcf_manta(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
	{
		vcf_get_sv_type(header, vcf, vcf_a->type);
		vcf_get_sv_END(header, vcf, &(vcf_a->ed));
		strcpy(vcf_a->ID, vcf->d.id);
		vcf_a->rid = vcf->rid;
		vcf_a->st = vcf->pos;
		return true;
	}
	static int VCF_ANALYSIS_cmp_sample_type_pos(const void*a_,const void*b_)
	{
		VCF_ANALYSIS *a = (VCF_ANALYSIS *)a_;
		VCF_ANALYSIS *b = (VCF_ANALYSIS *)b_;

		if(a->sample != b->sample)
			return a->sample > b->sample;
		int type_cmp = strcmp(a->type, b->type);
		if(type_cmp != 0)
			return (type_cmp > 0);
		if(a->rid != b->rid)
			return a->rid > b->rid;
		return a->st > b->st;

	}

	void vcf_manta_load(char *fn_in, VCF_A_L *vcf_l)
	{
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read
		bcf_hdr_t *header = vcf_r.header;
		char *GT[3]; char GT_1[3];char GT_2[3];char GT_3[3]; GT[0] = GT_1; GT[1] = GT_2; GT[2] = GT_3;
		int nsmpl = bcf_hdr_nsamples(header);
		while(VCF_next(&vcf_r))//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			VCF_ANALYSIS c_va;
			get_one_vcf_manta(c_r, header, &c_va);

			int N_GT = vcf_get_sv_GT(header, c_r, GT)/nsmpl;
			for(int i = 0; i < N_GT; i++)
			{
				if(GT[i][0] == 4 || GT[i][1] == 4)
				{
					c_va.sample = i;
					kv_push_2(VCF_ANALYSIS, vcf_l, c_va);
				}
			}
		}
		qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
		//close
		bcf_close(vcf_r.file);
	}

	int get_one_vcf_jlra(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
	{
		char sample[128];
		vcf_get_sample(header, vcf, sample);
		int sample_int = 100;
		if(strcmp(sample, "0") == 0)
			sample_int = 0;
		else if(strcmp(sample, "1") == 0)
			sample_int = 1;
		else if(strcmp(sample, "2") == 0)
			sample_int = 2;
		if(sample_int > 2)
			return false;
		vcf_a->sample = sample_int;
		vcf_get_sv_type(header, vcf, vcf_a->type);
		vcf_get_sv_END(header, vcf, &(vcf_a->ed));
		strcpy(vcf_a->ID, vcf->d.id);
		vcf_a->rid = vcf->rid;
		vcf_a->st = vcf->pos;
		return true;
	}
	void vcf_jlra_load(char *fn_in, VCF_A_L *vcf_l)
	{
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read
		bcf_hdr_t *header = vcf_r.header;
		while(VCF_next(&vcf_r))//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			VCF_ANALYSIS *c_va = NULL;
			kv_pushp_2(VCF_ANALYSIS, vcf_l, c_va);
			if(! get_one_vcf_jlra(c_r, header, c_va))
				vcf_l->n--;
		}
		qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
		//close
		bcf_close(vcf_r.file);
	}

	void vcf_jlra_dump(char *fn_in, char *fn_out)
	{
		VCF_A_L vcf_l = {0};
		FILE * out = xopen(fn_out, "w");
		vcf_jlra_load(fn_in, &vcf_l);
		for(unsigned int i = 0; i < vcf_l.n; i++)
			print_one_vcf_analysis(out, &(vcf_l.a[i]));
		fclose(out);
	}
	int get_one_vcf_nstd(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
	{
		char sample[128];
		vcf_get_sample(header, vcf, sample);
		int sample_int = 100;
		if(strcmp(sample, "HG00512") == 0)
			sample_int = 0;
		else if(strcmp(sample, "HG00513") == 0)
			sample_int = 1;
		else if(strcmp(sample, "HG00514") == 0)
			sample_int = 2;
		if(sample_int > 2)
			return false;
		vcf_a->sample = sample_int;
		vcf_get_sv_type(header, vcf, vcf_a->type);
		vcf_get_sv_END(header, vcf, &(vcf_a->ed));
		strcpy(vcf_a->ID, vcf->d.id);
		vcf_a->rid = vcf->rid;
		vcf_a->st = vcf->pos;
		return true;
	}

	void vcf_nstd_load(char *fn_in, VCF_A_L *vcf_l)
	{
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read
		bcf_hdr_t *header = vcf_r.header;
		while(VCF_next(&vcf_r))//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			VCF_ANALYSIS *c_va = NULL;
			kv_pushp_2(VCF_ANALYSIS, vcf_l, c_va);
			if(! get_one_vcf_nstd(c_r, header, c_va))
				vcf_l->n--;
		}
		qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
		bcf_close(vcf_r.file);	//close
	}

	void vcf_nstd_dump(char *fn_in, char *fn_out)
	{
		VCF_A_L vcf_l = {0};
		FILE * out = xopen(fn_out, "w");
		vcf_nstd_load(fn_in, &vcf_l);
		for(unsigned int i = 0; i < vcf_l.n; i++)
			print_one_vcf_analysis(out, &(vcf_l.a[i]));
		fclose(out);
	}

	void vcf_compare(char *fn_1, char * file_type1, char * fn_2, char *file_type2, char * fn_log)
	{
		//open log file
		FILE *flog = xopen(fn_log, "w");
		//open files
		VCF_A_L vcf_l1 = {0};
		VCF_A_L vcf_l2 = {0};
		for(int loop = 0; loop < 2; loop++)
		{
			char *fn = (loop == 0)?fn_1:fn_2;
			char *ft = (loop == 0)?file_type1:file_type2;
			VCF_A_L *vcf_l = (loop == 0)?(&vcf_l1):(&vcf_l2);
			if(strcmp(ft, "jlra") == 0)
				vcf_jlra_load(fn, vcf_l);
			if(strcmp(ft, "manta") == 0)
				vcf_manta_load(fn, vcf_l);
			if(strcmp(ft, "nstd") == 0)
				vcf_nstd_load(fn, vcf_l);
		}
		//for one sample
		int sample_ID = 0;
		unsigned int sample_st = 0, sample_ed = 0;
		unsigned int sample_st2 = 0, sample_ed2 = 0;
		while(1)
		{
			//get st/ed and sample ID for method1
			if(sample_st >= vcf_l1.n)//end of data
				break;
			sample_ID = vcf_l1.a[sample_st].sample;
			for(sample_ed = sample_st; sample_ed < vcf_l1.n && vcf_l1.a[sample_ed].sample == sample_ID; sample_ed++);
			//get st/ed for method2 for sample ID
			int old_st = sample_st2, old_ed = sample_ed2;
			for(; sample_st2 < vcf_l2.n && vcf_l2.a[sample_st2].sample != sample_ID; sample_st2++);
			for(sample_ed2 = sample_st2; sample_ed2 < vcf_l2.n && vcf_l2.a[sample_ed2].sample == sample_ID; sample_ed2++);
			if(sample_st2 == sample_ed2) {sample_st2 = old_st;sample_ed2 = old_ed;}//reset when no results

			//for one type
			char type_ID[32];
			int type_st = sample_st, type_ed = type_st;
			int type_st2 = sample_st2, type_ed2 = type_st2;
			while(1)
			{
				//get st/ed and type ID for method1
				if(type_st >= sample_ed)//end of data
					break;
				strcpy(type_ID, vcf_l1.a[type_st].type);
				for(type_ed = type_st; type_ed < sample_ed && (strcmp(vcf_l1.a[type_ed].type, type_ID) == 0); type_ed++);
				//get st/ed for method2 for type ID
				int old_st = type_st2, old_ed = type_ed2;
				for(; type_st2 < sample_ed2 && (strcmp(vcf_l2.a[type_st2].type, type_ID) != 0); type_st2++);
				for(type_ed2 = type_st2; type_ed2 < sample_ed2 &&  (strcmp(vcf_l2.a[type_ed2].type, type_ID) == 0); type_ed2++);
				if(type_st2 == type_ed2) {type_st2 = old_st;type_ed2 = old_ed;}//reset when no results
				//for one chromosome
				//for one type
				int chr_ID;
				int chr_st = type_st, chr_ed = chr_st;
				int chr_st2 = type_st2, chr_ed2 = chr_st2;
				while(1)
				{
					//get st/ed and type ID for method1
					if(chr_st >= type_ed)//end of data
						break;
					chr_ID = vcf_l1.a[chr_st].rid;
					for(chr_ed = chr_st; chr_ed < type_ed && vcf_l1.a[chr_ed].rid == chr_ID; chr_ed++);
					//get st/ed for method2 for type ID
					int old_st = chr_st2, old_ed = chr_ed2;
					for(; chr_st2 < type_ed2 && vcf_l2.a[chr_st2].rid != chr_ID; chr_st2++);
					for(chr_ed2 = chr_st2; chr_ed2 < type_ed2 && vcf_l2.a[chr_ed2].rid == chr_ID; chr_ed2++);
					if(chr_st2 == chr_ed2) {chr_st2 = old_st;chr_ed2 = old_ed;}//reset when no results
					//compare
					int overlap_90_percent, overlap_50_percent;
					int number_method1 = chr_ed - chr_st;
					int number_method2 = chr_ed2 - chr_st2;
					vcf_compare_sample_type_chromosome(vcf_l1.a + chr_st, number_method1,
							vcf_l2.a + chr_st2, number_method2, flog, &overlap_90_percent, &overlap_50_percent);
					//print results
					fprintf(stderr,
							"sample_ID: %d\t"
							"type_ID:%s\t"
							"chr_ID:%d\t"
							"number_method1:%d\t"
							"number_method2:%d\t"
							"overlap_90_percent:%d\t"
							"overlap_50_percent:%d\t"
							"SEN_90:%f\t"
							"SEN_90:%f\t"
							"\n",
							sample_ID,
							type_ID,
							chr_ID,
							number_method1,
							number_method2,
							overlap_90_percent,
							overlap_50_percent,
							((float)overlap_90_percent)/number_method1,
							((float)overlap_90_percent)/number_method1);
					//end
					chr_st = chr_ed;
					chr_st2 = chr_ed2;
				}
				//end
				type_st = type_ed;
				type_st2 = type_ed2;
			}
			//end
			sample_st = sample_ed;
			sample_st2 = sample_ed2;
		}
		fclose(flog);
	}
};

struct BED_MERGING_handler{
	std::vector<std::vector<int>> count;

	int level_1_persent;
	int level_5_persent;
	int level_50_persent;

	int get_level(int count, int sample_number){
		int level = 0;
		if(count == 0)						level = 0;//zero
		if(count*20 <= sample_number)		level = 1;//less then 5%;
		else 								level = 2;//other
		return level;
	}

	int bed_merging(int argc, char *argv[]){
		char * fa_fn_in = argv[1];//separate by ','
		char * bed_fn_to_be_merged_fn = argv[2];//separate by ','
		//
		int STEP_LEN = 20;

		//load the reference index
		FILE* try_open = xopen(fa_fn_in, "r");
		fclose(try_open);
		faidx_t * c_ref_idx = reference_index_load(fa_fn_in);
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(26, N_seq);
		//init the list
		for(int chr_id = 0; chr_id < N_seq; chr_id++){
			int chrLength = faidx_seq_len(c_ref_idx, faidx_iseq(c_ref_idx, chr_id));
			int window_size = (chrLength/STEP_LEN) + 1;
			count.emplace_back();
			count.back().resize(window_size);
			std::vector<int> &chr_count = count.back();
			for(uint j = 0; j < chr_count.size() ; j++)
				chr_count[j] = 0;
		}

		//load the BED files to merge
		std::vector<std::string> bed_fn_to_be_merged_v;
		load_string_list_from_file(bed_fn_to_be_merged_fn, bed_fn_to_be_merged_v);
		size_t sample_num = bed_fn_to_be_merged_v.size();//skip the header line
		std::vector<std::string> bed_in_sample;
		std::vector<std::string> item_value;
		std::string old_buff;
		int old_chr_ID = -1;
		std::map <int, int> SVLEN_counter;
		for(uint sampleID = 0; sampleID < sample_num; sampleID++){
			//load data:
			fprintf(stderr, "Current handle sample is %s\n", bed_fn_to_be_merged_v[sampleID].c_str());
			bed_in_sample.clear();
			load_string_list_from_file(bed_fn_to_be_merged_v[sampleID].c_str(), bed_in_sample);
			//for all record
			for(uint bed_record = 0; bed_record < bed_in_sample.size(); bed_record++){
				split_string(item_value, bed_in_sample[bed_record].c_str(), "\t");
				int record_chr_ID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), old_buff.c_str(), old_chr_ID);
				old_buff = item_value[0]; old_chr_ID = record_chr_ID;
				int record_pos = (int) atoi(item_value[1].c_str()); int record_pos_wid = record_pos/STEP_LEN;
				int record_end = (int) atoi(item_value[2].c_str()); int record_end_wid = record_end/STEP_LEN + 1;
				if(record_end < record_pos + 100)
					continue;
				for(int k = record_pos_wid; k < record_end_wid; k++)
					count[record_chr_ID][k] ++;
			}
		}

		//merging and outputting
		for(int chr_id = 0; chr_id < N_seq; chr_id++){
			const char * chrName = faidx_iseq(c_ref_idx, chr_id);
			std::vector<int> &chr_count = count[chr_id];
			int old_st = -1; int old_ed = -1; int old_level = -1;
			int total_window = 0; int total_count = 0;
			for(uint j = 0; j < chr_count.size() + 1; j++){
				if(j == chr_count.size()){
					if(old_st != -1 && old_ed >= old_st + 500)
						fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%f\n", chrName, old_st, old_ed, old_level, total_window, total_count, (float)total_count/total_window);
				}
				else if(chr_count[j] > 0){
					int count_level = get_level(chr_count[j], sample_num);
					int st_pos = j*STEP_LEN;
					int ed_pos = (j + 1)*STEP_LEN;
					//fprintf(stdout, "TEST::%s\t%d\t%d\t%d\n", chrName, st_pos, ed_pos, chr_count[j]);
					if(st_pos == old_ed && count_level == old_level){
						/*DO NOTHING*/
					}else{
						if(old_st != -1 && old_ed >= old_st + 500)
							fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%f\n", chrName, old_st, old_ed, old_level, total_window, total_count, (float)total_count/total_window);
						total_window = 0;
						total_count = 0;
						old_st = st_pos;
					}
					total_window ++;
					total_count += chr_count[j];
					old_ed = ed_pos;
					old_level = count_level;
				}
			}

		}
		return 0;
	}



};
struct FA_stat_handler{
	int N_count_ANA(const std::string &word){
		int N = 0;
		const char *s = word.c_str();
		int n = word.size();
		for(int i = 0; i < n; i++){
			if(s[i] == 'N')
				N++;
		}
		return N;
	}

	int fa_stat_full(int argc, char *argv[]){
		char * fa_fn_in = argv[1];//separate by ','
		int chrID = atoi(argv[2]);

		int REGION_LEN = atoi(argv[3]);
		int REGION_STEP_LEN = atoi(argv[4]);

		//FILE * bed_file = xopen(bed_fn_out, "w");
		//default parameters
		//int REGION_LEN = 500;
		//int REGION_STEP_LEN = 250;
		int MIN_K_size = 5;//4^5=
		int MAX_K_size = 200;
		int K_step_LEN = 1;

		//test file exist
		FILE* try_open = xopen(fa_fn_in, "r");
		fclose(try_open);
		faidx_t * c_ref_idx = reference_index_load(fa_fn_in);
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(26, N_seq);

		//header line
		fprintf(stdout, "%d\t%d\t%d\n", REGION_LEN, REGION_STEP_LEN, 0);
		std::set<std::string> kmerCountingSet;
		for(int i = chrID; i < chrID+1; i++){
			int len = 0;
			const char * chrName = faidx_iseq(c_ref_idx, i);
			char * c_reference = fai_fetch(c_ref_idx, chrName, &len);
			int beginPos = 0; int endPos = 0; int LENGTH = REGION_LEN;
			int chrLength = faidx_seq_len(c_ref_idx, chrName);
			std::string regionString;
			for(;beginPos < chrLength; beginPos += REGION_STEP_LEN){
				endPos = beginPos + LENGTH;
				endPos = MIN(endPos, chrLength);
				regionString.clear();
				for(int j = beginPos; j < endPos; j++){
					regionString += c_reference[j];
				}
				bool withloop = false;
				int kmerSizeF = 0;
				for(int kmerSize = MIN_K_size; kmerSize <= MAX_K_size; kmerSize += K_step_LEN){
					kmerCountingSet.clear();
					const std::string &seq = regionString;
					const unsigned seqLen = seq.size();
					// track all words from the read, including repetitive k
					unsigned kmer_number = seqLen - kmerSize + 1;
					for (unsigned j = 0; j < kmer_number; ++j) {
						const std::string word(seq.substr(j, kmerSize));
						if (N_count_ANA(word) > 0 )
							continue;
						std::set<std::string>::iterator it = kmerCountingSet.find(word);
						if(it == kmerCountingSet.end()){
							kmerCountingSet.emplace(word);
						}else{
							withloop = true;
							kmerSizeF = kmerSize;
							break;
						}
					}
					if(!withloop){//am.DBGTest(kmerSize)){
						break;
					}
				}
				//debug show results
				//if(REPEAT_DEBUG) {for( DBG_INFO_item_ * kmerIdx_i : id2info) kmerIdx_i->showData();}
				if(withloop){
					fprintf(stdout, "%d\t%d\t%d\n",i, beginPos, kmerSizeF);
				}
			}
			free(c_reference);
		}
		//fclose(bed_file);
		return 0;
	}

	int string_local_repeat(int argc, char *argv[]){
		char * string = argv[1];//separate by ','

		int MIN_K_size = 5;//4^5=
		int MAX_K_size = 200;
		int K_step_LEN = 1;
		std::set<std::string> kmerCountingSet;
		int kmerSizeF = 0;
		for(int kmerSize = MIN_K_size; kmerSize <= MAX_K_size; kmerSize += K_step_LEN){
			kmerCountingSet.clear();
			const std::string seq(string);
			const unsigned seqLen = seq.size();
			// track all words from the read, including repetitive k
			int kmer_number = seqLen - kmerSize + 1;
			int j = 0;
			for (j = 0; j < kmer_number; ++j) {
				const std::string word(seq.substr(j, kmerSize));
				if (N_count_ANA(word) > 0 )
					continue;
				std::set<std::string>::iterator it = kmerCountingSet.find(word);
				if(it == kmerCountingSet.end()){
					kmerCountingSet.emplace(word);
				}else{
					kmerSizeF = kmerSize;
					break;
				}
			}
			if(j == kmer_number)
				break;
		}
		//debug show results
		//if(REPEAT_DEBUG) {for( DBG_INFO_item_ * kmerIdx_i : id2info) kmerIdx_i->showData();}
		fprintf(stdout, "%d\n", kmerSizeF);
		return 0;
	}
};

struct Ref_region_handler{
	int get_ref_region(int argc, char *argv[]){
	  //load reference string;
	  faidx_t *fai = fai_load(argv[optind]);//load reference
	  int true_region_load_len = 0;
	  char *ref_seq = fai_fetch(fai, argv[optind + 1], &true_region_load_len);
	  fprintf(stdout, "%s\n", ref_seq);
	  fai_destroy(fai);
	  return 0;
	}

	static int ref_split(int argc, char *argv[])
	{
		char *fasta_file = argv[1];
		gzFile fp = xzopen(fasta_file, "r");
		kstream_t *_fp = ks_init(fp);

		kseq_t temp = {0};
		temp.f = _fp;
		while( kseq_read(&temp) >= 0)
		{
			char out_file_name[1024];
			strcpy(out_file_name, temp.name.s);
			strcat(out_file_name, ".fa");
			FILE *out_file = xopen(out_file_name, "w");
			printf(">%s %s\n", temp.name.s, temp.comment.s);
			fprintf(out_file, ">%s %s\n", temp.name.s, temp.comment.s);
			//print content
			for(unsigned int i = 0; i < temp.seq.l; i++)
			{
				fprintf(out_file, "%c", temp.seq.s[i]);
				if(i %70 == 69)
					fprintf(out_file, "\n");
			}
			fclose(out_file);
		}
		return 0;
	}

	void ref_dump(char *fasta_file, int maxChr)
	{
		gzFile fp = xzopen(fasta_file, "r");
		kstream_t *_fp = ks_init(fp);
		kseq_t temp = {0};
		temp.f = _fp;
		int chr = 0;
		while( kseq_read(&temp) >= 0 && chr++ < maxChr)
		{
			char out_file_name[1024];
			strcpy(out_file_name, temp.name.s);
			strcat(out_file_name, ".fna");
			FILE *out_file = xopen(out_file_name, "w");
			printf(">%s %s\n", temp.name.s, temp.comment.s);
			fprintf(out_file, ">%s %s\n", temp.name.s, temp.comment.s);
			//print content
			for(unsigned int i = 0; i < temp.seq.l; i++)
			{
				fprintf(out_file, "%c", temp.seq.s[i]);
				if(i %70 == 69)
					fprintf(out_file, "\n");
			}
			fclose(out_file);
		}
	}

	#define MAX_LINE_LENGTH 1000
	void getSVRegion(char * region_fn, std::vector<R_region>&regionList, std::vector<std::string> &id2Name){
		//get sv region:
		char *temp = new char[MAX_LINE_LENGTH];//1M
		std::map<std::string, uint32_t> name2ID;

		std::ifstream regionListFile(region_fn);
		while(true){
			regionListFile.getline(temp, MAX_LINE_LENGTH);
			if(*temp == 0)
				break;
			//get chrID
			char *token = strtok(temp, "\t");
			std::string strID(token); int intID = 0;
			std::map<std::string, uint32_t>::iterator it = name2ID.find(strID);
			if(it!=name2ID.end()){
				intID = it->second;
			}
			else{
				intID = name2ID.size();
				name2ID[strID] = name2ID.size();
				id2Name.emplace_back(strID);
			}

			//get refSt
			token = strtok(NULL, "\t");
			uint32_t regionSt = strtoul(token, NULL, 10);
			//get refEd
			token = strtok(NULL, "\t");
			uint32_t regionEnd = strtoul(token, NULL, 10);
			//
			regionList.emplace_back();
			R_region & r = regionList.back();
			r.chr_ID = intID; r.st_pos = regionSt; r.ed_pos = regionEnd;
		}
		regionListFile.close();
		delete [] temp;
	}

	int dump_ref_by_region(int argc, char *argv[])
	{
		char *fasta_file = argv[1];
		char * region_fn = argv[2];
		//load regions
		std::vector<R_region> regionList;
		std::vector<std::string> id2Name;
		getSVRegion(region_fn, regionList, id2Name);

		//open fasta file
	//	gzFile fp = xzopen(fasta_file, "r");
	//	kstream_t *_fp = ks_init(fp);

		//open fai file
		faidx_t * fai = fai_load(fasta_file);

		//load regions
		for(R_region &r :regionList){
			char reg[1024]; int load_len = 0;
			sprintf(reg, "%s:%d-%d", id2Name[r.chr_ID].c_str(), r.st_pos, r.ed_pos);
			//char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);
			char *ref = fai_fetch(fai, reg, &load_len);

			printf(">%s_%d_%d len_%d_%d_%d\n",  id2Name[r.chr_ID].c_str(), r.st_pos, r.ed_pos, load_len, load_len/70, load_len % 70);
			for(int i = 0; i < load_len; i++){
				printf("%c", ref[i]);
				if(i %70 == 69)
					printf("\n");
			}
			if((load_len - 1) % 70 != 69)
				printf("\n");
			free(ref);
		}

		fai_destroy(fai);
		return 0;
	}

	int getReverseStr_simple(int argc, char *argv[]){
		char * str = argv[1];
		int len = strlen(str);
		getReverseStr_char(str, len);
		fprintf(stderr, "%s\n", str);
		return 0;
	}
};

struct KMER_COUNTING_handler{

#define KMER_COUNT_LEN 20
void simple_kmer_counter(char *fn){
	fprintf(stderr, "version 1.00\n");
	extern uint64_t kmerMask[33];
	FILE *f_read = xopen(fn, "r");
	char buff_[1000];
	uint8_t buff_bin[1000];
	uint8_t buff_bin_rev[1000];
	std::map<uint64_t, int> kmer_set;
	std::map<uint64_t, int>::iterator kmer_set_it;
	std::map<uint64_t, int> global_kmer_set;

	uint64_t MASK = kmerMask[KMER_COUNT_LEN];
	int read_number = 0;
	while(fgets(buff_, 1000, f_read)){
		read_number++;
		if(read_number % 100 == 0)
			fprintf(stderr, "%d\r", read_number);
		//find \t
		int start_idx = 0;
		while(buff_[start_idx++] != '\t');
		char * str_p = buff_ + start_idx;
		 int string_n = 0;
		for(;str_p[string_n] != '\n';string_n++)
		{
			switch(str_p[string_n]){
			case 'A': buff_bin[string_n] = 0; break;
			case 'C': buff_bin[string_n] = 1; break;
			case 'G': buff_bin[string_n] = 2; break;
			case 'T': buff_bin[string_n] = 3; break;
			default : xassert(0, "");
			}
		}
		//reverse string:
		for(int i = 0; i < string_n; i++){
			buff_bin_rev[string_n - i - 1] = 3 - buff_bin[i];
		}
		uint64_t kmer     = bit2_nextKmer_init(buff_bin, KMER_COUNT_LEN);
		uint64_t kmer_rev = bit2_nextKmer_init(buff_bin_rev, KMER_COUNT_LEN);
		int kmer_number = string_n - KMER_COUNT_LEN + 1;
		kmer_set.clear();
		for(int i = 0; i < kmer_number; i++){
			kmer     = bit2_nextKmerMASK( buff_bin     + i, kmer, KMER_COUNT_LEN);
			kmer_set_it = kmer_set.find(kmer);		if(kmer_set_it!=kmer_set.end()){kmer_set_it->second ++;	}	else{kmer_set[kmer] = 1;	}
		}
		for(int i = 0; i < kmer_number; i++){
			kmer_rev = bit2_nextKmerMASK( buff_bin_rev + i, kmer_rev, KMER_COUNT_LEN);
			kmer_set_it = kmer_set.find(kmer_rev);	if(kmer_set_it!=kmer_set.end()){kmer_set_it->second ++;	}	else{kmer_set[kmer_rev] = 1;}
		}
		kmer     = bit2_nextKmer_init(buff_bin, KMER_COUNT_LEN);
		kmer_rev = bit2_nextKmer_init(buff_bin_rev, KMER_COUNT_LEN);

		for(kmer_set_it = kmer_set.begin(); kmer_set_it != kmer_set.end(); kmer_set_it++){
			if(kmer_set_it->second >= 3){
				auto g_it = global_kmer_set.find(kmer_set_it->first);
				if(g_it!=global_kmer_set.end()){g_it->second += kmer_set_it->second;}	else{global_kmer_set[kmer_set_it->first] = kmer_set_it->second;	}
			}
		}
	}

	for(auto it = global_kmer_set.begin(); it != global_kmer_set.end(); it++){
		char kmerStr[100];
		uint64_t kmer = it->first;
		for(int i = 0; i < KMER_COUNT_LEN; i++){
			kmerStr[i] = "ACGT"[(kmer >> ((KMER_COUNT_LEN - 1 - i) * 2)) & 0x3];
		}
		kmerStr[KMER_COUNT_LEN] = 0;
		fprintf(stdout, "%d\t%s\n", it->second, kmerStr);
	}

	fclose(f_read);
}

int read_ACGT_analysis(int argc, char *argv[]){
	const char * input_bam_fn = argv[1]; const char * ref_fn = argv[2];
	Bam_file c_b;
	memset(&c_b, 0, sizeof(Bam_file));
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
	bam_hdr_t* hdr = c_b._hdr;

	bam1_t b = {0};//BAM record for the first read in a pair

	char * seq_buff = (char *)xmalloc(10000);
	int * analysis_cnt = (int *)xcalloc(10000, sizeof(int));
	int read_ID = 0;
	while (sam_read1(c_b._hfp, hdr, &b) >= 0){

		const int read_len = b.core.l_qseq;
		get_bam_seq(0, read_len, seq_buff, &b);//store in binary format
		int A_count = 0;
		int C_count = 0;
		int G_count = 0;
		int T_count = 0;
		int N_count = 0;
		for(int i = 0; i < read_len; i++){
			switch(seq_buff[i])
			{
			case 'A': case 'a': A_count++; break;
			case 'C': case 'c': C_count++; break;
			case 'G': case 'g': G_count++; break;
			case 'T': case 't': T_count++; break;
			default: N_count++; break;
			}
		}

		if(N_count > 25){
			fprintf(stderr, "[N_count %d @ %d]%s\n", N_count, read_ID ,seq_buff );
		}

//		if(A_count >= 75 || C_count >= 75 || G_count >= 75 || T_count >= 75){
//			fprintf(stdout, "[%d %d %d %d]@ %d %s\n", A_count, C_count, G_count, T_count, read_ID, seq_buff );
//		}

		if((G_count >= 75 && C_count < 20) || (C_count >= 75 && G_count < 20)){
			fprintf(stdout, "[%d %d %d %d]@ ID %d @ pos %d:%d %s\n", A_count, C_count, G_count, T_count, read_ID,b.core.tid, b.core.pos, seq_buff );
		}

		A_count = A_count / 16;
		C_count = C_count / 16;
		G_count = G_count / 16;
		T_count = T_count / 16;

		int count_number = 0;
		xassert(A_count < 10, ""); count_number*= 10; count_number += A_count;
		xassert(C_count < 10, ""); count_number*= 10; count_number += C_count;
		xassert(G_count < 10, ""); count_number*= 10; count_number += G_count;
		xassert(T_count < 10, ""); count_number*= 10; count_number += T_count;

		analysis_cnt[count_number] ++;
		if(read_ID % 1000000 == 0){
			for(int i = 0; i < 10000; i++){
				if(analysis_cnt[i] > 0){
					fprintf(stderr, "[%d %d]\n", i, analysis_cnt[i] );
				}
			}
		}
		read_ID++;
	}

	bam_file_close(&c_b);

	return 0;
}


};

struct SV_REGION_handler{

	//gcSV tools vntr_analysis ref.fa single_sample_vcf_list.txt VNTR_region.bed
	//1kgp_sample_info.txt

	#define MIN_GCA_sv_len_GIAB 50
	//get the SV in the original GIAB vcf file(mixed with all types of variation) and store it into output file
	//file type: GIAB CCS DATA: deepvariant_.GRCh38_15kb_37X_SequelII.vcf
	//Only get SV when : the ref_len >= 50 base pair or the Alt_len >= 50 base pair
	void vcf_GIAB_getSV(char *fn_in, char * fn_out)
	{
		bool is_compression = false;//BCF or VCF
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read
		BCF_FILE vcf_w;//vcf for write
		VCF_open_write(&vcf_w, fn_out, is_compression);
		bcf_hdr_write(vcf_w.file, vcf_r.header);
		char *c_sample = (char *)malloc(1000);
		char *c_sv_type = (char *)malloc(1000);
		int all_sample = false;
		int all_type = false;
		int all_chrom = false;
		char sample_name[10] = "ALL";
		char sv_type[10] = "ALL";
		char CHROM_ID[10] = "ALL";
		if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)	all_sample = true;
		if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)			all_type = true;
		if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)		all_chrom = true;
		while(VCF_next(&vcf_r))//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			//unpack the vcf data to get the alt string
			bcf_unpack(c_r, BCF_UN_STR);
			//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
			if(all_sample == false)
			{
				vcf_get_sample(vcf_r.header, c_r, c_sample);
				if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
					continue;
			}
			if(all_type == false)
			{
				vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
				if(strcmp(c_sv_type, sv_type) != 0)
					continue;
			}
			if(all_chrom == false)
			{
				if(c_r->rid != strtol(CHROM_ID,0,10))
					continue;
			}
			//check the alleles
			uint32_t refLen = strlen(c_r->d.allele[0]);
			uint32_t maxAlleleLen = 0;
			for(uint32_t i = 1; i < c_r->n_allele; i++){
				uint32_t allele_len = strlen(c_r->d.allele[i]);
				maxAlleleLen = MAX(maxAlleleLen, allele_len);
			}
			if(refLen >= MIN_GCA_sv_len_GIAB || maxAlleleLen >= MIN_GCA_sv_len_GIAB){
				fprintf(stderr, "%d\t%d\t%d\t%d\t\n", c_r->rid, c_r->pos, refLen, maxAlleleLen);
				bcf_write(vcf_w.file, vcf_r.header, c_r);
			}

			//vcf_get_genotype(vcf_r.header, c_r);
		}
		//close
		bcf_close(vcf_r.file);
		bcf_close(vcf_w.file);
	}

	//get the SV in the original GIAB vcf file(mixed with all types of variation) and store it into output file
	//file type: GIAB CCS DATA: deepvariant_.GRCh38_15kb_37X_SequelII.vcf
	//Only get SV when : the ref_len >= 50 base pair or the Alt_len >= 50 base pair
	void vcf_SURVIVOR_getSV(char *fn_in, char *sv_type)
	{
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read

		char *c_sample = (char *)malloc(1000);
		char *c_sv_type = (char *)malloc(1000);
		int all_sample = false;
		int all_type = false;
		int all_chrom = false;
		char sample_name[10] = "ALL";
		char CHROM_ID[10] = "ALL";
		bcf_hdr_t *header = vcf_r.header;

		std::set<std::string> SV_TYPE_SET;

		if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)	all_sample = true;
		if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)			all_type = true;
		if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)		all_chrom = true;

		do//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			//unpack the vcf data to get the alt string
			bcf_unpack(c_r, BCF_UN_STR);
			//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
			if(all_sample == false)
			{
				vcf_get_sample(vcf_r.header, c_r, c_sample);
				if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
					continue;
			}
			if(all_type == false)
			{
				vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
				std::string str_SV_TYPE(c_sv_type);
				SV_TYPE_SET.emplace(str_SV_TYPE);
				if(strcmp(c_sv_type, sv_type) != 0)
					continue;
			}
			if(all_chrom == false)
			{
				if(c_r->rid != strtol(CHROM_ID,0,10))
					continue;
			}
			//get the start position:
			int32_t stPos = c_r->pos;
			//get the end position:
			int32_t edPos = 0;
			vcf_get_sv_END(vcf_r.header, c_r, &edPos);
			const char *chrID = bcf_hdr_id2name(header, c_r->rid);
			int32_t SV_LEN = edPos - stPos;
			if(SV_LEN < 0 || SV_LEN > 10000)
				SV_LEN = 0;
			printf("%s\t%d\t%d\t\n",chrID , stPos, SV_LEN);// chrID+st+len

		}while(VCF_next(&vcf_r));

		for(auto & t :SV_TYPE_SET)
			fprintf(stderr, "%s\t",t.c_str());

		//close
		bcf_close(vcf_r.file);
	}
#define MAX_LINE_LENGTH 1000
	void SVRegionCombine(char* regionListFN, int SV_EDGE_LENGTH){
		//step 1: load region, 其中每一行的c_r->pos - 1000 作为开始 position； c_r->pos + refLen + 1000 作为结束 position
		char *temp = new char[MAX_LINE_LENGTH];//1M
		std::vector<RefRegion>regionList;
		std::ifstream regionListFile(regionListFN);
		std::map<std::string, uint32_t> name2ID;
		std::vector<std::string> id2Name;

		while(true){
			regionListFile.getline(temp, MAX_LINE_LENGTH);
			if(*temp == 0)
				break;
			//get chrID
			char *token = strtok(temp, "\t");
			std::string strID(token); int intID = 0;
			std::map<std::string, uint32_t>::iterator it = name2ID.find(strID);
			if(it!=name2ID.end()){
				intID = it->second;
			}
			else{
				intID = name2ID.size();
				name2ID[strID] = name2ID.size();
				id2Name.emplace_back(strID);
			}

			//get refSt
			token = strtok(NULL, "\t");
			uint32_t refSt = strtoul(token, NULL, 10);
			//get refEd
			token = strtok(NULL, "\t");
			uint32_t refLen = strtoul(token, NULL, 10);
			regionList.emplace_back(intID, refSt - SV_EDGE_LENGTH, refSt + refLen + SV_EDGE_LENGTH);
		}
		regionListFile.close();
		delete [] temp;
		std::sort(regionList.begin(), regionList.end(), RefRegion::cmp_by_pos);
		//merge region
		auto r_ed = regionList.end();//for each sve
		int total_output_ITEM = 0;
		for(auto r = regionList.begin(); r < r_ed;)	{
			auto r_try = r + 1;
			for(; r_try < r_ed && r->region_overlap(*r_try); r_try++)//use sve as main, than try to combine
				r->Combine(*r_try, true);
			printf("%s\t%d\t%d\t\n", id2Name[r->chr_ID].c_str(), r->st_pos, r->ed_pos);
			total_output_ITEM++;
			r = r_try;
		}
		fprintf(stderr, "total_output_ITEM: %d\n", total_output_ITEM);
	}

	#define MAX_read_LEN 256
	void GIAB_SV_region_full_test(char *input_bam_fn, char *ref_fn, char * region_fn, int readLen)
	{
		xassert(readLen < MAX_read_LEN, "MAX read length: 250 bp");
		fprintf(stderr, "\n\n V1.02\n\n");
		//read bam/cram file:
		Bam_file c_b;
		bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
		bam_hdr_t* hdr = c_b._hdr;

		//get sv region:
		char *temp = new char[MAX_LINE_LENGTH];//1M
		std::vector<R_region>regionList;
		std::ifstream regionListFile(region_fn);
		while(true){
			regionListFile.getline(temp, MAX_LINE_LENGTH);
			if(*temp == 0 || *temp == 'N')
				break;
			//get chrID
			char *token = strtok(temp, "\t");
			const char *strID = token;
			int chrID = bam_name2id(hdr, strID) ;

			//get refSt
			token = strtok(NULL, "\t");
			uint32_t regionSt = strtoul(token, NULL, 10);
			//get refEd
			token = strtok(NULL, "\t");
			uint32_t regionEnd = strtoul(token, NULL, 10);
			//
			regionList.emplace_back();
			R_region & r = regionList.back();
			r.chr_ID = chrID; r.st_pos = regionSt; r.ed_pos = regionEnd;
		}
		regionListFile.close();
		delete [] temp;

		char SA_TAG[3] = "SA";
		for(R_region & r: regionList){
			//reset region
			resetRegion_ID(&c_b, &r);
			//reset calculators
			uint32_t readNum = 0;
			//mapQ, region 0~60, if mapQ > 60; MAPQ[61]++
			uint32_t MAPQCount[62] = {0};
			uint32_t MAPQ_less_than30 = 0;
			//insert size: when it is [0 ~1023]: 0~15 store in insertSize[0], 16~31 in insertSize[1].....
			uint32_t insertSizeCount[64] = {0};
			uint32_t insertSizeOver1K = 0;//the read of insert size > 1023
			uint64_t totalInsert_SIZE = 0;
			uint32_t insertSizeOver100K = 0;//the read of insert size > 10000
			//mate unmapped:
			uint32_t mate_unmapped_count = 0;
			// mate_different_chromsome
			uint32_t mate_different_chromsome_count = 0;

			uint32_t insertSizeFF = 0;//both forward:
			uint32_t insertSizeRR = 0;//both reverse:
			uint32_t insertSizeRF = 0;//the read of smaller position is reverse, the other is forward
			uint32_t insertSizeFR = 0;//the read of smaller position is forward, the other is reverse

			//SA signal
			uint32_t SA_read_len[MAX_read_LEN] = {0};
			//SNP and INDEL
			uint32_t SNP_INDEL_len[MAX_read_LEN] = {0};
			uint32_t SNP_INDEL_total = 0;
			//soft/hard clip
			uint32_t clip_len_left[MAX_read_LEN] = {0};
			uint32_t clip_len_right[MAX_read_LEN] = {0};
			uint32_t total_clip_len = 0;
			//flags
			uint32_t flag_count[MAX_uint16_t] = {0};

			// analysis the reads
			while (bam_next(&c_b)) {
				bam1_t *b = &(c_b._brec);
				readNum++;
				fprintf(stderr, "%d\n", readNum);
				//mapq:
				uint8_t mapq = b->core.qual;
				if(mapq > 60)	MAPQCount[61]++;
				else{			MAPQCount[mapq]++;	if(mapq < 30)	MAPQ_less_than30 ++;	}
				bool abnormal_insert_size = false;
				//mate read:
				if(bam_is_mate_unmapped(b)) 	{abnormal_insert_size = true; mate_unmapped_count ++;}
				if(b->core.tid != b->core.mtid)	{abnormal_insert_size = true; mate_different_chromsome_count ++;}
				//insert size
				if(!abnormal_insert_size){
					int32_t insertSizeOri = b->core.isize;
					uint32_t insertSize = ABS(insertSizeOri);
					if(insertSize <= 1023) insertSizeCount[insertSize/16]++;
					else insertSizeOver1K++;

					if(insertSize <= 100000) { totalInsert_SIZE += insertSize;}
					else insertSizeOver100K++;

				}
				//read pair orientation
				bool direction = bam_is_fwd_strand(b);
				bool mateDirection = bam_is_mate_fwd_strand(b);
				if(b->core.pos > b->core.mpos) std::swap(direction, mateDirection);

				if		(direction == FORWARD && mateDirection == FORWARD) insertSizeFF++;
				else if (direction == FORWARD && mateDirection == REVERSE) insertSizeFR++;
				else if (direction == REVERSE && mateDirection == FORWARD) insertSizeRF++;
				else if (direction == REVERSE && mateDirection == REVERSE) insertSizeRR++;

				//SA signal:
				const char * SA_String = bam_get_string_tag(b, SA_TAG);
				if(SA_String != NULL){
					uint32_t SA_len = strlen(SA_String);
					SA_read_len[SA_len]++;
				}

				//SNP and INDEL
				int INDEL_len = 0;
				int NM_len = 0;
				bam_get_INDEL_NM(b,&INDEL_len, &NM_len);

				SNP_INDEL_len[NM_len]++;
				SNP_INDEL_total+= NM_len;

				//clip:
				int soft_left; int soft_right;
				bam_has_SH_cigar(b, &soft_left, &soft_right);
				clip_len_left[soft_left]++;
				clip_len_right[soft_right]++;
				total_clip_len += (soft_left + soft_right);

				//flag
				flag_count[b->core.flag]++;
			}

			// out put information for that region
			//title:
			printf("[%d:%d~%d]", r.chr_ID, r.st_pos, r.ed_pos);
			printf("\t[RN:%d]", readNum);
			float readDepth = (float)readLen * readNum / (r.ed_pos - r.st_pos);
			printf("\t[DP:%f]", readDepth);
			if(readNum == 0) readNum = 1;
			printf("\t[MUM:%d, %f%%]", mate_unmapped_count,(float)mate_unmapped_count*100/readNum );
			printf("\t[MAPQ<30:%f%%]", (float)MAPQ_less_than30*100/readNum);
			printf("\t[AVG.ISIZE(<100000):%f]", (float)totalInsert_SIZE/(readNum - insertSizeOver100K));
			printf("\t[>1023:%f%%]", (float)insertSizeOver1K*100/readNum);

			printf("\t[AVG.SNP_INDEL:%f]", (float)SNP_INDEL_total/(readNum));
			printf("\t[AVG.CLIP:%f]", (float)total_clip_len/(readNum));
			printf("\t[AVG.FR:%f%%]", (float)insertSizeFR*100/(readNum));

			//mapq:
			printf("\tMAPQ: ");
			for(int i = 0; i < 61; i++)
				if(MAPQCount[i] > 0)
					printf("[%d:%d] ", i, MAPQCount[i]);
			if(MAPQCount[61] > 0)
				printf("[> 60:%d] ", MAPQCount[61]);
			//mate_unmapped_count:

			//mate_different_chromsome_count:
			printf("\t[MDC:%d]", mate_different_chromsome_count);
			//insert size:
			printf("\tIS: ");
			for(int i = 0; i < 64; i++)
				if(insertSizeCount[i] > 0)
					printf("[%d:%d] ", i, insertSizeCount[i]);
			if(insertSizeOver1K > 0)
				printf("[> 1023:%d] ", insertSizeOver1K);
			//orientation
			printf("\t[FF:%d]", insertSizeFF);
			printf("\t[FR:%d]", insertSizeFR);
			printf("\t[RF:%d]", insertSizeRF);
			printf("\t[RR:%d]", insertSizeRR);
			// SA len
			printf("\tSA_L: ");
			for(int i = 0; i < MAX_read_LEN; i++)
				if(SA_read_len[i] > 0)
					printf("[%d:%d] ", i, SA_read_len[i]);
			//SNP INDEL
			printf("\tSNP_INDEL: ");
			for(int i = 0; i < MAX_read_LEN; i++)
				if(SNP_INDEL_len[i] > 0)
					printf("[%d:%d] ", i, SNP_INDEL_len[i]);
			//clip_len_left
			printf("\tCL: ");
			for(int i = 0; i < MAX_read_LEN; i++)
				if(clip_len_left[i] > 0)
					printf("[%d:%d] ", i, clip_len_left[i]);
			//clip_len_right
			printf("\tCR: ");
			for(int i = 0; i < MAX_read_LEN; i++)
				if(clip_len_right[i] > 0)
					printf("[%d:%d] ", i, clip_len_right[i]);
			//flags
			printf("\tFLAG: ");
			for(int i = 0; i < MAX_uint16_t; i++)
				if(flag_count[i] > 0)
					printf("[%d:%d] ", i, flag_count[i]);
			//END

			printf("\n");
		}
		//close file
		bam_file_close(&c_b);
	}
	void vcf_sample(char *fn_in, char * fn_out, char * sample_name, char * sv_type, char * CHROM_ID)
	{
		bool is_compression = false;//BCF or VCF
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, fn_in);//open for read
		BCF_FILE vcf_w;//vcf for write
		VCF_open_write(&vcf_w, fn_out, is_compression);
		bcf_hdr_write(vcf_w.file, vcf_r.header);
		char *c_sample = (char *)malloc(1000);
		char *c_sv_type = (char *)malloc(1000);
		int all_sample = false;
		int all_type = false;
		int all_chrom = false;
		if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)
			all_sample = true;
		if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)
			all_type = true;
		if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)
			all_chrom = true;
		while(VCF_next(&vcf_r))//read one
		{
			bcf1_t *c_r = &( vcf_r.r);
			//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
			if(all_sample == false)
			{
				vcf_get_sample(vcf_r.header, c_r, c_sample);
				if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
					continue;
			}
			if(all_type == false)
			{
				vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
				if(strcmp(c_sv_type, sv_type) != 0)
					continue;
			}
			if(all_chrom == false)
			{
				if(c_r->rid != strtol(CHROM_ID,0,10))
					continue;
			}
			if(c_r->rid > 24)//discard decoy sequence
				continue;
			bcf_write(vcf_w.file, vcf_r.header, c_r);
			//vcf_get_genotype(vcf_r.header, c_r);
		}
		//close
		bcf_close(vcf_r.file);
		bcf_close(vcf_w.file);
	}

	int vcf_dump(int argc, char *argv[]){
		char *fn_in = argv[1];
		char * fn_out = argv[2];
		char * sample_name = argv[3];
		char * sv_type = argv[4];
		char * CHROM_ID = argv[5];
		vcf_sample(fn_in, fn_out, sample_name, sv_type, CHROM_ID);
		return 0;
	}

	struct SV_INFO_vcf_add_alt_string{
		std::vector<std::string> info;
		int sample_ID;
		int SV_len;
		int chrID;
		int pos;
		int sample_numer;//todo::
		int max_same_sv_len_sample_numer;//todo::

		std::string ALT_str;
		std::string REF_str;
		int true_sample_pos;

		void store_basic(std::vector<std::string> &item_value, faidx_t *c_ref_idx){
			for(int j = 0; j < 9; j ++)
				info.emplace_back(item_value[j]);
			chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			pos = atoi(item_value[1].c_str());
		}
	};

	struct SAMPLE_INFO_vcf_add_alt_string{
		std::string sample_name;
		std::vector<std::vector<int>> SV_IDs;
	};

	int vcf_add_alt_string(int argc, char *argv[]){
		//char *joint_vcf_fn = argv[1];
		faidx_t *c_ref_idx = reference_index_load("/media/fenghe/Data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa");
		const char *joint_vcf_fn = "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/filltag_S5.vcf";
		std::vector<SV_INFO_vcf_add_alt_string> SV_info;
		std::vector<SAMPLE_INFO_vcf_add_alt_string> sample_info;
		//load joint data
		{
			std::vector<std::string> joint_vcf_data;
			load_string_list_from_file_MAX_line(joint_vcf_fn, joint_vcf_data, 5000);
			size_t line_num = joint_vcf_data.size();//skip the header line
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_2;
			std::map <int, int> SVLEN_counter;
			for(uint64_t i = 0; i < line_num; i++){
				//show header
				const char * line_data = joint_vcf_data[i].c_str();
				if(line_data[0] == '#' && line_data[1] == '#'){
					printf("%s\n", joint_vcf_data[i].c_str());
				}
				else if(line_data[0] == '#'){
					split_string(item_value, joint_vcf_data[i].c_str(), "\t");//skip the header line
					for(int j = 0; j < 9; j ++)
						printf("%s\t", item_value[j].c_str());
					printf("SAMPLE\n");
					//store sample name
					for(uint j = 9; j < item_value.size(); j ++){
						sample_info.emplace_back();
						sample_info.back().sample_name = item_value[j];
						sample_info.back().SV_IDs.resize(24);
					}
				}
				else{
					split_string(item_value, joint_vcf_data[i].c_str(), "\t");//skip the header line
					SV_info.emplace_back();
					SV_info.back().store_basic(item_value,c_ref_idx);
					//get the sample with max same SV length number
					SVLEN_counter.clear();
					for(uint j = 9; j < item_value.size(); j++){
						split_string(item_value_2, item_value[j].c_str(), ":");//skip the header line
						if(strcmp(item_value_2[0].c_str(), "./.") != 0){
							int SVLEN = atoi(item_value_2[1].c_str());
							std::map<int, int >::iterator it = SVLEN_counter.find(SVLEN);
							if(it != SVLEN_counter.end()){
								it->second++;
							}else
								SVLEN_counter[SVLEN] = 1;
						}
					}
					int max_SV_len = -1;
					int max_count = -1;
					for(std::map<int, int >::iterator it = SVLEN_counter.begin(); it != SVLEN_counter.end(); it++){
						if(max_count < it->second){
							max_count = it->second;
							max_SV_len = it->first;
						}
					}
					SV_info.back().SV_len = max_SV_len;
					for(uint j = 9; j < item_value.size(); j ++){
						split_string(item_value_2, item_value[j].c_str(), ":");//skip the header line
						if(strcmp(item_value_2[0].c_str(), "./.") != 0){
							int SVLEN = atoi(item_value_2[1].c_str());
							if(SVLEN == max_SV_len){
								SV_info.back().sample_ID = j - 9;
								sample_info[SV_info.back().sample_ID].SV_IDs[SV_info.back().chrID].emplace_back(SV_info.size() - 1);
								break;
							}
						}
					}
				}
			}
		}

		//recover from single sample vcf data:
		{
			//load single sample vcf
			for(uint i = 0; i < sample_info.size(); i++){
				//if(i > 30) break;
				char vcf_fn[1024];
				sprintf(vcf_fn, "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/sorted_vcf/%s_sv_spa.sort.vcf.gz", sample_info[i].sample_name.c_str());
				fprintf(stderr, "Current handle %s\n", vcf_fn);
				BCF_FILE input_vcf;

				VCF_open_read(&input_vcf, vcf_fn);
				bcf_hdr_t *vcf_header = input_vcf.header;
				int recover_SV_num = 0;
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
					int SV_POS = c_r->pos;
					int SV_length = 0;
					vcf_get_sv_LENGTH(input_vcf.header, c_r, &SV_length);
					//search
					std::vector<int> &SV_list = sample_info[i].SV_IDs[chrID];
					//search if this SV is needed
					for(int SV_ID: SV_list){
						if(ABS_U(SV_info[SV_ID].pos,SV_POS) < 1000){
							if(ABS(SV_length) == ABS(SV_info[SV_ID].SV_len)){
								//store the ALT string; //store the REF string
								SV_info[SV_ID].REF_str.append(c_r->d.allele[0]);
								SV_info[SV_ID].ALT_str.append(c_r->d.allele[1]);
								SV_info[SV_ID].true_sample_pos = SV_POS;
								recover_SV_num++;
							}
						}
					}
				}while(VCF_next(&input_vcf));//read one
				fprintf(stderr, "\n recover_SV_num %d\n", recover_SV_num);
				bcf_close(input_vcf.file);
			}

			//output the final result
			{
				//std::vector<SV_INFO_vcf_add_alt_string> SV_info;
				int SV_number = SV_info.size();
				for(int i = 0; i < SV_number; i++){
					SV_INFO_vcf_add_alt_string & sv = SV_info[i];
					std::vector<std::string> & info = sv.info;
					bool STR_reset = (!sv.ALT_str.empty());
					int store_pos = (STR_reset)?(sv.true_sample_pos + 1):(sv.pos);
					printf("%s\t%d\t%s\t", info[0].c_str(), store_pos, info[2].c_str());//chrID, pos, name
					if(STR_reset){
						printf("%s\t%s\t", sv.REF_str.c_str(), sv.ALT_str.c_str());//chrID
					}else
						printf("%s\t%s\t", info[3].c_str(), info[4].c_str());//chrID
					// sample_info[sv.sample_ID].sample_name.c_str()
					printf("%s\t%s\t%s\tGT\t0/1\n", info[5].c_str(), info[6].c_str(), info[7].c_str() );//chrID, pos, name
				}
			}

		}
		fprintf(stderr, "\nEND\n");
		return 0;
	}

};

struct LIFT_OVER_HANDLER{
	struct ASS_INFO{
		ASS_INFO(){}
		ASS_INFO(const char * name_){strcpy(name, name_); }
		uint32_t length = 0;
		char name[100];
		void print(FILE * o, char endl){
			fprintf(o, "[%s %d]\t%c", name, length, endl);
		}
	};

	struct LIFTOVER{
		LIFTOVER(
				uint32_t tid_ass_,  uint32_t pos_ass_st_,  uint32_t pos_ass_ed_,
				uint32_t tid_hs37_, uint32_t pos_hs37_st_, uint32_t pos_hs37_ed_,
				uint32_t length_, bool direction_):
					tid_ass(tid_ass_),   pos_ass_st(pos_ass_st_),   pos_ass_ed(pos_ass_ed_),
					tid_hs37(tid_hs37_), pos_hs37_st(pos_hs37_st_), pos_hs37_ed(pos_hs37_ed_),
					length(length_), direction(direction_){}
		LIFTOVER():
			tid_ass(0),  pos_ass_st(0),  pos_ass_ed(0),
			tid_hs37(0), pos_hs37_st(0), pos_hs37_ed(0),
			length(0), direction(0){}

		uint32_t tid_ass;
		uint32_t pos_ass_st;
		uint32_t pos_ass_ed;

		uint32_t tid_hs37;
		uint32_t pos_hs37_st;
		uint32_t pos_hs37_ed;

		uint32_t length;

		bool direction;

		void print(FILE * o, std::vector<ASS_INFO> & assemble_name_str_list){
			fprintf(o,
					"\t\t[ASS:INFO %s %d ]"
					" [ %d %d %d ]"
					"[ %d %d %d ] "
					"%d %c\n",
					assemble_name_str_list[tid_ass].name, assemble_name_str_list[tid_ass].length,
					tid_ass,  pos_ass_st,  pos_ass_ed,
					tid_hs37, pos_hs37_st, pos_hs37_ed,
					length, (direction==FORWARD)?'+':'-');
		}

		static inline int cmp_by_ref_pos(const LIFTOVER &a, const LIFTOVER &b){
			if(a.tid_hs37 != b.tid_hs37)
				return a.tid_hs37 < b.tid_hs37;
			else
				return a.pos_hs37_st < b.pos_hs37_st;
		}

		static inline int cmp_by_ass_pos(const LIFTOVER &a, const LIFTOVER &b){//todo::
			if(a.tid_ass != b.tid_ass)
				return a.tid_ass < b.tid_ass;
			else
				return a.pos_ass_st < b.pos_ass_st;
		}

	};

	struct liftoverRegion_direction{

		liftoverRegion_direction(
				uint32_t tid_ass_,  uint32_t pos_ass_st_,  uint32_t pos_ass_ed_,
				bool direction_):
					tid_ass(tid_ass_),   pos_ass_st(pos_ass_st_),   pos_ass_ed(pos_ass_ed_),
					direction(direction_){}
		liftoverRegion_direction(){}

		uint32_t tid_ass = 0;
		uint32_t pos_ass_st  = 0;
		uint32_t pos_ass_ed = 0;

		bool direction = 0;

		static inline int cmp_by_ass_pos(const liftoverRegion_direction &a, const liftoverRegion_direction &b){//todo::
			if(a.tid_ass != b.tid_ass)
				return a.tid_ass < b.tid_ass;
			else
				return a.pos_ass_st < b.pos_ass_st;
		}
	};

	struct LIFTOVER_index{
		std::map<std::string, uint32_t> assembly_name_2_ID;
		std::vector<ASS_INFO> id2Name;
		std::vector<LIFTOVER> liftIndex;

		std::vector<liftoverRegion_direction> LD;

		void dump(char * fn){
			FILE * f = xopen(fn, "wb");
			uint64_t id2Name_size = id2Name.size();
			fwrite(&id2Name_size, sizeof(uint64_t), 1, f);
			fwrite(&(id2Name[0]), sizeof(ASS_INFO), id2Name_size, f);

			uint64_t liftIndex_size = liftIndex.size();
			fwrite(&liftIndex_size, sizeof(uint64_t), 1, f);
			fwrite(&(liftIndex[0]), sizeof(LIFTOVER), liftIndex_size, f);

			//LD
			uint64_t LD_size = LD.size();
			fwrite(&LD_size, sizeof(uint64_t), 1, f);
			fwrite(&(LD[0]), sizeof(liftoverRegion_direction), LD_size, f);

			fclose(f);
		}

		void load(char * fn){
			FILE * f = xopen(fn, "rb");

			uint64_t id2Name_size = 0;
			xread(&id2Name_size, sizeof(uint64_t), 1, f);
			id2Name.clear();
			id2Name.resize(id2Name_size);
			xread(&(id2Name[0]), sizeof(ASS_INFO), id2Name_size, f);

			uint64_t liftIndex_size = 0;
			xread(&liftIndex_size, sizeof(uint64_t), 1, f);
			liftIndex.clear();
			liftIndex.resize(liftIndex_size);
			xread(&(liftIndex[0]), sizeof(LIFTOVER), liftIndex_size, f);

			//LD
			uint64_t LD_size = 0;
			xread(&LD_size, sizeof(uint64_t), 1, f);
			LD.clear();
			LD.resize(liftIndex_size);
			xread(&(LD[0]), sizeof(liftoverRegion_direction), LD_size, f);

			//debug code:
			//for(auto l: liftIndex)	{l.print();}
			fclose(f);
		}

	};

	void liftoverBuildingIndex(char *input_bam_fn, char *liftOverData_fn){
		fprintf(stderr, "\n\n V1.09\n\n");
			//read bam/cram file:
		Bam_file c_b;
		bam_file_open(input_bam_fn, NULL, NULL, &c_b);
		bam_hdr_t* hdr = c_b._hdr;
		uint64_t readNum = 0;

		//print header
		for(int i = 0; i < hdr->n_targets; i++)
			fprintf(stderr, "%d %s %d \n", i, hdr->target_name[i], hdr->target_len[i]);

		bam1_t b = {0};//BAM record for the first read in a pair
		// analysis the reads'
		//reset region
		//resetRegion_ID(&c_b, &r);
		int sam_rst1 = 0;

		LIFTOVER_index idx;
		std::map<std::string, uint32_t> & assembly_name_2_ID = idx.assembly_name_2_ID;
		std::vector<ASS_INFO> &id2Name = idx.id2Name;
		std::vector<LIFTOVER> &liftIndex = idx.liftIndex;
		std::vector<liftoverRegion_direction> &LD =  idx.LD;

		while (1){
			//load SAM 1 & 2
			readNum++; // debug code //if(readNum == 10) break;
			sam_rst1 = sam_read1(c_b._hfp, hdr, &b);
			if(sam_rst1 < 0) break;
			//load the read name
			std::string strID((char *)bam_qname(&b)); int assemblyID = 0;
			std::map<std::string, uint32_t>::iterator it = assembly_name_2_ID.find(strID);
			if(it!=assembly_name_2_ID.end()){assemblyID = it->second;}
			else{
				assemblyID = assembly_name_2_ID.size();
				assembly_name_2_ID[strID] = assembly_name_2_ID.size();
				id2Name.emplace_back(strID.c_str());
			}
			if(!bam_is_secondary(&b) && !bam_is_supplementary(&b) && !bam_is_duplicate(&b))
				id2Name[assemblyID].length = b.core.l_qseq;
			//get cigar
			bool direction1 = bam_is_fwd_strand(&b);

			uint32_t tid = b.core.tid;
			uint64_t st_pos_in_ref = b.core.pos;
			uint64_t st_pos_in_assembly = 0;
			uint32_t* bam_cigar = bam_get_cigar(&b);

			for (uint32_t i = 0; i < b.core.n_cigar; ++i)
			{
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				switch (type)
				{
				case CIGAR_MATCH:
				case CIGAR_SEQ_MATCH:
					liftIndex.emplace_back(
						assemblyID, st_pos_in_assembly, st_pos_in_assembly + length,
						tid, st_pos_in_ref, st_pos_in_ref + length,
						length, direction1);
					//liftIndex.back().print(stderr);
					st_pos_in_ref += length; st_pos_in_assembly += length;
					break;
				case CIGAR_INSERT:
					st_pos_in_assembly += length;
					break;
				case CIGAR_DELETE:
					st_pos_in_ref += length;
					break;
				case CIGAR_SKIP:
					break;
				case CIGAR_SOFT_CLIP:
				case CIGAR_HARD_CLIP:
					st_pos_in_assembly += length;
					break;
				case CIGAR_PAD:
					break;
				case CIGAR_SEQ_MISMATCH:
					break;
				default:
					break;
				}
			}
			int cigar_st_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
			int cigar_st_length = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
			uint32_t region_st_in_ass = (cigar_st_type == CIGAR_SOFT_CLIP || cigar_st_type == CIGAR_HARD_CLIP)?cigar_st_length:0;
			LD.emplace_back(assemblyID, region_st_in_ass, liftIndex.back().pos_ass_ed,direction1);
		}
		//close file
		bam_file_close(&c_b);

		//modify LD
		for(auto &l : LD){
			if(l.direction == REVERSE){
				uint32_t ass_len = id2Name[l.tid_ass].length;
				l.pos_ass_st = ass_len - l.pos_ass_st;
				l.pos_ass_ed = ass_len - l.pos_ass_ed;
				std::swap(l.pos_ass_st, l.pos_ass_ed);
			}
		}
		//sort by pos
		std::sort(LD.begin(), LD.end(), liftoverRegion_direction::cmp_by_ass_pos);

		//modify the position for REVERSE records
	//	for(auto &l : liftIndex){
	//		if(l.direction == REVERSE){
	//			uint32_t ass_len = id2Name[l.tid_ass].length;
	//			l.pos_ass_st = ass_len - l.pos_ass_st;
	//			l.pos_ass_ed = ass_len - l.pos_ass_ed;
	//			std::swap(l.pos_ass_st, l.pos_ass_ed);
	//		}
	//	}

		//store lift over data
		idx.dump(liftOverData_fn);
	}

	void liftoverSearchRegion(char *liftOverData_fn){
			//read bam/cram file:
		LIFTOVER_index idx;
		idx.load(liftOverData_fn);

		std::vector<ASS_INFO> &id2Name = idx.id2Name;
		std::vector<LIFTOVER> &liftIndex = idx.liftIndex;

		uint64_t ass_total_len = 0;
		for(auto & a: id2Name){	a.print(stderr, '\n');	ass_total_len += a.length;}
		fprintf(stderr, "[total lenghth: %ld]\n", ass_total_len);
		std::sort(liftIndex.begin(), liftIndex.end(), LIFTOVER::cmp_by_ref_pos);
		//for(auto & a: liftIndex){	a.print(stderr, id2Name);}
		std::vector<RefRegion> r_list;

		r_list.emplace_back(1, 1013668, 1014819);
		r_list.emplace_back(7, 699616, 700715);
		r_list.emplace_back(12, 129674126, 129675212);
		r_list.emplace_back(20, 62757546, 62759452);
		r_list.emplace_back(2, 91600041, 91601285);
		r_list.emplace_back(7, 103549585, 103550702);
		r_list.emplace_back(19, 8772442, 8773733);
		r_list.emplace_back(7, 104895073, 104896074);
		r_list.emplace_back(21, 11141346, 11142347);
		r_list.emplace_back(10, 71752046, 71753047);
		r_list.emplace_back(18, 77712216, 77713220);
		r_list.emplace_back(10, 42673500, 42674501);
		for(auto & r : r_list){
			int chrID_ref = r.chr_ID - 1;
			int st_pos_ref = r.st_pos;
			int ed_pos_ref = r.ed_pos;

			auto l_bg = liftIndex.begin();
			auto l_ed = liftIndex.end();
			auto l_c_bg = l_bg;
			for(; l_c_bg < l_ed; l_c_bg++){
				if(l_c_bg->tid_hs37 != chrID_ref)			continue;
				if(l_c_bg->pos_hs37_st < st_pos_ref)			continue;
				else										break;
			}
			auto l_c_ed = l_c_bg;
			for(; l_c_ed < l_ed; l_c_ed++){
				if(l_c_ed->tid_hs37 != chrID_ref)				break;
				if(l_c_ed->pos_hs37_st < ed_pos_ref)				continue;
				else											break;
			}

			LIFTOVER new_l_bg = l_c_bg[-1];

			fprintf(stderr, "Region in hs37d5:\t");
			r.print(stderr);

			int offset_bg = st_pos_ref - new_l_bg.pos_hs37_st;
			new_l_bg.pos_hs37_st += offset_bg;
			new_l_bg.pos_ass_st += offset_bg;

			LIFTOVER new_l_ed = l_c_ed[-1];
			int offset_ed = ed_pos_ref - new_l_ed.pos_hs37_st;
			new_l_ed.pos_hs37_st += offset_ed;
			new_l_ed.pos_ass_st += offset_ed;

			if(new_l_bg.direction == REVERSE){
				//fprintf(stderr, "\t\t before reverse BG:\t");
				//new_l_bg.print(stderr, id2Name);
				uint32_t length = id2Name[new_l_bg.tid_ass].length;
				new_l_bg.pos_ass_st = length - new_l_bg.pos_ass_st;
			}

			if(new_l_ed.direction == REVERSE){
				//fprintf(stderr, "\t\t before reverse ED:\t");
				//new_l_bg.print(stderr, id2Name);
				//new_l_ed.print(stderr, id2Name);
				std::swap(new_l_bg, new_l_ed);
				uint32_t length = id2Name[new_l_ed.tid_ass].length;
				new_l_ed.pos_ass_st = length - new_l_ed.pos_ass_st;
			}

			if(new_l_bg.direction == REVERSE && new_l_ed.direction == REVERSE){		std::swap(new_l_bg, new_l_ed);	}

			fprintf(stderr, "\t\t region liftover begin:\t");
			new_l_bg.print(stderr, id2Name);
			fprintf(stderr, "\t\t region liftover end:\t");
			new_l_ed.print(stderr, id2Name);
			fprintf(stderr, "\t\t final region: [%s:%d-%d]\n", id2Name[new_l_bg.tid_ass].name, new_l_bg.pos_ass_st, new_l_ed.pos_ass_st);
		}
	}

	void liftoverRead(char *liftOverData_fn, char * liftoverBam_fn){
		//load lift over file

		LIFTOVER_index idx;
		idx.load(liftOverData_fn);
		std::vector<ASS_INFO> &id2Name = idx.id2Name;
		std::vector<LIFTOVER> &liftIndex = idx.liftIndex;
		std::vector<liftoverRegion_direction> &LD = idx.LD;

		//lift over file check/sort and building index
		uint64_t ass_total_len = 0;
		for(auto & a: id2Name){	a.print(stderr, '\n');	ass_total_len += a.length;}
		fprintf(stderr, "[total lenghth: %ld]\n", ass_total_len);
		std::sort(liftIndex.begin(), liftIndex.end(), LIFTOVER::cmp_by_ass_pos); //sort by assembly position
		//for(auto & a: liftIndex){	a.print(stderr, id2Name);}
		//build simple index for the lift index
		std::vector<uint32_t> ass_ID_2_st_pos_in_liftover;
		ass_ID_2_st_pos_in_liftover.emplace_back(0);
		uint32_t old_ass_ID = MAX_uint32_t;
		//may be duplication in some position
		for(uint32_t i = 0; i < liftIndex.size(); i++){
			auto & a = liftIndex[i];
			if(old_ass_ID != a.tid_ass){
				old_ass_ID = a.tid_ass;
				uint32_t old_index = ass_ID_2_st_pos_in_liftover.back();
				ass_ID_2_st_pos_in_liftover.resize(old_ass_ID + 1, old_index);
				ass_ID_2_st_pos_in_liftover[ass_ID_2_st_pos_in_liftover.size() - 1] = i;
			}
		}
		ass_ID_2_st_pos_in_liftover.emplace_back(liftIndex.size());

		//pos in read ----> pos in assmbly ----- > pos in hs37d5
		//read bam/cram file:
		//read bam/cram file:
		Bam_file c_b;
		bam_file_open(liftoverBam_fn, NULL, NULL, &c_b);
		bam_hdr_t* hdr = c_b._hdr;
		uint64_t readNum = 0;
		bam1_t b; int sam_read_rst = 0;

		uint32_t cigar_buff[1000];
		uint32_t cigar_len = 0;

		while (1){
			readNum++; // debug code //if(readNum == 10) break;
			sam_read_rst = sam_read1(c_b._hfp, hdr, &b);
			if(sam_read_rst < 0) break;

			int st_pos_read_in_ass = b.core.pos;
			int ass_tid_in_read = b.core.tid;
			//get end pos:
			uint32_t* bam_cigar = bam_get_cigar(&b);
			int ed_pos_read_in_ass = st_pos_read_in_ass + bam_cigar2rlen( b.core.n_cigar, bam_cigar);

			//search direction in LD list
			bool direction_in_ass = true;
			for(liftoverRegion_direction &l: LD ){
				if(l.tid_ass != ass_tid_in_read)				continue;
				if(l.pos_ass_st < st_pos_read_in_ass)			continue;
				else{direction_in_ass = l.direction; break;}
			}

			//load cigar
			//reverse CIGAR when needed
			cigar_len = b.core.n_cigar;
			bam_cigar = bam_get_cigar(&b);

			//copy the cigar
			if(direction_in_ass == FORWARD)
				for (uint32_t i = 0; i < cigar_len; ++i)
					cigar_buff[i] = bam_cigar[i];
			else{
				for (uint32_t i = 0; i < cigar_len; ++i)
					cigar_buff[i] = bam_cigar[cigar_len - i - 1];
				uint32_t ass_len = id2Name[ass_tid_in_read].length;
				st_pos_read_in_ass = ass_len - st_pos_read_in_ass;
				ed_pos_read_in_ass = ass_len - ed_pos_read_in_ass;
				std::swap(st_pos_read_in_ass, ed_pos_read_in_ass);
			}

			//get beginning lift over block
			//search the lift over index
			auto l_bg = liftIndex.begin() + ass_ID_2_st_pos_in_liftover[ass_tid_in_read];
			auto l_ed = liftIndex.end() + ass_ID_2_st_pos_in_liftover[ass_tid_in_read + 1];
			//for start position
			auto l_c_bg = l_bg;
			for(; l_c_bg < l_ed; l_c_bg++)
				if(l_c_bg->direction != direction_in_ass || l_c_bg->pos_ass_st < st_pos_read_in_ass)		continue;
				else										break;
		}
	}
};

struct COMPLEX_SV_region{
private:
	struct SV_basic_infomation{
		int SV_CHR_ID;
		int pos_bgn;
		int pos_end;
		int is_CPX;
		SV_basic_infomation(int SV_CHR_ID, int pos_begin, int end_begin){
			this->SV_CHR_ID = SV_CHR_ID;
			this->pos_bgn = pos_begin;
			this->pos_end = end_begin;
			this->is_CPX = false;
		}
	};

	int CPX_distance;
	int min_SV_LEN;
	std::vector<std::vector<SV_basic_infomation>> HAP1;
	std::vector<std::vector<SV_basic_infomation>> HAP2;
	bcf_hdr_t *header;

	bool near(SV_basic_infomation & a, SV_basic_infomation & b){
		if(		ABS(a.pos_bgn - b.pos_bgn) <= CPX_distance ||
				ABS(a.pos_end - b.pos_end) <= CPX_distance ||
				ABS(a.pos_bgn - b.pos_end) <= CPX_distance ||
				ABS(a.pos_end - b.pos_bgn) <= CPX_distance ){
			return true;
		}
		return false;
	}

	void CPX_SV_detection(std::vector<std::vector<SV_basic_infomation>> & HAP1, const char * INFO){
		for(int chrID = 0; chrID <= 23; chrID++){
			std::vector<SV_basic_infomation> & h1 = HAP1[chrID];
			//load HAP1
			const char *chrID_str = bcf_hdr_id2name(header, chrID);
			for(uint i = 0; i < h1.size(); i++){
				for(uint j = i+1; j < h1.size(); j++){
					if(near(h1[i], h1[j])){
						h1[i].is_CPX = true;
						h1[j].is_CPX = true;
					}
				}
			}
			for(uint i = 0; i < h1.size(); i++){
				if(h1[i].is_CPX){
					int bg_POS = MAX(0, h1[i].pos_bgn - 1);
					int ed_POS = h1[i].pos_end + 1;
					fprintf(stdout, "%s\t%d\t%d\t%s\tTRUE_POS\n", chrID_str, bg_POS , ed_POS, INFO);
					{
						int bg_500W = (bg_POS/500) *500;
						int ed_500W = (ed_POS/500) *500;
						for(int i_500W = bg_500W; i_500W <= ed_500W; i_500W+=500){
							fprintf(stdout, "%s\t%d\t%d\t%s\t500W\n", chrID_str, i_500W , i_500W + 500, INFO);
						}
					}
				}
			}
		}
	}
public:
	int complex_sv_region_detection(int argc, char *argv[]){
		char * vcf_fn_in = argv[1];//separate by ','
		CPX_distance = atoi(argv[2]);
		min_SV_LEN = atoi(argv[3]);
		//load vcf files
		//open vcf files
		BCF_FILE vcf_read;//vcf for read
		VCF_open_read(&vcf_read, vcf_fn_in);//open for read
		char *c_sv_type = (char *)malloc(1000);
		int SV_CHR_ID;
		int SV_POS;
		int SV_length;
		int SV_END;
		char *GT[1]; char GT_1[3]; GT[0] = GT_1;

		HAP1.resize(24);
		HAP2.resize(24);
		header = vcf_read.header;

		do{
			bcf1_t *c_r = &( vcf_read.r);
			if(c_r->d.flt != NULL)
				*c_r->d.flt = 0;
			//unpack the vcf data to get the Filter
			bcf_unpack(c_r, BCF_UN_INFO);
			//vcf filters
			if(c_r->d.flt != NULL && *c_r->d.flt != 0)//filter: PASS
				continue;

			SV_CHR_ID = c_r->rid;
			SV_POS = c_r->pos;
			//vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
			vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
			vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

			//if(SV_POS < 893700) continue;
			//skip * SVs:
			if(c_r->d.allele[1][0] == '*')
				continue;

			if(SV_length < min_SV_LEN)
				continue;
			if(strcmp(c_sv_type, "INS") == 0)
				SV_END = SV_POS + 1;
			else
				SV_END = SV_POS + SV_length;

			vcf_get_sv_GT(header, c_r, GT);
			//store data:
			if(GT[0][0] == 4)
				HAP1[SV_CHR_ID].emplace_back(SV_CHR_ID, SV_POS, SV_END);
			if(GT[0][1] == 5)
				HAP2[SV_CHR_ID].emplace_back(SV_CHR_ID, SV_POS, SV_END);
		}
		while(VCF_next(&vcf_read));//read one
		bcf_close(vcf_read.file);
		CPX_SV_detection(HAP1, "HAP1");
		CPX_SV_detection(HAP2, "HAP2");
		return 0;
	}

	int complex_sv_region_detection_500W(int argc, char *argv[]){
		char * vcf_fn_in = argv[1];//separate by ','
		CPX_distance = atoi(argv[2]);
		min_SV_LEN = atoi(argv[3]);
		//load vcf files
		//open vcf files
		BCF_FILE vcf_read;//vcf for read
		VCF_open_read(&vcf_read, vcf_fn_in);//open for read
		char *c_sv_type = (char *)malloc(1000);
		int SV_CHR_ID;
		int SV_POS;
		int SV_length;
		int SV_END;
		char *GT[1]; char GT_1[3]; GT[0] = GT_1;

		HAP1.resize(24);
		HAP2.resize(24);
		header = vcf_read.header;

		do{
			bcf1_t *c_r = &( vcf_read.r);
			if(c_r->d.flt != NULL)
				*c_r->d.flt = 0;
			//unpack the vcf data to get the Filter
			bcf_unpack(c_r, BCF_UN_INFO);
			//vcf filters
			if(c_r->d.flt != NULL && *c_r->d.flt != 0)//filter: PASS
				continue;

			SV_CHR_ID = c_r->rid;
			SV_POS = c_r->pos;
			//vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
			vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
			vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

			//if(SV_POS < 893700) continue;
			//skip * SVs:
			if(c_r->d.allele[1][0] == '*')
				continue;

			if(SV_length < min_SV_LEN)
				continue;
			if(strcmp(c_sv_type, "INS") == 0)
				SV_END = SV_POS + 1;
			else
				SV_END = SV_POS + SV_length;

			vcf_get_sv_GT(header, c_r, GT);
			//store data:
			if(GT[0][0] == 4)
				HAP1[SV_CHR_ID].emplace_back(SV_CHR_ID, SV_POS, SV_END);
			if(GT[0][1] == 5)
				HAP2[SV_CHR_ID].emplace_back(SV_CHR_ID, SV_POS, SV_END);
		}
		while(VCF_next(&vcf_read));//read one
		bcf_close(vcf_read.file);
		CPX_SV_detection(HAP1, "HAP1");
		CPX_SV_detection(HAP2, "HAP2");
		return 0;
	}

	void get_region_string(std::vector<std::string> &item_value, std::vector<std::string> & tmp, std::string & region_str){
		region_str.clear();
		split_string(tmp, item_value[0].c_str() + 11, ":");
		region_str += tmp[0] + "_";
		std::string POS = tmp[1];
		split_string(tmp, POS.c_str(), "-");
		region_str += tmp[0] + "_" + tmp[1];
	}

	void get_region_for_loading_ref(std::vector<std::string> &item_value, std::vector<std::string> & tmp, std::string & region_str){
		region_str.clear();
		split_string(tmp, item_value[0].c_str() + 11, ":");
		region_str += tmp[0] + ":";
		std::string POS = tmp[1];
		split_string(tmp, POS.c_str(), "-");
		region_str += tmp[0] + "-" + tmp[1];
	}

	int contig_file_split(int argc, char *argv[]){
		const char * FC_Contig_L = argv[1];
		const char * work_dir = argv[2];
		const char * ref_fn = argv[3];

		faidx_t *fai = fai_load(ref_fn);//load reference

		//load all the FC_Contig_L
		std::vector<std::string> contig_FN_l;
		load_string_list_from_file(FC_Contig_L, contig_FN_l);
		std::vector<std::string> contigs_in_f;
		std::vector<std::string> item_value;
		std::vector<std::string> item_value1;
		std::unordered_map<std::string ,  FILE*> region_index;
		for(uint sample_idx = 0; sample_idx < contig_FN_l.size(); sample_idx++ ){
			fprintf(stderr, "Current handle file is: %s\n", contig_FN_l[sample_idx].c_str());
			split_string(item_value, contig_FN_l[sample_idx].c_str(), "/");
			std::string sample_name = item_value[item_value.size() -2];
			contigs_in_f.clear();
			load_string_list_from_file(contig_FN_l[sample_idx].c_str(), contigs_in_f);
			for(uint c_idx = 0; c_idx < contigs_in_f.size(); c_idx++ ){
				//for each contig:
				split_string(item_value, contigs_in_f[c_idx].c_str(), "\t");
				if(item_value.size() < 5)
					continue;
				std::string region;
				get_region_string(item_value, item_value1, region);
				//store contig strings
				split_string(item_value1, item_value[4].c_str(), "=");
				std::string contig_str = item_value1[1];
				//store the contig_ID
				split_string(item_value1, item_value[2].c_str(), ":");
				std::string contig_ID = item_value1[9];
				std::string store_string = ">" + sample_name + "_" + region + "_" + contig_ID + "\n" + contig_str + "\n";
				//store:
				std::unordered_map<std::string , FILE*>::iterator it = region_index.find(region);
				if(it == region_index.end()){
					//open a new file:
					char region_fn[1024];
					sprintf(region_fn, "%s/%s", work_dir, region.c_str());
					FILE* try_open = xopen(region_fn, "w");
					region_index[region] = try_open;
					//store data
					//store the reference
					int true_region_load_len = 0;
					std::string ref_region;
					get_region_for_loading_ref(item_value, item_value1, ref_region);
					char *ref_seq_char = fai_fetch(fai, ref_region.c_str(), &true_region_load_len);
					std::string ref_seq_char_string(ref_seq_char);
					if(ref_seq_char != NULL) free(ref_seq_char);
					std::string store_string_ref = ">ref" + region + "\n" + ref_seq_char_string + "\n";
					fwrite(&(store_string_ref[0]), 1, store_string_ref.size(), try_open);
					fwrite(&(store_string[0]), 1, store_string.size(), try_open);
				}
				else{
					//store data
					fwrite(&(store_string[0]), 1, store_string.size(), it->second);
				}
			}
		}
		//close all files
		std::unordered_map<std::string , FILE*>::iterator it = region_index.begin();
		for(;it!=region_index.end();it++){
			fclose(it->second);
		}
		return 0;
	}

	int vcf_2_contig(int argc, char *argv[]){
		char *ref_fn_in = argv[1];//ref
		char *vcf_fn_in = argv[2]; //vcf
		char *ref_region = argv[3]; //region

		int DEL_number = 0;
		int INS_number = 0;
		//load reference string;
		faidx_t *fai = fai_load(ref_fn_in);//load reference
		int true_region_load_len = 0;
		char *ref_seq_char = fai_fetch(fai, ref_region, &true_region_load_len);
		std::vector<std::string> sp1;
		std::vector<std::string> sp2;
		std::string ref_seq = ref_seq_char;
		char *c_sv_type = (char *)malloc(1000);
		split_string(sp1, ref_region, ":");
		split_string(sp2, sp1[1].c_str(), "-");
		int region_bgn = atoi(sp2[0].c_str()) - 1;
		//int region_end = atoi(sp2[1].c_str());
		//todo::
		fprintf(stderr, "%s\n", ref_seq_char);
		fai_destroy(fai);

		//load vcf files
		BCF_FILE vcf_read;//vcf for read
		VCF_open_read(&vcf_read, vcf_fn_in);//open for read
		int SV_POS;

		header = vcf_read.header;
		int indel_offset = 0;
		do{
		 bcf1_t *c_r = &( vcf_read.r);
		 if(c_r->d.flt != NULL)
			 *c_r->d.flt = 0;
		 //unpack the vcf data to get the Filter
		 bcf_unpack(c_r, BCF_UN_INFO);
		 //vcf filters
		 if(c_r->d.flt != NULL && *c_r->d.flt != 0)//filter: PASS
			 continue;
		 if(c_r->d.allele == NULL)
			 continue;
		 //SV_CHR_ID = c_r->rid;
		 SV_POS = c_r->pos;
		 //skip * SVs:
		 if(c_r->d.allele[1][0] == '*')
			 continue;
		 std::string REF_s = c_r->d.allele[0];
		 std::string ALT_s = c_r->d.allele[1];
		 int sv_type_tag_n = vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);
		 if(sv_type_tag_n > 0){
			 //count the region number and SV number
			 if(REF_s.size() > ALT_s.size()) DEL_number ++;
			 else							 INS_number ++;
			 //todo:: other tags:
//			 /RM_clsfam
			 //REMAP
			 //RM_repeat
			 {


			 }
		 }
		 //modify;
		 int sv_pos_local = SV_POS - region_bgn + indel_offset;
		 if(sv_pos_local < 0)
			 continue;
		 //check:
		 if(REF_s[0] != ref_seq[sv_pos_local]){
			 fprintf(stderr, "string ERROR at SV position %d\n", SV_POS);
			 continue;
		 }
		 indel_offset += ALT_s.size() - REF_s.size();
		 //modify
		 if(sv_pos_local + REF_s.size()    > ref_seq.size() ||  sv_pos_local < 0) continue;
		 ref_seq.erase(sv_pos_local, REF_s.size());
		 if(sv_pos_local >= (int)ref_seq.size() ||  sv_pos_local < 0) continue;
		 ref_seq.insert(ref_seq.begin() + sv_pos_local, ALT_s.begin(), ALT_s.end());
		}
		while(VCF_next(&vcf_read));//read one
		bcf_close(vcf_read.file);

		fprintf(stdout, "%s %d %d \n", ref_seq.c_str(), INS_number, DEL_number);
		//output the final results
		return 0;
	}
};


#endif /* SVCALLING_CORE_ANALYSIS_HPP_ */
