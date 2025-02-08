/*
 * Contig_aligner.hpp
 *
 *  Created on: 2024年2月4日
 *      Author: fenghe
 */

#ifndef SV_SPA_CORE_CONTIG_ALIGNER_HPP_
#define SV_SPA_CORE_CONTIG_ALIGNER_HPP_

#include<vector>
#include<algorithm>
#include "../clib/bam_file.h"
#include "../SVcalling_core/NovaSVRst.hpp"
extern "C"{
	#include "../clib/utils.h"
	extern int vcf_write_line(htsFile *fp, kstring_t *line);
	#include "../kswlib/kalloc.h"
	#include "../kswlib/ksw2.h"
	#include "../clib/desc.h"
	#include "../clib/vcf_lib.h"
}

struct SUM_variations{
	int ref_idx;
	int read_idx;

	int total_var_number;
	int total_var_size;
	char ref_alt_char;//

	int cigar_index_bg;
	int cigar_index_ed;

	void clear(){
		memset(this, 0, sizeof(SUM_variations));
		ref_idx = -1;
	}
	void store_signal(int ref_idx_, int read_idx_, int c_length, char ref_alt_char_, int cigar_idx){
		if(ref_idx == -1) {
			ref_idx = ref_idx_;
			read_idx = read_idx_;
			ref_alt_char = ref_alt_char_;
			cigar_index_bg = cigar_idx;
		}
		total_var_number ++;
		total_var_size += c_length;
		cigar_index_ed = cigar_idx;
	}

	bool minSV_filter(int minSVlen){
		if(total_var_number > 1 && ref_idx != -1 && ((total_var_size >= minSVlen) || (-total_var_size >= minSVlen)))
			return true;
		return false;
	}

};

void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);

struct SV_REGION_TGS{
	int ref_begin;
	int ref_end;
	int contig_begin;
	int contig_end;
	int cigar_idx;
	void show(){
		fprintf(stderr, "SV_REGION_TGS: ref_begin %d , ref_end %d , contig_begin %d , contig_end %d \n", ref_begin, ref_end, contig_begin, contig_end);
	}
};

struct Contig_String_aligner{

public:
	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
		setRandomHeadTail();
	}

	void setRef(uint8_t * ref_string, int ref_len, 	int chr_ID_, int global_ref_pos_)
	{
		tseq = ref_string;
		tlen = ref_len;
		chr_ID = chr_ID_;
		global_ref_pos = global_ref_pos_;
		xassert(tseq != NULL, "tseq not empty!");
	}

	void printf_ref(FILE * output){
		for(int i = 0; i < tlen; i++){
			fprintf(output, "%c", "ACGTNN"[tseq[i]]);
		}
		fprintf(output, "\n");
	}

	void get_ref_info(int & chr_id, int & pos, int & len){
		chr_id = chr_ID;
		pos = global_ref_pos;
		len = tlen;
	}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}



	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos, int additionEndLen){
		xassert(tseq != NULL, "tseq not empty1!");
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		ref_end_pos += additionEndLen;
		last_ref_st_pos = ref_st_pos;
		last_ref_end_pos = ref_end_pos;

		if(false){
			//debug code:
			uint32_t debug_qlen = qlen;
			uint8_t* debug_qseq = qseq;
			for(uint32_t i = 0; i < debug_qlen; i++)
				fprintf(stderr, "%c", "ACGTNNN"[ debug_qseq[i]]);
			fprintf(stderr, "\n");
			uint32_t debug_tlen = ref_end_pos - ref_st_pos;
			uint8_t* debug_tseq =  tseq + ref_st_pos;
			for(uint32_t i = 0; i < debug_tlen; i++)
				fprintf(stderr, "%c", "ACGTNNN"[ debug_tseq[i]]);
			fprintf(stderr, "\n");
		}
//fprintf(stderr, "S52 , km %ld, qlen  %ld, qseq %ld, ref_end_pos - ref_st_pos  %ld , tseq + ref_st_pos  %ld, ref_st_pos  %ld, 5  %ld, mata_D  %ld, gap_open_D  %ld, gap_ex_D  %ld, gap_open2_D  %ld, gap_ex2_D  %ld, bandwith  %ld, zdrop_D  %ld, -1  %ld, flag  %ld, &ez  %ld, \n",				   km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}
#define ALN_ADDITION_LOAD_MAX 12000
	void align_non_splice_super_long(std::vector<uint8_t> &bin_contig, std::vector<uint8_t> &bin_ref, int &addition_load){
		if(addition_load > ALN_ADDITION_LOAD_MAX){
			addition_load = ALN_ADDITION_LOAD_MAX;
		}

		combine_contig.clear();
		combine_contig.insert(combine_contig.end(), random_head.begin(), random_head.begin() + addition_load/2);
		combine_contig.insert(combine_contig.end(), bin_contig.begin(), bin_contig.end());
		combine_contig.insert(combine_contig.end(), random_tail.begin(), random_tail.begin() + addition_load/2);

		combine_ref.clear();
		combine_ref.insert(combine_ref.end(), random_head.begin(), random_head.begin() + addition_load/2);
		combine_ref.insert(combine_ref.end(), bin_ref.begin(), bin_ref.end());
		combine_ref.insert(combine_ref.end(), random_tail.begin(), random_tail.begin() + addition_load/2);

		qseq = &(combine_contig[0]);
		qlen = combine_contig.size();

		tseq = &(combine_ref[0]);
		tlen = combine_ref.size();

		last_ref_st_pos = 0;
		last_ref_end_pos = tlen;
		bandwith = 10000; zdrop_D = 10000;
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	bool align_non_splice_super_long_adjust_cigar(int addition_load){
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0)
			return false;
		int shift_size = addition_load/2;
		if((bam_cigar[0] >> BAM_CIGAR_SHIFT) >= shift_size ){
			bam_cigar[0] -= (shift_size << BAM_CIGAR_SHIFT);
		}
		else
			return false;
		if((bam_cigar[cigar_len - 1] >> BAM_CIGAR_SHIFT) >= shift_size ){
			bam_cigar[cigar_len - 1] -= (shift_size << BAM_CIGAR_SHIFT);
		}else
			return false;
		return true;
	}


	void align_non_splice_HARD_LOCAL(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		xassert(tseq != NULL, "tseq not empty1!");
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		last_ref_st_pos = ref_st_pos;
		last_ref_end_pos = ref_end_pos;

		if(0){
			//debug code:
			uint32_t debug_qlen = qlen;
			uint8_t* debug_qseq = qseq;
			for(uint32_t i = 0; i < debug_qlen; i++)
				fprintf(stderr, "%c", "ACGTNNN"[ debug_qseq[i]]);
			fprintf(stderr, "\n");
			uint32_t debug_tlen = ref_end_pos - ref_st_pos;
			uint8_t* debug_tseq =  tseq + ref_st_pos;
			for(uint32_t i = 0; i < debug_tlen; i++)
				fprintf(stderr, "%c", "ACGTNNN"[ debug_tseq[i]]);
			fprintf(stderr, "\n");
		}
//fprintf(stderr, "S52 , km %ld, qlen  %ld, qseq %ld, ref_end_pos - ref_st_pos  %ld , tseq + ref_st_pos  %ld, ref_st_pos  %ld, 5  %ld, mata_D  %ld, gap_open_D  %ld, gap_ex_D  %ld, gap_open2_D  %ld, gap_ex2_D  %ld, bandwith  %ld, zdrop_D  %ld, -1  %ld, flag  %ld, &ez  %ld, \n",				   km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	int adjustCIGAR(){
		uint32_t n_cigar = ez.n_cigar;
		int adj_size = cigar_adjust(&n_cigar, ez.cigar, false, 15);
		ez.n_cigar = n_cigar;
		return adj_size;
	}

	void printf_alignment_detail(FILE * output, int suggest_st_pos, uint16_t * contig_depth,
			int ref_pos_enough_match_base){
		log_output = output;
		printCIGAR(output);
		int cigar_len = ez.n_cigar;
		if(cigar_len == 0) return;
		if(true) printContigSeq(suggest_st_pos);//print contig sequence
		int contig_coverage = 0;
		if(true){ contig_coverage = print_X_E_sequence(suggest_st_pos, ref_pos_enough_match_base); } //print X/= sequence
		if(true && contig_depth != NULL){ print_coverage(suggest_st_pos, contig_coverage, contig_depth); } 		//print coverage sequence
		print_SV_canditate(suggest_st_pos);
		fprintf(output, "\n");
	}

	int select_suggest_sv_length(int SV_length, std::vector<int> &suggest_SV_length){
		int best_suggest = 0; int min_dis = 9999999;
		for(int sug: suggest_SV_length){
			float diff_ratio = (float)SV_length/sug - 1;
			diff_ratio = ABS(diff_ratio);
			int dis = ABS_U(SV_length, sug);

			if(diff_ratio > 0.12 && dis > 8)
				continue;

			if(dis < min_dis){
				min_dis = dis;
				best_suggest = sug;
			}
		}
		return best_suggest;
	}//

	int getMiddleMatchLen(int cigar_index_bg, int cigar_index_ed, uint32_t* bam_cigar){
		int matchBaseTotal = 0;
		for(int cigar_ID = cigar_index_bg; cigar_ID< cigar_index_ed; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			if(c_type == 0){//match
				matchBaseTotal += c_length;
			}
		}
		fprintf(stderr, "Middle length of SV is %5d,\n ", matchBaseTotal);
		return matchBaseTotal;
	}

	void get_canditate_SVs(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> &SVs,
			int minSVlen, int ref_st_pos, std::vector<int> &suggest_SV_length, int region_ref_global_position){
		uint64_t SVs_size_ori = SVs.size();
		if(print_log) fprintf(stderr, "get_canditate_SVs BG\n");
		uint32_t* bam_cigar = ez.cigar;
		int n_cigar = ez.n_cigar;
		int ref_ed_pos = last_ref_end_pos;
		int ref_index = ref_st_pos;
		int read_index = 0;
		uint8_t * ref = NULL;
		uint8_t * alt = NULL;
		int ABS_ref_distance_to_region_end;
		//for SUM deletions
		SUM_variations del_sum; del_sum.clear();
		SUM_variations ins_sum; ins_sum.clear();
		SUM_variations cpx_sum; cpx_sum.clear();
		//for SUM insertions
		int suggest_sv_len = 0;
		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:	ref_index += c_length; read_index += c_length; break;//M
			case 1:	//insertion
				//if(c_length < minSVlen) break;
				ABS_ref_distance_to_region_end = ABS_U(ref_index, ref_ed_pos);
				if(ABS_ref_distance_to_region_end < 20) break;//
				ref = tseq + ref_index; ref_buff.resize(2);//ref
				ref_buff[0] = "ACGTNNN"[ref[0]]; ref_buff[1] = 0;
				alt = qseq + read_index; alt_buff.resize(c_length + 1);//read
				for(int i = 0; i < c_length; i++) alt_buff[i] = "ACGTNNN"[alt[i]];
				alt_buff[c_length] = 0;
				alt_buff.insert(alt_buff.begin(), ref_buff[0]);
				//sum insertions
				if(cigar_ID!=0 && cigar_ID!=n_cigar-1){
					ins_sum.store_signal(ref_index + global_ref_pos, read_index, c_length, "ACGTNNN"[qseq[read_index]], cigar_ID);
					cpx_sum.store_signal(ref_index + global_ref_pos, read_index, c_length, "ACGTNNN"[qseq[read_index]], cigar_ID);
					if(print_log) fprintf(stderr, "ins_sum ADD %d\n", c_length);
				}
				//suggest SV length:
				suggest_sv_len = select_suggest_sv_length(c_length, suggest_SV_length);
				if(print_log) fprintf(stderr, "INS:suggest_sv_len %d\n", suggest_sv_len);
				if(suggest_sv_len != 0 && c_length >= minSVlen && cigar_ID!=0 && cigar_ID!=n_cigar-1)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, ref_index + global_ref_pos, "INS", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, cigar_ID, cigar_ID, ref_st_pos, region_ref_global_position);
				read_index += c_length;
				break;//I, int chr_ID; int
			case 2:	 //Deletions
				//if(c_length < minSVlen) break;
				ref = tseq + ref_index; ref_buff.resize(c_length + 1);//ref
				for(int i = 0; i < c_length; i++) ref_buff[i] = "ACGTNNN"[ref[i]];
				ref_buff[c_length] = 0;
				alt = qseq + read_index; alt_buff.resize(2);//read
				alt_buff[0] = "ACGT"[alt[0]]; alt_buff[1] = 0;
				ref_buff.insert(ref_buff.begin(), alt_buff[0]);
				//sum deletion:
				if(cigar_ID!=0 && cigar_ID!=n_cigar-1){
					del_sum.store_signal(ref_index + global_ref_pos, read_index, c_length, "ACGTNNN"[qseq[read_index]], cigar_ID);
					cpx_sum.store_signal(ref_index + global_ref_pos, read_index, -c_length, "ACGTNNN"[qseq[read_index]], cigar_ID);
					if(print_log) fprintf(stderr, "del_sum ADD %d\n", c_length);
				}
				suggest_sv_len = select_suggest_sv_length(-c_length, suggest_SV_length);
				if(print_log) fprintf(stderr, "DEL:suggest_sv_len %d\n", suggest_sv_len);
				if(suggest_sv_len != 0 && c_length >= minSVlen && cigar_ID!=0 && cigar_ID!=n_cigar-1)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, ref_index + global_ref_pos, "DEL", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len,  qseq, qlen, bam_cigar, n_cigar, cigar_ID, cigar_ID, ref_st_pos, region_ref_global_position);
				ref_index += c_length;
				break;
			case 3:	ref_index += c_length; read_index += c_length; break;//N, print N
			case 4:	ref_index += c_length; break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}

		{
			//sum deletions:
			SUM_variations &c_sum = del_sum;
			if(c_sum.minSV_filter(minSVlen)){
				ref = tseq + c_sum.ref_idx - global_ref_pos; ref_buff.resize(c_sum.total_var_size + 1);//ref
				for(int i = 0; i < c_sum.total_var_size; i++) ref_buff[i] = "ACGTNNN"[ref[i]];
				ref_buff[c_sum.total_var_size] = 0;
				alt_buff.resize(2);	alt_buff[0] = 'N'; alt_buff[1] = 0; //read

				int suggest_sv_len = select_suggest_sv_length(-c_sum.total_var_size, suggest_SV_length);
				if(print_log) fprintf(stderr, "DEL:SUM suggest_sv_len %d SV_len %d\n", suggest_sv_len, -c_sum.total_var_size);

				int matchBaseTotal = getMiddleMatchLen(c_sum.cigar_index_bg, c_sum.cigar_index_ed, bam_cigar);
				if(suggest_sv_len != 0 && matchBaseTotal < 30)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "DEL", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
			}
		}

		{
			//sum insertion:
			SUM_variations &c_sum = ins_sum;
			if(c_sum.minSV_filter(minSVlen)){
				ref_buff.resize(2); ref_buff[0] = 'N'; ref_buff[1] = 0;
				alt = qseq + c_sum.read_idx; alt_buff.resize(c_sum.total_var_size + 1);//ref
				for(int i = 0; i < c_sum.total_var_size; i++) alt_buff[i] = "ACGTNNN"[alt[i]];
				alt_buff[c_sum.total_var_size] = 0;

				int matchBaseTotal = getMiddleMatchLen(c_sum.cigar_index_bg, c_sum.cigar_index_ed, bam_cigar);

				int suggest_sv_len = select_suggest_sv_length(c_sum.total_var_size, suggest_SV_length);
				if(print_log) fprintf(stderr, "INS:SUM suggest_sv_len %d SV_len %d\n", suggest_sv_len, c_sum.total_var_size);
				if(suggest_sv_len != 0 && matchBaseTotal < 30)
					NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "INS", &(ref_buff[0]), &(alt_buff[0]),
							suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
			}
		}

		//only consider CPX when no other SVs are called
		if(SVs_size_ori == SVs.size()){
			//sum complex vars:
			SUM_variations &c_sum = cpx_sum;
			if(c_sum.minSV_filter(minSVlen) && c_sum.total_var_number <= 4){
				//{std::string A = "<NA_CPX_REF>";A.push_back(0); ref_buff.resize(A.size()); memcpy(&(ref_buff[0]), &(A[0]), A.size());}
				//{std::string A = "<NA_CPX_ALT>";A.push_back(0); alt_buff.resize(A.size()); memcpy(&(alt_buff[0]), &(A[0]), A.size());}
				if(c_sum.total_var_size < 0){

					ref = tseq + c_sum.ref_idx - global_ref_pos; ref_buff.resize(-c_sum.total_var_size + 2);//ref
					for(int i = 0; i < -c_sum.total_var_size + 1; i++) ref_buff[i] = 'N';
					ref_buff[-c_sum.total_var_size + 1] = 0;
					alt_buff.resize(2);	alt_buff[0] = 'N'; alt_buff[1] = 0; //read

					int matchBaseTotal = getMiddleMatchLen(c_sum.cigar_index_bg, c_sum.cigar_index_ed, bam_cigar);

					int suggest_sv_len = select_suggest_sv_length(c_sum.total_var_size, suggest_SV_length);
					if(print_log) fprintf(stderr, "CPX_DEL:: suggest_sv_len %d SV_len %d\n", suggest_sv_len, c_sum.total_var_size);
					if(suggest_sv_len != 0 && matchBaseTotal < 30)
						NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "DEL", &(ref_buff[0]), &(alt_buff[0]),
								suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
				}
				else if(c_sum.total_var_size > 0){

					ref_buff.resize(2); ref_buff[0] = 'N'; ref_buff[1] = 0;
					alt_buff.resize(c_sum.total_var_size + 2);//ref
					for(int i = 0; i < c_sum.total_var_size + 1; i++) alt_buff[i] = 'N';
					alt_buff[c_sum.total_var_size + 1] = 0;
					int matchBaseTotal = getMiddleMatchLen(c_sum.cigar_index_bg, c_sum.cigar_index_ed, bam_cigar);
					int suggest_sv_len = select_suggest_sv_length(c_sum.total_var_size, suggest_SV_length);
					if(print_log) fprintf(stderr, "CPX_INS:: suggest_sv_len %d SV_len %d\n", suggest_sv_len, c_sum.total_var_size);
					if(suggest_sv_len != 0 && matchBaseTotal < 30)
						NOVA_SV_FINAL_RST_item::add_to_vector(SVs, chr_ID, c_sum.ref_idx, "INS", &(ref_buff[0]), &(alt_buff[0]),
								suggest_sv_len, qseq, qlen, bam_cigar, n_cigar, c_sum.cigar_index_bg, c_sum.cigar_index_ed, ref_st_pos, region_ref_global_position);
				}
			}
		}
	}

	int get_canditate_SVs_TGS(bool print_log, int minSVlen, int ref_st_pos, std::vector<SV_REGION_TGS> &sv_r, bool output_small_var,
			std::string &full_ref_str, std::string &full_contig){
		int truSV_num = 0;
		if(print_log) fprintf(stderr, "get_canditate_SVs BG\n");
		uint32_t* bam_cigar = ez.cigar;
		int n_cigar = ez.n_cigar;
		int ref_ed_pos = last_ref_end_pos;
		int ref_index = ref_st_pos;
		int contig_index = 0;
		int ABS_ref_distance_to_region_end;
		const char * ref_p = full_ref_str.c_str();
		const char * contig_p = full_contig.c_str();
		//for SUM insertions
		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int c_length =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:
			if(output_small_var){
				for(int i = 0; i < c_length; i++){
					if(ref_p[ref_index + i] != contig_p[contig_index + i]){
						sv_r.emplace_back();
						sv_r.back().ref_begin = ref_index + i; sv_r.back().ref_end = ref_index + i;
						sv_r.back().contig_begin = contig_index + i; sv_r.back().contig_end = contig_index + i;
						sv_r.back().cigar_idx = cigar_ID;
					}
				}
			}
			ref_index += c_length; contig_index += c_length; break;//M
			case 1:	//insertion
				ABS_ref_distance_to_region_end = ABS_U(ref_index, ref_ed_pos);
				if(ABS_ref_distance_to_region_end < 20) break;//
				//suggest SV length:
				if((c_length >= minSVlen || output_small_var) && cigar_ID!=0 && cigar_ID!=n_cigar-1){
					if(c_length >= minSVlen) truSV_num++;
					sv_r.emplace_back();
					sv_r.back().ref_begin = ref_index; sv_r.back().ref_end = ref_index;
					sv_r.back().contig_begin = contig_index; sv_r.back().contig_end = contig_index + c_length;
					sv_r.back().cigar_idx = cigar_ID;
				}
				contig_index += c_length;
				break;//I, int chr_ID; int
			case 2:	 //Deletions
				if((c_length >= minSVlen || output_small_var) && cigar_ID!=0 && cigar_ID!=n_cigar-1){
					if(c_length >= minSVlen) truSV_num++;
					sv_r.emplace_back();
					sv_r.back().ref_begin = ref_index; sv_r.back().ref_end = ref_index + c_length;
					sv_r.back().contig_begin = contig_index; sv_r.back().contig_end = contig_index;
					sv_r.back().cigar_idx = cigar_ID;
				}
				ref_index += c_length;
				break;
			case 3:	ref_index += c_length; contig_index += c_length; break;//N, print N
			case 4:	ref_index += c_length; break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}
		return truSV_num;
	}

	inline uint8_t getTseq(int i){ return tseq[i];}

	uint8_t* get_ref(){return tseq;}
	int get_tlen(){return tlen;}
	void set_tseq(uint8_t *tseq_){this->tseq = tseq_;}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

	int get_cigar(uint32_t** bam_cigar){
		*bam_cigar = ez.cigar;
		return ez.n_cigar;
	}

	//part7: candidate SVs
	int last_ref_st_pos;  //the reference start and end position of last alignment
	int last_ref_end_pos;

	void printCIGAR(FILE * log_output){
		fprintf(log_output, "\n");
		int cigar_len = ez.n_cigar;
		uint32_t* bam_cigar = ez.cigar;
		if(cigar_len == 0){
			for(int i = 0; i < (int)qlen; i++) {fprintf(log_output, "%c", "ACGTNNN"[qseq[i]]);}
			fprintf(log_output, "\n ???NO alignment result\n");
		}else{
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int type = 1 + (bam_cigar[cigar_ID] & BAM_CIGAR_MASK);
				char c_type = segment_type_to_cigar_code(type);
				int length = (bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT);
				fprintf(log_output,"%d%c",length,  c_type);
			}
			fprintf(log_output,"\n");
		}
	}

	//part1 : basic informations
	int chr_ID;
	int global_ref_pos;
private:

	//part2: reference and target
	uint8_t *tseq; int tlen;
	uint8_t *qseq; uint32_t qlen; //query temp: used for analysis
	//part3: alignment result
	ksw_extz_t ez;
	//part4: buffs
	void *km;
	//part5: options
	int8_t mata_D[25]; int8_t match_D; int8_t mismatch_D;
	int8_t gap_open_D; int8_t gap_ex_D; int8_t gap_open2_D; int8_t gap_ex2_D;
	int zdrop_D; int bandwith; //for DNA zdrop = 400, 200 for RNA
	int flag;
	//part6: logs
	FILE * log_output;
	std::vector<char> ref_buff;
	std::vector<char> alt_buff;

	//for long SV comparing
	//part2: reference and target
	std::vector<uint8_t> random_head;
	std::vector<uint8_t> random_tail;

	std::vector<uint8_t> combine_ref;
	std::vector<uint8_t> combine_contig;

	void setRandomHeadTail(){
		//part3: alignment result
		if(false)fprintf(stderr, "CA: random head string generate\n");
		uint64_t random_seed = 10086;
		for(int i = 0; i < ALN_ADDITION_LOAD_MAX/2 + 100; i++){
			simple_rand(&random_seed);
			uint8_t c = random_seed % 4;
			random_head.emplace_back(c);
			if(false)fprintf(stderr, "%c", "ACGTN"[c]);
		}
		if(false)fprintf(stderr, "\n");
		for(int i = 0; i < ALN_ADDITION_LOAD_MAX/2 + 100; i++){
			simple_rand(&random_seed);
			uint8_t c = random_seed % 4;
			random_tail.emplace_back(c);
			if(false)fprintf(stderr, "%c", "ACGTN"[c]);
		}
		if(false)fprintf(stderr, "\n");
	}

	void copy_option(){
//M1
			match_D = 2;
			mismatch_D= 10;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 0;
			zdrop_D= gap_open2_D + 600; //for DNA zdrop = 400, 200 for RNA
			bandwith = 600;//zdrop_D;
			flag = 0;
//M2:
			if(true){
				match_D = 1;
				mismatch_D= 15;
				gap_open_D= 15;
				gap_ex_D= 3;
				gap_open2_D= 50;
				gap_ex2_D= 1;
				zdrop_D= gap_open2_D + 600; //for DNA zdrop = 400, 200 for RNA
				bandwith = 600;//zdrop_D;
				flag = 0;
			}

			//
			fprintf(stderr,
					"SET ALN PARA: "
					"match_D %d, mismatch_D %d, gap_open_D %d, gap_ex_D %d,"
					" gap_open2_D %d, gap_ex2_D %d, zdrop_D %d,"
					" bandwith %d, flag %d\n",
					match_D, mismatch_D, gap_open_D, gap_ex_D,
					gap_open2_D,	gap_ex2_D, zdrop_D, bandwith, flag);
	}

	void ksw_gen_mat_D(){
		int8_t l,k,m;
		for (l = k = 0; l < 4; ++l) {
			for (m = 0; m < 4; ++m) { mata_D[k] = l == m ? match_D : -(mismatch_D);	/* weight_match : -weight_mismatch */ k++; }
			mata_D[k] = 0; // ambiguous base
			k++;
		}
		for (m = 0; m < 5; ++m) { mata_D[k] = 0; k++; }
	}

	void printContigSeq(int suggest_st_pos){
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		int output_index = 0;
		int seq_i = 0;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", "ACGTNNN"[qseq[seq_i]]); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
	}

	int print_X_E_sequence(int suggest_st_pos, int ref_pos_enough_match_base){
		int contig_coverage = 0;
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++)
				{
					if(qseq[seq_i] == tseq[output_index]){
						contig_coverage++;
						if(output_index < ref_pos_enough_match_base) fprintf(log_output, "M");//1 == M; 2 == X; 0 == -;
						else fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
					}
					else{
						fprintf(log_output, "X");//1 == M; 2 == X; 0 == -;
					}
					output_index++;
				}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\n");
		return contig_coverage;
	}

public:
	int print_X_E(FILE * output, int suggest_st_pos){
		log_output = output;
		return print_X_E_sequence(suggest_st_pos, 0);
	}

	void show_score(){
		fprintf(stderr, "Score is %d\n", ez.score);
	}

	int check_match_base(int suggest_st_pos){
		int contig_coverage = 0;
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) { output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++){
					if(qseq[seq_i] == tseq[output_index])
						contig_coverage++;
					output_index++;
				}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {output_index++;}	break;//S, print -
			default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		return contig_coverage;
	}
	int get_contig_NM(int suggest_st_pos, int &match_size){
		int NM = 0;
		match_size = 0;
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		output_index += suggest_st_pos;
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:
				for(int i = 0; i < cigar_len; i++, seq_i++)	{
					if(qseq[seq_i] != tseq[output_index]) NM++;
					output_index++;
				}
				match_size += cigar_len; break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	output_index += cigar_len;	break;//D, print '-'
			case 3:	output_index += cigar_len;	seq_i += cigar_len;break;//N, print N
			case 4:	output_index += cigar_len;	seq_i += cigar_len;	break;//S, print -
			}
		}
		return NM;
	}
private:
	void print_SV_canditate(int suggest_st_pos){
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		bool have_long_SV = false;
		int minSVlen = 50;

		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	break;//M
			case 1:	case 2: case 3:	case 4:
				if(cigar_len >= minSVlen) have_long_SV = true;
				break;//I, print nothing
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}

		int ref_index = suggest_st_pos;
		if(have_long_SV){
			fprintf(log_output, "Long SV detected:");
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
				case 0:	ref_index += cigar_len; break;//M
				case 1:	if(cigar_len >= minSVlen) {fprintf(log_output, "Insertion: "); fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } seq_i += cigar_len; break;//I, int chr_ID; int
				case 2:	if(cigar_len >= minSVlen) {fprintf(log_output, "Deletion: ");  fprintf(log_output, "Start %d:%d; length %d: ", chr_ID, ref_index + global_ref_pos, cigar_len);  } ref_index += cigar_len; break;//D, print '-'
				case 3:	ref_index += cigar_len; break;//N, print N
				case 4:	ref_index += cigar_len; break;//S, print -
				default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
				}
			}
		}
		fprintf(log_output, "\n");
	}

	void print_coverage(int suggest_st_pos, int contig_coverage, uint16_t * contig_depth){
		int output_index = 0;
		int seq_i = 0;
		uint32_t* bam_cigar = ez.cigar;
		int cigar_len = ez.n_cigar;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log_output, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "%c", contig_depth[seq_i] + '#'); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log_output, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "0"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log_output, "-"); output_index++;}	break;//S, print -
			default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log_output, "\t");
		fprintf(log_output, "\nCigar sequence: ");
		for(int i = 0; i < cigar_len; i++){
			fprintf(log_output, "%d%c", (bam_cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(bam_cigar[i] & BAM_CIGAR_MASK)]);
		}
		fprintf(log_output, "contig_coverage: [%d bp]\t", contig_coverage);

		fprintf(log_output, "\t");
		fprintf(log_output, "\n");
	}
};




#endif /* SV_SPA_CORE_CONTIG_ALIGNER_HPP_ */
