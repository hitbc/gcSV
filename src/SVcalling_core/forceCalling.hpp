/*
 * forceCalling.hpp
 *
 *  Created on: 2024年12月30日
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_FORCECALLING_HPP_
#define SVCALLING_CORE_FORCECALLING_HPP_

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include "cpp_lib/RefRegion.hpp"
#include "cpp_lib/get_option_cpp.hpp"
#include "cpp_lib/cpp_utils.hpp"
#include <math.h>
extern "C"
{
#include "../clib/utils.h"
#include "../clib/bam_file.h"
#include "../clib/vcf_lib.h"
#include "htslib/hfile.h"
}
//#-----------------------------
void showGT(int GT);

struct Haplotype_read_aligner{
public:
	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
	}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}

	int get_score_reach_end_of_read(){ return MAX(0, ez.mqe);}

	//when align_to_contig == true, aligned to contig, otherwise, aligned to refernece
	void align_genotyping(uint8_t *qseq, uint32_t qlen,  uint8_t *tseq, uint32_t tlen){
		ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	int gap_penalty(int gap_len){
		int penalty1 = gap_open_D + gap_len*gap_ex_D;
		int penalty2 = gap_open2_D + gap_len*gap_ex2_D;
		return MIN(penalty1, penalty2);
	}

	int getScoreByMismatch(int search_length, int mismatch_num){
		return ((search_length - mismatch_num) * match_D) - (mismatch_num * mismatch_D);
	}

	int getMaxScore(int search_length){
		return ((search_length) * match_D);
	}

	int adjustCIGAR(){
		uint32_t n_cigar = ez.n_cigar;
		int adj_size = cigar_adjust(&n_cigar, ez.cigar, false, 15);
		ez.n_cigar = n_cigar;
		return adj_size;
	}

	void printf_alignment_detail(FILE * output, int suggest_st_pos){
		log_output = output;
		int cigar_len = ez.n_cigar;
		//fprintf(output," CIGAR number: %d ", cigar_len);
		if(cigar_len == 0) return;
		fprintf(output," CIGAR: ", cigar_len);
		uint32_t* bam_cigar =  ez.cigar;
		 for (int i = 0; i < cigar_len; ++i){
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			char c_type = segment_type_to_cigar_code(type);
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			fprintf(output,"%d%c",length,  c_type);
		}
		fprintf(output," ");
	}

	bool aln_with_indel(){
		int cigar_len = ez.n_cigar;
		//fprintf(output," CIGAR number: %d ", cigar_len);
		if(cigar_len == 0) return true;
		uint32_t* bam_cigar =  ez.cigar;
		 for (int i = 0; i < cigar_len; ++i){
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			if(type == CIGAR_DELETE || type == CIGAR_INSERT || type == CIGAR_SOFT_CLIP){
				return true;
			}
		}
		return false;
	}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

	int getGenotypingMinScore(int used_read_length){
		int min_score = (used_read_length -80) * match_D;
		int global_min_match_score = 50 * match_D;
		return MAX(global_min_match_score, min_score);
	}

	int getGenotypingMinScoreFC(int used_read_length){
		return (used_read_length) * match_D * 0.9;
	}

private:
	//part3: alignment result
	ksw_extz_t ez;
	//part4: buffs
	void *km;
	//part5: options
	int8_t mata_D[25]; int8_t match_D; int8_t mismatch_D;
	int8_t gap_open_D; int8_t gap_ex_D; int8_t gap_open2_D; int8_t gap_ex2_D;
	uint16_t zdrop_D; int bandwith;
	int flag;
	//part6: logs
	FILE * log_output;

	void copy_option(){
			match_D = 2;
			mismatch_D= 6;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 1;
			zdrop_D= gap_open2_D + 200;
			bandwith = 200;//zdrop_D;
			flag = 0;
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
};

struct read_info_GT{
	std::string read_name;
	int read_pos;
	int read_ID;
	void set(const char * rn, int pos, int read_ID){
		this->read_pos = pos;
		this->read_ID = read_ID;
		this->read_name.append(rn);
	}
};

struct ALN_result{
	int overlap_mode;
	int read_score_p1;
	int read_score_p2;
	int aln_method_id1;
	int aln_method_id2;
	int aln_with_indel1;
	int aln_with_indel2;
	int read_score_MAX;

	int readID;

	void clear(){
		readID = -1;
		overlap_mode = -1;
		read_score_MAX = -1;
		aln_method_id1 = -1;
		aln_method_id2 = -1;
		read_score_p1 = -1;
		read_score_p2 = -1;
		aln_with_indel1 = -1;
		aln_with_indel2 = -1;
	}

	void show_result(int hap_ID){
		fprintf(stderr, "%d:%d:%d:%d[%d:%d]:%d[%d:%d]:%d\t",
				readID, hap_ID, overlap_mode,
				read_score_p1, aln_method_id1, aln_with_indel1,
				read_score_p2, aln_method_id2, aln_with_indel2, read_score_MAX);
	}
};

struct Haplotype_pairing{
	int hap_idx1;
	int hap_idx2;
	int total_read_with_indel;
	int hap_read_count_i;
	int hap_read_count_j;
	int hap_read_count_UN;
	int total_score;
	void set(int hap_idx1, int hap_idx2, int total_read_with_indel, int hap_read_count_i, int hap_read_count_j, int hap_read_count_UN, int total_score){
		this->hap_idx1 = hap_idx1;
		this->hap_idx2 = hap_idx2;
		this->total_read_with_indel = total_read_with_indel;
		this->hap_read_count_i = hap_read_count_i;
		this->hap_read_count_j = hap_read_count_j;
		this->hap_read_count_UN = hap_read_count_UN;
		this->total_score = total_score;
	}
	void show(int total_sig_read_n){
		fprintf(stderr, "%d:[%d:%d]:[%d:%d:%d]:%d:",total_read_with_indel, hap_idx1,hap_idx2,hap_read_count_i, hap_read_count_j,hap_read_count_UN, total_score);
		if(total_sig_read_n != 0)	fprintf(stderr, "%d\t", total_score/(total_sig_read_n));
		else						fprintf(stderr, "0\t");
	}

};

struct Haplotype_handler{
	int chr_ID;
	int sv_pos;
	int sv_end;
	int SV_length;
	std::string REF;
	std::string ALT;

	int ref_global_position;
	int contig_global_position;

	int contig_pos_in_ref_when_aln_to_bp1;
	int contig_pos_in_ref_when_aln_to_bp2;

	std::vector<uint8_t> contig_bin;

	std::vector<ALN_result> aln_r;

	int max_score_read_count;

	void set(int chr_ID_, int st_pos_, std::string & REF, std::string & ALT, int SV_length){
		this->chr_ID = chr_ID_;
		this->sv_pos = st_pos_;
		this->ALT = ALT;
		this->REF = REF;
		this->sv_end = st_pos_ + REF.size() - 1;
		this->SV_length = SV_length;
	}

	int get_contig_golbal_position_core(bool is_bp1){
		if(is_bp1 == false)	return contig_global_position - SV_length;
		else				return contig_global_position;
	}

	void get_contig_global_position(){
		contig_pos_in_ref_when_aln_to_bp1 = get_contig_golbal_position_core(true);
		contig_pos_in_ref_when_aln_to_bp2 = get_contig_golbal_position_core(false);
		//fprintf(stderr, "SUG:contig:bp1 %d; SUG:contig:bp2 %d\t", contig_pos_bp1, contig_pos_bp2);
	}

	void store_bin_from_char_append(const char * contig_seq, int contig_seq_len, std::vector<uint8_t> &bin_contig){
		//store bin contig
		for (int i = 0; i < contig_seq_len; ++i)
			bin_contig.emplace_back(charToDna5n[(uint8_t)contig_seq[i]]);
	}

	void store_contig(char * ref, int ref_len, int ref_position){
		contig_bin.clear();
		store_bin_from_char_append(ref, sv_pos + 1 - ref_position + 1, contig_bin);
		store_bin_from_char_append(ALT.c_str() + 1, ALT.size() - 1, contig_bin);
		int p2_pos = sv_end + 1 - ref_position + 1;
		int p2_len = ref_position + ref_len - sv_end - 1 - 1;
		store_bin_from_char_append(ref + p2_pos, p2_len, contig_bin);
		contig_global_position = ref_position;
		get_contig_global_position();
	}
};

struct GT_LOCAL_handler{
	int region_support_number[3][3]; //[region number] * [support read type number]
	int chr_ID;

	Haplotype_read_aligner ga;
	int normal_read_length;
	Haplotype_read_aligner * hra;
	Bam_file *bam_f;
	faidx_t *c_ref_idx;

	int global_read_count = 0;
	uint64_t global_qual_sum = 0;

	void init(Bam_file *c_b, int normal_read_length_, faidx_t *c_ref_idx){
		ga.init();
		this->normal_read_length = normal_read_length_;
		bam_f = c_b;
		hra = &ga;
		this->c_ref_idx = c_ref_idx;
		global_read_count = 0;
		global_qual_sum = 0;
	}

	int GT_main(bool print_log, int chrID, int ref_bgn_pos, int ref_end_pos, std::vector<Haplotype_handler> &hap_small_SV_l){
		FC_support_read_counting(print_log, bam_f, hra, normal_read_length, chrID, ref_bgn_pos, ref_end_pos, hap_small_SV_l);
		int NGS_suggest_GT = getGT_from_support_read(region_support_number[2][0], region_support_number[2][1], region_support_number[2][2]);
		fprintf(stderr,"Region analysis: for region(1+2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[2][0],
				region_support_number[2][1], region_support_number[2][2]);
		fprintf(stderr, "NGS_suggest_GT is: ");
		showGT(NGS_suggest_GT);
		fprintf(stderr, " \n");
		return NGS_suggest_GT;
	}

	int getMAX_score(int normal_read_len){
		return hra->getMaxScore(normal_read_len);
	}
	int get_alignment_core(bool printLog, int & CIGAR_Method_id, int & aln_with_indel, Haplotype_read_aligner * hra, bam1_t *br, uint8_t *qseq, std::vector<uint8_t> &hap_bin, int read_in_contig_st_pos, int minMismatchQUAL){
		CIGAR_Method_id = 1;
		aln_with_indel = false;
		int qlen = br->core.l_qseq;
		uint8_t *qual_str = bam_get_qual(br);
		int read_in_contig_ed_pos = read_in_contig_st_pos + qlen;
		//realigned for contig regions
		if(read_in_contig_st_pos < 0)
			return -1;
		if((int)hap_bin.size() < read_in_contig_ed_pos)
			return -1;
		int tlen = qlen;
		uint8_t* tseq = &(hap_bin[0]) + read_in_contig_st_pos;
		int search_len = MIN(qlen, tlen);
		if(printLog){
			fprintf(stderr, "Alignment (read_in_contig_st_pos %5d): ", read_in_contig_st_pos);
			for(int i = 0; i < search_len; i++)
				fprintf(stderr, "%c", "ACGT"[tseq[i]]);
		}
		int wrong_base = 0;
		int wrong_base_HIGH = 0;
		for(int i = 0; i < search_len && wrong_base_HIGH < 6; i++){
			if(tseq[i] != qseq[i]){
				if(qual_str[i] > minMismatchQUAL)
					wrong_base_HIGH ++;
				wrong_base++;
			}
		}
		if(printLog && wrong_base_HIGH < 6) fprintf(stderr,"(Simple search used: [%d-%d, len: %d], mismatch %d) \t", 0,	search_len, search_len, wrong_base_HIGH);
		if(wrong_base_HIGH < 6){
			int score = hra->getScoreByMismatch(search_len, wrong_base);
			return score;
		}
		//tseq -= 6; tlen += 6;
		hra->align_genotyping(qseq, qlen, tseq, tlen);
		hra->adjustCIGAR();
		int score = hra->get_score_reach_end_of_read();

		aln_with_indel = hra->aln_with_indel();

		if(printLog){
			//fprintf(stderr,"( Detail search used Contig:");
			hra->printf_alignment_detail(stderr, read_in_contig_st_pos);
			fprintf(stderr,")\tscore:%d\t", score);
		}
		CIGAR_Method_id = 2;
		return hra->get_score_reach_end_of_read();
	}

	int getGT_from_support_read(int read_sup, int read_not_sup, int read_unknown){
	    //signal read number adjust for insertions and deletions
		int GT;
			if(read_sup > read_not_sup * 4)     GT = 3;//1/1
			else if(read_sup*3 < read_not_sup)  GT = 1;//0/0
			else                                GT = 2;//0/1
			if((read_sup + read_not_sup) * 3 < read_unknown || (read_sup < 4 && read_not_sup < 4))
												GT = 0;// ./.
		return GT;
	}

	//return:
	//0: not overlap with both
	//1: overlap with bp1
	//2: overlap with bp2
	//3: overlap with both
	int FCread_overlap_breakpoint(int min_overlap_len, int true_read_pos_bg, int true_read_pos_ed, int breakpoint1, int breakpoint2){
		bool read_overlap_with_breakpoint1 = (true_read_pos_bg <= breakpoint1 && true_read_pos_ed > breakpoint1 && ((breakpoint1 - true_read_pos_bg) > min_overlap_len));
		bool read_overlap_with_breakpoint2 = (true_read_pos_bg <= breakpoint2 && true_read_pos_ed > breakpoint2 && ((true_read_pos_ed - breakpoint2) > min_overlap_len));
		return read_overlap_with_breakpoint1 + read_overlap_with_breakpoint2 * 2;
	}

	//get the read position of the first base and the last base in the reference
	void get_true_read_pos(bam1_t *br, int * true_read_pos_bg, int *true_read_pos_ed){
		int begin_pos = br->core.pos + 1;//change to 1-base
		int end_pos = br->core.pos + 1;//change to 1-base

		uint32_t* bam_cigar = bam_get_cigar(br);
		uint32_t n_cigar = br->core.n_cigar;
		for (uint i = 0; i < n_cigar; ++i)
		{
			int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			switch (c_type){
			case CIGAR_MATCH: case CIGAR_SEQ_MATCH: end_pos += c_size; break;
			case CIGAR_INSERT:	break; //do nothing
			case CIGAR_DELETE:	end_pos += c_size; break;
			case CIGAR_SOFT_CLIP:	case CIGAR_HARD_CLIP:
				if(i == 0) 	begin_pos -= c_size; else end_pos += c_size;
				break;
			default:	break;
			}
		}
		*true_read_pos_bg = begin_pos;
		*true_read_pos_ed = end_pos;

		if(false){
			fprintf(stderr,"\t Ori Cigar: ");
			for (unsigned int i = 0; i < n_cigar; ++i){
				int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
				char c_type = segment_type_to_cigar_code(type);
				int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				fprintf(stderr,"%d%c",length,  c_type);
			}
		}
	}

	void FC_support_read_counting(bool print_log,  Bam_file *bam_f, Haplotype_read_aligner * hra, int normal_read_length,
			int chrID, int ref_bgn_pos, int ref_end_pos, std::vector<Haplotype_handler> &hap_small_SV_l){
		uint8_t qseq_buff[1024];
		for(int region_ID = 0; region_ID < 3; region_ID++)
			for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
				region_support_number[region_ID][read_genotype_type] = 0;

		R_region region;
		region.chr_ID =chrID;
		region.st_pos = ref_bgn_pos;
		region.ed_pos = ref_end_pos;

		resetRegion_ID(bam_f, &region);	//reset region
		//reference check:
		bool show_aln_detail = false;
		for(uint i = 0; i < hap_small_SV_l.size(); i++)
			hap_small_SV_l[i].max_score_read_count = 0;

		int read_id = -1;
		std::vector<read_info_GT> rinfo;
		int read_len = -1;
		while (bam_next(bam_f)) {
			read_id++;
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))       	continue;
			if(br->core.qual < 10)       	continue;

			{
				int total_qual = 0;
				read_len = br->core.l_qseq;
					uint8_t *qual_str = bam_get_qual(br);
				for(int i = 0; i < read_len; i++)
					total_qual += qual_str[i];
				fprintf(stderr, "total_qual is %d\n", total_qual);

				global_read_count++;
				global_qual_sum += total_qual;
				if(global_read_count > 100 && total_qual * global_read_count < global_qual_sum * 0.8){
					continue;
				}
			}
			//load reads
			int true_read_pos_bg; int true_read_pos_ed; int min_overlap_len;
			{
				read_len = br->core.l_qseq;
				min_overlap_len = 0.15*read_len;
				get_bam_seq_bin(0, read_len, qseq_buff, br);
				get_true_read_pos(br, &true_read_pos_bg, &true_read_pos_ed);
				if(show_aln_detail)// show reads
				{
					fprintf(stderr, "\n\n@Handle New read: Read name: %s ALN POS: %d \t", bam_qname(br), br->core.pos);
					fprintf(stderr, "true_read_pos: [%d, %d]\t", true_read_pos_bg, true_read_pos_ed);
					fprintf(stderr, "@Read length %d, @Read string: ", read_len);
					for(int i = 0; i < read_len; i++)
						fprintf(stderr, "%c", "ACGT"[qseq_buff[i]]);
					fprintf(stderr, "\n");
				}
			}
			bool is_overlap_with_BP = false;
			//for each haplotype:
			//store the result:
			rinfo.emplace_back();
			rinfo.back().set((char*)bam_qname(br), br->core.pos, read_id);

			for(uint i = 0; i < hap_small_SV_l.size(); i++){
				Haplotype_handler & hh  = hap_small_SV_l[i];
				hh.aln_r.emplace_back();
				ALN_result & alr = hh.aln_r.back();
				alr.clear();
				bool is_ref = (i ==  hap_small_SV_l.size() - 1);
				int overlap_mode = 0;
				if(is_ref && is_overlap_with_BP)	  {	overlap_mode = 1;}
				else if(is_ref && !is_overlap_with_BP){	overlap_mode = 0;}
				else								  { overlap_mode = FCread_overlap_breakpoint(min_overlap_len, true_read_pos_bg, true_read_pos_ed, hh.sv_pos, hh.sv_end);}

				if(show_aln_detail)fprintf(stderr, "\n@Read name: %s @[SV_POS: %d-%d] Overlap_mode %d;\n", bam_qname(br), hh.sv_pos, hh.sv_end, overlap_mode);
				if(overlap_mode == 0)
					continue;
				else
					is_overlap_with_BP = true;
				//realigned for contig regions
				int read_in_contig_st_pos_bp1_rbg = (true_read_pos_bg) - hh.contig_pos_in_ref_when_aln_to_bp1;
				int read_in_contig_st_pos_bp2_rbg = (true_read_pos_bg) - hh.contig_pos_in_ref_when_aln_to_bp2;

				int minMismatchQUAL = 12;
				alr.overlap_mode = overlap_mode;
				alr.readID = read_id;
				if(overlap_mode == 1){
					int read_contig_score = get_alignment_core(show_aln_detail, alr.aln_method_id1, alr.aln_with_indel1, hra, br, qseq_buff, hh.contig_bin, read_in_contig_st_pos_bp1_rbg, minMismatchQUAL);
					alr.read_score_MAX = alr.read_score_p1 = read_contig_score;
				}else if(overlap_mode == 2){
					int read_contig_score = get_alignment_core(show_aln_detail, alr.aln_method_id2, alr.aln_with_indel2, hra, br, qseq_buff, hh.contig_bin, read_in_contig_st_pos_bp2_rbg, minMismatchQUAL);
					alr.read_score_MAX = alr.read_score_p2 = read_contig_score;
				}
				else if(overlap_mode == 3){
					int read_contig_score1 = get_alignment_core(show_aln_detail, alr.aln_method_id1, alr.aln_with_indel1, hra, br, qseq_buff, hh.contig_bin, read_in_contig_st_pos_bp1_rbg, minMismatchQUAL);
					int read_contig_score2 = get_alignment_core(show_aln_detail, alr.aln_method_id2, alr.aln_with_indel2, hra, br, qseq_buff, hh.contig_bin, read_in_contig_st_pos_bp2_rbg, minMismatchQUAL);
					alr.read_score_MAX = MAX(read_contig_score1, read_contig_score2);
					alr.read_score_p1 = read_contig_score1;
					alr.read_score_p2 = read_contig_score2;
				}
			}
			if(!is_overlap_with_BP){
				for(int i = (int)hap_small_SV_l.size() - 1; i >= 0;  i--){
					hap_small_SV_l[i].aln_r.pop_back();
				}
				rinfo.pop_back();
			}
		}

		int total_sig_read_n = hap_small_SV_l[0].aln_r.size();
		//logs:
		//(1)
		int total_hap_n = hap_small_SV_l.size();
		for(int read_idx = 0; read_idx < total_sig_read_n; read_idx++){
			int max_score = -1;
			int max_idx = -1;
			int max_with_indel = -1;
			//if(print_log)fprintf(stderr, "\n");
			bool is_uniq_max = true;
			for(int hap_idx = (int)hap_small_SV_l.size() - 1; hap_idx >= 0;  hap_idx--){
				ALN_result & alr = hap_small_SV_l[hap_idx].aln_r[read_idx];
				if(alr.read_score_MAX > max_score){
					max_score = alr.read_score_MAX;
					max_idx = hap_idx;
					max_with_indel = (alr.aln_with_indel1 && alr.aln_with_indel2);
				}
				else if(alr.read_score_MAX == max_score){
					is_uniq_max = false;
				}
				if(print_log)
					alr.show_result(hap_idx);
			}
			if(is_uniq_max && max_idx != -1)
				hap_small_SV_l[max_idx].max_score_read_count ++;
			if(print_log){
				fprintf(stderr, "max_score: %d @[%d:%d]\t", max_score, max_idx, max_with_indel);
				fprintf(stderr, "@Read name: %s @POS %d ID %d \n", rinfo[read_idx].read_name.c_str(),  rinfo[read_idx].read_pos, rinfo[read_idx].read_ID );
			}
		}

		Haplotype_handler & ref_hap = hap_small_SV_l.back();
		std::vector<ALN_result> & rr = ref_hap.aln_r;

		std::vector<Haplotype_pairing> hp;
		int max_score_hp_idx = -1;
		int max_score_hp_score = -1;
		bool is_same_max_score = false;
		//fit chech:
		for(int i = 0; i < total_hap_n; i++){
			std::vector<ALN_result> & alr_i = hap_small_SV_l[i].aln_r;
			for(int j = i; j < total_hap_n; j++){
				std::vector<ALN_result> & alr_j = hap_small_SV_l[j].aln_r;
				int total_read_with_indel = 0;
				int hap_read_count_i = 0;
				int hap_read_count_j = 0;
				int hap_read_count_UN = 0;
				int total_score = 0;
				for(int read_idx = 0; read_idx < total_sig_read_n; read_idx++){
					ALN_result & ri = alr_i[read_idx]; ALN_result & rj = alr_j[read_idx];
					if(ri.readID == -1)	ri = rr[read_idx];
					if(rj.readID == -1) rj = rr[read_idx];
					bool aln_with_indel_i = ri.aln_with_indel1 && ri.aln_with_indel2;
					bool aln_with_indel_j = rj.aln_with_indel1 && rj.aln_with_indel2;
					if(aln_with_indel_i && aln_with_indel_j) total_read_with_indel++;
					if(ri.read_score_MAX > rj.read_score_MAX + 10)		hap_read_count_i ++;
					else if(rj.read_score_MAX > ri.read_score_MAX + 10)	hap_read_count_j ++;
					else												hap_read_count_UN ++;
					total_score += MAX(rj.read_score_MAX, ri.read_score_MAX);
				}
				if(total_score >= max_score_hp_score){
					if(total_score == max_score_hp_score)	is_same_max_score = true;
					else									is_same_max_score = false;
					max_score_hp_score = total_score;
					max_score_hp_idx = hp.size();
				}
				hp.emplace_back();
				hp.back().set(i, j, total_read_with_indel, hap_read_count_i, hap_read_count_j, hap_read_count_UN, total_score);
				hp.back().show(total_sig_read_n);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "The max pair is: ");
		hp[max_score_hp_idx].show(total_sig_read_n);
		fprintf(stderr, "is_same_max_score: %d\n", is_same_max_score);

		int pass_filter = 0;
		//filter1: UNKNOWN-INDEL check
		if(pass_filter == 0 && hp[max_score_hp_idx].total_read_with_indel > 1){
			pass_filter = 1;
		}

		//filter2: HOM-REF check:
		{
			if(pass_filter == 0 && max_score_hp_idx == hp.size() - 1){
				pass_filter = 2;
			}
			//one of hap is REF
			if(pass_filter == 0 && hp[max_score_hp_idx].hap_idx2 == total_hap_n - 1){
				if(hp[max_score_hp_idx].hap_read_count_i*3 <= hp[max_score_hp_idx].hap_read_count_j){
					pass_filter = 3;
				}
			}
			//similar ref check
			if(pass_filter == 0 && hp.back().total_read_with_indel < 3){
				if(hp[max_score_hp_idx].total_score < hp.back().total_score + 12 * total_sig_read_n){
					pass_filter = 4;
				}
			}
			//similar ref check
			if(pass_filter == 0 && hp[max_score_hp_idx].total_score + 35 * total_sig_read_n < read_len * total_sig_read_n){
				if(hp[max_score_hp_idx].total_score < hp.back().total_score + 12 * total_sig_read_n){
					pass_filter = 5;
				}
			}
			//similar ref check
			if(pass_filter == 0 && hp[max_score_hp_idx].hap_read_count_i < 3){
				pass_filter = 6;
			}
			//similar ref check
			if(pass_filter == 5 && (hp[max_score_hp_idx].hap_read_count_i +  hp[max_score_hp_idx].hap_read_count_j)*4 < hp[max_score_hp_idx].hap_read_count_UN){
				pass_filter = 7;
			}
		}
		fprintf(stderr, "The final pass_filter is %d\n", pass_filter);


		bool with_sig = false;
		for(uint i = 0; i < hap_small_SV_l.size() - 1; i++){
			if(hap_small_SV_l[i].max_score_read_count >= 5){
				with_sig = true;
			}
		}
		//with_sig = (max_alt_read_n > 3);
		if(with_sig){
			fprintf(stderr, "\nWith SIG\t");
		}
		else
			fprintf(stderr, "\nNo SIG\t");

		if(true){
			fprintf(stderr, "max_score_read_count\t");
			fprintf(stderr, "\n");
			for(uint i = 0; i < hap_small_SV_l.size(); i++){
				fprintf(stderr, "%d:%d\t", i, hap_small_SV_l[i].max_score_read_count);
			}
			fprintf(stderr, "\n");
			//fprintf(stderr, "unknown_read_n %d\tmax_ref_read_n %d\t max_alt_read_n %d\tmax_SIM_read_n %d\t\n",unknown_read_n, max_ref_read_n, max_alt_read_n, max_SIM_read_n);
			for(uint i = 0; i < hap_small_SV_l.size(); i++){
				fprintf(stderr, "%d:%d\t",hap_small_SV_l[i].SV_length, hap_small_SV_l[i].sv_pos);
			}
			fprintf(stderr, "\n");

		}

		//fprintf(stderr,"Region analysis: for region1: ALT:REF:UNKNOW [%d, %d, %d]\t", region_support_number[0][0], region_support_number[0][1], region_support_number[0][2]);
		//fprintf(stderr,"Region analysis: for region2: ALT:REF:UNKNOW [%d, %d, %d]\t", region_support_number[1][0], region_support_number[1][1], region_support_number[1][2]);
		//for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)region_support_number[2][read_genotype_type] = region_support_number[0][read_genotype_type] + region_support_number[1][read_genotype_type];
	}
};

struct FC_handler{

	//data for vcf files
	bcf_hdr_t *vcf_header = NULL;
	//data for bam files
	Bam_file c_b;
	bam_hdr_t* bam_header = NULL;

	void init_bam(char * input_bam_fn, char * ref_fn){
		bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
		bam_header = c_b._hdr;
	}

	int ana_block_length;
	int ana_block_step_size;

	void genotypingDRP_long_INS(int read_chrID,  int SV_POS, int SV_END, int &DRP_support_SV, int &DRP_support_ref ){
		int max_isize = 1000;
		int normal_read_length = 150;
		DRP_support_SV = 0;
		DRP_support_ref = 0;
		RefRegion withInRef_p1(read_chrID, SV_POS - (max_isize - normal_read_length), SV_POS);
		RefRegion withInRef_p2(read_chrID, SV_POS, SV_POS + (max_isize - normal_read_length - normal_read_length));

		//only consider reads in the first region:
		R_region bam_load_region;
		bam_load_region.chr_ID = read_chrID;
		bam_load_region.st_pos = SV_POS - max_isize;
		bam_load_region.ed_pos = SV_POS + normal_read_length;

		resetRegion_ID(&c_b, &bam_load_region);
		while (bam_next(&c_b)) {
			bam1_t *br = &(c_b._brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))		continue;
			if(br->core.qual < 10)			continue;
			//get iSIZE
			int isize = br->core.isize;
			if(br->core.mtid == -1 || isize < 0)
				continue;
			if(isize <= 0 || isize > 1000000){
				DRP_support_SV ++;
			}
			else{
				int pos = br->core.pos + br->core.l_qseq;
				int mpos = br->core.mpos;

				if(withInRef_p1.pos_within_region(pos) && withInRef_p2.pos_within_region(mpos)){//for read pairs that support the reference
					DRP_support_ref++;
				}
			}
		}
	}
	void genotypingDRP_long_DEL(int read_chrID,  int SV_POS, int SV_END, int &DRP_support_SV, int &DRP_support_ref){
		int max_isize = 1000;
		int normal_read_length = 150;

		RefRegion withInRef_p1(read_chrID, SV_POS - (max_isize - normal_read_length), SV_POS);
		RefRegion withInSV_p2(read_chrID, SV_END, SV_END + (max_isize - normal_read_length - normal_read_length));
		RefRegion withInRef_p2(read_chrID, SV_POS, SV_POS + (max_isize - normal_read_length - normal_read_length));

		//only consider reads in the first region:
		R_region bam_load_region;
		bam_load_region.chr_ID = read_chrID;
		bam_load_region.st_pos = SV_POS - max_isize;
		bam_load_region.ed_pos = SV_POS + normal_read_length;

		int64_t total_isize_deletion = 0;
		int64_t total_isize_reference = 0;

		resetRegion_ID(&c_b, &bam_load_region);
		while (bam_next(&c_b)) {
			bam1_t *br = &(c_b._brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))		continue;
			if(br->core.qual < 10)			continue;
			//get iSIZE
			int isize = br->core.isize;
			if(isize <= 0 || isize > 1000000) continue;

			int pos = br->core.pos + br->core.l_qseq;
			int mpos = br->core.mpos;

			if(	withInRef_p1.pos_within_region(pos) && withInSV_p2.pos_within_region(mpos)){//for read pairs that support the deletion
				DRP_support_SV++;
				total_isize_deletion += isize;
			}

			if(withInRef_p1.pos_within_region(pos) && withInRef_p2.pos_within_region(mpos)){//for read pairs that support the reference
				DRP_support_ref++;
				total_isize_reference += isize;
			}
		}
	}

	GT_LOCAL_handler gth;
	faidx_t *c_ref_idx;

	struct merging_SV_info{
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
		std::vector<SAMPLE_info> supportSAMPLE_ID;

		void store_basic(int pos, std::string &ALT_str, std::string &REF_str, int SV_len){
			this->pos = pos;
			this->SV_len = SV_len;
			this->ALT_str = ALT_str;
			this->REF_str = REF_str;
		}

		void load(std::vector<std::string> &item_value, std::vector<std::string> &buff){
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

		void show(){
			fprintf(stderr, "Original: %d\t%d\t%d\t%s\t%s\n", chrID, pos, SV_len, REF_str.c_str(), ALT_str.c_str());
		}

		static inline int cmp_by_pos(const merging_SV_info &a, const merging_SV_info &b){
			if(a.pos != b.pos)
				return a.pos < b.pos;
			return a.SV_len < b.SV_len;
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

	};
	int forceCalling_run(int argc, char *argv[]){
		bool print_log = true;
		//parameters
		char * ref_fn = argv[1];
		char * vcf_fn_in = argv[2];//separate by ','
		char * input_bam_fn = argv[3];
		char * local_repeat_fn = argv[4];
		int sampleID = atoi(argv[5]);

		std::vector<merging_SV_info> SV_info;
		//load vcf binary
		{
			std::vector<std::string> vcf_record_l;
			load_string_list_from_file_MAX_line("/home/fenghe/E", vcf_record_l, 50000);
			std::vector<std::string> item_value;
			std::vector<std::string> item_value_buff;
			for(std::string &vcf_r :vcf_record_l){
				split_string(item_value, vcf_r.c_str(), "\t");
				SV_info.emplace_back();
				SV_info.back().load(item_value, item_value_buff);
			}
		}
		//open bam file
		init_bam(input_bam_fn, ref_fn);
		std::vector<Haplotype_handler> hap_small_SV_l;
		std::vector<Haplotype_handler> hap_long_del_l;
		std::vector<Haplotype_handler> hap_long_ins_l;
		std::vector<Haplotype_handler> hap_CPX_l;
		c_ref_idx = reference_index_load(ref_fn);
		gth.init(&c_b, 150, c_ref_idx);
		//for each SVs
		int clusterSV_number = 0;
		int Addition_load_size = 300;
		for(uint SV_ID = 0; SV_ID < SV_info.size(); SV_ID += clusterSV_number)
		{
			//load the SV cluster
			bool already_withSV = false;
			clusterSV_number = 0;
			int clusterSampleNumber = 0;
			hap_small_SV_l.clear();
			hap_long_del_l.clear();
			hap_long_ins_l.clear();
			int max_pos = SV_info[SV_ID].pos;
			for(uint j = SV_ID; j < SV_info.size(); j++){
				bool is_nearby = (SV_info[j].pos - max_pos) < 30;
				max_pos = MAX(max_pos, SV_info[j].pos);
				if(is_nearby){
					clusterSV_number ++;
					if(SV_info[j].SV_len > 800){
						hap_long_ins_l.emplace_back();
						hap_long_ins_l.back().set(SV_info[j].chrID, SV_info[j].pos, SV_info[j].REF_str, SV_info[j].ALT_str, SV_info[j].SV_len);
					}
					else if(SV_info[j].SV_len < -700){
						hap_long_del_l.emplace_back();
						hap_long_del_l.back().set(SV_info[j].chrID, SV_info[j].pos, SV_info[j].REF_str, SV_info[j].ALT_str, SV_info[j].SV_len);
					}else if(SV_info[j].REF_str.size() > 30 && SV_info[j].ALT_str.size() > 30){
						hap_CPX_l.emplace_back();
						hap_CPX_l.back().set(SV_info[j].chrID, SV_info[j].pos, SV_info[j].REF_str, SV_info[j].ALT_str, SV_info[j].SV_len);
					}
					else{
						hap_small_SV_l.emplace_back();
						hap_small_SV_l.back().set(SV_info[j].chrID, SV_info[j].pos, SV_info[j].REF_str, SV_info[j].ALT_str, SV_info[j].SV_len);
					}
					clusterSampleNumber += SV_info[j].supportSAMPLE_ID.size();
					for(merging_SV_info::SAMPLE_info s : SV_info[j].supportSAMPLE_ID){
						if(s.ID == sampleID){
							already_withSV = true;
							//output the SVs
							//todo::
							SV_info[j].show();
						}
					}
				}else
					break;
			}
			if(SV_info[SV_ID].pos < 23504400) continue;
			fprintf(stderr, "clusterSV_number %d clusterSampleNumber %d hap_small_SV_l.size() %ld @ %d\n", clusterSV_number, clusterSampleNumber, hap_small_SV_l.size(), SV_info[SV_ID].pos);
			//load original SV info
			std::string BLANK = "N";
			if(already_withSV == false && clusterSampleNumber != 1){//begin the force calling
			//if(clusterSampleNumber != 1){//begin the force calling
				//S1: store the haplotype list
				if(!hap_small_SV_l.empty() && hap_small_SV_l.size() < 20)
				{
					int min_bgn = MAX_int32t;
					int max_end = 0;
					for(Haplotype_handler &hh : hap_small_SV_l){
						min_bgn = MIN(min_bgn, hh.sv_pos);
						max_end = MAX(max_end, hh.sv_end);
					}
					int chrID = SV_info[SV_ID].chrID;
					int addition_load_ref_left = Addition_load_size;
					int addition_load_ref_right = Addition_load_size;
					int ref_bgn_pos = min_bgn - addition_load_ref_left + 1 + 1;
					int ref_end_pos = max_end + addition_load_ref_right + 1;

					//load reference
					char reg[1024]; int load_len = 0;
					sprintf(reg, "%s:%d-%d", faidx_iseq(c_ref_idx, chrID), ref_bgn_pos, ref_end_pos);
					char *ref = fai_fetch(c_ref_idx, reg, &load_len);
					//store the reference as a hap
					hap_small_SV_l.emplace_back();
					hap_small_SV_l.back().set(chrID, ref_bgn_pos, BLANK, BLANK, 0);

					for(Haplotype_handler &hh : hap_small_SV_l){
						hh.store_contig(ref, load_len, ref_bgn_pos);
					}

					//end
					free(ref); ref = NULL;
					//show all haplotype
					if(print_log){
						for(Haplotype_handler &hh : hap_small_SV_l){
							for(uint i = 0; i < hh.contig_bin.size(); i++){
								fprintf(stderr, "%c", "ACGTNNNN"[hh.contig_bin[i]]);
							}
							fprintf(stderr, "\n");
						}
					}
					gth.GT_main(print_log, chrID, ref_bgn_pos, ref_end_pos, hap_small_SV_l);
				}

				//GT all SVs:
				{
//					fprintf(stderr, "Running FC: clusterSV_number %d clusterSampleNumber %d @ %d\n", clusterSV_number, clusterSampleNumber, SV_info[SV_ID].pos);
//					//for DRP reads
//					int DRP_support_SV = 0;
//					int DRP_support_ref = 0;
//					bool using_DRP_for_GT = false;
//					if(SV_END - SV_POS > 700){//long deletions
//						using_DRP_for_GT = true;
//						genotypingDRP_long_DEL(read_chrID, SV_POS, SV_END, DRP_support_SV, DRP_support_ref);
//					}
//					if(SV_length > 800){//long insertion
//						using_DRP_for_GT = true;
//						genotypingDRP_long_INS(read_chrID, SV_POS, SV_END, DRP_support_SV, DRP_support_ref);
//					}
				}
			}
			//break;

		}
		//close files
		bam_file_close(&c_b);
		return 0;
	}
};


#endif /* SVCALLING_CORE_FORCECALLING_HPP_ */
