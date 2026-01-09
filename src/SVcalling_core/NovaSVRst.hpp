/*
 * NovaSVRst.hpp
 *
 *  Created on: 2021-8-20
 *      Author: fenghe
 */

#ifndef NOVASVGENERATEVCF_NOVASVRST_HPP_
#define NOVASVGENERATEVCF_NOVASVRST_HPP_

#include <string>
#include <vector>
#include <cstring>

#include "../SVcalling_core/RefHandler.hpp"

extern "C"{
	#include "../clib/desc.h"
	#include "../clib/utils.h"
	#include "../clib/bam_file.h"
	#include "../kswlib/kalloc.h"
	#include "../kswlib/ksw2.h"
}

struct Genotyping_read_aligner{
public:
	//as input
	void init(){
		km = km_init();
		memset(&ez, 0, sizeof(ksw_extz_t));
		//mapping options
		copy_option();
		ksw_gen_mat_D();
	}

	void setRef(uint8_t * ori_ref_string, int ori_ref_len_, const char * contig_string, int contig_len)
	{
		bin_contig.resize(contig_len);
		//change char contig into bin contig
		for (int i = 0; i < contig_len; ++i)
			bin_contig[i] = charToDna5n[(uint8_t)contig_string[i]];
		ori_ref = ori_ref_string;
		ori_ref_len = ori_ref_len_;
	}

	void setRef_M2(uint8_t * ori_ref_string, int ori_ref_len_, uint8_t * contig_string, int contig_len)
	{
		bin_contig.resize(contig_len);
		//change char contig into bin contig
		for (int i = 0; i < contig_len; ++i)
			bin_contig[i] = contig_string[i];
		ori_ref = ori_ref_string;
		ori_ref_len = ori_ref_len_;
	}

	void destory(){
		if(ez.cigar != NULL)
			free(ez.cigar);
		km_destroy(km);
	}

	void align_non_splice(uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		qseq = qseq_;
		qlen = qlen_;
		if(ref_end_pos > (int)tlen)	ref_end_pos = tlen;
		//simple check:
		if(0){
			//debug code:
			uint32_t debug_qlen = qlen;
			uint8_t* debug_qseq = qseq;
			for(uint32_t i = 0; i < debug_qlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_qseq[i]]);
			fprintf(stderr, "\n");
			uint32_t debug_tlen = ref_end_pos - ref_st_pos;
			uint8_t* debug_tseq =  tseq + ref_st_pos;
			for(uint32_t i = 0; i < debug_tlen; i++)
				fprintf(stderr, "%c", "ACGT"[ debug_tseq[i]]);
			fprintf(stderr, "\n");
		}
		ksw_extd2_sse(km, qlen, qseq, ref_end_pos - ref_st_pos, tseq + ref_st_pos, 5, mata_D, gap_open_D, gap_ex_D, gap_open2_D, gap_ex2_D, bandwith, zdrop_D, -1, flag, &ez);
	}

	int get_score_reach_end_of_read(){ return MAX(0, ez.mqe);}

	//when align_to_contig == true, aligned to contig, otherwise, aligned to refernece
	void align_genotyping(bool align_to_contig, uint8_t *qseq_, uint32_t qlen_, int ref_st_pos, int ref_end_pos){
		if(align_to_contig){
			tseq = &(bin_contig[0]);
			tlen = bin_contig.size();
		}else{
			tseq = ori_ref;
			tlen = ori_ref_len;
		}
		align_non_splice(qseq_, qlen_, ref_st_pos, ref_end_pos);
	}

	int gap_penalty(int gap_len){
		int penalty1 = gap_open_D + gap_len*gap_ex_D;
		int penalty2 = gap_open2_D + gap_len*gap_ex2_D;
		return MIN(penalty1, penalty2);
	}
	int getScoreByCigar_with_skip_region(bam1_t *br, int read_skip_left, int read_skip_right,
			uint8_t * qseq_buff, RefHandler *refHandler, int MIN_misMAtchQuality){
		int read_left_boundary = read_skip_left;
		int read_right_boundary = br->core.l_qseq - read_skip_right;

		int score = 0;
		uint32_t* bam_cigar = bam_get_cigar(br);
		uint8_t* bam_quality = bam_get_qual(br);
		int match_base, mis_match_base;
		int q_seq_idx = 0;
		int t_seq_idx = 0;
		uint8_t *qseq_str = qseq_buff;
		uint8_t *tseq_str = refHandler->getRefStr(br->core.pos);
		for (uint i = 0; i < br->core.n_cigar; ++i)
		{
			int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			switch (c_type){
			case CIGAR_MATCH:
			case CIGAR_SEQ_MATCH:
				match_base = 0;
				mis_match_base = 0;
				for(int i = 0; i < c_size; i++, q_seq_idx++, t_seq_idx++){
					if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
						continue;
					if(qseq_str[q_seq_idx] != tseq_str[t_seq_idx] && bam_quality[q_seq_idx] > MIN_misMAtchQuality){
						mis_match_base++;
					}
					else
						match_base++;
				}
				score += ((match_base * match_D) - (mis_match_base * mismatch_D)); break;
			case CIGAR_INSERT:	case CIGAR_SOFT_CLIP:	case CIGAR_HARD_CLIP:
				mis_match_base = 0;
				for(int i = 0; i < c_size; i++, q_seq_idx++){
					if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
						continue;
					mis_match_base++;
				}
				score -= gap_penalty(mis_match_base); break;
			case CIGAR_DELETE:
				t_seq_idx += c_size;
				if(q_seq_idx < read_left_boundary || q_seq_idx >= read_right_boundary)
					break;
				score -= gap_penalty(c_size); break;
			default:
				break;
			}
		}
		score = MAX(0, score);
		return score;
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
		fprintf(output," CIGAR number: %d ", cigar_len);
		if(cigar_len == 0) return;
		fprintf(output," \tCIGAR: ", cigar_len);
		uint32_t* bam_cigar =  ez.cigar;
		 for (int i = 0; i < cigar_len; ++i){
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			char c_type = segment_type_to_cigar_code(type);
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			fprintf(output,"%d%c",length,  c_type);
		}
		fprintf(output," ");
		if(0){
			print_X_E_sequence(suggest_st_pos); //print X/= sequence
			fprintf(output, "\n");
		}
	}

	inline uint8_t getTseq(int i){ return tseq[i];}

	void setZdrop(uint16_t zdrop_D_, int bandwith_){
		zdrop_D = zdrop_D_;
		bandwith = bandwith_;
	}

	uint8_t * get_bin_contig(){ return &(bin_contig[0]); }
	uint8_t * get_ori_ref(){ return ori_ref; }

	int getGenotypingMinScore(int used_read_length){
		int min_score = (used_read_length -80) * match_D;
		int global_min_match_score = 50 * match_D;
		return MAX(global_min_match_score, min_score);
	}

	int getGenotypingMinScoreFC(int used_read_length){
		return (used_read_length) * match_D * 0.9;
	}

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
	uint16_t zdrop_D; int bandwith;
	int flag;
	//part6: logs
	FILE * log_output;
	std::vector<uint8_t> bin_contig;

	uint8_t* ori_ref; int ori_ref_len;

	void copy_option(){
			match_D = 2;
			mismatch_D= 6;
			gap_open_D= 24;
			gap_ex_D= 2;
			gap_open2_D= 32;
			gap_ex2_D= 1;
			zdrop_D= gap_open2_D + 30;
			bandwith = 30;//zdrop_D;
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

	int print_X_E_sequence(int suggest_st_pos){
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
			case 0://M
			case 7://=
			case 8://X
				for(int i = 0; i < cigar_len; i++, seq_i++)
				{
					if(qseq[seq_i] == tseq[output_index]){
						contig_coverage++;
						fprintf(log_output, "=");//1 == M; 2 == X; 0 == -;
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
};

bool bam_aligned_analysis(bam1_t *b, int *clip_left, int *clip_right, int *gap_mismatch_inside);

struct NOVA_SV_TGS_INFO{
	int global_region_ID = 0;
	int supp_read_n;
	int not_supp_read_n;
	int unknown_read_n;
	int haplotype_ID;
	int total_hap_number_in_region;
	int SV_in_HAP_ID;
	int total_SV_number_in_HAP;
	int SV_hap_position;

	void setINFO(
			int global_region_ID, int supp_read_n, int not_supp_read_n, int unknown_read_n,
			int haplotype_ID, int total_hap_number_in_region, int SV_in_HAP_ID, int total_SV_number_in_HAP, int SV_hap_position){
		this->global_region_ID = global_region_ID;
		this->supp_read_n = supp_read_n;
		this->not_supp_read_n = not_supp_read_n;
		this->unknown_read_n = unknown_read_n;
		this->haplotype_ID = haplotype_ID;
		this->total_hap_number_in_region = total_hap_number_in_region;
		this->SV_in_HAP_ID = SV_in_HAP_ID;
		this->total_SV_number_in_HAP = total_SV_number_in_HAP;
		this->SV_hap_position = SV_hap_position;
	}
};

#define SOMATIC_UNKNOWN 0
#define SOMATIC_TUMOR_UNIQ 1
#define SOMATIC_NORMAL_UNIQ 2
#define SOMATIC_TUMOR_NORMAL_SHARE 3

struct TUMOR_FLAG{
	TUMOR_FLAG(){
		is_SOMATIC = false;
		flag = SOMATIC_UNKNOWN;
	}
	bool is_SOMATIC;
	int flag;
	void set(int flag){
		is_SOMATIC = true;
		this->flag = flag;
	}

	void get_TUMOR_FLAG_STR(std::string & s){
		//CONTIG POSITION:
		switch (flag) {
			case SOMATIC_UNKNOWN: s="UNKNOWN"; break;
			case SOMATIC_TUMOR_UNIQ: s="TUMOR_UNIQ"; break;
			case SOMATIC_NORMAL_UNIQ: s="NORMAL_UNIQ"; break;
			case SOMATIC_TUMOR_NORMAL_SHARE: s="TUMOR_NORMAL_SHARE"; break;
			default:
				break;
		}
	}
};

int calGT_LRS_CORE(int not_supp_read_n, int supp_read_n);
void get_SV_quality_score_core(int sup_ref_number, int sup_alt_number, int sup_both_number, bool SVLEN_is_Imprecise, int final_GT, int &QUAL_INT, int &GQ_INT);
void getGT_string(int final_genotype, std::string & GT_STR);

struct NOVA_SV_FINAL_RST_item{

public:
	//new a blank nodes
	NOVA_SV_FINAL_RST_item(){
		chr_ID = 0; st_pos = 0; final_genotype = 0;
		suggenst_SV_length = 0;	SV_length = 0;	endPos = 0;
		will_be_output_to_vcf = false;
		dis_SV_len_suggset_length = 0;
		cigar_idx_bg = 0;
		cigar_idx_ed = 0;
		contig_st_pos_in_region = 0;
		region_ref_global_position = 0;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		is_BND = false;
		BND_BP_is_polishing = false;
		SVLEN_is_Imprecise = false;
		force_ALU = false;
		isVNTRcallerRst = false;
		QUAL_presetVNTR = 0;
		passVNTR_depthFilter = true;
		is_LRS_SV = false;
	}

	static inline int cmp_by_position(const NOVA_SV_FINAL_RST_item &a, const NOVA_SV_FINAL_RST_item &b){
		if(a.chr_ID != b.chr_ID)		return a.chr_ID < b.chr_ID;
		else if(a.st_pos != b.st_pos)	return a.st_pos < b.st_pos;
		else 							return a.isVNTRcallerRst < b.isVNTRcallerRst;
	}

	static inline int cmp_by_position_TGS(const NOVA_SV_FINAL_RST_item &a, const NOVA_SV_FINAL_RST_item &b){
		if(a.chr_ID != b.chr_ID)		return a.chr_ID < b.chr_ID;
		else if(a.st_pos != b.st_pos)	return a.st_pos < b.st_pos;
		else if(a.LRS_INFO.global_region_ID != b.LRS_INFO.global_region_ID)	return a.LRS_INFO.global_region_ID < b.LRS_INFO.global_region_ID;
		else 							return a.LRS_INFO.haplotype_ID < b.LRS_INFO.haplotype_ID;
	}

	bool SV_overlap(NOVA_SV_FINAL_RST_item &B){
		return (chr_ID == B.chr_ID && st_pos - 50  < B.endPos && endPos > B.st_pos - 50);
	}

	static inline int cmp_by_position_INV(const NOVA_SV_FINAL_RST_item &a, const NOVA_SV_FINAL_RST_item &b){
		if(a.chr_ID != b.chr_ID)
			return a.chr_ID < b.chr_ID;
		int MIN_A_POS = MIN(a.st_pos,a.endPos);
		int MIN_B_POS = MIN(b.st_pos,b.endPos);
		return MIN_A_POS < MIN_B_POS;
	}

	//Store INS/DEL/INV
	NOVA_SV_FINAL_RST_item(int chr_ID_, int st_pos_, const char * SV_type_,const char * ref_,const char * alt_,
			int suggest_sv_len, uint8_t * contig_, int contig_len, uint32_t *cigar_, int cigar_num, int cigar_idx_bg_,
			int cigar_idx_ed_, int contig_st_pos_in_ref_, int region_ref_global_position_){
		//basic informations
		chr_ID = chr_ID_; st_pos = st_pos_; final_genotype = 0;
		ref.clear(); alt.clear(); SV_type.clear();
		SV_type.append(SV_type_);
		SV_name.append("TODO_NO_NAME");
		ref.append(ref_);
		alt.append(alt_);
		if(SV_type.compare("INV") == 0){
			SV_length = suggest_sv_len;
			endPos = st_pos + SV_length;
		}else{
			SV_length = alt.size() - ref.size();
			endPos = st_pos + ref.size();
		}

		//additional informations
		suggenst_SV_length = suggest_sv_len;
		will_be_output_to_vcf = false;
		dis_SV_len_suggset_length = ABS_U(SV_length, suggenst_SV_length);
		contig.resize(contig_len);
		for(int i = 0; i < contig_len; i++) contig[i] = "ACGTN"[contig_[i]];
		cigar.clear();
		for(int i = 0; i < cigar_num; i++)
			cigar.emplace_back(cigar_[i]);
		cigar_idx_bg = cigar_idx_bg_;
		cigar_idx_ed = cigar_idx_ed_;
		contig_st_pos_in_region = contig_st_pos_in_ref_;
		region_ref_global_position = region_ref_global_position_;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		region_support_number[2][0] = 0;
		region_support_number[2][1] = 0;
		region_support_number[2][2] = 0;
		is_BND = false;
		isVNTRcallerRst = false;
		QUAL_presetVNTR = 0;
		BND_BP_is_polishing = false;
		SVLEN_is_Imprecise = false;
		force_ALU = false;
		passVNTR_depthFilter = true;
		is_LRS_SV = false;
	}

	static void add_to_vector(std::vector<NOVA_SV_FINAL_RST_item> & v, int chr_ID_, int st_pos_,
			const char * SV_type_, const char * ref_,const char * alt_,
			int suggest_sv_len, uint8_t * contig_, int contig_len,
			uint32_t *cigar, int cigar_num, int cigar_idx_bg_,
			int cigar_idx_ed_, int contig_st_pos_in_region, int region_ref_global_position_){
		//xassert(strlen(ref_) != strlen(alt_), "");
		v.emplace_back(chr_ID_, st_pos_, SV_type_, ref_, alt_,
				suggest_sv_len,  contig_, contig_len, cigar, cigar_num,
				cigar_idx_bg_, cigar_idx_ed_, contig_st_pos_in_region, region_ref_global_position_);
	}

	void setGenotype_directly(int final_genotype){
		this->final_genotype = final_genotype;
	}

	//final_genotype = 0: 0/0; final_genotype = 1: 0/1; final_genotype = 2: 1/1
	//final_genotype = 3: 0|0; final_genotype = 4: 0|1; final_genotype = 5: 1|0
	//final_genotype = 6: 1|1;
	//if is_first == true: 0/1 --> 0|1
	//if is_first == false: 0/1 --> 1|0
	void randomly_phasing_Genotype(int random_num){
		if(final_genotype == 0)
			final_genotype = 3;
		else if(final_genotype == 1){
			bool is_first = true;
			//simple local phasing
			if(LRS_INFO.total_hap_number_in_region > 1){
				is_first = (LRS_INFO.haplotype_ID == 0);
			}else{
				is_first = random_num%2;
			}
			final_genotype = (is_first)?4:5;
		}else if(final_genotype == 2){
			final_genotype = 6;
		}else{
			xassert(0, "ERROR: randomly_phasing_Genotype");
		}
	}

	void setGenotype(bool isHOM){
		if(isHOM)
			final_genotype = 2;
		else
			final_genotype = 1;
	}

	//store BND
	NOVA_SV_FINAL_RST_item(int chr_ID_, int st_pos_, int end_pos, const char * ref_,const char * alt_, const char * BND_type,
			int region_support_numnber_read, int region_support_numnber_ref, int region_support_numnber_unknown,
			int genotype, bool BND_BP_is_polishing_){
		//basic informations
		chr_ID = chr_ID_; st_pos = st_pos_;
		final_genotype = genotype;
		ref.clear(); alt.clear(); SV_type.clear();
		SV_type.append("BND");
		SV_name.append("LOG_SV(NOT_FINAL_RESULTS)");
		ref.append(ref_);
		alt.append(alt_);
		bnd_type.append(BND_type);
		SV_length = ABS_U(st_pos,end_pos);
		this->endPos = end_pos;

		//additional informations
		suggenst_SV_length = 0;
		will_be_output_to_vcf = true;
		dis_SV_len_suggset_length = 0;
		cigar_idx_bg = 0;
		cigar_idx_ed = 0;
		contig_st_pos_in_region = 0;
		region_ref_global_position = 0;
		region_is_overlap = false;
		deletion_pass_DR_filter = true;
		is_duplication_sv = false;
		region_support_number[2][0] = region_support_numnber_read;
		region_support_number[2][1] = region_support_numnber_ref;
		region_support_number[2][2] = region_support_numnber_unknown;
		is_BND = true;
		QUAL_presetVNTR = 0;
		this->BND_BP_is_polishing = BND_BP_is_polishing_;
		SVLEN_is_Imprecise = false;
		force_ALU = false;
		isVNTRcallerRst = false;
		passVNTR_depthFilter = true;
		is_LRS_SV = false;
	}

	bool is_BND_TRA(){
		return (is_BND && (bnd_type.compare("TRA") == 0));
	}

	//return 0 if is not BND; 1 for INV_1; 2 for INV_2
	int get_BND_TYPE_INV(){
		if(is_BND && (bnd_type.compare("INV_1") == 0)) return 1;
		if(is_BND && (bnd_type.compare("INV_2") == 0)) return 2;
		return 0;
	}

	static void store_BND(std::vector<NOVA_SV_FINAL_RST_item> & v, int chr_ID_, int st_pos_, int end_pos,
			const char * ref_,const char * alt_, const char * BND_type,
			int region_support_numnber_read, int region_support_numnber_ref, int region_support_numnber_unknown,
			int genotype, bool BND_BP_is_polishing){
		v.emplace_back(chr_ID_, st_pos_, end_pos, ref_, alt_, BND_type,
				region_support_numnber_read,  region_support_numnber_ref, region_support_numnber_unknown, genotype, BND_BP_is_polishing);
	}

#define MIN_QUAL 20
#define MAX_QUAL 500
//	int static get_quality_score(int sup_ref_number, int sup_alt_number, int final_GT){
//		float RR_log_P = -sup_alt_number*5-5.3;
//		float AA_log_P = -sup_ref_number*5;
//		float RA_log_P = -0.3*(sup_ref_number + sup_alt_number)-5;
//		float qual_float = 0;
//		if(final_GT == 0){//0/0
//			qual_float = (RR_log_P - (MAX(AA_log_P,RA_log_P)))*10;
//		}else if(final_GT == 1){//0/1
//			qual_float = (RA_log_P - (MAX(RR_log_P,AA_log_P)))*10;
//		}else{
//			qual_float = (AA_log_P - (MAX(RR_log_P,RA_log_P)))*10;
//		}
//
//		return MAX(MIN_QUAL, qual_float);
//	}

	void get_SV_quality_score(int sup_ref_number, int sup_alt_number, int sup_both_number,int final_GT, int &QUAL_INT, int &GQ_INT){
		get_SV_quality_score_core(sup_ref_number, sup_alt_number, sup_both_number, SVLEN_is_Imprecise, final_GT, QUAL_INT, GQ_INT);
	}

	void calGT_LRS(){
		final_genotype = calGT_LRS_CORE(LRS_INFO.not_supp_read_n, LRS_INFO.supp_read_n);
	}

	bool LRS_need_reGT_by_SRS(){
		return (LRS_INFO.not_supp_read_n == 0 && LRS_INFO.supp_read_n < 4);
		//return (TGS_INFO.supp_read_n < 15);
	}

	bool LRS_need_filter_by_SRS(){
		//return (TGS_INFO.not_supp_read_n == 0 && TGS_INFO.supp_read_n < 4);
		return (LRS_INFO.supp_read_n == 1 && LRS_INFO.not_supp_read_n >= 1);
	}

	bool contig_not_support_by_NGS(){
		//region_support_number[2][1] : not support contig
		//region_support_number[2][0] : support contig
		return (region_support_number[2][1] >= 20*region_support_number[2][0]);
	}


	bool writeVCF_final(kstring_t *s, bam_hdr_t * header, int *SV_global_ID){
		//skip SVs
		if(final_genotype == 0 || final_genotype == 3 || is_duplication_sv){ return false; }
		//DEL
		s->l = 0;
		//show results
		sprintf(s->s + s->l, "%s\t%d\t", header->target_name[chr_ID], (st_pos < 1)?1:st_pos); s->l += strlen(s->s + s->l);
		//generate Names
		//format: GLOBALID_CHRID_TYPE_ST_ED_LEN
		int cur_SV_ID = 0;
		if(SV_global_ID != NULL){
			(*SV_global_ID) += 1;
			cur_SV_ID = (*SV_global_ID);
		}
		int NAME_SV_length = ABS(SV_length);
		if(is_BND) NAME_SV_length = 0;

		sprintf(s->s + s->l, "%d_%s_%s_%d_%d_%d\t", (cur_SV_ID), header->target_name[chr_ID], SV_type.c_str(), st_pos,
				(is_BND)?0:endPos, NAME_SV_length);
		s->l += strlen(s->s + s->l);
		//REF and ALT
		if(!is_BND && SV_length < -50000){sprintf(s->s + s->l, "%c\t<DEL>\t", ref.c_str()[0]); s->l += strlen(s->s + s->l);}
		else{				 			 sprintf(s->s + s->l, "%s\t%s\t", 	 ref.c_str(), alt.c_str()); s->l += strlen(s->s + s->l); }
		//QUAL score
		int QUAL_INT = 0;
		int GQ_INT = 0;		//simple GT for TGS reads
		if(is_LRS_SV){
			get_SV_quality_score( LRS_INFO.not_supp_read_n, LRS_INFO.supp_read_n, LRS_INFO.unknown_read_n, final_genotype, QUAL_INT, GQ_INT);
		}
		else{
			get_SV_quality_score(region_support_number[2][1],  region_support_number[2][0], region_support_number[2][2], final_genotype, QUAL_INT, GQ_INT);
			if(isVNTRcallerRst){
				QUAL_INT = GQ_INT = QUAL_presetVNTR;
			}
		}

		sprintf(s->s + s->l, "%d\t",QUAL_INT); s->l += strlen(s->s + s->l);
		//FILTERS//IMPRECISE
		if(SV_length < 50 && SV_length > -50)			sprintf(s->s + s->l, "%s\t", "PASS");
		else if(is_duplication_sv)						sprintf(s->s + s->l, "%s\t","Duplication_sv");
		else if(final_genotype == 0)					sprintf(s->s + s->l, "%s\t","LOW_DEPTH");
		else if(!deletion_pass_DR_filter)				sprintf(s->s + s->l, "%s\t","DR_not_support");
		else if(isVNTRcallerRst && QUAL_presetVNTR <=5)	sprintf(s->s + s->l, "%s\t","LOW_QUAL");
		else											sprintf(s->s + s->l, "%s\t","PASS");
		s->l += strlen(s->s + s->l);
		//INFO:type and length
		if(is_BND){
			if(!BND_BP_is_polishing){ sprintf(s->s + s->l, "IMPRECISE;"); s->l += strlen(s->s + s->l); }
			sprintf(s->s + s->l, "SVTYPE=%s;BND_TYPE=%s" , SV_type.c_str(), bnd_type.c_str()); s->l += strlen(s->s + s->l);
		}
		else{
			sprintf(s->s + s->l, "SVTYPE=%s;" , SV_type.c_str()); s->l += strlen(s->s + s->l);
			sprintf(s->s + s->l, "END=%d;SVLEN=%d" , endPos, SV_length); s->l += strlen(s->s + s->l);
			if(SVLEN_is_Imprecise){
				sprintf(s->s + s->l, ";IMPRECISE"); s->l += strlen(s->s + s->l);
			}
			if(isVNTRcallerRst){
				sprintf(s->s + s->l, ";VNTR"); s->l += strlen(s->s + s->l);
				if(!passVNTR_depthFilter){
					sprintf(s->s + s->l, ";VNTR_DepthLow"); s->l += strlen(s->s + s->l);
				}
			}
		}
		//output the contig when needed:
		if(!cigar.empty()){
			//CONTIG POSITION:
			sprintf(s->s + s->l, ";CONTIG_POSITION=%s:%d", header->target_name[chr_ID], contig_st_pos_in_region + region_ref_global_position); s->l += strlen(s->s + s->l);
			//CIGAR
			sprintf(s->s + s->l, ";CIGAR=" ); s->l += strlen(s->s + s->l);
			for (uint i = 0; i < cigar.size(); ++i){
				int type = (int)(1 + (cigar[i] & BAM_CIGAR_MASK));
				char c_type = segment_type_to_cigar_code(type);
				int length = (cigar[i] >> BAM_CIGAR_SHIFT);
				sprintf(s->s + s->l, "%d%c", length, c_type); s->l += strlen(s->s + s->l);
			}
			//CONTIG:
			sprintf(s->s + s->l, ";CONTIG=%s", contig.c_str()); s->l += strlen(s->s + s->l);
		}

		if(tumor_flag.is_SOMATIC){
			std::string tumor_flag_str; tumor_flag.get_TUMOR_FLAG_STR(tumor_flag_str);
			sprintf(s->s + s->l, ";SOMATIC=%s", tumor_flag_str.c_str()); s->l += strlen(s->s + s->l);
		}
		//end of INFO:
		sprintf(s->s + s->l, "\t" ); s->l += strlen(s->s + s->l);

		//IMPRECISE flag:
		sprintf(s->s + s->l, "GT:GQ:LR:SR:SL:IF\t"); s->l += strlen(s->s + s->l);
		//SAMPLE::GT
		std::string GT_STR;
		getGT_string(final_genotype, GT_STR);
		sprintf(s->s + s->l, "%s:", GT_STR.c_str());
		s->l += strlen(s->s + s->l);
		//sprintf(s->s + s->l, "%s:", (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1")); s->l += strlen(s->s + s->l);
		//SAMPLE::GQ
		sprintf(s->s + s->l, "%d:", GQ_INT); s->l += strlen(s->s + s->l);
		//FORMAT
		//fprintf(vcf_output, "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Reads number that supported for the reference, allele and unknown when genotyping\">\n");
		//fprintf(vcf_output, "##FORMAT=<ID=IF,Number=.,Type=Integer,Description=\"RegionID,Hap idx in Region,SV idx in HAP.\">\n");
		//fprintf(vcf_output, "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"[BY short reads]Reads number that supported for the reference, allele and unknown when genotyping\">\n");
		//fprintf(vcf_output, "##FORMAT=<ID=LR,Number=.,Type=Integer,Description=\"[BY long reads]Reads number that supported for the reference, allele and unknown when genotyping\">\n");
		if(!is_LRS_SV)
			LRS_INFO.setINFO(-1, -1, -1, -1, -1, -1, -1, -1, 0);
		else if(region_support_number[2][1] == 0 && region_support_number[2][0] == 0 && region_support_number[2][2] == 0){
			region_support_number[2][1] = -1;
			region_support_number[2][0] = -1;
			region_support_number[2][2] = -1;
		}

		//SAMPLE::LR
		sprintf(s->s + s->l, "%d,%d,%d:", LRS_INFO.not_supp_read_n,  LRS_INFO.supp_read_n, LRS_INFO.unknown_read_n); s->l += strlen(s->s + s->l);
		//SAMPLE::SR
		sprintf(s->s + s->l, "%d,%d,%d:", region_support_number[2][1],  region_support_number[2][0], region_support_number[2][2]); s->l += strlen(s->s + s->l);
		//SAMPLE::SL
		sprintf(s->s + s->l, "%d:", suggenst_SV_length); s->l += strlen(s->s + s->l);
		//SAMPLE::IF
		sprintf(s->s + s->l, "%d,%d,%d,%d,%d", LRS_INFO.global_region_ID,  LRS_INFO.haplotype_ID, LRS_INFO.total_hap_number_in_region, LRS_INFO.SV_in_HAP_ID, LRS_INFO.total_SV_number_in_HAP ); s->l += strlen(s->s + s->l);
		sprintf(s->s + s->l, "\n"); s->l += strlen(s->s + s->l);

		return true;
	}

	void printContigSeq(int suggest_st_pos, FILE * log){
		uint32_t* bam_cigar = &(cigar[0]);
		int cigar_len = cigar.size();
		int output_index = 0;
		const char * qseq = contig.c_str();
		int seq_i = 0;
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0://M
			case 7://=
			case 8://X
				for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "%c", qseq[seq_i]); output_index++;}	break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "-"); output_index++;}	break;//S, print -
			default: fprintf(log, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log, "\n");
	}

	void print_REF_sequence(int t_len, uint8_t * tseq, FILE * log){
		for(int i = 0;i < t_len; i++){
			fprintf(log, "%c","ACGTN"[tseq[i]]);
		}
		fprintf(log, "\n");
	}

	void printSimple_M2(FILE * log){
		fprintf(log, "%d\t%d\t%s\t", chr_ID, st_pos, SV_name.c_str());
		fprintf(log, "%s\t%s\t", 	 ref.c_str(), alt.c_str());
		fprintf(log, ".\t%s\tSVTYPE=%s;END=%d;SVLEN=%d\t" "GT:DP:SG\t" "%s:", (final_genotype == 0)?"LOW_DEPTH":"PASS" ,SV_type.c_str(), endPos, SV_length , (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));
		fprintf(log, "%d,%d:%d", region_support_number[2][0],  region_support_number[2][1], suggenst_SV_length);
		fprintf(log, "\n");
	}

	void printSimple(FILE * log){
		fprintf(log, "%d\t%d\t%s\t", chr_ID, st_pos, SV_name.c_str());
		fprintf(log, "%s\t%s\t", 	 ref.c_str(), alt.c_str());
		fprintf(log, ".\t%s\tSVTYPE=%s;END=%d;SVLEN=%d\t" "GT:DP:SG\t" "%s:", (final_genotype == 0)?"LOW_DEPTH":"PASS" ,SV_type.c_str(), endPos, SV_length , (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));

		fprintf(log, "%d,%d:%d", region_support_number[2][0],  region_support_number[2][1], suggenst_SV_length);
		fprintf(log, "\nother info: contig: %s; \ncigar_idx: [%d, %d]; contig in ref: %d, ref in global %d \t", contig.c_str(),
				cigar_idx_bg, cigar_idx_ed, contig_st_pos_in_region, region_ref_global_position);

		for(uint32_t c: cigar)
			fprintf(log, "%d%c", (c >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(c & BAM_CIGAR_MASK)]);
		fprintf(log, "\n");
		printContigSeq(contig_st_pos_in_region, log);
		fprintf(log, "\n");
	}

	void print_X_E_sequence(int suggest_st_pos, uint8_t * tseq, FILE * log){
		int output_index = 0;
		int seq_i = 0;
		const char * qseq = contig.c_str();
		uint32_t* bam_cigar = &(cigar[0]);
		int cigar_len = cigar.size();
		for(int i = 0; i < suggest_st_pos; i++) {fprintf(log, "-");  output_index++;}
		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0://M
			case 7://=
			case 8://X
				for(int i = 0; i < cigar_len; i++, seq_i++){ fprintf(log, "%c", (qseq[seq_i] == "ACGT"[tseq[output_index]])?'=':'X'); output_index++;}
				break;//M
			case 1:	seq_i += cigar_len;	break;//I, print nothing
			case 2:	for(int i = 0; i < cigar_len; i++) {fprintf(log, "-");  output_index++;}	break;//D, print '-'
			case 3:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "N"); output_index++;}	break;//N, print N
			case 4:	for(int i = 0; i < cigar_len; i++, seq_i++) {fprintf(log, "-"); output_index++;}	break;//S, print -
			default: fprintf(log, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		fprintf(log, "\n");
	}

	void print(FILE * log,  bam_hdr_t * header, RefHandler *refHandler){
		fprintf(log, "%s\t%d\t%s\t", header->target_name[chr_ID], st_pos, SV_name.c_str());
		fprintf(log, "%s\t%s\t", 	 ref.c_str(), alt.c_str());
		fprintf(log, ".\t%s\tSVTYPE=%s;END=%d;SVLEN=%d\t" "GT:DP:SG\t" "%s:", (final_genotype == 0)?"LOW_DEPTH":"PASS" ,
				SV_type.c_str(), endPos, SV_length , (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));

		fprintf(log, "%d,%d:%d", region_support_number[2][0],  region_support_number[2][1], suggenst_SV_length);
		fprintf(log, "\nother info: contig: %s; \ncigar_idx: [%d, %d]; contig in ref: %d, ref in global %d \t",
				contig.c_str(), cigar_idx_bg, cigar_idx_ed, contig_st_pos_in_region, region_ref_global_position);

		for(uint32_t c: cigar)
			fprintf(log, "%d%c", (c >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(c & BAM_CIGAR_MASK)]);
		fprintf(log, "\n");

		print_REF_sequence(contig_st_pos_in_region, refHandler->getRefStr(region_ref_global_position), log);

		print_REF_sequence(1000, refHandler->getRefStr(region_ref_global_position), log);
		printContigSeq(contig_st_pos_in_region, log);
		print_X_E_sequence(contig_st_pos_in_region, refHandler->getRefStr(region_ref_global_position), log);

		fprintf(log, "\n");
	}

	static void resultFilter(std::vector<NOVA_SV_FINAL_RST_item> & rst_l, int MIN_sv_len){
		int smallest_SV_dis = MAX_int32t;
		int rst_num = rst_l.size();
		//SV length filter
		for(int i = 0; i < rst_num; i++)
			if(rst_l[i].suggenst_SV_length <= -MIN_sv_len || rst_l[i].suggenst_SV_length >= MIN_sv_len)
				{ smallest_SV_dis = MIN(smallest_SV_dis, rst_l[i].dis_SV_len_suggset_length); }
		if(smallest_SV_dis == MAX_int32t) return;
		std::map<int, int> duplication_count;
		for(int i = 0; i < rst_num; i++)
			if((rst_l[i].suggenst_SV_length <= -40 || rst_l[i].suggenst_SV_length >= 40) && rst_l[i].dis_SV_len_suggset_length <= smallest_SV_dis){
				rst_l[i].will_be_output_to_vcf = true;
				std::map<int, int>::iterator find_rst = duplication_count.find(rst_l[i].suggenst_SV_length);
				if(find_rst == duplication_count.end())
					duplication_count[rst_l[i].suggenst_SV_length] = 1;
				else
					find_rst->second++;
			}
		//remove duplications
		for(auto &dup: duplication_count){
			if(dup.second > 1){
				int c_suggenst_SV_length = dup.first;
				int c_index = 0;
				for(int i = 0; i < rst_num; i++){
					if(rst_l[i].will_be_output_to_vcf && rst_l[i].suggenst_SV_length == c_suggenst_SV_length){
						if(c_index > 0)
							rst_l[i].will_be_output_to_vcf = false;
						c_index++;
					}
				}
			}
		}
		//remove other wrong results
		for(int i = 0; i < rst_num; i++)
			if(rst_l[i].will_be_output_to_vcf && ((rst_l[i].SV_length < 0 && rst_l[i].suggenst_SV_length > 0) || (rst_l[i].SV_length > 0 && rst_l[i].suggenst_SV_length < 0)))
				rst_l[i].will_be_output_to_vcf = false;
	}

	//return:
	//0: not overlap with both
	//1: overlap with bp1
	//2: overlap with bp2
	//3: overlap with both

	int read_overlap_breakpoint(bam1_t *br, int region_ID, bool with_supp, int normal_read_length, int * left_clip_output){

		int breakpoint1 = st_pos;
		int breakpoint2 = endPos;
		*left_clip_output = 0;
		//check whether reads overlap with breakpoint
		int read_st_pos = -1;
		if(br->core.n_cigar <= 1)
			read_st_pos = br->core.pos;
		else{
			int clip_left, clip_right, gap_mismatch_inside; //the results of analysis alignment
			bam_aligned_analysis(br, &clip_left, &clip_right, &gap_mismatch_inside);
			read_st_pos = br->core.pos - clip_left;
			*left_clip_output = clip_left;
		}
		int min_overlap_len = 0.3*normal_read_length;
		bool read_overlap_with_breakpoint1 = (read_st_pos <= breakpoint1 && read_st_pos + normal_read_length > breakpoint1 && ((breakpoint1 - read_st_pos) > min_overlap_len));
		bool read_overlap_with_breakpoint2 = (read_st_pos <= breakpoint2 && read_st_pos + normal_read_length > breakpoint2 && ((read_st_pos + normal_read_length - breakpoint2) > min_overlap_len));

		if(0){
			if(read_overlap_with_breakpoint1)	fprintf(stderr, "overlap bp1\t"); else	fprintf(stderr, "NOTover bp1\t");
			if(read_overlap_with_breakpoint2)	fprintf(stderr, "overlap bp2\t"); else  fprintf(stderr, "NOTover bp2\t");
		}
		int overlap_mode = 0;

		if(with_supp){
			if(region_ID == 0 && read_overlap_with_breakpoint1)
				overlap_mode = read_overlap_with_breakpoint1;
			if(region_ID == 1 && read_overlap_with_breakpoint2 )
				overlap_mode = read_overlap_with_breakpoint2 * 2;
		}else
			overlap_mode = read_overlap_with_breakpoint1 + read_overlap_with_breakpoint2 * 2;

		return overlap_mode;
	}

	int get_contig_golbal_position_core(bool is_bp1){
		int used_cigar_num = cigar_idx_bg;
		if(is_bp1 == false)
			used_cigar_num = cigar_idx_ed + 1;
		//suggest contig_region_st for BP1
		int contig_region_st = 0;
		int small_indel_count = 0;

		for(int c_idx = 0; c_idx < used_cigar_num; c_idx++){
			int c_size = cigar[c_idx] >> BAM_CIGAR_SHIFT;
			int c_type = cigar[c_idx] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0:
				//contig_region_st += c_size;
				if(c_size < 20){ small_indel_count -= c_size; } else small_indel_count = 0;
				break;//M
			case 1:
				contig_region_st += c_size;
				if(c_size < 20){ small_indel_count -= c_size; } else small_indel_count = 0;
				break;//I
			case 2:
				contig_region_st -= c_size;
				if(c_size < 20){ small_indel_count += c_size; } else small_indel_count = 0;
				break;//M or D
			default: break;
			}
		}
		contig_region_st -= small_indel_count;

		return region_ref_global_position + contig_st_pos_in_region - contig_region_st;
	}

	void get_contig_global_position(int *contig_pos_bp1, int *contig_pos_bp2){
		*contig_pos_bp1 = get_contig_golbal_position_core(true);
		*contig_pos_bp2 = get_contig_golbal_position_core(false);
		fprintf(stderr, "SUG:contig:bp1 %d; SUG:contig:bp2 %d\n", *contig_pos_bp1, *contig_pos_bp2);
	}

	int get_ori_alignment_score(Genotyping_read_aligner * gra, bam1_t *br, int true_read_pos_bg,int region_ref_global_position, uint8_t * qseq_buff,
			int read_skip_left, int read_skip_right, int minMismatchQUAL, bool printLog){
		int score = 0;
		if(true){//running re-alignment
			//re-alignment and get score for both ref-region
			int read_in_ref_st_pos = true_read_pos_bg - region_ref_global_position;
			int read_in_ref_ed_pos = read_in_ref_st_pos + br->core.l_qseq;

			int qlen = br->core.l_qseq - read_skip_left - read_skip_right;
			uint8_t * qseq = qseq_buff + read_skip_left;
			uint8_t *qual_str = bam_get_qual(br) + read_skip_left;
			int tar_len = read_in_ref_ed_pos - read_in_ref_st_pos;
			uint8_t* tar_str = gra->get_ori_ref() + read_in_ref_st_pos + read_skip_left;
			int search_len = MIN(qlen, tar_len);
			if(printLog){
				fprintf(stderr, "Alignment to reference (true_read_pos_bg %5d, region_ref_global_position %5d, read_in_ref_st_pos %5d): ", true_read_pos_bg, region_ref_global_position, read_in_ref_st_pos);
				for(int i = 0; i < search_len; i++)
					fprintf(stderr, "%c", "ACGT"[tar_str[i]]);
			}
			int wrong_base = 0;
			for(int i = 0; i < search_len && wrong_base < 6; i++){
				if(tar_str[i] != qseq[i] && qual_str[i] > minMismatchQUAL)
					wrong_base ++;
			}
			if(printLog)
				fprintf(stderr,"(Ref: Simple search used, read region: read_in_ref_st_pos %d [%d-%d, len: %d], mismatch %d) \t", read_in_ref_st_pos,
					read_skip_left, search_len + read_skip_left, search_len, wrong_base);
			if(wrong_base < 6){
				int score = gra->getScoreByMismatch(search_len, wrong_base);
				return score;
			}
			gra->align_genotyping(false, qseq_buff + read_skip_left, qlen,
					read_in_ref_st_pos + read_skip_left, read_in_ref_ed_pos - read_skip_left);
			gra->adjustCIGAR();
			if(printLog){
				fprintf(stderr,"( Detail search used Ref:");
				gra->printf_alignment_detail(stderr, read_in_ref_st_pos);
				fprintf(stderr,")\t");
			}

			return gra->get_score_reach_end_of_read();
		}
		{//direct get score from cigar:
			//int score2 =  gra->getScoreByCigar_with_skip_region(br, read_skip_left, read_skip_right, qseq_buff, refHandler, minMismatchQUAL);
			//score = MAX(score, score2);
		}
		return score;
	}

	int get_contig_alignment_score_core(Genotyping_read_aligner * gra,
			bam1_t *br, uint8_t * qseq_buff, int read_in_contig_st_pos, int *read_skip_left_, int *read_skip_right_,
			int minMismatchQUAL){
		//realigned for contig regions
		int read_skip_left = 0;
		int read_skip_right = 0;
		if(read_in_contig_st_pos < 0){
			read_skip_left -= read_in_contig_st_pos;
			read_in_contig_st_pos = 0;
		}
		int read_in_contig_ed_pos = read_in_contig_st_pos + br->core.l_qseq;
		if(read_in_contig_ed_pos > (int)contig.size()){
			read_skip_right = read_in_contig_st_pos + br->core.l_qseq - contig.size();
			read_in_contig_ed_pos = contig.size();
		}

		*read_skip_left_ = read_skip_left;
		*read_skip_right_ = read_skip_right;

		int qlen = br->core.l_qseq - read_skip_left - read_skip_right;
		uint8_t * qseq = qseq_buff + read_skip_left;
		uint8_t *qual_str = bam_get_qual(br) + read_skip_left;
		int tar_len = read_in_contig_ed_pos - read_in_contig_st_pos;
		uint8_t* tar_str = gra->get_bin_contig() + read_in_contig_st_pos;
		int search_len = MIN(qlen, tar_len);
		int wrong_base = 0;
		for(int i = 0; i < search_len && wrong_base < 6; i++){
			if(tar_str[i] != qseq[i] && qual_str[i] > minMismatchQUAL)
				wrong_base ++;
		}
		if(false) fprintf(stderr,"(Contig: Simple search used, read region: [%d-%d, len: %d], mismatch %d) \t", read_skip_left,
				search_len + read_skip_left, search_len, wrong_base);
		if(wrong_base < 6){
			int score = gra->getScoreByMismatch(search_len, wrong_base);
			return score;
		}
		gra->align_genotyping(true, qseq, qlen, read_in_contig_st_pos - read_skip_left, read_in_contig_ed_pos - read_skip_right);
		gra->adjustCIGAR();
		if(false){
			fprintf(stderr,"( Detail search used Contig:");
			gra->printf_alignment_detail(stderr, read_in_contig_st_pos);
			fprintf(stderr,")\t");
		}
		return gra->get_score_reach_end_of_read();
	}


	int get_contig_alignment_score_middle(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff,
			int read_in_contig_st_pos1, int read_in_contig_st_pos2, int *read_skip_left_, int *true_used_read_in_contig,
			int *read_skip_right_, int minMismatchQUAL){
		int score = 0;
		if(read_in_contig_st_pos1 == read_in_contig_st_pos2){
			score = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos1, read_skip_left_,
					read_skip_right_, minMismatchQUAL);
			*true_used_read_in_contig = read_in_contig_st_pos1;
		}else{
			int read_contig_score1 = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos1,
					read_skip_left_, read_skip_right_, minMismatchQUAL);
			int read_skip_left_bp2, read_skip_right_bp2;
			int read_contig_score2 = get_contig_alignment_score_core(gra, br, qseq_buff, read_in_contig_st_pos2,
					&read_skip_left_bp2, &read_skip_right_bp2, minMismatchQUAL);
			if(read_contig_score1 >= read_contig_score2){
				score = read_contig_score1;
				*true_used_read_in_contig = read_in_contig_st_pos1;
			}else{
				*read_skip_left_ = read_skip_left_bp2; *read_skip_right_ = read_skip_right_bp2;
				score = read_contig_score2;
				*true_used_read_in_contig = read_in_contig_st_pos2;
			}
		}
		if(false)
			fprintf(stderr, "\tread_skip_left %d, read_skip_right % d\t",  *read_skip_left_, *read_skip_right_);
		return score;
	}


	int get_contig_alignment_scoreM2(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff,
			int read_in_contig_st_pos1, int read_in_contig_st_pos2, int *read_skip_left_, int *true_used_read_in_contig,
			int *read_skip_right_, int minMismatchQUAL){
		return get_contig_alignment_score_middle(gra, br, qseq_buff, read_in_contig_st_pos1, read_in_contig_st_pos2, read_skip_left_, true_used_read_in_contig, read_skip_right_, minMismatchQUAL);
	}

	int get_contig_alignment_score(Genotyping_read_aligner * gra, bam1_t *br, uint8_t * qseq_buff,
			int read_in_contig_st_pos1, int read_in_contig_st_pos2, int *read_skip_left_,
			int *read_skip_right_, int minMismatchQUAL){
		int true_used_read_in_contig;
		return get_contig_alignment_score_middle(gra, br, qseq_buff, read_in_contig_st_pos1, read_in_contig_st_pos2, read_skip_left_,
				&true_used_read_in_contig, read_skip_right_, minMismatchQUAL);
	}

	//get the read position of the first base and the last base in the reference
	void get_true_read_pos(bam1_t *br, int * true_read_pos_bg, int *true_read_pos_ed){
		int begin_pos = br->core.pos;
		int end_pos = br->core.pos;

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

	int get_length(){ return SV_length; }
	//this step running before genotyping
	//return false only when distance of two SVs is over 300
	bool duplication_SV_filter(NOVA_SV_FINAL_RST_item & to_compare){
		//position is nearby:
		int ABS_POS = ABS_U(st_pos, to_compare.st_pos);
		if(chr_ID != to_compare.chr_ID || ABS_POS > 300) return false;
		//filter:
		if(is_duplication_sv || to_compare.is_duplication_sv) return true;
		//SV type is same and SV length is similar
		float length_rate = (float)SV_length/to_compare.SV_length;
		if(length_rate < 0.9 || length_rate > 1.1) return true;

		//skip short SVs
		//if(SV_length < 50 && SV_length > -50)		return true;
		//if(to_compare.SV_length < 50 && to_compare.SV_length > -50)		return true;
		//skip BND SVs
		if(is_BND || to_compare.is_BND) return true;

		if(SV_length < 0){//for deletions:
			if(SV_length < to_compare.SV_length)	to_compare.is_duplication_sv = true;
			else									is_duplication_sv = true;
		}else{//insertions
			if(SV_length > to_compare.SV_length)	to_compare.is_duplication_sv = true;
			else									is_duplication_sv = true;
		}
		return true;
	}

	//return false only when distance of two SVs is over 300
	bool is_same_INV(NOVA_SV_FINAL_RST_item & to_compare){
		int MIN_POS = MIN(this->st_pos, this->endPos);
		int MIN_POS_CMP = MIN(to_compare.st_pos, to_compare.endPos);
		//position is nearby:
		int ABS_POS = ABS_U(MIN_POS, MIN_POS_CMP);

		int MAX_POS = MAX(this->st_pos, this->endPos);
		int MAX_POS_CMP = MAX(to_compare.st_pos, to_compare.endPos);
		int MAX_ABS_POS = ABS_U(MAX_POS, MAX_POS_CMP);
		if(chr_ID != to_compare.chr_ID || ABS_POS > 200 || MAX_ABS_POS > 200) return false;
		//SV type is same and SV length is similar
		float length_rate = (float)SV_length/to_compare.SV_length;
		if(length_rate < 0.7 || length_rate > 1.5) return false;
		//skip short SVs
		//if(SV_length < 50 && SV_length > -50)		return true;
		//if(to_compare.SV_length < 50 && to_compare.SV_length > -50)		return true;
		//skip BND SVs
		if(is_BND || to_compare.is_BND) return true;
		return true;
	}

	bool is_BND_IMPRECISE(){
		return is_BND && (BND_BP_is_polishing == false);
	}

	void get_INV_BG_ED_POS(int &BG, int &ED){
		BG = st_pos;
		ED = endPos;
		if(ED < BG) std::swap(BG, ED);
	}

	void INV_ADD_DEPTH(int64_t *region_support_number_to_add){
		region_support_number_to_add[0] += region_support_number[2][0];
		region_support_number_to_add[1] += region_support_number[2][1];
		region_support_number_to_add[2] += region_support_number[2][2];
	}

	void INV_SET_GT(int *region_support_number_to_store, int GT){
		final_genotype = GT;
		region_support_number[2][0] = region_support_number_to_store[0];
		region_support_number[2][1] = region_support_number_to_store[1];
		region_support_number[2][2] = region_support_number_to_store[2];
	}

	int get_chr_ID(){ return chr_ID; }

	bool SV_is_duplicated(){ return is_duplication_sv; }
	void SV_set_duplicated(){ is_duplication_sv = true; }

	float long_deletion_read_depth_filter(Bam_file *bam_f, int normal_read_length, int global_read_depth){

		//Analyzer
		int64_t total_read_number = 0;

		R_region bam_load_region;
		bam_load_region.chr_ID = chr_ID;
		bam_load_region.st_pos = st_pos + normal_read_length;
		bam_load_region.ed_pos = st_pos - SV_length - normal_read_length;

		resetRegion_ID(bam_f, &bam_load_region);
		while (bam_next(bam_f)) {
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))		continue;
			if(br->core.qual < 5)			continue;
			total_read_number ++;
		}
		float AVG_read_depth = (float)total_read_number*normal_read_length/(bam_load_region.ed_pos - bam_load_region.st_pos + 1);

		fprintf(stderr, "long_deletion_read_depth_filter: AVG read depth:  %f total_read_number %ld, global depth %d\n",
				AVG_read_depth, total_read_number, global_read_depth);
		return AVG_read_depth/global_read_depth;

	}

	int long_deletion_DR_filter(Bam_file *bam_f, int normal_read_length, int max_isize){
		RefRegion withInR1(chr_ID, st_pos - (max_isize - normal_read_length), st_pos);
		RefRegion withInR2(chr_ID, endPos, endPos + (max_isize - normal_read_length - normal_read_length));
		RefRegion withInDel(chr_ID, st_pos, st_pos + (max_isize - normal_read_length - normal_read_length));

		//only consider reads in the first region:
		R_region bam_load_region;
		bam_load_region.chr_ID = chr_ID;
		bam_load_region.st_pos = st_pos - max_isize;
		bam_load_region.ed_pos = st_pos + normal_read_length;

		//Analyzer
		int read_num_support_deletion = 0;
		int64_t total_isize_deletion = 0;
		int read_num_support_reference = 0;
		int64_t total_isize_reference = 0;


		resetRegion_ID(bam_f, &bam_load_region);
		while (bam_next(bam_f)) {
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))		continue;

			//get iSIZE
			int isize = br->core.isize;
			if(isize <= 0 || isize > 1000000) continue;

			int pos = br->core.pos + br->core.l_qseq;
			int mpos = br->core.mpos;

			if(	withInR1.pos_within_region(pos) && withInR2.pos_within_region(mpos)){//for read pairs that support the deletion
				read_num_support_deletion++;
				total_isize_deletion += isize;
			}

			if(withInR1.pos_within_region(pos) && withInDel.pos_within_region(mpos)){//for read pairs that support the reference
				read_num_support_reference++;
				total_isize_reference += isize;
			}
		}

		fprintf(stderr, "Read_num_support_deletion %d %f\n",
				read_num_support_deletion, (float)total_isize_deletion/read_num_support_deletion);
		fprintf(stderr, "Read_num_support_reference(IN DEL) %d %f\n",
				read_num_support_reference, (float)total_isize_reference/read_num_support_reference);

		if((read_num_support_deletion >= read_num_support_reference * 2))			return 2;
		else if((read_num_support_deletion * 3 >= read_num_support_reference))		return 1;
		else																		return 0;
	}

	void support_read_counting(bool print_log,  Bam_file *bam_f, Genotyping_read_aligner * gra,
			RefHandler *refHandler, RefRegion &main, RefRegion &supp, int normal_read_length){
		if(print_log){
			fprintf(stderr, "Region 1: "); main.show();
			if(!region_is_overlap) {fprintf(stderr, "Region 2: "); supp.show();}
			fprintf(stderr, "SV breakpoint1 %d SV; breakpoint2 %d\n", st_pos, endPos);
			fprintf(stderr, "The contig is %s %ld\n", contig.c_str(), contig.size());
			fprintf(stderr, "The Ref is:\n");
			uint8_t *ref = refHandler->getRefStr(region_ref_global_position);
			for(int i = 0; i < 2000; i++){
				fprintf(stderr, "%c", "ACGTNNNN"[ref[i]]);
			}
			fprintf(stderr, "\n");
		}

		uint8_t qseq_buff[1024]; int left_clip = 0;
		if(print_log)fprintf(stderr, "SV breakpoint1 %d SV; breakpoint2 %d\n", st_pos, endPos);

		//suggest contig_region_st for BP1/2
		int the_contig_position_in_ref_when_align_to_bp1_IN_SV; int the_contig_position_in_ref_when_align_to_bp2_IN_SV;
		get_contig_global_position(&the_contig_position_in_ref_when_align_to_bp1_IN_SV, &the_contig_position_in_ref_when_align_to_bp2_IN_SV);

		for(int region_ID = 0; region_ID < 3; region_ID++)
			for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
				region_support_number[region_ID][read_genotype_type] = 0;

		for(int region_ID = 0; region_ID < 2; region_ID++){ //mode == 0 for region1, mode == 1 for region2
			if(region_is_overlap == true && region_ID == 1) continue;

			R_region region;
			region.chr_ID = main.chr_ID;
			region.st_pos = ((region_ID == 0)?main.st_pos:supp.st_pos) + 1;
			region.ed_pos = ((region_ID == 0)?main.ed_pos:supp.ed_pos) + 1;

			gra->setRef(refHandler->getRefStr(region_ref_global_position), 100000, contig.c_str(), contig.size());

			if(print_log)fprintf(stderr, "Loading reads for region %d\n", region_ID + 1);

			resetRegion_ID(bam_f, &region);	//reset region
			//reference check:
			while (bam_next(bam_f)) {
				bam1_t *br = &(bam_f->_brec);
				if(bam_is_secondary(br))		continue;
				if(bam_is_supplementary(br))	continue;
				if(bam_is_duplicate(br))       continue;
				int overlap_mode = read_overlap_breakpoint(br, region_ID, !region_is_overlap, normal_read_length, &left_clip);
				if(overlap_mode == 0) {
					if(false) fprintf(stderr, "@Read name: %s Overlap_mode %d; SKIP\n", bam_qname(br), overlap_mode);
					continue;
				}

				get_bam_seq_bin(0, br->core.l_qseq, qseq_buff, br);
				//realigned for contig regions
				int read_contig_score = 0;

				int true_read_pos_bg; int true_read_pos_ed;
				get_true_read_pos(br, &true_read_pos_bg, &true_read_pos_ed);

				if(print_log) fprintf(stderr, "\ttrue_read_pos: [%d, %d]\t", true_read_pos_bg, true_read_pos_ed);

				int read_in_contig_st_pos_bp1_rbg = (true_read_pos_bg) - the_contig_position_in_ref_when_align_to_bp1_IN_SV;
				int read_in_contig_st_pos_bp2_rbg = (true_read_pos_bg) - the_contig_position_in_ref_when_align_to_bp2_IN_SV;

				int read_in_contig_st_pos_bp1_red = (true_read_pos_ed - br->core.l_qseq) - the_contig_position_in_ref_when_align_to_bp1_IN_SV;
				int read_in_contig_st_pos_bp2_red = (true_read_pos_ed - br->core.l_qseq) - the_contig_position_in_ref_when_align_to_bp2_IN_SV;

				int read_skip_left, read_skip_right;
				int minMismatchQUAL = 12;
				if(overlap_mode == 1){
					read_contig_score = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red, &read_skip_left, &read_skip_right, minMismatchQUAL);
				}else if(overlap_mode == 2){
					read_contig_score = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red, &read_skip_left, &read_skip_right, minMismatchQUAL);
				}else if(overlap_mode == 3){
					int read_contig_score1 = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red, &read_skip_left, &read_skip_right, minMismatchQUAL);
					int read_skip_left_bp2, read_skip_right_bp2;
					int read_contig_score2 = get_contig_alignment_score(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red, &read_skip_left_bp2, &read_skip_right_bp2, minMismatchQUAL);
					if(read_contig_score1 >= read_contig_score2){
						read_contig_score = read_contig_score1;
					}else{
						read_skip_left = read_skip_left_bp2; read_skip_right = read_skip_right_bp2;
						read_contig_score = read_contig_score2;
					}
				}

				int read_ori_score =
						get_ori_alignment_score(gra, br, true_read_pos_bg, region_ref_global_position,
								qseq_buff, read_skip_left, read_skip_right, minMismatchQUAL, print_log);

				if(print_log) {
					fprintf(stderr,"Read_contig_score %d, Read_ori_score is %d,diff is %d ", read_contig_score, read_ori_score, read_contig_score - read_ori_score);
					fprintf(stderr, "@Read name: %s POS: %d Overlap_mode %d;\t", bam_qname(br), br->core.pos, overlap_mode);
					if(true){
						fprintf(stderr, "@Read string: ");
						for(int i = 0; i < br->core.l_qseq; i++)
							fprintf(stderr, "%c", "ACGT"[qseq_buff[i]]);
					}

					fprintf(stderr, "\n");
				}
				//analysis:
				int minDiff = 4;
				int min_score = gra->getGenotypingMinScore(br->core.l_qseq - read_skip_left - read_skip_right);
				if(read_contig_score > read_ori_score + minDiff && read_contig_score > min_score)		region_support_number[region_ID][0]++;
				else if(read_contig_score + minDiff < read_ori_score && read_ori_score> min_score)	region_support_number[region_ID][1]++;
				else																		region_support_number[region_ID][2]++;
			}
		}
		fprintf(stderr,"Region analysis: for region1: ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[0][0], region_support_number[0][1], region_support_number[0][2]);
		fprintf(stderr,"Region analysis: for region2: ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[1][0], region_support_number[1][1], region_support_number[1][2]);
		for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
				region_support_number[2][read_genotype_type] = region_support_number[0][read_genotype_type] + region_support_number[1][read_genotype_type];

	}

	int getGT_from_support_read(int read_sup, int read_not_sup, int read_unknown){
	    //signal read number adjust for insertions and deletions
		int GT;
		if(SV_length > 0){//INS
			if(read_sup > read_not_sup * 4)     GT = 2;//1/1
			else if(read_sup*3 < read_not_sup)  GT = 0;//0/0
			else                                    GT = 1;//0/1
			if((read_sup + read_not_sup) * 3 < read_unknown)
				GT = 0;//0/0
			if(read_sup < 2)
				GT = 0;//0/0
		}else{//DEL
			if(read_sup > read_not_sup * 4)       GT = 2;//1/1
			else if(read_sup * 3 < read_not_sup)  GT = 0;//0/0
			else                                    GT = 1;//0/1
			if((read_sup + read_not_sup) * 3 < read_unknown)
				GT = 0;//0/0
			if(read_sup < 2)
				GT = 0;//0/0
		}
		return GT;
	}

	int genotyping_LRS_Using_SRS(bool print_log, int normal_read_length, Bam_file *bam_f,
			Genotyping_read_aligner * gra, RefHandler *refHandler){

		//step1: show break point region:
		//load reference:
		region_is_overlap = false;
		//load reference in break point 1:
		int edge_len = normal_read_length;
		RefRegion main(chr_ID, st_pos - 10, st_pos + edge_len);
		RefRegion supp(chr_ID, endPos - 10, endPos + edge_len);
		if(main.region_overlap(supp)){ main.Combine(supp, true); region_is_overlap = true; }

		support_read_counting(print_log, bam_f, gra, refHandler, main, supp, normal_read_length);

		int NGS_suggest_GT = getGT_from_support_read(region_support_number[2][0], region_support_number[2][1], region_support_number[2][2]);

		fprintf(stderr,"Region analysis: for region(1+2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[2][0],
				region_support_number[2][1], region_support_number[2][2]);
		fprintf(stderr, "NGS_suggest_GT is : %s", (NGS_suggest_GT == 0)?"0/0":((NGS_suggest_GT == 1)?"0/1":"1/1"));
		return NGS_suggest_GT;
	}

	void genotyping(int normal_read_length, int global_read_depth, int max_isize, Bam_file *bam_f,
			Genotyping_read_aligner * gra, RefHandler *refHandler){
		if(is_BND) return;
		bool print_log = false;
		//step1: show break point region:
		//load reference:
		region_is_overlap = false;
		//load reference in break point 1:
		int edge_len = normal_read_length;
		RefRegion main(chr_ID, st_pos - 10, st_pos + edge_len);
		RefRegion supp(chr_ID, endPos - 10, endPos + edge_len);
		if(main.region_overlap(supp)){ main.Combine(supp, true); region_is_overlap = true; }

		support_read_counting(print_log, bam_f, gra, refHandler, main, supp, normal_read_length);

		int final_genotype_old = final_genotype;
		if(force_ALU)
			region_support_number[2][1] /=2;

		final_genotype = getGT_from_support_read(region_support_number[2][0], region_support_number[2][1], region_support_number[2][2]);

		if(isVNTRcallerRst){
			//restore the final geno-type when nearly all reads support both
			if((region_support_number[2][0] + region_support_number[2][1]) * 3 < region_support_number[2][2])
				final_genotype = final_genotype_old;
		}

		fprintf(stderr,"Region analysis: for region(1+2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[2][0],
				region_support_number[2][1], region_support_number[2][2]);
		fprintf(stderr, "final_genotype : %s", (final_genotype == 0)?"0/0":((final_genotype == 1)?"0/1":"1/1"));

		//DR filter:
		deletion_pass_DR_filter = true;
		if(final_genotype > 0 && SV_length < -400){
			fprintf(stderr, "Running long deletion DR filter\n");
			bool pass_dr_filter = long_deletion_DR_filter(bam_f, normal_read_length, max_isize);
			if( 0 == pass_dr_filter){
				deletion_pass_DR_filter = false;
				fprintf(stderr, "Deletion failed the DR filter\n");
			}
		}
		else if(final_genotype == 0 && SV_length < -700){
			fprintf(stderr, "Running long deletion DR filter\n");
			bool pass_dr_filter = long_deletion_DR_filter(bam_f, normal_read_length, max_isize);
			if(pass_dr_filter > 0){
				final_genotype = pass_dr_filter;
				fprintf(stderr, "Long deletion genotype reset using DR signals\n");
			}
		}

		if(final_genotype > 0 && SV_length < -50000){
			final_genotype = 0;
		}

		//read depth filter:
		if(final_genotype > 0 && SV_length < -1000){
			float pass_RD_filter = long_deletion_read_depth_filter(bam_f, normal_read_length, global_read_depth);
			if(pass_RD_filter > 0.9){
				final_genotype = 0;
				fprintf(stderr, "Long deletion failed the read depth filter: pass_RD_filter %f\n", pass_RD_filter);
			}else{
				fprintf(stderr, "Long deletion PASS the read depth filter: pass_RD_filter %f\n", pass_RD_filter);
			}
		}
	}

	static void write_to_vcf_header(FILE *vcf_output, bam_hdr_t * bam_header){
		//BASIC PART
		time_t c_time; struct tm*p;
		time(&c_time);
		p = gmtime(&c_time);
		//date
		fprintf(vcf_output, "##fileformat=VCFv4.1\n");
		fprintf(vcf_output, "##fileDate=%d%d%d\n",1990+p->tm_year, 1 + p->tm_mon, p->tm_mday);
		fprintf(vcf_output, "##source=%s%s\n",PACKAGE_NAME, PACKAGE_VERSION);//software
		///INFO PART
		//fprintf(vcf_output, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"sample_id from dbVar submission; every call must have SAMPLE\">\n");
		fprintf(vcf_output, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
		fprintf(vcf_output, "##INFO=<ID=BND_TYPE,Number=1,Type=String,Description=\"Type of BND variant\">\n");
		fprintf(vcf_output, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
		fprintf(vcf_output, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n");
		fprintf(vcf_output, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise breakpoints or SV length\">\n");
		fprintf(vcf_output, "##INFO=<ID=VNTR,Number=0,Type=Flag,Description=\"SV is VNTR\">\n");
		fprintf(vcf_output, "##INFO=<ID=VNTR_DepthLow,Number=0,Type=Flag,Description=\"VNTR is Low Depth\">\n");
		fprintf(vcf_output, "##INFO=<ID=CONTIG,Number=1,Type=String,Description=\"the contig string\">\n");
		fprintf(vcf_output, "##INFO=<ID=CIGAR,Number=1,Type=String,Description=\"the CIGAR of contig string\">\n");
		fprintf(vcf_output, "##INFO=<ID=CONTIG_POSITION,Number=1,Type=String,Description=\"the contig alignment position in the reference\">\n");
		fprintf(vcf_output, "##INFO=<ID=SOMATIC,Number=1,Type=String,Description=\"Somatic flag, one of TUMOR_UNIQ, NORMAL_UNIQ, TUMOR_NORMAL_SHARE or UNKNOWN\">\n");

		///FORMAT PART
		fprintf(vcf_output, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		fprintf(vcf_output, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");

		fprintf(vcf_output, "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"[BY short reads]Reads number that supported for the reference, allele and unknown when genotyping.(set -1 when data is NA)\">\n");
		fprintf(vcf_output, "##FORMAT=<ID=LR,Number=.,Type=Integer,Description=\"[BY long reads]Reads number that supported for the reference, allele and unknown when genotyping.(set -1 when data is NA)\">\n");

		fprintf(vcf_output, "##FORMAT=<ID=SL,Number=.,Type=Integer,Description=\"Suggested SV length by the assembler.\">\n");
		fprintf(vcf_output, "##FORMAT=<ID=IF,Number=.,Type=Integer,Description=\"(1)RegionID, (2)Hap idx in Region, (3)total hap number in region, (4)SV idx in HAP, (5)total SV number in hap (set -1 when data is NA).\">\n");
		///FILTER PART
		fprintf(vcf_output, "##FILTER=<ID=LOW_DEPTH,Description=\"Variant with genotype [0/0]\">\n");
		fprintf(vcf_output, "##FILTER=<ID=lt50bp,Description=\"Supported variant but smaller than 50bp\">\n");
		fprintf(vcf_output, "##FILTER=<ID=DR_not_support,Description=\"Long deletion (>300) without enough Discordant read pair signal support\">\n");
		fprintf(vcf_output, "##FILTER=<ID=Duplication_sv,Description=\"When it is true, this SV is similar with others, and will be removed\">\n");
		fprintf(vcf_output, "##FILTER=<ID=LOW_QUAL,Description=\"LOW QUALITY\">\n");

		//ALT PART
		fprintf(vcf_output, "##ALT=<ID=DEL,Description=\"Deletion\">\n");
		fprintf(vcf_output, "##ALT=<ID=INS,Description=\"Insertion\">\n");
		fprintf(vcf_output, "##ALT=<ID=SNP,Description=\"SNP\">\n");
		fprintf(vcf_output, "##ALT=<ID=INDEL,Description=\"Small Insertion or Deletion\">\n");
		fprintf(vcf_output, "##ALT=<ID=TRA,Description=\"Trans-location\">\n");
		fprintf(vcf_output, "##ALT=<ID=INV,Description=\"Inversion\">\n");

		//CONTIG PART
		for(int i = 0; i < bam_header->n_targets; i++)
			fprintf(vcf_output, "##contig=<ID=%s,length=%d>\n", bam_header->target_name[i], bam_header->target_len[i]);
		fprintf(vcf_output, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample\n");
	}

	bool will_be_output_to_vcf;

	void get_sv_ST_EN_in_contig(int &SV_contig_bg, int &SV_contig_ed, int &SV_TYPE_int){
		int contig_index = 0;
		int cigar_ID = 0;
		for(;cigar_ID < cigar_idx_bg; cigar_ID++){
			int c_length =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = 	cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0://M
			case 7://=
			case 8://X
			case 1:
			case 3:
				contig_index += c_length; break;//M or I
			case 2:
			case 4:
				break; //Deletions
			default:
				fprintf(stderr, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}
		SV_contig_bg = contig_index;

		for(;cigar_ID <= cigar_idx_ed; cigar_ID++){
			int c_length =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int c_type = 	cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(c_type){
			case 0://M
			case 7://=
			case 8://X
			case 1:
			case 3:
				contig_index += c_length; break;//M or I
			case 2:	case 4:	break; //Deletions
			default: fprintf(stderr, "ERROR CIGAR  %d %d ", c_type, c_length);
			}
		}
		SV_contig_ed = contig_index;
		SV_TYPE_int = get_SV_TYPE_int();
	}

	int get_SV_TYPE_int(){
		if(SV_type.compare("INS") == 0){ return 0;}
		else if(SV_type.compare("DEL") == 0){ return 1;}
		else { return 2;}
	}

	void set_SVLEN_is_Imprecise(){
		SVLEN_is_Imprecise = true;
	}

	void set_force_ALU(){
		force_ALU = true;
	}

	void setVNTR_QUAL(float QUAL_presetVNTR){
		isVNTRcallerRst = true;;
		this->QUAL_presetVNTR = QUAL_presetVNTR;
	}
	bool is_duplication_sv;//when it is true, this SV is similar with others, and will be removed
	bool passVNTR_depthFilter;
	bool isVNTRcallerRst;
	int QUAL_presetVNTR;
	//basic
	int chr_ID;
	int st_pos;
	int SV_length;

	void set_LRS_INFO(int global_region_ID, int supp_read_n, int not_supp_read_n, int unknown_read_n, int haplotype_ID,
			int total_hap_number, int SV_in_HAP_ID, int total_SV_number_in_HAP, int contig_position){
		is_LRS_SV = true;
		LRS_INFO.setINFO(global_region_ID, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, total_hap_number, SV_in_HAP_ID, total_SV_number_in_HAP, contig_position);//info for TGS SV
	}

	bool is_same_var(NOVA_SV_FINAL_RST_item & b){
		if(st_pos == b.st_pos && ref.compare(b.ref) == 0 && alt.compare(b.alt) == 0)
			return true;
		else
			return false;
	}

	bool is_similar_var(NOVA_SV_FINAL_RST_item & b){
		if(alt.size() == 1){//del
			if(ABS_U(st_pos, b.st_pos) > 15)			return false;
			if(ABS_U(ref.size(), b.ref.size()) > 5)		return false;
			if(ABS_U(alt.size(), b.alt.size()) != 0)	return false;
			//if(false == string_homopolymer_same(ref, b.ref))return false;
			//if(false == string_homopolymer_same(alt, b.alt))return false;
		}else{//INS
			if(ABS_U(st_pos, b.st_pos) > 0)				return false;
			if(ABS_U(ref.size(), b.ref.size()) != 0)	return false;
			if(ABS_U(alt.size(), b.alt.size()) > 10)	return false;
			//if(false == string_homopolymer_same(ref, b.ref))return false;
			//if(false == string_homopolymer_same(alt, b.alt))return false;
		}
		return true;
	}

	bool string_homopolymer_same(std::string &a, std::string &b){
		char old_same_char = -1;
		uint a_i = 0, b_i = 0;
		while(true){
			if(a_i < a.size() && a[a_i] == old_same_char) { a_i++; continue; }
			if(b_i < b.size() && b[b_i] == old_same_char) { b_i++; continue; }
			if(a_i < a.size() && b_i < b.size() ){
				if(a[a_i] != b[b_i]){return false;}	else{ old_same_char = a[a_i]; a_i++; b_i++; }
			}else break;
		}
		if(a_i == a.size() && b_i == b.size()) return true;
		else 								   return false;
	}

	bool is_same_homopolymer_var(NOVA_SV_FINAL_RST_item & b){
		if(ABS_U(st_pos, b.st_pos) > 20)				return false;
		if(ABS_U(ref.size(), b.ref.size()) > 5) 		return false;
		if(ABS_U(alt.size(), b.alt.size()) > 5)			return false;
		if(false == string_homopolymer_same(ref, b.ref))return false;
		if(false == string_homopolymer_same(alt, b.alt))return false;
		return true;
	}

	bool is_same_supp_read_LRS(NOVA_SV_FINAL_RST_item & b){
		int A = ABS_U(LRS_INFO.supp_read_n, b.LRS_INFO.supp_read_n);
		int B = ABS_U(LRS_INFO.not_supp_read_n, b.LRS_INFO.not_supp_read_n);
		int C = LRS_INFO.supp_read_n + b.LRS_INFO.supp_read_n + LRS_INFO.not_supp_read_n + b.LRS_INFO.not_supp_read_n;
		bool is_Similar = ((A + B) / (float) C) < 1;

		if(is_Similar && LRS_INFO.global_region_ID != b.LRS_INFO.global_region_ID)
			return true;
		else
			return false;
	}

	bool is_same_global_region_ID_TGS(NOVA_SV_FINAL_RST_item & b){
		if(LRS_INFO.haplotype_ID != b.LRS_INFO.haplotype_ID && LRS_INFO.global_region_ID == b.LRS_INFO.global_region_ID)
			return true;
		else
			return false;
	}

	void combine_supp_read_TGS(NOVA_SV_FINAL_RST_item & b){
		LRS_INFO.supp_read_n += b.LRS_INFO.supp_read_n;
		LRS_INFO.not_supp_read_n -= b.LRS_INFO.not_supp_read_n;
		LRS_INFO.not_supp_read_n = MAX(LRS_INFO.not_supp_read_n, 0);
	}

	void set_supp_read_NGS(int region_support_numnber_ref, int region_support_numnber_read , int region_support_numnber_unknown){
		region_support_number[2][0] = region_support_numnber_read;
		region_support_number[2][1] = region_support_numnber_ref;
		region_support_number[2][2] = region_support_numnber_unknown;
	}
	NOVA_SV_TGS_INFO LRS_INFO;//info for TGS SV
	TUMOR_FLAG tumor_flag;
	std::string ref;
	std::string alt;

private:
	std::string SV_type;
	std::string SV_name;

	//variation information
	int endPos;
	std::string bnd_type;

	///supplementary informations///

	//suggest SV length, those data are generated from the analysis of positions of assembled reads
	int suggenst_SV_length;//the suggested SV length given by the assembly problems
	int dis_SV_len_suggset_length;//distance between SV length and suggest SV length

	//assembly results
	std::string contig;//the assembly results this is derived from
	std::vector<uint32_t> cigar; //the CIGAR of contig aligned to the reference
	int cigar_idx_bg;//the cigar index of this SV
	int cigar_idx_ed;//the cigar index of this SV
	int region_ref_global_position;//the global position of reference region, a reference region is the region collected signals and calling SVs in the SV calling steps
	int contig_st_pos_in_region;//the position of contig in the reference region

	//genotyping results
	bool region_is_overlap;//region around bp1 is overlapped with region around bp2
	//region include: [region1, region2 and region(1+2)]; read type include: [support ALT, support REF, and (not) support both]
	int region_support_number[3][3]; //[region number] * [support read type number]
	//final_genotype = 0: 0/0; final_genotype = 1: 0/1; final_genotype = 2: 1/1
	//final_genotype = 3: 0|0; final_genotype = 4: 0|1; final_genotype = 5: 1|0
	//final_genotype = 6: 1|1;
	int final_genotype;
	bool deletion_pass_DR_filter;
	//
	bool is_LRS_SV;

	//
	bool is_BND;
	bool BND_BP_is_polishing;

	//
	bool SVLEN_is_Imprecise;
	bool force_ALU;
};

#endif /* NOVASVGENERATEVCF_NOVASVRST_HPP_ */
