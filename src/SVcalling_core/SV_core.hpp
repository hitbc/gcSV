/*
 * SveHandler.hpp
 *
 *  Created on: 2020-4-29
 *      Author: fenghe
 */

#ifndef SRC_SIGNAL_SVEHANDLER_HPP_
#define SRC_SIGNAL_SVEHANDLER_HPP_

#include "../cpp_lib/statistics/StatsManager.hpp"
#include "../cpp_lib/Assembler/GreedyAssembler.hpp"
#include "../SVcalling_core/RefHandler.hpp"
#include "../cpp_lib/Assembler/DBGAssembler.hpp"
#include "../SVcalling_core/bam_stat.hpp"
#include "../SVcalling_core/NovaSVRst.hpp"
#include "../SVcalling_core/ReadHandler.hpp"
#include "../SVcalling_core/Contig_polishing.hpp"

extern "C"{
	extern int vcf_write_line(htsFile *fp, kstring_t *line);
	#include "../kswlib/kalloc.h"
	#include "../kswlib/ksw2.h"
	#include "../clib/desc.h"
	#include "../clib/vcf_lib.h"
}

#include <array>
#include <map>
#include <queue>

#define PRESET_UNKNOWN 0
#define PRESET_ONT_Q20 1
#define PRESET_CCS 2
#define PRESET_ASM 3
#define PRESET_ERR 4

#define OUT_MODE_UNSET 0
#define OUT_MODE_PURE_STR 1
#define OUT_MODE_VCF 2

#define MAX_WRONG_BASE 8

struct SIGNAL_PARAMETER{
	//for DR signal
	int insert_size_min; //0.01%
	int insert_size_max;//0.99%
	int insert_size_combine;    //0.80%
	int insert_region_len;
	int MaxReadLen;
	double read_depth;

	int max_del_dup_length;//the max length of deletion or duplication, when deletion or duplication longer than this value, variations will be treated as BND

	int SVE_MIN_SOLID_SCORE;
	int SVE_MIN_READ_NUM;

	int SVE_combine_min_score_step1;

	void init( BAM_STATUS *bs){
		insert_size_max = bs->maxInsertLen;
		insert_size_min = bs->minInsertLen; if(insert_size_min < 0) insert_size_min = 0;
		insert_size_combine = bs->maxInsertLen;
		MaxReadLen = bs->analysis_read_length;
		insert_region_len = insert_size_max - MaxReadLen*2;
		read_depth = bs->ave_read_depth;

		SVE_MIN_SOLID_SCORE = bs->ave_read_depth * 0.08;// 40X ----> 3; 60X -----> 4; 29X -------> 2;
		SVE_MIN_READ_NUM = bs->ave_read_depth * 0.1;//40X ---> 4; 60X ------> 6; 29X ------> 2;
		SVE_MIN_SOLID_SCORE = MAX(2, SVE_MIN_SOLID_SCORE);//at lease 2
		SVE_MIN_READ_NUM = MAX(3, SVE_MIN_READ_NUM);//at least 3

		SVE_combine_min_score_step1 = bs->ave_read_depth * 0.4; //66-> 26; 44 -> 17; 30 -> 12
		SVE_combine_min_score_step1 = MIN(SVE_combine_min_score_step1, 29);
		SVE_combine_min_score_step1 = MAX(SVE_combine_min_score_step1, 16);

		max_del_dup_length = 50000;
	}
};

struct ORI_REF_Depth_Item{
	int16_t ACGTD_num[6];//4 for a deletion ,5 for insertion
	int16_t ACGTD_num_tmp[6];//4 for a deletion ,5 for insertion at next base
	int16_t total_depth;
	int16_t cur_ab_block;
	uint8_t ref_base;
	uint8_t max_base;
	void set_base(uint8_t base, int ab_block, int depth){
		if(cur_ab_block == ab_block){
			ACGTD_num_tmp[base] = MAX(ACGTD_num_tmp[base], depth);
		}else{
			cur_ab_block = ab_block;
			ACGTD_num[base] += ACGTD_num_tmp[base];
			ACGTD_num_tmp[base] = depth;
		}
	}

	void add_tmp_data(){
		int max_base_count = 0;
		for(int base = 0; base < 6; base++){
			ACGTD_num[base] += ACGTD_num_tmp[base];
			total_depth += ACGTD_num[base];
			if(ACGTD_num[base] > max_base_count){
				max_base = base; max_base_count = ACGTD_num[base];
			}
		}
		//ref_base = ref_base_;
	}

	int event_info(){
		if(total_depth == 0)
			return 1;
		else if(max_base != ref_base)//3~8
			return 3 + max_base;
		else if(ACGTD_num[max_base] != total_depth){
			return 2;
		}
		return 0;
	}

	void print(FILE * output){
		fprintf(output, "[ Global depth: [ref %c, MAX %c]"
				" [depth: ",
				"ACGT-I"[ref_base],
				(ref_base == max_base)?'=':"ACGT-I"[max_base]);
		for(uint8_t base = 0; base < 6; base++){
			int base_depth = ACGTD_num[base];
			if(base_depth != 0) fprintf(output, "%c:%d ", "ACGT-I"[base], base_depth);
		}
		fprintf(output, " ]]\t");
	}

};

struct read_depth_counter{

private:
	struct read_depth_counter_item{
		int read_idx_bg;
		int read_idx_ed;
		int passed_read_num;
	};

	std::vector<read_depth_counter_item> depth_counter;
	int dc_region_st;
	int dc_shift_offset;//5
	int dc_total_filter_block_num;
	int max_high_read_counter_in_block;

public:
	void init(int region_length, int average_depth, int normal_read_length){
		dc_shift_offset = 5;
		dc_total_filter_block_num = (region_length >> dc_shift_offset);
		max_high_read_counter_in_block = average_depth* (0x1 << dc_shift_offset) / normal_read_length;
		depth_counter.resize(dc_total_filter_block_num);
	}
	read_depth_counter_item * get_rsf_by_pos(int pos){
		int st_idx = pos - dc_region_st;
		//left and right has additional blocks to store additional read
		int sf_idx = (st_idx >> dc_shift_offset);
		if(sf_idx < 0 || sf_idx >= dc_total_filter_block_num){
			return NULL;
		}
		return &(depth_counter[sf_idx]);
	}
	void add_read_depth_item(int read_idx, int pos){
		read_depth_counter_item * rdc = get_rsf_by_pos(pos);
		if(rdc == NULL)	 return;
		rdc->read_idx_bg = MIN(read_idx, rdc->read_idx_bg);
		rdc->read_idx_ed = MAX(read_idx, rdc->read_idx_ed);
	}
	void clear_read_depth_list(){
		for(read_depth_counter_item & rdc : depth_counter){
			rdc.read_idx_ed = 0;
			rdc.read_idx_bg = MAX_int32t;
			rdc.passed_read_num = 0;
		}
	}
	void set_dc_st(int dc_region_st_){
		dc_region_st = dc_region_st_;
	}

	float get_ave_depth(int pos_st, int pos_ed, float *max_depth){
		*max_depth = 0;
		read_depth_counter_item * rdc_st = get_rsf_by_pos(pos_st);
		read_depth_counter_item * rdc_ed = get_rsf_by_pos(pos_ed);
		if(rdc_st == NULL && rdc_ed == NULL) return 0;
		if(rdc_st == NULL) rdc_st = rdc_ed;
		if(rdc_ed == NULL) rdc_ed = rdc_st;
		int total_depth_read_num = 0;
		int max_depth_read_number = 0;
		for(read_depth_counter_item * c_rdc = rdc_st; c_rdc <= rdc_ed; c_rdc ++){
			if(c_rdc->read_idx_ed - c_rdc->read_idx_bg < 0)
				continue;
			int current_block_read_number = (c_rdc->read_idx_ed - c_rdc->read_idx_bg + 1);
			total_depth_read_num += current_block_read_number;
			max_depth_read_number = MAX(max_depth_read_number, current_block_read_number);
		}
		*max_depth = (float)max_depth_read_number/max_high_read_counter_in_block;
		return (float)total_depth_read_num/(rdc_ed - rdc_st + 1)/max_high_read_counter_in_block;
	}

	bool is_high_coverage(int pos){
		read_depth_counter_item * rdc = get_rsf_by_pos(pos);
		if(rdc == NULL)	 return false;
		return (rdc->read_idx_ed - rdc->read_idx_bg + 1) > (max_high_read_counter_in_block * 2);
	}
};

struct BND_ASS_read{
	BND_ASS_read(int read_in_ref_offset_, int read_list_index_, bool store_in_reverse){
		read_in_ref_offset = read_in_ref_offset_;
		read_list_index = read_list_index_;
		this->store_in_reverse = store_in_reverse;
	}
	int read_list_index;
	int read_in_ref_offset;
	bool store_in_reverse;
	void print(){
		fprintf(stderr, "read_list_index %d, read_in_ref_offset %d\t%c\t", read_list_index, read_in_ref_offset, store_in_reverse?'R':'F');
	}
	static inline int cmp_by_pos(const BND_ASS_read &a, const BND_ASS_read &b){
		return a.read_in_ref_offset < b.read_in_ref_offset;
	}
};

struct BND_ASS_Block{
	//input read info
	std::vector<BND_ASS_read> read_list;
	//input read data
	std::vector<std::string> reads;
	//output: read result
	std::vector<AssemblyContig> contigs;
	uint8_t *seq;
	void init(){
		seq = (uint8_t *)xcalloc(500,1);
	}

	void clear(){
		read_list.clear();
		reads.clear();
		contigs.clear();
	}

	void run_assembly(MainAssemblyHandler *am){
		std::swap(am->reads, reads);
		am->assembley();
		std::swap(am->contigs, contigs);
		std::swap(am->reads, reads);//swap back the read list
	}

	void storeReadCore(bam1_t *br){
		 int soft_left; int soft_right; int gap_mismatch_inside;
		bam_aligned_analysis(br, &soft_left, &soft_right, &gap_mismatch_inside);

	    bool mate_forward = bam_is_mate_fwd_strand(br);
	    bool read_forward = bam_is_fwd_strand(br);
	    bool store_in_reverse = false;
	    if(mate_forward == read_forward){
	    	store_in_reverse = true;
	    }
		//if(soft_left == 0 && soft_right == 0 && mate_forward != read_forward) return;
	    if(soft_left == 0 && soft_right == 0) return;
		uint8_t *qual = bam_get_seq(br);
		int read_len = br->core.l_qseq;
		if(get_bam_seq_bin(0, read_len, seq, br) == false) return;//get string failed, return
		int min_base_qual = 6;
		reads.emplace_back(); std::string &s_store = reads.back(); //string to store: string pointer type
		s_store.resize(read_len+1);
			int last_N_base = -1;
			for(int i = 0; i < read_len; i++){
				//char store_char = "ACGTN"[seq[i]];
				char store_char = store_in_reverse?"TGCANNN"[seq[read_len - i - 1]]:"ACGTNNN"[seq[i]];
				s_store[i] = store_char;//store binary to acgt
				if(qual[i] < min_base_qual){
					s_store[i] = 'N';
					if(i - last_N_base < 5)
						for(int j = last_N_base+1;j<i;j++) s_store[j] = 'N';
					last_N_base = i;
				}
			}
			s_store[read_len] = 0;
			if(false){
				fprintf(stderr, "XXX%s \n%d %d %d %d %d %d ", s_store.c_str(), soft_right,(soft_left),read_len, read_forward, mate_forward, store_in_reverse);
				fprintf(stderr, "%d %d \n", br->core.pos, br->core.mpos);
			}
		//store
		int read_in_ref_offset = br->core.pos - (store_in_reverse?(-soft_left):(soft_left));
		read_list.emplace_back(read_in_ref_offset, reads.size() - 1,store_in_reverse);
	}

	void show_all_reads(){
		fprintf(stderr, "BS show_all_reads \n\n");
		if(read_list.empty()) {
			fprintf(stderr, "No reads \n\n");
			return;
		}
		//std::sort(read_list.begin(), read_list.end(), BND_ASS_read::cmp_by_pos);
		int bg_pos = read_list[0].read_in_ref_offset;
		for(BND_ASS_read & r:read_list){
			for(int i = 0; i < r.read_in_ref_offset - bg_pos; i++)
				fprintf(stderr, "-");
			fprintf(stderr, "%s @ pos %d\t", reads[r.read_list_index].c_str(), r.read_in_ref_offset - bg_pos);
			r.print();
			fprintf(stderr, "\n");
		}
	}
};

struct TRANS_INS_info{
	AssemblyContig * contig_p ;
	bool is_R_not_L ;
	bool source_is_forward;

	int source_tid;
	int source_BP_pos;
	int source_contig_bg;
	int source_contig_ed;

	int target_tid;
	int target_BP_pos;
	int target_contig_bg;
	int target_contig_ed;

	void set(AssemblyContig * contig_p, bool source_is_forward, bool is_R_not_L, int BP_in_contig,
			int source_tid, int source_BP_pos, int target_tid, int target_BP_pos){
		this->contig_p = contig_p;
		this->source_is_forward = source_is_forward;
		this->is_R_not_L = is_R_not_L;
		this->source_tid = source_tid;
		this->source_BP_pos = source_BP_pos;
		this->target_tid = target_tid;
		this->target_BP_pos = target_BP_pos;

		int contig_len = contig_p->seq.size();
		if(is_R_not_L){ //right: S+T
			this->source_contig_bg = 0;
			this->source_contig_ed = BP_in_contig;
			this->target_contig_bg = BP_in_contig + 1;
			this->target_contig_ed = contig_len;
		}else{//left: T+S
			this->target_contig_bg = 0;
			this->target_contig_ed = BP_in_contig;
			this->source_contig_bg = BP_in_contig + 1;
			this->source_contig_ed = contig_len;
		}
	}
	void print(){
		//basic:
		fprintf(stderr, "\nContig: \n%s\n", contig_p->seq.c_str());
		fprintf(stderr, "Mode: %c source_is_forward %c\n", "LR"[is_R_not_L], "RF"[source_is_forward]);
		fprintf(stderr, "source BP: [%d:%d] @ contig: [%d~%d]\n",source_tid, source_BP_pos, source_contig_bg, source_contig_ed);
		const char * con_seq =  &(contig_p->seq[0]);
		for(int i = source_contig_bg; i < source_contig_ed; i++)
			fprintf(stderr, "%c", con_seq[i]);
		fprintf(stderr, "\n");

		fprintf(stderr, "target BP: [%d:%d] @ contig: [%d~%d]\n",target_tid, target_BP_pos, target_contig_bg, target_contig_ed);
		for(int i = target_contig_bg; i < target_contig_ed; i++)
			fprintf(stderr, "%c", con_seq[i]);
		fprintf(stderr, "\n");
	}
};

struct Tran_ins_suggest_pos{
	int ref_chr_ID;
	int SUG_CON_POS_in_ref;
	int contig_type;
	int contig_is_revered_to_ref;
	int contig_ID;
	int repeat;

	void print(){
		fprintf(stderr, "[%d~%d] %c %c @ contig_ID %d repeat %d\t", ref_chr_ID, SUG_CON_POS_in_ref, "NLR"[contig_type], "FR"[contig_is_revered_to_ref],contig_ID, repeat);
	}

	Tran_ins_suggest_pos(	int ref_chr_ID,	int SUG_CON_POS_in_ref,	int contig_type,int contig_is_revered_to_ref,int contig_ID){
		this->ref_chr_ID = ref_chr_ID;
		this->SUG_CON_POS_in_ref = SUG_CON_POS_in_ref;
		this->contig_type = contig_type;
		this->contig_is_revered_to_ref =contig_is_revered_to_ref;
		this->contig_ID = contig_ID;
		repeat = 1;
	}

	static inline int cmp_by_ref_position(const Tran_ins_suggest_pos &a, const Tran_ins_suggest_pos &b){
		if(a.ref_chr_ID == b.ref_chr_ID)	return a.SUG_CON_POS_in_ref < b.SUG_CON_POS_in_ref;
		else								return a.ref_chr_ID < b.ref_chr_ID;
	}

	bool is_same(Tran_ins_suggest_pos & c){
		return (this->ref_chr_ID == c.ref_chr_ID && this->SUG_CON_POS_in_ref == c.SUG_CON_POS_in_ref && this->contig_ID == c.contig_ID);
	}

	bool is_paired(Tran_ins_suggest_pos & c){
		if(this->ref_chr_ID != c.ref_chr_ID) return false;
		if(this->contig_is_revered_to_ref != c.contig_is_revered_to_ref) return false;
		if(this->contig_type == c.contig_type) return false;
		int64_t pos_diff = ABS_U(this->SUG_CON_POS_in_ref, c.SUG_CON_POS_in_ref);
		if(pos_diff < 20) return false;
		if(pos_diff > 10000) return false;
		return true;
	}
};
void debug_code_(std::vector<NOVA_SV_FINAL_RST_item> & result_sv_l, const char * vcf_fn, RefHandler *ref);

#define TYPE_NUM 4

#define LRS_SIG_INS 0
#define LRS_SIG_DEL 1
#define LRS_SIG_INV 2
#define LRS_SIG_TRA 3
#define LRS_SIG_CLIP_LEFT 4
#define LRS_SIG_CLIP_RIGHT 5
#define LRS_SIG_X 6

struct LRS_SIG{
	int8_t type;
	int chrID;
	int64_t POS;
	int64_t END;
	//is used
	bool isUsed;
	int support_read_Num;

	//return true when similar
	//return false when NOT
	bool similar_SIG_Normal(LRS_SIG &to_compare, int sig_nearby){
		//position is nearby:
		int ABS_P1 = ABS_U(POS, to_compare.POS);
		int ABS_P2 = ABS_U(END, to_compare.END);
		int ABS_P3 = ABS_U(POS, to_compare.END);
		int ABS_P4 = ABS_U(END, to_compare.POS);
		//R1 nearby R2:
		if(chrID == to_compare.chrID && (ABS_P1 < sig_nearby || ABS_P2 < sig_nearby || ABS_P3 < sig_nearby || ABS_P4 < sig_nearby))
			return true;
		//R1 including R2
		if(chrID == to_compare.chrID && (POS <= to_compare.POS && END >= to_compare.END))
			return true;
		return false;
	}

	bool similar_SIG_SMALL(LRS_SIG &to_compare, int sig_nearby){
		//position is nearby:
		int ABS_P1 = ABS_U(POS, to_compare.POS);
		int ABS_P2 = ABS_U(END, to_compare.END);
		//R1 nearby R2:
		if(chrID == to_compare.chrID && (ABS_P1 < sig_nearby && ABS_P2 < sig_nearby))
			return true;
		return false;
	}

	bool region_overlap(LRS_SIG &to_compare, int overlap_edge){
		RefRegion A(chrID, POS - overlap_edge, END + overlap_edge);
		RefRegion B(chrID, to_compare.POS - overlap_edge, to_compare.END + overlap_edge);
		if(A.region_overlap(B))
			return true;
		return false;
	}

	bool similar_SIG_CLIP(LRS_SIG &to_compare){
		//position is nearby:
		int ABS_POS = ABS_U(POS, to_compare.POS);
		if(chrID == to_compare.chrID && type == to_compare.type && ABS_POS < 1000) return true;
		return false;
	}

	bool SIG_CLIP_pairing(LRS_SIG &to_compare){
		//position is nearby:
		int ABS_POS = ABS_U(POS, to_compare.POS);
		if(chrID == to_compare.chrID && type != to_compare.type && ABS_POS < 11000) return true;
		return false;
	}

	void store(int8_t type, int chrID, int64_t POS, int64_t END){
		this->type = type;
		this->chrID = chrID;
		this->POS = POS;
		this->END = END;
		//INIT:
		isUsed = false;
		support_read_Num = 1;
	}

	void show_simple(){
		fprintf(stderr, "LRS SIG: \t");
		fprintf(stderr, "[%d:%ld-%ld];[type %d; isUsed:%d, support_read_Num:%d]\n", chrID, POS, END, type, isUsed, support_read_Num);
	}

	static inline int cmp_by_ref_pos(const LRS_SIG &a, const LRS_SIG &b){
		if(a.chrID != b.chrID)
			return a.chrID < b.chrID;
		return a.POS < b.POS;
	}

	static inline int cmp_by_ref_pos_ignore_type(const LRS_SIG &a, const LRS_SIG &b){
		if(a.chrID != b.chrID)
			return a.chrID < b.chrID;
		return a.POS < b.POS;
	}
};

struct LRS_SA_INFO{
	int sa_idx;
	int SA_chrID;
	int SA_POS;
	std::vector<path_segment> sa_cigar_l;
	int SA_MAPQ;
	int SA_number_of_mismatch;
	int sa_in_read_position_bg = 0;
	int sa_in_read_position_ed = 0;
	int sa_direction;

	void char_cigar_to_path(const char *string_cigar, std::vector<path_segment> & cigar_l)
	{
		cigar_l.clear();
		int length = 0;
		while((*string_cigar) != 0)
		{
			char c = *string_cigar;
			string_cigar++;
			if(c <= '9' && c >= '0')
				length = (length * 10) + c - '0';
			else{
				int type = cigar_code_to_segment_type(c);
				xassert(type != CIGAR_NONE, "Wrong Cigar!\n");
				cigar_l.emplace_back();
				cigar_l.back().type = type;
				cigar_l.back().length = length;
				length = 0;
			}
		}
	}

	void store(int sa_idx, std::vector<std::string> &item_tmp, faidx_t * c_ref_idx, int total_read_len){
		this->sa_idx = sa_idx;
		//chr5,47154842,+,87624S77370M374I1749543S,1,12794
		SA_chrID =  faidx_get_chrID(c_ref_idx, item_tmp[0].c_str(), NULL, 0);
		SA_POS = (int)atoi(item_tmp[1].c_str());
		sa_direction = ((item_tmp[2][0] == '+')?1:0);
		char_cigar_to_path(item_tmp[3].c_str(), sa_cigar_l);
		SA_MAPQ = (int)atoi(item_tmp[4].c_str());
		SA_number_of_mismatch = (int)atoi(item_tmp[5].c_str());
		//analysis cigar:
		sa_in_read_position_bg = 0;
		if(!sa_cigar_l.empty() && (sa_cigar_l[0].type == CIGAR_SOFT_CLIP || sa_cigar_l[0].type == CIGAR_HARD_CLIP)){
			sa_in_read_position_bg = sa_cigar_l[0].length;
		}
		sa_in_read_position_ed = total_read_len;
		if(!sa_cigar_l.empty() && (sa_cigar_l.back().type == CIGAR_SOFT_CLIP || sa_cigar_l.back().type == CIGAR_HARD_CLIP)){
			sa_in_read_position_ed = total_read_len - sa_cigar_l.back().length;
		}
	}

	void store_main(int SA_chrID, int SA_POS, path_segment *main_cigar_st, path_segment *main_cigar_ed, int total_read_len){
		//chr5,47154842,+,87624S77370M374I1749543S,1,12794
		this->SA_chrID = SA_chrID;
		this->SA_POS = SA_POS;
		sa_direction = 1;
		for(path_segment * ps = main_cigar_st; ps < main_cigar_ed; ps++){
			sa_cigar_l.emplace_back();
			sa_cigar_l.back().length = ps->length;
			sa_cigar_l.back().type = ps->type;
		}
		SA_MAPQ = 999;
		//analysis cigar:
		sa_in_read_position_bg = 0;
		if(!sa_cigar_l.empty() && (sa_cigar_l[0].type == CIGAR_SOFT_CLIP || sa_cigar_l[0].type == CIGAR_HARD_CLIP)){
			sa_in_read_position_bg = sa_cigar_l[0].length;
		}
		sa_in_read_position_ed = total_read_len;
		if(!sa_cigar_l.empty() && (sa_cigar_l.back().type == CIGAR_SOFT_CLIP || sa_cigar_l.back().type == CIGAR_HARD_CLIP)){
			sa_in_read_position_ed = total_read_len - sa_cigar_l.back().length;
		}
	}

	void show(){
		fprintf(stderr,
				"sa_idx %d SA_chrID %d, SA_POS %d, SA_MAPQ %d SA_number_of_mismatch %d "
				"sa_in_read_position: [%d-%d, len %d,diff %d], sa_direction %d"
				"\n"
				,
				sa_idx, SA_chrID, SA_POS, SA_MAPQ, SA_number_of_mismatch, sa_in_read_position_bg, sa_in_read_position_ed,
				sa_in_read_position_ed - sa_in_read_position_bg,
				SA_POS - sa_in_read_position_bg,
				 sa_direction );
	}

    static inline int cmp_by_r_pos(const LRS_SA_INFO &a, const LRS_SA_INFO &b){
        //var basic
        return a.sa_in_read_position_bg < b.sa_in_read_position_bg;
    }

};

struct LRS_SA_IN_READ{
	int read_ID;
	int chrID;
	int pos;
	bool is_primary;
	int total_read_len;
	std::vector<LRS_SA_INFO> sa_l;
	LRS_SA_INFO sa_primary;
	void store_basic(int read_ID, int chrID, int64_t pos, int bam_is_primary, int total_read_len, path_segment *main_cigar_st, path_segment *main_cigar_ed){
		this->read_ID = read_ID;
		this->chrID = chrID;
		this->pos = pos;
		this->is_primary = bam_is_primary;
		this->total_read_len = total_read_len;
		sa_l.clear();
		sa_l.emplace_back();
		sa_l.back().store_main(chrID, pos, main_cigar_st, main_cigar_ed, total_read_len);
		if(bam_is_primary)
			sa_primary.store_main(chrID, pos, main_cigar_st, main_cigar_ed, total_read_len);
	}

	void store_SA_full(std::vector<std::string> & item_tmp, faidx_t * c_ref_idx){
		std::vector<std::string> item_tmp_1;
		for(uint sa_idx = 0; sa_idx < item_tmp.size(); sa_idx ++){
			split_string(item_tmp_1, item_tmp[sa_idx].c_str(), ",");
			sa_l.emplace_back();
			sa_l.back().store(sa_idx,item_tmp_1, c_ref_idx, total_read_len);
			if( sa_idx == 0 && false == is_primary){
				//when the current record is not primary, select the "FIRST" SA as primary
				sa_primary.store(sa_idx,item_tmp_1, c_ref_idx, total_read_len);
			}
		}
	}

	//sort
	void sort_by_read_position(){
		std::sort(sa_l.begin(), sa_l.end(), LRS_SA_INFO::cmp_by_r_pos);
	}

	void show(){
		fprintf(stderr, "READ with SA read_ID %d,chrID  %d, MAIN_POS %d,is_primary %d total_read_len %d\n ", read_ID,chrID, pos,is_primary, total_read_len);
		for(uint i = 0; i < sa_l.size(); i++){
			LRS_SA_INFO & sa = sa_l[i];
			sa.show();
		}
	}

};
struct LRS_SA_Handler{
	std::vector<std::string> item_tmp;

	std::vector<LRS_SA_IN_READ> sa_read_l;
	void store(int read_ID, int chrID, int64_t pos, bool bam_is_primary,
			path_segment *main_cigar_st, path_segment *main_cigar_ed,
			const char * SA_string, int total_read_len, faidx_t * c_ref_idx){
		//fprintf(stderr, "SA_string %s \n", SA_string);
		split_string(item_tmp, SA_string, ";");

		sa_read_l.emplace_back();
		sa_read_l.back().store_basic(read_ID, chrID, pos, bam_is_primary, total_read_len, main_cigar_st, main_cigar_ed);
		sa_read_l.back().store_SA_full(item_tmp, c_ref_idx);
		//sort by the position of reads
		sa_read_l.back().sort_by_read_position();
		//show all data
		if(false) sa_read_l.back().show();
	}
};
struct LRS_SIG_Handler{

    std::vector<LRS_SIG> cigar_sig;
    std::vector<LRS_SIG> clip_sig;
    std::vector<LRS_SIG> small_sig;
    std::vector<LRS_SIG> FINAL_sig;
    LRS_SA_Handler sa_handler;
    void clear(){
        cigar_sig.clear();
        clip_sig.clear();
        small_sig.clear();
        FINAL_sig.clear();
    }

    void show(){
		fprintf(stderr, "CIGAR_sig show:\n");   for(LRS_SIG &s : cigar_sig) s.show_simple();
		fprintf(stderr, "CLIP_sig show:\n");    for(LRS_SIG &s : clip_sig)  s.show_simple();
		fprintf(stderr, "SMALL_sig show:\n");   for(LRS_SIG &s : small_sig) s.show_simple();
    }
};

struct FC_LRS_REGIONS_HANDLER{
	#define HOMO_REF_MIN_SCORE 0.95

	bool region_is_loaded;
	std::vector<R_region> region_l;

	FC_LRS_REGIONS_HANDLER(){
		region_is_loaded = false;
	}

	void load_region_list(const char * region_bed_in, faidx_t * c_ref_idx){
		region_is_loaded = true;
		region_l.clear();
		//load region list[BED file]
		std::vector<std::string> load_line;
		std::vector<std::string> item_value;
		load_string_list_from_file(region_bed_in, load_line);
		size_t line_num = load_line.size();//skip the header line
		//skip the head line
		for(uint64_t i = 0; i < line_num; i++){
			split_string(item_value, load_line[i].c_str(), "\t");//skip the header line
			if(item_value.size() < 3) continue;
			region_l.emplace_back();
			int chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			if(chrID <= 23 && chrID >= 0){
				region_l.back().chr_ID =chrID;
				region_l.back().st_pos = atoi(item_value[1].c_str());
				region_l.back().ed_pos = atoi(item_value[2].c_str());
			}
		}
	}

	void store_region_to_sig(RefRegion* rr, LRS_SIG_Handler &lrs_sig_h){
		lrs_sig_h.clear();
		for(R_region & r :region_l){
			if(r.chr_ID == rr->chr_ID && r.st_pos >= rr->st_pos && r.ed_pos <= rr->ed_pos){
				lrs_sig_h.FINAL_sig.emplace_back();
				lrs_sig_h.FINAL_sig.back().store(LRS_SIG_X, r.chr_ID, r.st_pos, r.ed_pos);
				lrs_sig_h.FINAL_sig.back().support_read_Num = 1000;//BIG enough number
			}
		}
	}
};

struct SUPPORT_READ_LIST_ITEM{
	int read_id;
	int share_target_num;
	bool is_best_target;
	SUPPORT_READ_LIST_ITEM(int read_id){
		this->read_id = read_id;
		share_target_num = 1;
		is_best_target = true;
	}
	void show(){ fprintf(stderr, "%d:%d\t", read_id, share_target_num); }
	static inline int cmp_by_readid(const SUPPORT_READ_LIST_ITEM &a, const SUPPORT_READ_LIST_ITEM &b){ return a.read_id < b.read_id; }
};

struct LONG_BRANCH_CHECK{
	int read_id;
	int target_id;
	int query_total_node;
	int targt_total_node;
	bool query_is_absorb;
	int begin_UID;
	void set(int read_id, int target_id, int query_total_node, int targt_total_node, bool query_is_absorb, int begin_UID){
		this->read_id = read_id;
		this->target_id = target_id;
		this->query_total_node = query_total_node;
		this->targt_total_node = targt_total_node;
		this->query_is_absorb = query_is_absorb;
		this->begin_UID = begin_UID;
	}

	bool is_similar(LONG_BRANCH_CHECK & l){
		if(begin_UID == l.begin_UID && query_total_node == l.query_total_node && targt_total_node == l.targt_total_node)
			return true;
		return false;
	}

	void show(){
		fprintf(stderr, "Read_id %d, target_id %d, query_total_node %d, targt_total_node %d query_is_absorb %d, begin_UID %d\n",
				read_id, target_id, query_total_node, targt_total_node, query_is_absorb, begin_UID);
	}

	static inline int cmp_by_uid(const LONG_BRANCH_CHECK &a, const LONG_BRANCH_CHECK &b){
		if(a.query_is_absorb != b.query_is_absorb)
			return a.query_is_absorb < b.query_is_absorb;
		if(a.begin_UID != b.begin_UID)
			return a.begin_UID < b.begin_UID;
		if(a.query_total_node != b.query_total_node)
			return a.query_total_node < b.query_total_node;
		if(a.targt_total_node != b.targt_total_node)
			return a.targt_total_node < b.targt_total_node;
		return a.read_id < b.read_id;
	}
};

struct SPECIAL_DETECT_BRANCH{
	int query_total_node;
	int targt_total_node;
	int begin_UID;

	void set(int query_total_node, int targt_total_node, int begin_UID){
		this->query_total_node = query_total_node;
		this->targt_total_node = targt_total_node;
		this->begin_UID = begin_UID;
	}
	bool same_as(int query_total_node, int targt_total_node, int begin_UID){
		return (begin_UID == this->begin_UID &&
				((query_total_node == this->query_total_node && targt_total_node == this->targt_total_node) ||
				( targt_total_node == this->query_total_node && query_total_node == this->targt_total_node) )
		);
	}
};

struct SPECIAL_BRANCH_HANDLER{
	std::vector<LONG_BRANCH_CHECK> long_branch_check_list;
	std::vector <SPECIAL_DETECT_BRANCH> special_branch_check_list;
	bool need_reclassify;
	void clear(){
		long_branch_check_list.clear();
		need_reclassify = false;
	}

	void special_branch_check(bool print_log, int minDepth_LRS);
};

//SDP index
struct TARGET_path_and_index{
	std::vector<ASS_TARGET_ITEM> target_link_path;//target path, a link table
	std::vector<std::vector<int>> target_index;//index
	int target_total_len = 0;

	std::vector<SUPPORT_READ_LIST_ITEM> support_read_l;//support read list:

	int the_bg_node_id;
	//match result:
	int max_score_final = -1;
	int max_score_idx = -1;
	std::vector<SDP_MATCH_ITEM> match_list;
	//SDP link result:
	std::vector<int> SDP_link_path;
	//
	std::string contig;

	int target_id;

	int clip_target_mode;

	void get_the_contig_string(bool print_log, std::vector<UNITIG_INFO_item> &unitig_l,
			std::vector<std::string> & read_string);
	void show_the_final_result(bool print_log, int cur_target_ID, std::vector<UNITIG_INFO_item> &unitig_l);
	void show_read_list();
	void show_the_target_path();
	void show_total_length(std::vector<UNITIG_INFO_item> &unitig_l);
	void build(std::vector<READ_PATH_COMBINE_ITEM> & target_seq, int read_id, std::vector<UNITIG_INFO_item> & unitig_l,
			int target_id, int clip_target_mode);
	void search_match(std::vector<READ_PATH_COMBINE_ITEM> &query, std::vector<SDP_MATCH_ITEM> &SDP_match_list);
	void SDP(bool print_log, std::vector<READ_PATH_COMBINE_ITEM> &query, bool LRS_full_begin, bool LRS_full_end, int read_length);
	void showSDP_result(READ_PATH_COMBINE_ITEM * q_p, int read_id, int cur_target_ID);
	bool true_branch_check_is_same_run_length(bool print_log, std::string &a, std::string &b );
	bool node_is_branch(bool print_log,	bool is_clip_read, SDP_MATCH_ITEM * r, int read_id,
			int target_bg_id, int target_ed_id, READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed, std::vector<UNITIG_INFO_item> & unitig_l,
			std::vector<LONG_BRANCH_CHECK> &long_branch_check_list, std::vector <SPECIAL_DETECT_BRANCH> &special_branch_check_list,
			int & new_long_branch_check_node);
	void target_length_count(SDP_MATCH_ITEM * r, int target_bg_id, int target_ed_id);
	void query_length_count(SDP_MATCH_ITEM * r, READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed);
	void update_one_region(bool print_log,  int read_id, int target_bg_id, int target_ed_id,
			READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed, std::vector<UNITIG_INFO_item> & unitig_l);
	bool SDP_update_and_branch_condition_check(
			bool print_log, int MIN_score_final,
			std::vector<READ_PATH_COMBINE_ITEM> &query, int read_id, bool LRS_full_begin, bool LRS_full_end, bool is_clip_read, //read info
			std::vector<UNITIG_INFO_item> & unitig_l,//background
			std::vector<LONG_BRANCH_CHECK> &long_branch_check_list, std::vector <SPECIAL_DETECT_BRANCH> &special_branch_check_list);
};

struct CONTIG_SELECT{
	int support_read_number;
	int clip_read_number;

	int clip_left_read_number;
	int clip_right_read_number;
	int full_cover_read_n;
	int target_ID;
	int contig_cmp_len;
	int contig_len;
	bool is_final_path;
	int not_support_read_number;

	void set(int support_read_number, int contig_ID,
			int clip_read_number, int clip_read_left_n, int clip_read_right_n, int full_cover_read_n,
			int not_support_read_number, int contig_cmp_len, int contig_len){
		this->support_read_number = support_read_number;
		this->target_ID = contig_ID;
		this->clip_read_number = clip_read_number;
		this->clip_left_read_number = clip_read_left_n;
		this->clip_right_read_number = clip_read_right_n;
		this->full_cover_read_n = full_cover_read_n;
		this->not_support_read_number = not_support_read_number;
		this->contig_cmp_len = contig_cmp_len;
		this->contig_len = contig_len;
		if(support_read_number >= 3) 	is_final_path = true;
		else 							is_final_path = false;
	}
	static inline int cmp_by_supp_read(const CONTIG_SELECT &a, const CONTIG_SELECT &b){
		if(a.is_final_path != b.is_final_path)
			return a.is_final_path > b.is_final_path;
		return a.support_read_number > b.support_read_number;
	}

	static inline int cmp_by_final(const CONTIG_SELECT &a, const CONTIG_SELECT &b){
		return a.is_final_path > b.is_final_path;
	}

	static inline int cmp_by_contig_len(const CONTIG_SELECT &a, const CONTIG_SELECT &b){
		return a.contig_cmp_len > b.contig_cmp_len;
	}

	void show(bool with_endle){
		fprintf(stderr, "CONTIG_SELECT: full_or_clip_support_read_number %d [full cover %d], clip_read_number %d [left %d, right %d],"
				" target_ID %d contig_cmp_len %d contig_len %d is_final_path %d",
				support_read_number,full_cover_read_n, clip_read_number, clip_left_read_number, clip_right_read_number,
				target_ID, contig_cmp_len, contig_len, is_final_path);
		if(with_endle) fprintf(stderr, "\n");
	}

	bool is_long_range_contig()		{ return (clip_read_number*2 > support_read_number); }
	bool is_full_cover_contig(int MIN_READ_N)		{ return (full_cover_read_n >= MIN_READ_N || (clip_left_read_number >= MIN_READ_N && clip_right_read_number >= MIN_READ_N));	}
	bool is_left_clip_contig(int MIN_READ_N)		{ 	return (clip_left_read_number >= MIN_READ_N); }
	bool is_right_clip_contig(int MIN_READ_N)		{	return (clip_right_read_number >= MIN_READ_N);}
	bool is_single_support_contig()	{ return (support_read_number == 1); }
};

struct CONTIG_SELECTER{
	CONTIG_SELECTER(){
		use_contig_num = 0;
	}
	std::vector<CONTIG_SELECT> select_list;
	int use_contig_num;

	void contig_basic_count(bool print_log, std::vector<TARGET_path_and_index> &target_path_idx_l, std::vector<READ_MATCH_NUM> & rml);
	bool select_contig(bool print_log, int preset_platform);
	void get_haplotype_support_read_n(int haplotype_ID, int &supp_read_n, int &not_supp_read_n){
		supp_read_n = select_list[haplotype_ID].support_read_number;
		not_supp_read_n = 0;
		if(use_contig_num > 1) not_supp_read_n = (haplotype_ID == 0)?select_list[1].support_read_number:select_list[0].support_read_number;
	}
};

struct LONG_SV_CONTIG_REF_MATCH{
	int contig_position;
	int ref_position;
	int position_diff;
	void set(int contig_position, int ref_position){
		this->contig_position = contig_position;
		this->ref_position = ref_position;
		position_diff = ref_position - contig_position;
	}
	static inline int cmp_by_position_diff(const LONG_SV_CONTIG_REF_MATCH &a, const LONG_SV_CONTIG_REF_MATCH &b){
			return a.position_diff < b.position_diff;
	}
};

struct LONG_SV_MATCH_CLUSTER{
	int count;
	int position_diff;

	int min_ref_position = MAX_int32t;
	int max_ref_position = -1000000;
	int min_contig_position = MAX_int32t;
	int max_contig_position = -1000000;

	void set(int count, int position_diff, 	int min_ref_position, int max_ref_position, int min_contig_position, int max_contig_position){
		this->count = count;
		this->position_diff = position_diff;
		this->min_ref_position = min_ref_position;
		this->max_ref_position = max_ref_position;
		this->min_contig_position = min_contig_position;
		this->max_contig_position = max_contig_position;
	}
	static inline int cmp_by_diff_position(const LONG_SV_MATCH_CLUSTER &a, const LONG_SV_MATCH_CLUSTER &b){
		return a.position_diff < b.position_diff;
	}
	static inline int cmp_by_ref_position(const LONG_SV_MATCH_CLUSTER &a, const LONG_SV_MATCH_CLUSTER &b){
		return a.min_ref_position < b.min_ref_position;
	}
	static inline int cmp_by_contig_position(const LONG_SV_MATCH_CLUSTER &a, const LONG_SV_MATCH_CLUSTER &b){
		return a.min_contig_position < b.min_contig_position;
	}

	int get_ave_ref_pos(){ return (min_ref_position + max_ref_position) /2 ; }
	int get_ave_con_pos(){ return (min_contig_position + max_contig_position) /2 ; }

	void show(){
		fprintf(stderr, "position_diff %d, count %d min_ref_position %d, max_ref_position %d,"
				" min_contig_position %d, max_contig_position %d\n", position_diff, count,
				min_ref_position, max_ref_position, min_contig_position, max_contig_position);
	}
};

struct Clip_contig_combine_Handler{
	struct CONTIG_KMER_MATCH_ITEM{
		int c2_pos;	int c1_pos;
		int position_diff;
		void set(int c1_pos, int c2_pos){
			this->c2_pos = c2_pos;
			this->c1_pos = c1_pos;
			position_diff = c1_pos - c2_pos;
		}
		static inline int cmp_by_position_diff(const CONTIG_KMER_MATCH_ITEM &a, const CONTIG_KMER_MATCH_ITEM &b){
			if(a.position_diff != b.position_diff)
				return a.position_diff < b.position_diff;
			return a.c1_pos < b.c1_pos;
		}
	};

	struct KMER_MATCH_MEM{
		int count;
		int position_diff;
		int c_R_pos_bg;
		//score
		int max_score;
		int max_previous_node;

		void store_final(int position_diff, int count, int c1_pos_bg){
			this->count = count;
			this->position_diff = position_diff;
			this->c_R_pos_bg = c1_pos_bg;
		}

		int get_c_L_pos(){ return c_R_pos_bg - position_diff; }
		int get_c_R_end(){ return c_R_pos_bg + count; }
		int get_c_L_end(){ return get_c_L_pos() + count; }

		static inline int cmp_by_diff_position(const KMER_MATCH_MEM &a, const KMER_MATCH_MEM &b){
			return a.position_diff < b.position_diff;
		}
		static inline int cmp_by_contig_position(const KMER_MATCH_MEM &a, const KMER_MATCH_MEM &b){
			return a.c_R_pos_bg < b.c_R_pos_bg;
		}

		void show(int id){
			fprintf(stderr, "[id:%d]:[C:%d:S:%d]:[PRE%d]:C1(CLIP@right)[%d~%d]:C2(CLIP@left)[%d~%d]:[P_diff:%d]-->\n",
					id, count, max_score, max_previous_node,
					c_R_pos_bg, c_R_pos_bg + count, get_c_L_pos(), get_c_L_pos() + count, position_diff);
		}
	};
	struct Index_Item{
		int position;
		int contig_id;
		Index_Item(int position, int contig_id){
			this->position = position;
			this->contig_id = contig_id;
		}
	};

	struct SDP_ANALYSIS_ITEM{
		int bg_SDP_cR;	int ed_SDP_cR; int CLIP_AT_right_DIS; int len_match_C_AT_right;
		double full_cover_rate_C_AT_right; int CLIP_AT_right_contig_len;
		int bg_SDP_cL;	int ed_SDP_cL; int CLIP_AT_left_DIS;  int len_match_C_AT_left;
		double full_cover_rate_C_AT_left;  int CLIP_AT_left_contig_len;
		bool condition1_absorbed;
		bool condition2_linked;
		int clip_AT_right_idx; int clip_AT_left_idx;
		int max_score;

		void ana(std::vector<KMER_MATCH_MEM> & match_MEM_l, std::vector<int> &SDP_link_path,
				int CLIP_AT_right_contig_len, int CLIP_AT_left_contig_len,
				int clip_AT_right_idx, int clip_AT_left_idx){
			this->CLIP_AT_right_contig_len = CLIP_AT_right_contig_len;
			this->CLIP_AT_left_contig_len = CLIP_AT_left_contig_len;
			bg_SDP_cR = match_MEM_l[SDP_link_path[0]].c_R_pos_bg;
			ed_SDP_cR = match_MEM_l[SDP_link_path.back()].get_c_R_end();
			bg_SDP_cL = match_MEM_l[SDP_link_path[0]].get_c_L_pos();
			ed_SDP_cL = match_MEM_l[SDP_link_path.back()].get_c_L_end();
			max_score = match_MEM_l[SDP_link_path.back()].max_score;
			this->clip_AT_right_idx = clip_AT_right_idx;
			this->clip_AT_left_idx = clip_AT_left_idx;
			//analysis:
			//CLIP@RIGHT:
			CLIP_AT_right_DIS = CLIP_AT_right_contig_len - ed_SDP_cR;
			CLIP_AT_left_DIS = bg_SDP_cL;
			len_match_C_AT_right = ed_SDP_cR - bg_SDP_cR;
			len_match_C_AT_left = ed_SDP_cL - bg_SDP_cL;
			full_cover_rate_C_AT_right = ((double)len_match_C_AT_right)/CLIP_AT_right_contig_len ;
			full_cover_rate_C_AT_left  = ((double)len_match_C_AT_left)/CLIP_AT_left_contig_len ;

			condition1_absorbed = ((full_cover_rate_C_AT_right > 0.95) || (full_cover_rate_C_AT_left > 0.95));
			condition2_linked   = (((double)CLIP_AT_right_DIS/len_match_C_AT_right < 0.1) && ((double)CLIP_AT_left_DIS/len_match_C_AT_left < 0.1));
		}

		void show(){
			fprintf(stderr,
					"CLIP@RIGHT [IDX:%d]:[%d-%d]:IN[0~%d]:[DIS:%d]:[MATCH:%d]:[COV:%f]\n"
					"CLIP@LEFT [IDX:%d]:[%d-%d] :IN[0~%d]:[DIS:%d]:[MATCH:%d]:[COV:%f]\n"
					"max_score %d condition1_absorbed %d, condition2_linked %d \n\n",
					clip_AT_right_idx, bg_SDP_cR, ed_SDP_cR, CLIP_AT_right_contig_len, CLIP_AT_right_DIS,len_match_C_AT_right, full_cover_rate_C_AT_right,
					clip_AT_left_idx,  bg_SDP_cL, ed_SDP_cL, CLIP_AT_left_contig_len,  CLIP_AT_left_DIS, len_match_C_AT_left,  full_cover_rate_C_AT_left,
					max_score, condition1_absorbed, condition2_linked);
		}
	};

	struct COMBINE_ITEM{
		std::string combine_s;
		int target_ID_part1;
		int target_ID_part2;
	};
	int left_contig_n;
	int right_contig_n;
	int total_contig_n;

	bool clip_contig_with_both_end;
	bool clip_contig_over2;

	std::vector<std::vector<CONTIG_KMER_MATCH_ITEM>> contig_match_l;
	std::vector<KMER_MATCH_MEM> match_MEM;
	std::vector<int> SDP_link_path;
	int REF_IDX_KMER_LEN;
	uint64_t MASK;
	bool show_log = false;
	std::vector<SDP_ANALYSIS_ITEM> SDP_analysis_l;

	std::unordered_map<int32_t, std::vector<Index_Item>> kmer_index;
	std::vector<uint8_t> bin_contig;
	int CLIP_AT_right_contigs_n;

	std::vector<std::string > CLIP_AT_right_contigs;
	std::vector<std::string > CLIP_AT_left_contigs;
	std::vector<std::string > full_contigs;

	std::vector<int> CLIP_AT_right_contigs_target_ID;
	std::vector<int> CLIP_AT_left_contigs_target_ID;
	std::vector<int> full_contigs_target_ID;

	//result:
	std::vector<COMBINE_ITEM> combine_contigs;
	void init();
	void clear();
	void add_right_contig_n(){ right_contig_n ++; }
	void add_left_contig_n(){ left_contig_n ++; }
	void add_total_contig_n(){ total_contig_n ++; }

	void add_CLIP_AT_right_contig(std::string & s, int target_ID){ CLIP_AT_right_contigs.emplace_back(s); 	CLIP_AT_right_contigs_target_ID.emplace_back(target_ID);}
	void add_CLIP_AT_left_contig(std::string & s, int target_ID) { CLIP_AT_left_contigs.emplace_back(s); 	CLIP_AT_left_contigs_target_ID.emplace_back(target_ID);}
	void add_full_contig(std::string & s, int target_ID) 		 { full_contigs.emplace_back(s);  			full_contigs_target_ID.emplace_back(target_ID); }

	bool is_contig_need_combine(){
		clip_contig_with_both_end = (left_contig_n > 0 && right_contig_n > 0);
		clip_contig_over2 = (total_contig_n >= 2 && left_contig_n + right_contig_n > 0);
		return (clip_contig_with_both_end || clip_contig_over2);
	}

	void build_idx_for_right_contig(std::vector<std::string > &right_contig);
	void get_contig_kmer_match(std::string &right_contig);
	void load_MEMs(std::vector<CONTIG_KMER_MATCH_ITEM> & cur_c_r_match);
	int sdp_score(KMER_MATCH_MEM & pre_try_match, KMER_MATCH_MEM & cur_handle_match);
	void SDP_run(std::vector<SDP_ANALYSIS_ITEM> &SDP_analysis_l,
			int CLIP_AT_right_contig_len, int CLIP_AT_left_contig_len,
			int clip_AT_right_idx, int clip_AT_left_idx);
	void classify_contig_ref_match(std::vector<std::string > &CLIP_AT_right_contigs,
			std::vector<std::string > &CLIP_AT_left_contigs, int clip_AT_left_idx);
	void combine();
};

//regenerate refenrece region for long or complex SVs after the LRS contig generation
struct Contig_align_region_Handler{
	Contig_align_region_Handler(){
		//output:
		suggest_sv_len = 0;
		using_local_reference = false;
		sv_handle_is_needed = false;
		ref_handler = NULL;
		REF_IDX_KMER_LEN = 0;
		MASK = 0;
		show_log = false;
		max_suggest_sv_len = 0;
		min_suggest_sv_len = 0;
		additional_load = 0;
		min_ref_position = 0;
		max_ref_position = 0;
		min_contig_position = 0;
		max_contig_position = 0;
		cluster_count = 0;

		is_left_clip_contig = false;
		is_right_clip_contig = false;//flag
		local_ref_region = NULL;
	}
	//output:
	RefRegion r_p1;
	RefRegion r_p2;
	int suggest_sv_len;

	bool using_local_reference;
	bool sv_handle_is_needed;

	bool is_left_clip_contig;
	bool is_right_clip_contig;//flag
	RefRegion *local_ref_region;

	RefHandler *ref_handler;

	int REF_IDX_KMER_LEN;
	uint64_t MASK;
	std::vector<LONG_SV_CONTIG_REF_MATCH> CONTIG_REF_MATCH_list;
	std::vector<LONG_SV_MATCH_CLUSTER> match_clusters;

	bool show_log = false;
	int max_suggest_sv_len;
	int min_suggest_sv_len;
	int additional_load;

	std::string load_ref;
	std::vector<uint8_t> bin_ref;

	void init(RefHandler *ref_){
		this->ref_handler = ref_;
		REF_IDX_KMER_LEN = ref_handler->REF_IDX_KMER_LEN;
		MASK = ref_handler->REF_IDX_KMER_MASK;
		show_log = true;
		max_suggest_sv_len = 50000;
		min_suggest_sv_len = 1000;
		additional_load = 800;

		load_ref.clear();
		bin_ref.clear();
	}

	int min_ref_position;
	int max_ref_position;
	int min_contig_position;
	int max_contig_position;
	int cluster_count;

	void get_contig_ref_match(uint8_t * buff_bin, int kmer_number);
	void init_cluster_position();
	void set_cluster_position(LONG_SV_CONTIG_REF_MATCH & m);
	void classify_contig_ref_match();
	bool generate_sv_region(int contig_length);
	void store_string(uint8_t * s, int l);
	void load_ref_string();
	bool get_the_long_range_ref_region(std::vector<uint8_t> &bin_contig);
	bool load_alignment_reference(bool is_long_range_contig, bool is_left_clip_contig, bool is_right_clip_contig,//flag
		RefRegion &local_ref_region, //input ori_reference
		std::vector<uint8_t> &ori_bin_reference, std::string & ori_ref_str,//input ori_reference
		std::vector<uint8_t> &bin_contig, //input ori_contig
		RefRegion &rst_r_p1, RefRegion &rst_r_p2,
		std::string & rst_ref_str, std::vector<uint8_t> &rst_bin_ref,
		std::string &full_ref_string//the full size reference string
		);
};

//ALU force handler
struct ALU_force_handler{
private:
	struct ALU_kmer{
		std::string s;
		int p;
		int dir;

		ALU_kmer(const char * s, int p, int dir){
			this->s.append(s);
			this->p = p;
			this->dir = dir;
		}
	};
	std::string ALU_full_string[2];
	std::vector<uint8_t> ALU_full_string_bin[2];
	int ALU_FULL_LEN;

	std::vector<ALU_kmer> ALU_kmer_list;
	int ALU_kmer_list_size;

	//INIT ALIGNER
	Contig_String_aligner aligner[2];

public:
	void init();
	bool ALU_INS_check(uint8_t *q_seq, int q_len , int ALU_DIR);
	int kmer_check(std::string & contig, int &alu_kmer_i);
	int get_suggest_BP_in_contig(int alu_kmer_i, int search_pos, int & ALU_DIR);
};

//store the SNPs into list
struct SNP_SIG_HANDLER{
	struct READ_SNP_L{
		int read_id; int st_pos; int ed_pos;
		std::vector<int> p;	std::vector<uint8_t> c;
		void show(){
			fprintf(stderr, "read_id:st_pos-ed_pos %d:%d-%d\n", read_id, st_pos, ed_pos);
			for(uint i = 0; i < p.size(); i ++)
				fprintf(stderr, "%d %c \n", p[i], "ACGT"[c[i]]);
		}
	};

	struct SNP_Read_Handler{
		std::vector<READ_SNP_L> l;
		bool with_data;
		//int total_influence_base;
		void clear(){ l.clear(); with_data = false; }
		//add signals read
		void add_sig_read(int read_id, int st_pos){
			with_data = true;
			l.emplace_back();
			l.back().read_id = read_id;
			l.back().st_pos = st_pos;
		}
		void add_sig_read_end(int end_pos){ l.back().ed_pos = end_pos; }
		void add_sig_snp(int position, uint8_t c){ l.back().p.emplace_back(position); l.back().c.emplace_back(c); }
		void show_all_data(){ for(uint i = 0; i < l.size(); i ++) l[i].show(); }
	};
	struct Link_Handler{
		struct Link_Diff_ITEM{
			int diff; int link; int n_match;
			Link_Diff_ITEM(){ diff = 0; link = 0; n_match = 0;}
		};
		std::vector<std::vector<Link_Diff_ITEM>> link_diff;
		void init(int all_read_number){
			link_diff.clear();
			for(int i = 0; i < all_read_number; i++){
				link_diff.emplace_back();
				link_diff.back().resize(all_read_number);
			}
		}

		int get_link_score(int node_i,int node_j){
			Link_Diff_ITEM & ldi = link_diff[node_i][node_j];
			return ldi.link * 3 - ldi.diff *4 + ldi.n_match;
		}

		void store_link(int node_i, int node_j, uint8_t c_i, uint8_t c_j){
			if(c_i != c_j)				link_diff[node_i][node_j].diff ++;
			else if(c_i != 5)			link_diff[node_i][node_j].link ++;
			else						link_diff[node_i][node_j].n_match ++;
		}

		int get_node_n(){ return link_diff.size(); }
		void show_all_link(){
			//show all the links
			for(int i = -1; i < (int)link_diff.size(); i++){
				fprintf(stderr, "%d:\t", i);
				for(uint j = 0; j < link_diff.size(); j++){
					if(i == -1) fprintf(stderr, "[%4d %4d %4d]\t",j,j,j);
					else 		fprintf(stderr, "[%4d %4d %4d]\t",link_diff[i][j].link,link_diff[i][j].diff,
							link_diff[i][j].n_match);
				}
				fprintf(stderr, "\n");
			}
		}
	};
	struct SNP_Position_Handler{
		struct ANA_ITEM{ int read_id; uint8_t c; ANA_ITEM(int read_id, uint8_t c){  this->read_id = read_id; this->c = c; } };
		std::unordered_map<int, std::vector<ANA_ITEM>> ana_l;
		std::map<int, uint8_t> read_id_char_map;
		int all_read_number;
		void store_all_position_SNP(std::vector<int> &used_read_list, SNP_Read_Handler & snp_read_h);
		void store_2_link_format(std::vector<int> &used_read_list, Link_Handler &link_handler, SNP_Read_Handler & snp_read_h);
	};
	struct Cluster_Handler{
		struct NODE_ITEM{ int read_id; int score; NODE_ITEM(int read_id){ this->read_id = read_id; this->score = 0; }};
		struct Cluster_ITEM{		//cluster all reads into group:
			std::vector<NODE_ITEM> nl;
			int clu_sum_score;
			int clu_low_score_node;
			void show(){
				fprintf(stderr, "Cluster: clu_sum_score %d; clu_low_score_node %d\n", clu_sum_score, clu_low_score_node);
				for(uint j = 0; j < nl.size(); j++) fprintf(stderr, "[read_id:%d score:%d]", nl[j].read_id, nl[j].score);
				fprintf(stderr, "\n");
			}
		};
		std::vector<Cluster_ITEM> clu_l;
		struct Read_Cluster_New_Add{int cid; int c_score_sum; };
		std::vector<Read_Cluster_New_Add> add_clusters;
		std::vector<int> read_final_cluster_id;

		void init(int read_number);
		void add_cluster_detection(int rid_add, Link_Handler &link_handler);
		void add_new_read_to_cluster(int rid_add, Link_Handler &link_handler);
		void calculate_sum_score_for_all_cluster();
		void run(Link_Handler &link_handler );
	};
	SNP_Read_Handler SNP_read_h;
	SNP_Position_Handler SNP_p_h;
	Cluster_Handler cluster_handler;
	Link_Handler link_handler;

	void clear(){ SNP_read_h.clear(); }
	bool is_with_data(){ return SNP_read_h.with_data; }
	void analysis(std::vector<int> &used_read_list);
	bool is_reject_read_target_pair(int read_id, std::vector<SUPPORT_READ_LIST_ITEM> &support_read_l);
};

struct SRS_REF_STRING_COMBINE_Handler{
	std::vector<uint8_t> s;
	int BP_region_st[2];int BP_region_ed[2];
	int s1_ref_st;	int s2_ref_st;
	bool s1_forward; bool s2_forward;
	int main_pos; int supp_pos;
	bool is_reverse;
	bool with_contig_supp;
	int get_len(){ if(with_contig_supp) return supp_pos- main_pos; else return 0;}
	void print_ref();
	void print_results();
	void store_ref( std::vector<uint8_t> & s1, bool s1_forward, int s1_ref_st, std::vector<uint8_t> & s2, bool s2_forward, int s2_ref_st, bool is_reverse );
	void check_and_store(int BP1, int BP2);
	bool check_and_store_TRAN_INS(int BP1, int BP2, bool is_right_part, int &source_BP_pos, int &target_BP_pos);
	void chece_and_store_M2(int BP1, int BP2);
};

struct SV_CALLING_Handler{
	struct LRS_SMALL_SIG_HANDLER{
		std::queue<int> position_queue;
		std::queue<int> influence_base_queue;
		int total_influence_base;
		void clear();
		void store_sig(bool print_log, int chrID,  int POS, int influence_base, std::vector<LRS_SIG> &SMALL_sig);
	};

	SV_CALLING_Handler(){
		//part 1: reference, it is only a reference for the refHandler
		ref_handler = NULL;
		Contig_align_region_Handler long_SV_Hanler;//todo::
		kmer_profile_handler = NULL;
		rlr = NULL;
		//part 2: bam files, the SVE handler keep and maintained the BAM files
		SRS_BAM_F = NULL;
		LRS_BAM_F = NULL;

		memset(&SRS_read, 0, sizeof(SRS_read));
		memset(&LRS_read, 0, sizeof(LRS_read));

		memset(&LRS_Tumor_read, 0, sizeof(LRS_Tumor_read));
		memset(&LRS_Normal_read, 0, sizeof(LRS_Normal_read));
		memset(&fc_lrs_handler, 0, sizeof(FC_LRS_REGIONS_HANDLER));

		cur_header = NULL;
		read_counter = 0;
		//buffs for unmapped signals handler
		UMQueryBuff = NULL;//buff used in handleUMSignal
		//buff for depth counter
		memset(&rdc, 0, sizeof(read_depth_counter));
		possibility_r = NULL;
		//part 5: parameters
		MIN_sv_len = 0;
		force_calling_ALU = false;
		random_phasing = false;
		output_small_var = false;
		memset(&sig_para, 0, sizeof(SIGNAL_PARAMETER));
		//part 6: data for break point distributions, those data used for calculating break point distribution for each signals
		MIN_ACCEPT_POSSIBILITY_DR = 0;
		MIN_ACCEPT_POSSIBILITY_SH = 0;
		MIN_ACCEPT_POSSIBILITY_UM = 0;
		START_OFFSET_UM = 0;
		//part 7.0: buffs for SVE signals
		SVE_L sve[SIG::LEN][SV::LEN]; //SVE signals list in the first combining loop

		memset(&SRS_sveVNTR, 0, sizeof(SVE_L));
		//SVE_L final; //SVE signals list in the second combining loop
		//part 7: assembly manager
		memset(&ass_block, 0, sizeof(Ass_Block));
		am = NULL;
		memset(&ca, 0, sizeof(Contig_String_aligner));
		memset(&somatic_same_ca, 0, sizeof(Contig_String_aligner));

		region_ref_global_position = 0;
		region_addition_load = 0;
		//part 8: kmer counter
		dbgHandler = NULL;

		vcf_w = NULL;
		memset(&vcfBuffStr, 0, sizeof(kstring_t));
		region_ID_global = 0;
		//genotyping buffs
		memset(&ga, 0, sizeof(Genotyping_read_aligner));
		//BND handlers
		 ;
		memset(&bnd_ass_block, 0, sizeof(BND_ASS_Block));
		memset(&ca_bnd, 0, sizeof(Contig_String_aligner));
		memset(&ca_re_locate, 0, sizeof(Contig_String_aligner));
		//SV output:
		global_SV_ID = 0;

		//TL reads
		TL_read_number = 0;

		//LRS SV-calling
		memset(&special_branch_handler, 0, sizeof(SPECIAL_BRANCH_HANDLER));
		memset(&read_path_handler, 0, sizeof(READ_PATH_HANDLER));
		memset(&contig_polsihing_handler, 0, sizeof(CONTIG_POLISHING_HANDLER));
		memset(&contig_selecter, 0, sizeof(CONTIG_SELECTER));

		alu_handler = NULL;

		LRS_Tumor_BAM_F = NULL;
		LRS_Normal_BAM_F = NULL;
		OUT_MODE_int = 0;
		FC_BED_F = NULL;
		INNER_TEST_FLAG = NULL;
		preset_platform = 0;
	}
public:
	FC_LRS_REGIONS_HANDLER fc_lrs_handler;
private:
	//part 1: reference, it is only a reference for the refHandler
	RefHandler *ref_handler;
	Contig_align_region_Handler contig_align_region_Handler;
	Clip_contig_combine_Handler clip_contig_combine_Handler;
	KMER_ERROR_PROFILE_HANDLER *kmer_profile_handler;
	RefLocalRepeat *rlr;
	std::vector<RefLocalRepeatItemR> rlr_result;
	//part 2: bam files, the SVE handler keep and maintained the BAM files
	char * SRS_BAM_F;
	char * LRS_BAM_F;
	char * LRS_Tumor_BAM_F;
	char * LRS_Normal_BAM_F;
	char * FC_BED_F;
	char * INNER_TEST_FLAG;

	BAM_handler SRS_read;
	BAM_handler LRS_read;

	BAM_handler LRS_Tumor_read;
	BAM_handler LRS_Normal_read;

	bam_hdr_t * cur_header;
	int read_counter;
	//buffs for unmapped signals handler
	uint8_t *UMQueryBuff;//buff used in handleUMSignal
	//buff for depth counter
	read_depth_counter rdc;
	std::vector<int> mismatch_position;
	//part 3: BUFFs for local alignment:
	//part 4: buff for signals combination
	float *possibility_r;
	std::vector<SVE> cmb_store_tmp;
	std::vector<int> cmb_try_list;
	//part 5: parameters
	int MIN_sv_len;
	bool force_calling_ALU;//when set true: ALU ins will be force-calling
	bool random_phasing;
	bool output_small_var;
	int preset_platform;
	int OUT_MODE_int;
	SIGNAL_PARAMETER sig_para;
	//part 6: data for break point distributions, those data used for calculating break point distribution for each signals
	std::vector<float> DR_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_DR = 0;
	std::vector<float> SH_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_SH = 0;
	std::vector<float> UM_stPos_distribution;

	float MIN_ACCEPT_POSSIBILITY_UM = 0;
	int START_OFFSET_UM = 0;
	//part 7.0: buffs for SVE signals
	SVE_L SRS_sve[SIG::LEN][SV::LEN]; //SVE signals list in the first combining loop
	SVE_L SRS_sveVNTR; //SVE signals for VNTR
	//part 7: assembly manager
	Ass_Block ass_block;//load reads and run assemblies
	//Ass_Block ass_block_r2;
	std::vector<uint8_t> bin_contig;
	MainAssemblyHandler *am;
	Contig_String_aligner ca;
	Contig_String_aligner somatic_same_ca;
	//for contig suggestion alignment positions
	std::map<int, int > suggest_st_pos_map;
	//std::vector<SUGGEST_POS_LIST_ITEM> SUGGEST_pos_list;
	std::vector<int> suggest_SV_length;
	int region_ref_global_position;
	//contig realignment
	int region_addition_load;

	//part 8: kmer counter
	MainDBGHandler1 *dbgHandler;

	std::vector<uint16_t> contig_depth;
	//AssemblyManager *am;
	//part 10: buffs used to VCF records
	FILE* vcf_w;
	kstring_t vcfBuffStr;
	std::vector<NOVA_SV_FINAL_RST_item> SV_result_TMP_buff;

	std::vector<NOVA_SV_FINAL_RST_item> result_SRS_INDEL;//large INS + DEL/ ALU
	std::vector<NOVA_SV_FINAL_RST_item> result_NGS_BND;//BND
	std::vector<NOVA_SV_FINAL_RST_item> result_BGS_INV;//INV

	std::vector<NOVA_SV_FINAL_RST_item> result_LRS_final;
	uint32_t region_ID_global = 0;
	//genotyping buffs
	Genotyping_read_aligner ga;

	//realignment bam handler
	std::vector<std::string> realn_sv_field_split;

	//BND handlers
	BND_ASS_Block bnd_ass_block;
	Contig_String_aligner ca_bnd;

	Contig_String_aligner ca_re_locate;//this aligner is used to re-locate the position of contig in trans reference
	//SV output:
	int global_SV_ID;

	//TL reads
	int TL_read_number;

	//LRS SV-calling
	SPECIAL_BRANCH_HANDLER special_branch_handler;
	std::vector<TARGET_path_and_index> target_path_idx_l;
	READ_PATH_HANDLER read_path_handler;
	std::map<int, int> updated_targets_by_current_read;
	CONTIG_POLISHING_HANDLER contig_polsihing_handler;
	CONTIG_SELECTER contig_selecter;

	//Somatic SV calling
	std::vector<TARGET_path_and_index> target_path_idx_l_Tumer;
	std::vector<TARGET_path_and_index> target_path_idx_l_Normal;
	CONTIG_SELECTER contig_selecter_Tumer;
	CONTIG_SELECTER contig_selecter_Normal;

	//ALU force handler
	ALU_force_handler *alu_handler;
	//
	SNP_SIG_HANDLER snp_sh;

private:
	/*************************************BASIC FUNCTIONs**********************************************/
	SV_CALLING_Handler (SV_CALLING_Handler &B);//SveHandler can`t be copied
	/*************************************SRS SIGNALS FUNCTIONs**********************************************/
	//S1: get signals from all reads, in this step, a read pair may be used as DR or SA or UM signals, for each read signals, a SVE will be stored //signals functions
	void SRS_set_min_accpet_possibility();
	bool is_mate_fwd(uint16_t flag) { return (!((flag & BAM_MATE_STRAND) != 0));}
	void SRS_GET_SIGNALS_FROM_READs();
	void SRS_load_realignment_read_from_file();
	void SRS_handleDRSignal(bam1_core_t *core, int middle_size);
	void SRS_storeMismatchSignals(bam1_t *br, READ_record &c_r);
	void SRS_storeClipSignals(bool isClipAtRight, uint32_t ori_pos, uint8_t read_mapq);
	void SRS_handleSASignal(READ_record& c_r, bam1_t *br);//handle a SA signal and store results
	void SRS_handleUMSignal(bam1_t *br);
	/*************************************SRS CLUSTERING FUNCTIONs**********************************************/
	//S2: combine signal from multiple reads. In this steps, signals from multiple reads will be clustered and joint together. At last, signals combined from many reads(SVEs) will form.
	void SRS_CLUSTERING_AND_COMBINE_signals();
	int  SRS_getTopPossibilityIdx(int r_min, int r_max, SVE_L & l, bool isR1, bool is_forward, float &max_poss, std::vector<float> &dis_bp_percent);//used in sve_combine_STEP1_DR
	void SRS_single_type_sve_combine(bool print_log, SVE_L & l, int min_score, SIG::T sigT, SV::T svt);
	void SRS_combine_duplication(SVE_L & l);
	void SRS_DR_SH_signal_combining(SVE_L &SH, SVE_L &DR);//S3: Joint Calling for single sv, multi-sample and multi-signal type //r1 and r2 of SH must be within r1 and r2 of DR
	/*************************************SRS ASSEMBLY FUNCTIONs**********************************************/
	void SRS_force_calling_MEI_ins(bool print_log, SVE &sv, std::vector<AssemblyContig> &contigs, int main_st_pos);
	void SRS_combine_repeat_tail_of_contigs(bool print_log, int ori_read_number, std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, std::vector<AssemblyContig> &contigs);\
	void SRS_get_BP_from_tran_based_ins(bool print_log, std::vector<Tran_ins_suggest_pos> &magic_pos_set, std::vector<AssemblyContig> &contigs, std::vector<TRANS_INS_info> &bp_l);
	void SRS_handle_tran_based_ins(bool print_log, int ori_read_number, std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, std::vector<AssemblyContig> &contigs);
	bool SRS_read_cover_repeat_filter(bool print_log, AssemblyContig &contig, uint32_t SV_size_ori,
			std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, int ori_read_number);
	/*************************************SRS ASSEMBLY FUNCTIONs**********************************************/
	bool SRS_SV_region_depth_filter(SVE & c_sve); //return whether pass assembly filter
	bool SRS_assembly_load_read(bool print_log, SVE &sv, RefRegion main, RefRegion supp);
	void LRS_assembly_load_read(bool print_log, RefRegion &main, Ass_Block & a_b, int min_LRS_read_length, BAM_handler *cur_LRS_read);
	void LRS_assembly_load_read_ASM(bool print_log, RefRegion &main, Ass_Block & a_b, int min_LRS_read_length, BAM_handler *cur_LRS_read);
	bool assembly_load_read(bool print_log, RefRegion r1, RefRegion r2, RefRegion main,
			RefRegion supp, int minQUAL, uint MAX_load_reads, bool using_LRS, bool using_SRS, int min_LRS_read_length, BAM_handler *cur_LRS_read);
	bool SRS_assembly_variations_MEI_AND_LONG_SV(SVE &sv, bool call_ALU_TL, bool &withRepeat);
	bool SRS_assembly_variations_INS_DEL_VNTR(SVE &sv, bool using_LRS, bool usingNGS);
	/*************************************SRS CALLING FUNCTIONs**********************************************/
	void SRS_get_suggention_alignment_position_list(AssemblyContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f );
	void SRS_getSuggestSVlength(AssemblyContig &contig);
	void SRS_alignment_and_get_var(bool print_log, int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len, int ref_region_length);
	/*************************************SRS INV/BND FUNCTIONs**********************************************/
	bool SRS_assembly_and_get_breakpoints_TRA(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am,
			int main_tid, int &main_pos, bool main_read_before_BP, bool main_read_forward,
			int supp_tid, int &supp_pos, bool supp_read_before_BP, bool supp_read_forward);
	int SRS_assembly_and_get_breakpoints_INV(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am, RefRegion &main_region,
			RefRegion &supp_region, int &main_BP, int &supp_BP);
	bool SRS_assembly_and_genetyping_BND(SVE &sv,
			bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward, bool main_ref_is_forward,
			bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward, bool supp_ref_is_forward,
			std::vector<NOVA_SV_FINAL_RST_item> &region_SVs, bam_hdr_t * header);//buff to store contig
	bool SRS_INV_genotyping_and_store(SVE &sv,
			bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward,
			bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward,
			bam_hdr_t * header, int *GT_result);//buff to store contig
	void SRS_print_signal_list(SIG::T sigT, SV::T svt);
	bool SRS_store_INV(int chr_ID, int position_bg,	int position_ed,
			std::vector<NOVA_SV_FINAL_RST_item> &region_SVs_TMP, int region_support_number[3], int final_genotype);
	//void BND and INVS
	void SRS_BND_WRITE();
	void SRS_show_signals();
	void SRS_REMOVE_DUPLICATED_INDEL();
	void SRS_INV_signal_combine(std::vector<SVE> &inv_sve_l);
	bool SRS_INV_ASSEMBLY(SVE &c_sve, int &ASS_SIGNAL_NUM);
	void SRS_INV_GT(SVE &c_sve, int ASS_SIGNAL_NUM);
	void SRS_INV_WRITE();
	/*************************************SRS GT FUNCTIONs**********************************************/
	void SRS_GENOTYPING_INS_DEL();
	/*************************************SRS MAIN PROCSEEs**********************************************/
	void SRS_SMALL_VAR_CALLING_PROCESS();
	void SRS_VNTR_VAR_CALLING_PROCESS();
	void SRS_INV_CALLING_AND_GT_PROCESS();
	void SRS_LARGE_DELETION_CALLING_PROCESS();
	void SRS_BND_CALLING_AND_GT_PROCESS();
	/*************************************LRS FUNCTIONs**********************************************/
	bool LRS_assembly_variations(bool print_log, BAM_handler *cur_LRS_read, LRS_SIG &sig, RefRegion &loadRef,std::string &full_ref_str, std::vector<uint8_t> &bin_ref);
	void LRS_SV_generating_Germline(bool print_log, LRS_SIG &sig, RefRegion &loadRef, std::string &full_ref_str, std::vector<uint8_t> &bin_ref,
			 std::vector<NOVA_SV_FINAL_RST_item> & SV_result_final, std::vector<TARGET_path_and_index> &target_path_idx_l, CONTIG_SELECTER &contig_selecter);
	void LRS_SV_generating_Somatic(bool print_log, LRS_SIG &sig,
			RefRegion &loadRef, std::string &full_ref_str, std::vector<uint8_t> &bin_ref,
			std::vector<NOVA_SV_FINAL_RST_item> &SV_result_final,
			std::vector<TARGET_path_and_index> &target_path_idx_l_Tumer, CONTIG_SELECTER &contig_selecter_Tumer,
			std::vector<TARGET_path_and_index> &target_path_idx_l_Normal,CONTIG_SELECTER &contig_selecter_Normal);
	void LRS_contig_candidate_generate(bool print_log, int minDepth_LRS, int KMER_len);
	bool LRS_contig_realignment(bool print_log, Contig_String_aligner *ca_p,//background
			std::string & full_contig, std::vector<uint8_t> &bin_contig,//input, contig
			RefRegion & ref_r1, RefRegion & ref_r2, //
			std::string & re_alignment_ref_str, std::vector<uint8_t> & re_alignment_bin_ref_string, //input, reference
			std::string & full_reference_string,
			int supp_read_n, int not_supp_read_n, int unknown_read_n, int haplotype_ID, int total_haplotype_number,  //info
			std::vector<NOVA_SV_FINAL_RST_item> & SV_output_buff,  bool is_replacement//output
			);//buff to store contig
	bool LRS_read_quality_check(bool print_log, bam1_t *br, uint8_t * seq_bin);
	void LRS_SIG_COLLECTION(bool print_log, BAM_handler &cur_read, LRS_SIG_Handler & sig_h,	float &average_read_depth);
	void LRS_SV_SIG_COMBINE(bool print_log,	LRS_SIG_Handler &sig_h,	float average_read_depth);
	void LRS_SV_CALLING_germline(bool print_log, std::vector<LRS_SIG> &FINAL_sig, float average_read_depth, bool output_pure_contig);
	void LRS_remove_duplications(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> &LRS_r_l);
	void LRS_combine_Homozygous_SVs(std::vector<NOVA_SV_FINAL_RST_item> &LRS_r_l);
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Hybird functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	void Hybrid_debug_code_load_SVs_from_vcf_f(std::vector<NOVA_SV_FINAL_RST_item> & result_sv_l, const char * vcf_fn, RefHandler *ref);
	void Hybrid_remove_similar_SVs_in_Target(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> & query,
			std::vector<NOVA_SV_FINAL_RST_item> & target);
	void Hybrid_regenotype_NGS_USING_TGS_data(bool print_log);

	void set_region_addition_load_short(){ 			  region_addition_load = 500;   ca.setZdrop(1200, 1200); }
	void set_region_addition_load_long(){  			  region_addition_load = 2000;  ca.setZdrop(2000, 2000); }
	void set_region_addition_load_extramly_long(){	  region_addition_load = 5000;  ca.setZdrop(5000, 5000); }
	void set_region_addition_load_super_super_long(){ region_addition_load = 50000; ca.setZdrop(50000, 50000);}

public:
	/*************************************PUBLIC FUNCTIONs**********************************************/
	void init(RefHandler *ref_, KMER_ERROR_PROFILE_HANDLER * kmer_profile_handler, RefLocalRepeat *rlr_,
			char * SRS_BAM_F, char * realigned_Filename, char * TL_read_Filename,
			char * LRS_BAM_F, char * LRS_Tumor_BAM_F, char * LRS_Normal_BAM_F, char * INNER_TEST_FLAG,
			BAM_STATUS *bs, char *outputFile, bool is_compression, int MIN_sv_len_,
			bool output_vcf_header, bool force_calling_ALU, bool random_phasing, bool output_small_var,
			int preset_platform, char * FC_BED_F, int OUT_MODE_int){
		this->random_phasing = random_phasing;
		this->output_small_var = output_small_var;
		this->preset_platform = preset_platform;
		this->OUT_MODE_int = OUT_MODE_int;
		this->kmer_profile_handler = kmer_profile_handler;

		ref_handler = ref_;
		rlr = rlr_;
		this->SRS_BAM_F = SRS_BAM_F;
		this->LRS_BAM_F = LRS_BAM_F;
		this->LRS_Tumor_BAM_F = LRS_Tumor_BAM_F;
		this->LRS_Normal_BAM_F = LRS_Normal_BAM_F;
		this->FC_BED_F = FC_BED_F;
		this->INNER_TEST_FLAG = INNER_TEST_FLAG;

		cur_header = NULL;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~germline~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//for NGS reads
		if(SRS_BAM_F != NULL){
			//S1: basic parameters
			sig_para.init(bs);
			//S2: open BAM file
			SRS_read.init(SRS_BAM_F, realigned_Filename, TL_read_Filename, ref_handler->get_refFileName());
			rdc.init(SEGMENT_LEN, bs->ave_read_depth, bs->analysis_read_length);
			//S3: get distribution
			bs->getBreakPoint_Distribution(DR_bp_distribution, SH_bp_distribution, UM_stPos_distribution, START_OFFSET_UM);
			SRS_set_min_accpet_possibility();
			//S4: allocated buffs for signaling
			UMQueryBuff = (uint8_t *)xcalloc(5000, sizeof(uint8_t));
			possibility_r = (float *)xcalloc(5000, sizeof(float));
			//ALU force
			alu_handler = new ALU_force_handler();
			alu_handler->init();
			this->force_calling_ALU = force_calling_ALU;
			cur_header = SRS_read.file._hdr;
		}
		if(LRS_BAM_F != NULL){
			LRS_read.init(LRS_BAM_F, NULL, NULL, ref_handler->get_refFileName());
			cur_header = LRS_read.file._hdr;
		}
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SOMATIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if(LRS_Tumor_BAM_F != NULL){
			LRS_Tumor_read.init(LRS_Tumor_BAM_F, NULL, NULL, ref_handler->get_refFileName());
			cur_header = LRS_Tumor_read.file._hdr;
		}

		if(LRS_Normal_BAM_F != NULL){
			LRS_Normal_read.init(LRS_Normal_BAM_F, NULL, NULL, ref_handler->get_refFileName());
			cur_header = LRS_Normal_read.file._hdr;
		}

		//anti registration for BAM header of REF handler
		xassert(cur_header != NULL, "Fatal ERROR: NO BAM\n");

		ref_handler->registHeader(cur_header);
		//S5: init assembler
		MIN_sv_len = MIN_sv_len_;
		am = new MainAssemblyHandler[1];
		dbgHandler = new MainDBGHandler1[1];
		dbgHandler->init();
		ca.init();
		somatic_same_ca.init();
		ga.init();
		//S6 open result vcf files
		//default: stdout
		if(strcmp(outputFile, "stdout") == 0) 	vcf_w = stdout;
		else								vcf_w = xopen(outputFile, "w");
		vcfBuffStr.s = (char *)xcalloc(1000000, sizeof(char));
		vcfBuffStr.l = 0; vcfBuffStr.m = 1000000;
		if(output_vcf_header && this->OUT_MODE_int != OUT_MODE_PURE_STR)
			NOVA_SV_FINAL_RST_item::write_to_vcf_header(vcf_w, cur_header);

		//realignment bam handler
		global_SV_ID = 0;

		//BND
		bnd_ass_block.init();
		ca_bnd.init();
		ca_bnd.setZdrop(1500, 1500);

		ca_re_locate.init();
		ca_re_locate.setZdrop(50000,50000);

		contig_polsihing_handler.init(kmer_profile_handler, &SRS_read.file);//todo::

		contig_align_region_Handler.init(ref_handler);
		clip_contig_combine_Handler.init();
	}

	void distory(){
		if(SRS_BAM_F != NULL){
			SRS_read.destroy();
			free(UMQueryBuff);
			free(possibility_r);
			delete(alu_handler);
		}

		if(LRS_BAM_F != NULL)		LRS_read.destroy();
		if(LRS_Tumor_BAM_F != NULL)LRS_Tumor_read.destroy();
		if(LRS_Normal_BAM_F != NULL)LRS_Normal_read.destroy();

		delete []am;
		ca.destory();
		somatic_same_ca.destory();
		ga.destory();
		fclose(vcf_w);
		free(vcfBuffStr.s);
	}

    void SV_calling_in_region_main_pipeline(){
    	bool is_somatic_calling = (LRS_Tumor_BAM_F != NULL);

    	bool is_FC_REGION_LRS = (true == fc_lrs_handler.region_is_loaded);
    	/*～～～～～～～～～～～～～～～～～～～～～～～～～～SOMATIC～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
    	if(is_somatic_calling){
    		//init results
    		result_LRS_final.clear();
    		bool print_log = false;
    		//clear SIG and RST:
    		LRS_SIG_Handler tumor_sig_h;
    		LRS_SIG_Handler normal_sig_h;
    		//load signals： tumor
    		float average_read_depth_tumer = 0;
    		LRS_SIG_COLLECTION(print_log, LRS_Tumor_read, tumor_sig_h, average_read_depth_tumer);
    		LRS_SV_SIG_COMBINE(print_log, tumor_sig_h, average_read_depth_tumer);
    		int min_support_read_number_tumer = average_read_depth_tumer*0.1;

    		float average_read_depth_normal = 0;
    		LRS_SIG_COLLECTION(print_log, LRS_Normal_read, normal_sig_h, average_read_depth_normal);
    		LRS_SV_SIG_COMBINE(print_log, normal_sig_h, average_read_depth_normal);
    		//int min_support_read_number_normal = average_read_depth_normal*0.1;

    		tumor_sig_h.show();
    		normal_sig_h.show();
    		std::vector<LRS_SIG> &l_T= tumor_sig_h.FINAL_sig;
    		//P3： calling
    		int signal_list_size = l_T.size();
    		RefRegion ref_r = *(ref_handler->get_cur_region());
    		for(int i = 0; i < signal_list_size; i++){
    			//skip when out of range
    			RefRegion r(l_T[i].chrID, l_T[i].POS, l_T[i].END);
    			if(!(r.region_overlap(ref_r)))
    				continue;
    			set_region_addition_load_short();
    			l_T[i].show_simple();
    			if(l_T[i].support_read_Num >= min_support_read_number_tumer){
    				RefRegion loadRef;	std::string full_ref_str;	std::vector<uint8_t> bin_ref;
    				//LRS_;
    				if(false == LRS_assembly_variations(print_log, &LRS_Tumor_read, l_T[i], loadRef, full_ref_str, bin_ref))
    					continue;
    				//store the results Tumer
    				std::swap(target_path_idx_l_Tumer, target_path_idx_l);
    				std::swap(contig_selecter_Tumer, contig_selecter);
    				if(LRS_Normal_BAM_F != NULL){
        				if(false == LRS_assembly_variations(print_log, &LRS_Normal_read, l_T[i], loadRef, full_ref_str, bin_ref))
        					continue;
        				//store the results Normal
        				std::swap(target_path_idx_l_Normal, target_path_idx_l);
        				std::swap(contig_selecter_Normal, contig_selecter);
    				}
    				//SOMATIC SV generation:: todo::?????
    				LRS_SV_generating_Somatic(print_log, l_T[i], loadRef, full_ref_str, bin_ref,
    						result_LRS_final,
							target_path_idx_l_Tumer,  contig_selecter_Tumer,
							target_path_idx_l_Normal, contig_selecter_Normal);
    			}
    		}

    		//combine CLIPING and MIDDLE, store the final results in the result_LRS_clip
			//remove duplication SVs in different reference regions
			LRS_remove_duplications(true, result_LRS_final);

			if(print_log) showSVList(result_LRS_final, "\n\n After remove_LRS_duplications:\n");

			//combine same SVs in different haplotype-s
			//0/1 + 0/1 ---> 1/1
			LRS_combine_Homozygous_SVs(result_LRS_final);

			if(print_log) showSVList(result_LRS_final, "\n\n After combine_Homozygous_SVs:\n");

			//set basic GT:
			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final)
				sv.calGT_LRS();

			fprintf(stderr, "genotyping_LRS_Using_SRS SKIP\n\n");

			if(random_phasing){
				for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final)
					sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
			}
			//output the final results
			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final){
				if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID))
					fprintf(vcf_w, "%s", vcfBuffStr.s);
			}
    	}

    	/*～～～～～～～～～～～～～～～～～～～～～～～～～～GERMLINE～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
    	//region force-calling
    	if(!is_somatic_calling && is_FC_REGION_LRS){
    		//init results
    		result_LRS_final.clear();
    		//load all regions
    		bool print_log = false;
    		float average_read_depth = 0;
    		LRS_SIG_Handler lrs_sig_h;
    		LRS_SIG_COLLECTION(print_log, LRS_read, lrs_sig_h, average_read_depth);
    		//FC sig storing:
    		fc_lrs_handler.store_region_to_sig(ref_handler->get_cur_region(), lrs_sig_h);
    		if(print_log && false) lrs_sig_h.show();
    		bool output_pure_contig = (OUT_MODE_int == OUT_MODE_PURE_STR);
    		LRS_SV_CALLING_germline(print_log, lrs_sig_h.FINAL_sig, average_read_depth, output_pure_contig);

    		if(OUT_MODE_int == OUT_MODE_VCF){
    			LRS_remove_duplications(true, result_LRS_final);
    			LRS_combine_Homozygous_SVs(result_LRS_final);
    			//output the final results
    			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final){
    				sv.calGT_LRS(); //set basic GT:
    				if(random_phasing) sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
    				if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID))
    					fprintf(vcf_w, "%s", vcfBuffStr.s);
    			}
    		}
    	}
    	//basic SV calling
    	if(!is_somatic_calling && !is_FC_REGION_LRS){
            /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART1: Signal Process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
            fprintf(stderr, "Signal process begin\n" );
            if(LRS_BAM_F != NULL){//handle ONT dataset
            	///return;
            	bool print_log = false;
            	result_LRS_final.clear();
        		//SIG list:
        		LRS_SIG_Handler lrs_sig_h;
        		float average_read_depth = 0;
        		LRS_SIG_COLLECTION(print_log, LRS_read, lrs_sig_h, average_read_depth);
        		if(print_log && false) lrs_sig_h.show();
        		LRS_SV_SIG_COMBINE(print_log, lrs_sig_h, average_read_depth);
        		LRS_SV_CALLING_germline(print_log, lrs_sig_h.FINAL_sig, average_read_depth, false);
        		//combine CLIPING and MIDDLE, store the final results in the result_LRS_clip
    			//remove duplication SVs in different reference regions
    			LRS_remove_duplications(true, result_LRS_final);
    			if(print_log) showSVList(result_LRS_final, "\n\n After remove_LRS_duplications:\n");
    			//combine same SVs in different haplotype-s
    			//0/1 + 0/1 ---> 1/1
    			LRS_combine_Homozygous_SVs(result_LRS_final);
    			if(print_log) showSVList(result_LRS_final, "\n\n After combine_Homozygous_SVs:\n");
    			//set basic GT:
    			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final)
    				sv.calGT_LRS();

    			//re-genotyping LRS result using NGS reads when data available
    			if(SRS_BAM_F != NULL && !result_LRS_final.empty() ){
    				fprintf(stderr, "genotypingTGS_Using_SRS begin:\n");
    				for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final){
    					if(ABS(sv.SV_length) < MIN_sv_len)//only re-genotype SVs NOT shorter than MIN_sv_len
    						continue;
    					showSV(sv, "\nLRS re-genotype using SRS, Genotyping for:");
    					int new_GT = sv.genotyping_LRS_Using_SRS(false, sig_para.MaxReadLen, &SRS_read.file, &ga, ref_handler);
    					if(sv.LRS_need_filter_by_SRS() && sv.contig_not_support_by_NGS())
    						sv.setGenotype_directly(0);//set to 0/0
    					else if(sv.LRS_need_reGT_by_SRS() && new_GT > 0 && new_GT <= 2)
    						sv.setGenotype_directly(new_GT);
    				}
    			}else fprintf(stderr, "genotyping_LRS_Using_SRS SKIP\n\n");
    			if(random_phasing){
    				for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final)
    					sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
    			}
    			//output the final results
    			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_LRS_final){
    				if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID))
    					fprintf(vcf_w, "%s", vcfBuffStr.s);
    			}
            }

            if(SRS_BAM_F != NULL){//debug code: loading NGS vcfs, mergin and output:// not used right now
                //handle NGS reads: clear result list:
    			result_SRS_INDEL.clear();//store the result for insertion, deletion, large tra-insertion
    			result_NGS_BND.clear(); //store the result for the INV and the BNDs
    			result_BGS_INV.clear();

    			for(int i = 0; i < SIG::LEN; i++)
    				for(int j = 0; j < SV::LEN; j++) SRS_sve[i][j].clear();
    			SRS_sveVNTR.clear();
    			SRS_read.clear();

    			//calling SVs using NGS reads
            	if(LRS_BAM_F != NULL && INNER_TEST_FLAG != NULL){//DEBUG mode: when with LRS dataset, run in simple mode ,30X, simple merge:
            		const char * vcf_fn = INNER_TEST_FLAG;
            		Hybrid_debug_code_load_SVs_from_vcf_f(result_SRS_INDEL, vcf_fn, ref_handler);
            		fprintf(stderr, "DEBUG_SRS using: %s\n", vcf_fn);
            	}
            	else
            	{  //pure NGS mode: truly handle NGS reads
                    //get signals from NGS CRAM
                    SRS_GET_SIGNALS_FROM_READs();//S2:load reads and get signals
                    //not used right now
                    if(false)SRS_load_realignment_read_from_file();//S3:try to load realignment reads into UNMAPPED read list
                    SRS_CLUSTERING_AND_COMBINE_signals();//S4: cluster signals and get SV signals
                    if(true)
                    	SRS_show_signals();            //log: print signals
                    /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART1: INDEL calling(assembly and genotyping) process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
                    //results is stored in "result_INDEL"
                    //call small INDELs
                    SRS_SMALL_VAR_CALLING_PROCESS();/* SMALL INDEL calling process*/
                    SRS_VNTR_VAR_CALLING_PROCESS();/* VNTR calling process*/
                    SRS_LARGE_DELETION_CALLING_PROCESS();/* LARGE DELETION calling process*/

                    SRS_REMOVE_DUPLICATED_INDEL();/*sort and removing duplication process(for INS and DEL)*/
                    SRS_GENOTYPING_INS_DEL();/*GENOTYPING and output for INS/DEL*/
                    /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART3: INV calling process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
        			SRS_INV_CALLING_AND_GT_PROCESS();
                    /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART2: BND calling(assembly and genotyping) process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
                    SRS_BND_CALLING_AND_GT_PROCESS();
            	}

                //remove SVs called by NGS data-set that is similar with LRS results
                if(!result_SRS_INDEL.empty() && LRS_BAM_F != NULL){
                	bool print_log = false;
                	Hybrid_regenotype_NGS_USING_TGS_data(print_log);
                }

            	if(random_phasing){
        			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_SRS_INDEL ){
        				sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
        			}
            	}

                //Output all the NGS results
                {
                	bool print_log = true;
        			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_SRS_INDEL ){
        				if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID)){
        					if(print_log) fprintf(stderr, "%s", vcfBuffStr.s);
        					fprintf(vcf_w, "%s", vcfBuffStr.s);
        				}
        			}
        			SRS_INV_WRITE();
                    SRS_BND_WRITE();//BND write
                }
            }
    	}
    }

public:
	void showSVList(std::vector<NOVA_SV_FINAL_RST_item> & sl, const char * title){
		fprintf(stderr, "%s\n", title);
		for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: sl){
			if(sv.writeVCF_final(&vcfBuffStr, cur_header, NULL))
				fprintf(stderr, "%s", vcfBuffStr.s);
		}
	}

	void showSV(NOVA_SV_FINAL_RST_item & sv, const char * title ){
		fprintf(stderr, "%s\n", title);
		if(sv.writeVCF_final(&vcfBuffStr, cur_header, NULL))
			fprintf(stderr, "%s", vcfBuffStr.s);
	}
};

bam_hdr_t* checkBamHeadersAreSAME(std::vector<SV_CALLING_Handler>& h);
void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
void ref_load_char_2_bin(bool cur_ref_is_forward, int load_len, std::vector<uint8_t> & to_bin_string, char* cur_load_ref);
void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
void store_bin_from_char(const char * contig_seq, int contig_seq_len, std::vector<uint8_t> &bin_contig);

#endif /* SRC_SIGNAL_SVEHANDLER_HPP_ */
