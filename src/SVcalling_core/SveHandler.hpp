/*
 * SveHandler.hpp
 *
 *  Created on: 2020年4月29日
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

#define CCS_SIG_INS 0
#define CCS_SIG_DEL 1
#define CCS_SIG_INV 2
#define CCS_SIG_TRA 3
#define CCS_SIG_CLIP_LEFT 4
#define CCS_SIG_CLIP_RIGHT 5

struct CCS_SIG{
	int8_t type;
	int chrID;
	int64_t POS;
	int END;
	//is used
	bool isUsed;
	int support_read_Num;

	//return true when similar
	//return false when NOT
	bool similar_SIG_Normal(CCS_SIG &to_compare){
		//position is nearby:
		int ABS_P1 = ABS_U(POS, to_compare.POS);
		int ABS_P2 = ABS_U(END, to_compare.END);
		int ABS_P3 = ABS_U(POS, to_compare.END);
		int ABS_P4 = ABS_U(END, to_compare.POS);
		if(chrID == to_compare.chrID && (ABS_P1 < 1000 || ABS_P2 < 1000 || ABS_P3 < 1000 || ABS_P4 < 1000)) return true;
		return false;
	}

	bool similar_SIG_CLIP(CCS_SIG &to_compare){
		//position is nearby:
		int ABS_POS = ABS_U(POS, to_compare.POS);
		if(chrID == to_compare.chrID && type == to_compare.type && ABS_POS < 200) return true;
		return false;
	}

	bool SIG_CLIP_pairing(CCS_SIG &to_compare){
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

	void show_core(char *TYPE_NAME){
		if(TYPE_NAME == NULL)
			fprintf(stderr, "TYPE:%d\t",type);
		else
			fprintf(stderr, "TYPE:%s\t",TYPE_NAME);
		fprintf(stderr, "[%d:%ld]\tEND:%d;[isUsed:%d, support_read_Num:%d]\n", chrID, POS, END, isUsed, support_read_Num);
	}

	void show(char TYPE_NAME[TYPE_NUM][16]){
		show_core(TYPE_NAME[type]);
	}

	void show_simple(){
		show_core(NULL);
	}

	static inline int cmp_by_ref_pos(const CCS_SIG &a, const CCS_SIG &b){
		if(a.type != b.type)
			return a.type < b.type;
		if(a.chrID != b.chrID)
			return a.chrID < b.chrID;
		return a.POS < b.POS;
	}

	static inline int cmp_by_ref_pos_ignore_type(const CCS_SIG &a, const CCS_SIG &b){
		if(a.chrID != b.chrID)
			return a.chrID < b.chrID;
		return a.POS < b.POS;
	}

};


struct SveHandler{
private:
	//part 1: reference, it is only a reference for the refHandler
	RefHandler *ref;
	RefLocalRepeat *rlr;
	std::vector<RefLocalRepeatItemR> rlr_result;

	//part 2: bam files, the SVE handler keep and maintained the BAM files
	char * NGS_bamFileName;
	char * CCS_bamFileName;
	char * ONT_bamFileName;

	BAM_handler NGS_read;
	BAM_handler CCS_read;
	BAM_handler ONT_read;

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
	SIGNAL_PARAMETER sig_para;
	//part 6: data for break point distributions, those data used for calculating break point distribution for each signals
	std::vector<float> DR_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_DR = 0;
	std::vector<float> SH_bp_distribution;
	float MIN_ACCEPT_POSSIBILITY_SH = 0;
	std::vector<float> UM_stPos_distribution;

	float MIN_ACCEPT_POSSIBILITY_TGS_CCS = 0;
	std::vector<float> UM_stPos_distribution_TGS_CCS;

	float MIN_ACCEPT_POSSIBILITY_UM = 0;
	int START_OFFSET_UM = 0;
	//part 7.0: buffs for SVE signals
	SVE_L sve[SIG::LEN][SV::LEN]; //SVE signals list in the first combining loop
	SVE_L sveVNTR; //SVE signals for VNTR
	//SVE_L final; //SVE signals list in the second combining loop

	//part 7: assembly manager
	Ass_Block ass_block;//load reads and run assemblies
	//Ass_Block ass_block_TGS;//load reads and run assemblies
	//Ass_Block ass_block_r2;
	std::vector<uint8_t> bin_contig;
	MainAssemblyHandler *am;
	Contig_String_aligner ca;
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

	std::vector<NOVA_SV_FINAL_RST_item> result_NGS_INDEL;//large INS + DEL/ ALU
	std::vector<NOVA_SV_FINAL_RST_item> result_NGS_BND;//BND
	std::vector<NOVA_SV_FINAL_RST_item> result_BGS_INV;//INV

	std::vector<NOVA_SV_FINAL_RST_item> result_TGS_middle;
	std::vector<NOVA_SV_FINAL_RST_item> result_TGS_clip;
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

	#define HUMAN_ALU_STRING "NNGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGNGGATCANGAGGTCAGGAGATCGAGACCATCCNGGCTAANANGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGNNGTGGCGGGCGCCTGTAGTCCCAGCTACTNGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
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
		void init(){
			//add string:
			ALU_full_string[0].append(HUMAN_ALU_STRING);
			ALU_FULL_LEN = ALU_full_string[0].size();

			const char * char_alu = ALU_full_string[0].c_str();
			ALU_full_string_bin[0].resize(ALU_FULL_LEN);
			ALU_full_string[1].resize(ALU_FULL_LEN);
			ALU_full_string_bin[1].resize(ALU_FULL_LEN);

			for(int i = 0; i < ALU_FULL_LEN; i++){
				uint8_t bin_c; char reverse_c; uint8_t reverse_bin_c;
				switch (char_alu[i]) {
					case 'A': bin_c = 0; reverse_c = 'T'; reverse_bin_c = 3; break;
					case 'C': bin_c = 1; reverse_c = 'G'; reverse_bin_c = 2; break;
					case 'G': bin_c = 2; reverse_c = 'C'; reverse_bin_c = 1; break;
					case 'T': bin_c = 3; reverse_c = 'A'; reverse_bin_c = 0; break;
					default:  bin_c = 4; reverse_c = 'N'; reverse_bin_c = 4; break;
				}
				ALU_full_string_bin[0][i] = bin_c;
				ALU_full_string[1][ALU_FULL_LEN - i - 1] = reverse_c;
				ALU_full_string_bin[1][ALU_FULL_LEN - i - 1] = reverse_bin_c;
			}
			// add kmers
			std::vector<int> kmer_pos = {28, 68, 109, 190, 213};
			std::vector<int> kmer_size = {20, 18, 16,  17,  16};

			for(int i = 0; i < (int)kmer_pos.size(); i++){
				ALU_kmer_list.emplace_back(ALU_full_string[0].substr(kmer_pos[i],							    kmer_size[i]).c_str(),  kmer_pos[i], 0);
				ALU_kmer_list.emplace_back(ALU_full_string[1].substr(ALU_FULL_LEN - kmer_pos[i] - kmer_size[i], kmer_size[i]).c_str(), -kmer_pos[i], 1);
			}

			ALU_kmer_list_size = ALU_kmer_list.size();

			//INIT CA
			aligner[0].init();
			aligner[0].setRef(&(ALU_full_string_bin[0][0]), ALU_full_string_bin[0].size(), 0, 0);
			aligner[1].init();
			aligner[1].setRef(&(ALU_full_string_bin[1][0]), ALU_full_string_bin[1].size(), 0, 0);
		}

		bool ALU_INS_check(uint8_t *q_seq, int q_len , int ALU_DIR){
			//first 200BP:
			if(ALU_DIR == 0){
				q_len = 200;
			}else{
				q_seq += 120;
				q_len = 200;
			}

			//N base check:
			int N_base = 0;
			for(int i = 0;i < q_len; i++ )
				if(q_seq[i] > 3)
					N_base++;
			fprintf(stderr, "ALU_INS_check BEGIN\n");

			Contig_String_aligner * cur_aligner = &(aligner[ALU_DIR]);
			cur_aligner->align_non_splice_HARD_LOCAL(q_seq, q_len, 0, ALU_FULL_LEN);
			int suggest_st_pos = 0;//cur_aligner->adjustCIGAR();
			cur_aligner->printf_ref(stderr);
			fprintf(stderr, "ALU_INS_check_detail\t");
			cur_aligner->printf_alignment_detail(stderr, suggest_st_pos, NULL, 0);
			int match_base = cur_aligner->check_match_base(suggest_st_pos);
			int UN_match_base = q_len - match_base - N_base - 10;
			bool pass_check = true;
			if(UN_match_base >= 10){
				pass_check = false;
			}
			fprintf(stderr, "match_base %d,UN_match_base %d, N_base %d %s\n",match_base,UN_match_base,N_base, (pass_check==true)?"PASS":"FAIL");
			fprintf(stderr, "ALU_INS_check END\n");
			return pass_check;
		}

		int kmer_check(std::string & contig, int &alu_kmer_i){
			alu_kmer_i = -1;
			for(int i = 0; i < ALU_kmer_list_size; i++){
				uint32_t search_pos = contig.find(ALU_kmer_list[i].s.c_str());
				if(search_pos != MAX_uint32_t){
					alu_kmer_i = i;
					return search_pos;
				}
			}
			return -1;
		}

		int get_suggest_BP_in_contig(int alu_kmer_i, int search_pos, int & ALU_DIR){
			ALU_DIR = ALU_kmer_list[alu_kmer_i].dir;
			int suggest_BP_in_contig = 0;
			if(ALU_DIR == 0){//type AA:
				suggest_BP_in_contig = search_pos - ALU_kmer_list[alu_kmer_i].p;
			}else{//type TT:
				suggest_BP_in_contig = search_pos - ALU_kmer_list[alu_kmer_i].p + ALU_kmer_list[alu_kmer_i].s.size();
			}
			return suggest_BP_in_contig;
		}

	};

	//ALU force handler
	ALU_force_handler *alu_handler;

public:
	/*************************************PUBLIC FUNCTIONs**********************************************/
	void set_min_accpet_possibility(){
		float MAX_POSSIBILITY_DR = 0;
		for(float poss:DR_bp_distribution)
			MAX_POSSIBILITY_DR = MAX(MAX_POSSIBILITY_DR, poss);
		MIN_ACCEPT_POSSIBILITY_DR = 2*MAX_POSSIBILITY_DR;

		float MAX_POSSIBILITY_SH = 0;
		for(float poss:SH_bp_distribution)
			MAX_POSSIBILITY_SH = MAX(MAX_POSSIBILITY_SH, poss);
		MIN_ACCEPT_POSSIBILITY_SH = 2*MAX_POSSIBILITY_SH;

		float MAX_POSSIBILITY_UM = 0;
		for(float poss:UM_stPos_distribution)
			MAX_POSSIBILITY_UM = MAX(MAX_POSSIBILITY_UM, poss);
		MIN_ACCEPT_POSSIBILITY_UM = 3*MAX_POSSIBILITY_UM;
	}

	void init(RefHandler *ref_, RefLocalRepeat *rlr_, char * NGS_BamFileName_, char * realigned_Filename, char * TL_read_Filename,
			char *CCS_BAM_F_, char *ONT_BAM_F_,
			BAM_STATUS *bs, char *outputFile, bool is_compression, int MIN_sv_len_, bool output_vcf_header, bool force_calling_ALU, bool random_phasing, bool output_small_var){
		this->random_phasing = random_phasing;
		this->output_small_var = output_small_var;
		ref = ref_;
		rlr = rlr_;
		NGS_bamFileName = NGS_BamFileName_;
		CCS_bamFileName = CCS_BAM_F_;
		ONT_bamFileName = ONT_BAM_F_;
		//for NGS reads
		if(NGS_bamFileName != NULL){
			//S1: basic parameters
			sig_para.init(bs);
			//S2: open BAM file
			NGS_read.init(NGS_bamFileName, realigned_Filename, TL_read_Filename, ref->get_refFileName());
			rdc.init(SEGMENT_LEN, bs->ave_read_depth, bs->analysis_read_length);
			//S3: get distribution
			bs->getBreakPoint_Distribution(DR_bp_distribution, SH_bp_distribution, UM_stPos_distribution, START_OFFSET_UM);
			set_min_accpet_possibility();
			//S4: allocated buffs for signaling
			UMQueryBuff = (uint8_t *)xcalloc(5000, sizeof(uint8_t));
			possibility_r = (float *)xcalloc(5000, sizeof(float));
			//ALU force
			alu_handler = new ALU_force_handler();
			alu_handler->init();
			this->force_calling_ALU = force_calling_ALU;
		}

		if(CCS_bamFileName != NULL){
			CCS_read.init(CCS_bamFileName, NULL, NULL, ref->get_refFileName());
		}

		if(ONT_bamFileName != NULL){
			ONT_read.init(ONT_bamFileName, NULL, NULL, ref->get_refFileName());
		}

		//anti registration for BAM header of REF handler
		cur_header = NULL;
		if(NGS_bamFileName != NULL) 	cur_header = NGS_read.file._hdr;
		else if(CCS_bamFileName != NULL)cur_header = CCS_read.file._hdr;
		else if(ONT_bamFileName != NULL)cur_header = ONT_read.file._hdr;

		xassert(cur_header != NULL, "Fatal ERROR: NO BAM\n");

		ref->registHeader(cur_header);
		//S5: init assembler
		MIN_sv_len = MIN_sv_len_;
		am = new MainAssemblyHandler[1];
		dbgHandler = new MainDBGHandler1[1];
		dbgHandler->init();
		ca.init();
		ga.init();
		//S6 open result vcf files
		//default: stdout
		if(strcmp(outputFile, "stdout") == 0) 	vcf_w = stdout;
		else								vcf_w = xopen(outputFile, "w");
		vcfBuffStr.s = (char *)xcalloc(1000000, sizeof(char));
		vcfBuffStr.l = 0; vcfBuffStr.m = 1000000;
		if(output_vcf_header)
			NOVA_SV_FINAL_RST_item::write_to_vcf_header(vcf_w, cur_header);

		//realignment bam handler
		global_SV_ID = 0;

		//BND
		bnd_ass_block.init();
		ca_bnd.init();
		ca_bnd.setZdrop(1500, 1500);

		ca_re_locate.init();
		ca_re_locate.setZdrop(50000,50000);

	}

	void distory(){
		if(NGS_bamFileName != NULL){
			NGS_read.destroy();
			free(UMQueryBuff);
			free(possibility_r);
			delete(alu_handler);
		}

		if(CCS_bamFileName != NULL)
			CCS_read.destroy();

		if(ONT_bamFileName != NULL)
			ONT_read.destroy();

		delete []am;
		ca.destory();
		ga.destory();
		fclose(vcf_w);
		free(vcfBuffStr.s);
	}

	void print_signal_list(SIG::T sigT, SV::T svt){
		SVE_L * cur_var_signals = &(sve[sigT][svt]);
		fprintf(stderr, "Print %s::%s signals, size %ld\n", SIG::STR[sigT].c_str(), SV::STR[svt].c_str(), cur_var_signals->size());
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
	}

	bool store_INV(int chr_ID, int position_bg, int position_ed, std::vector<NOVA_SV_FINAL_RST_item> &region_SVs_TMP, int region_support_number[3], int final_genotype){
		std::string ref_inv;
		std::string alt_inv;
		int INV_len = position_ed - position_bg;
		if(final_genotype == 0) { return false;}//length filter
		if(INV_len < 30 || INV_len > 50000){ return false;}//length filter
		if(INV_len <= 2000){
			int load_len;
			char* load_ref = ref->load_ref_by_region(chr_ID, position_bg, position_ed, &load_len);
			//set REF
			ref_inv.append(load_ref);
			//SET ALT
			alt_inv.resize(load_len + 1);
			char *cur_alt_string_p = &(alt_inv[0]);
			for(int i = 0; i < load_len; i++){
				switch( load_ref[load_len - i - 1]){
				case 'A': case 'a': cur_alt_string_p[i] = 'T'; break;
				case 'C': case 'c': cur_alt_string_p[i] = 'G'; break;
				case 'G': case 'g': cur_alt_string_p[i] = 'C'; break;
				case 'T': case 't': cur_alt_string_p[i] = 'A'; break;
				case 'n': case 'N': cur_alt_string_p[i] = 'N'; break;
				}
			}
			alt_inv[load_len] = 0;
			free(load_ref);
		}else{
			ref_inv.append("<INV_REF>");
			alt_inv.append("<INV_ALT>");
		}
		NOVA_SV_FINAL_RST_item::add_to_vector(region_SVs_TMP, chr_ID, position_bg, "INV", &(ref_inv[0]), &(alt_inv[0]), INV_len, NULL, 0, NULL, 0, 0, 0, 0, 0);
		region_SVs_TMP.back().INV_SET_GT(region_support_number, final_genotype);
		return true;
	}

	//void BND and INVS
	void NGS_BND_CALLING_AND_GT_PROCESS(){
		SVE_L * cur_var_signals = NULL;
		int signal_list_size;
		am->setRepeatMode();
		bool read_is_before_breakpoint_in_main;
		bool read_in_main_should_be_forward;
		bool read_is_before_breakpoint_in_supp;
		bool read_in_supp_should_be_forward;
		bool supp_ref_is_forward; bool main_ref_is_forward;
		fprintf(stderr, "BND calling process\n" );
		for(int mode = 0; mode < 2; mode++){
			if(mode == 0)
			{
				cur_var_signals = &(sve[SIG::DR][SV::TRA]);
				read_is_before_breakpoint_in_main = true;
				read_in_main_should_be_forward = true;
				read_is_before_breakpoint_in_supp = false;
				read_in_supp_should_be_forward = false;
				main_ref_is_forward = true;
				supp_ref_is_forward = true;
			}
			else if(mode == 1)
			{
				cur_var_signals = &(sve[SIG::DR][SV::TRA_INV]);
				read_is_before_breakpoint_in_main = false;
				read_in_main_should_be_forward = false;
				read_is_before_breakpoint_in_supp = true;
				read_in_supp_should_be_forward = true;
				main_ref_is_forward = true;
				supp_ref_is_forward = false;
			}
			if(!cur_var_signals->empty()){
				for(unsigned int i = 0; i < cur_var_signals->size(); i++)
					std::cerr << cur_var_signals[0][i];
				signal_list_size = cur_var_signals->size();
				for(int i = 0; i < signal_list_size; i++){
					if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
					assembly_and_genetyping_BND(cur_var_signals[0][i],
							read_is_before_breakpoint_in_main, read_in_main_should_be_forward, main_ref_is_forward,
							read_is_before_breakpoint_in_supp, read_in_supp_should_be_forward, supp_ref_is_forward,
							result_NGS_BND, NGS_read.file._hdr);
				}
			}
		}
	}

	void BND_WRITE(){
        std::vector<NOVA_SV_FINAL_RST_item> region_SVs_TMP;
        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_BND){
            if(sv.is_BND_TRA()){//only output precise results
                if(!sv.is_BND_IMPRECISE()){
                    fprintf(stderr, "OUTPUT BND::TRA\n: \n ");
                    if(sv.writeVCF_final(&vcfBuffStr, NGS_read.file._hdr, &global_SV_ID))
                        fprintf(vcf_w, "%s", vcfBuffStr.s);
                }
            }else{
                region_SVs_TMP.emplace_back();
                std::swap(sv, region_SVs_TMP.back());
            }
        }
        std::swap(region_SVs_TMP, result_NGS_BND);
	}

	//LARGE DELETION calling process
	void NGS_LARGE_DELETION_CALLING_PROCESS(){
		SVE_L * cur_var_signals = NULL;
		int signal_list_size;
		 //part 1 large deletions: (using DR/DEL signals)
		am->setNormalMode();
		fprintf(stderr, "Deletions calling process begin\n" );
		cur_var_signals = &(sve[SIG::DR][SV::DEL]);
		if(!cur_var_signals->empty()){
			for(unsigned int i = 0; i < cur_var_signals->size(); i++)
				std::cerr << cur_var_signals[0][i];
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
//              if(cur_var_signals[0][i].r1.st_pos == 108733125){
//                  std::cerr << "TTTT" <<cur_var_signals[0][i];
//              }
				set_region_addition_load_long();
				int total_region_length = cur_var_signals[0][i].r2.ed_pos - cur_var_signals[0][i].r1.st_pos;
				if(total_region_length > 1800){
					if(total_region_length < 5000)
						set_region_addition_load_extramly_long();
					else if(total_region_length < 50000)
						set_region_addition_load_super_super_long();
					else
						continue;
				}
				//calling the long deletions
				bool withRepeat;
				assembly_variations_ALU_AND_LONG_SV(cur_var_signals[0][i], false, withRepeat);
			}
		}
		fprintf(stderr, "Deletions calling process Done\n" );
	}

	void NGS_SMALL_VAR_CALLING_PROCESS(){
		//part 2: small variations: (using SH/INS signals)
		am->setRepeatMode();
		SVE_L * cur_var_signals = NULL;
		int signal_list_size;
		fprintf(stderr, "Small variations calling process begin\n" );
		cur_var_signals = &(sve[SIG::SH][SV::INS]);
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
		//std::sort(cur_var_signals.begin(), cur_var_signals.end(), SVE::cmp_by_position);//needed to be sorted?
		if(!cur_var_signals->empty()){
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				fprintf(stderr, "\nCurrent SVE ID: [%d]", i);
				if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
				set_region_addition_load_short();
				int cur_TL_read_num = 0;
				bool withRepeat = false;
				{
					int MIN_TL_read_number = (sig_para.read_depth *0.15); MIN_TL_read_number = MAX(4, MIN_TL_read_number);//at least 4 TL reads
					fprintf(stderr, "MIN_TL_read_number %d, TL_read_number %d\n", MIN_TL_read_number, cur_TL_read_num);
					//tran-ins and ALU calling process:
					if(cur_TL_read_num >= 0){
						char log_title[1024]; sprintf(log_title, "[INS_DEL_LONG_ALU]");
						double __cpu_time = cputime(); double __real_time = realtime();
						assembly_variations_ALU_AND_LONG_SV(cur_var_signals[0][i], true, withRepeat);
						fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
					}
				}

				{
					RefRegion &r_main = cur_var_signals[0][i].r1;
					rlr->search(r_main.chr_ID, r_main.st_pos, r_main.ed_pos, rlr_result);
					if(withRepeat){//trigger VNTR caller when repeat is found in assembly caller
						fprintf(stderr, "VNTR Handler is used in repeat region!\n");
						for(RefLocalRepeatItemR &rr : rlr_result)
							rr.show();
						//VNTR calling process
						char log_title[1024]; sprintf(log_title, "[INS_DEL_VNTR]");
						double __cpu_time = cputime(); double __real_time = realtime();
						assembly_variations_INS_DEL_VNTR(cur_var_signals[0][i],true,true);
						cur_TL_read_num = this->TL_read_number;
						fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
					}
					else{
						fprintf(stderr, "Skip VNTR caller\n");
					}
				}
			}
		}
		fprintf(stderr, "Small variations calling process Done\n" );
	}

	void NGS_VNTR_VAR_CALLING_PROCESS(){
		am->setRepeatMode();
		SVE_L * cur_var_signals = NULL;
		int signal_list_size;
		//part 2: small variations: (using SH/INS signals)
		fprintf(stderr, "VNTR variations calling process\n" );
		cur_var_signals = &(sveVNTR);
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
		//std::sort(cur_var_signals.begin(), cur_var_signals.end(), SVE::cmp_by_position);//needed to be sorted?
		if(!cur_var_signals->empty()){
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				fprintf(stderr, "\nCurrent SVE ID: [%d]", i);
				if(SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
				set_region_addition_load_short();
				{
					fprintf(stderr, "VNTR Handler is used in repeat region!\n");
					//VNTR calling process
					char log_title[1024]; sprintf(log_title, "[INS_DEL_VNTR]");
					double __cpu_time = cputime(); double __real_time = realtime();
					assembly_variations_INS_DEL_VNTR(cur_var_signals[0][i],true,true);
					fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
				}
			}
		}
		fprintf(stderr, "VNTR variations calling Done\n" );
	}

	void showNGS_signals(){
        //main type to call INS and DEL
        print_signal_list(SIG::DR, SV::DEL);
        print_signal_list(SIG::SH, SV::INS);
        //other type to call INS
        print_signal_list(SIG::DR, SV::INS);
        //types to call INV and BND
        print_signal_list(SIG::DR, SV::INV_1);
        print_signal_list(SIG::DR, SV::INV_2);
        print_signal_list(SIG::DR, SV::TRA);
        print_signal_list(SIG::DR, SV::TRA_INV);
	}

	void REMOVE_DUPLICATED_INDEL(){
        fprintf(stderr, "True SVs output\n" );
        std::sort(result_NGS_INDEL.begin(), result_NGS_INDEL.end(),
        		NOVA_SV_FINAL_RST_item::cmp_by_position);
        //part4: remove duplication
        int SV_number = result_NGS_INDEL.size();
        for(int c_SV_idx = 0; c_SV_idx < SV_number - 1; c_SV_idx++)
            for(int cmp_sv_idx = c_SV_idx + 1;
                    cmp_sv_idx < SV_number && result_NGS_INDEL[c_SV_idx].duplication_SV_filter(result_NGS_INDEL[cmp_sv_idx]);
                    cmp_sv_idx++);
	}

	void NGS_GENOTYPING_INDELS(){
        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_INDEL ){
            //running genotyping process:
            //if(sv.SV_is_duplicated()) continue;
            char log_title[1024]; sprintf(log_title, "[genotyping]");
            double __cpu_time = cputime(); double __real_time = realtime();
            if(sv.isVNTRcallerRst)
            	fprintf(stderr, "\nGenotyping for(VNTR): \n ");
            else
            	fprintf(stderr, "\nGenotyping for(short INDEL): \n ");
            sv.print(stderr, NGS_read.file._hdr, ref);
            sv.genotyping(sig_para.MaxReadLen,sig_para.read_depth, sig_para.insert_size_max, &NGS_read.file, &ga, ref);
            fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
        }
	}

	void INV_signal_combine(std::vector<SVE> &inv_sve_l){
		while(!sve[SIG::DR][SV::INV_1].empty() || !sve[SIG::DR][SV::INV_2].empty()){
			bool with_INV[2] = {0};
			SVE c_sve;
			if(!sve[SIG::DR][SV::INV_1].empty()){
				std::swap(c_sve, sve[SIG::DR][SV::INV_1].back());
				sve[SIG::DR][SV::INV_1].pop_back();
				with_INV[0] = true;
			}else{
				std::swap(c_sve, sve[SIG::DR][SV::INV_2].back());
				sve[SIG::DR][SV::INV_2].pop_back();
				with_INV[1] = true;
			}
			//try to combine
			for(int mode = 0; mode < 2; mode++){
				for(SVE & combine:sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2]){
					if(combine.r1.region_overlap(c_sve.r1) && combine.r2.region_overlap(c_sve.r2)){
						c_sve.r1.Combine(combine.r1, true);
						c_sve.r2.Combine(combine.r2, true);
						std::swap(combine, sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2].back());
						sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2].pop_back();
						with_INV[mode] = true;
					}
				}
			}
			inv_sve_l.emplace_back();
			std::swap(c_sve, inv_sve_l.back());
			fprintf(stderr, "INV LOGS: \n SIGNAL balance %d %d ; ", with_INV[0], with_INV[1]);
		}
	}

	bool INV_ASSEMBLY(SVE &c_sve, int &ASS_SIGNAL_NUM){
        if(SV_region_depth_filter(c_sve) == false)  return false;;
        int main_pos = c_sve.r1.getMiddle();
        int supp_pos = c_sve.r2.getMiddle();

        bool print_log = false;
        int edge_len = 400;
        RefRegion main_region(c_sve.r1.chr_ID, main_pos - edge_len,main_pos + edge_len);
        RefRegion supp_region(c_sve.r2.chr_ID, supp_pos - edge_len,supp_pos + edge_len);

        //length filter for INV
        int inv_length = ABS_U(main_pos, supp_pos);
        if(inv_length > 50000){//todo::
            fprintf(stderr, "Following INV is skiped because the length is over 50K\n");
            return false;
        }

        //assembly
        ASS_SIGNAL_NUM = assembly_and_get_breakpoints_INV(print_log, bnd_ass_block, am, main_region, supp_region, main_pos, supp_pos);
        c_sve.r1.st_pos = main_pos; c_sve.r1.st_pos = main_pos;
        c_sve.r2.st_pos = supp_pos; c_sve.r2.st_pos = supp_pos;

        fprintf(stderr, "INV: main_pos %d - supp_pos %d\n", main_pos, supp_pos);
        return true;
	}

	void INV_GT(SVE &c_sve, int ASS_SIGNAL_NUM){
		//genotyping
    	int main_pos = c_sve.r1.st_pos;
    	int supp_pos = c_sve.r2.st_pos;

		int GT_result[3][4];
        INV_genotyping_and_store(c_sve, true, true, true, true, NGS_read.file._hdr, GT_result[0]);
        std::swap(c_sve.r1, c_sve.r2);
        INV_genotyping_and_store(c_sve, false, false, false, false, NGS_read.file._hdr, GT_result[1]);
        GT_result[2][0] = GT_result[0][0] + GT_result[1][0];
        GT_result[2][1] = GT_result[0][1] + GT_result[1][1];
        GT_result[2][2] = GT_result[0][2] + GT_result[1][2];
        for(int i = 0; i < 3; i++){
            int *region_support_number = GT_result[i];
            if(region_support_number[0] > region_support_number[1] * 4)     region_support_number[3] = 2;//1/1
            else if(region_support_number[0]*4 < region_support_number[1])  region_support_number[3] = 0;//0/0
            else                                    region_support_number[3] = 1;//0/1
            if((region_support_number[0] + region_support_number[1]) * 3 < region_support_number[2])
                region_support_number[3] = 0;//0/0
            if(region_support_number[0] < 2)
                region_support_number[3] = 0;//0/0
        }

        //filter:
        if(ASS_SIGNAL_NUM != 0 && store_INV(c_sve.r1.chr_ID, main_pos, supp_pos, result_BGS_INV, GT_result[2], GT_result[2][3])){
            fprintf(stderr, "Read support balance %d:%d:%d %d:%d:%d ", GT_result[0][0], GT_result[0][1], GT_result[0][2], GT_result[1][0], GT_result[1][1], GT_result[1][2]);
            fprintf(stderr, "GT balance %d %d ", GT_result[0][3], GT_result[1][3]);
            fprintf(stderr, "ASS BALANCE: ASS_SIGNAL_NUM %d ; ", ASS_SIGNAL_NUM);
            fprintf(stderr, "\n");
            for(NOVA_SV_FINAL_RST_item &indel:result_NGS_INDEL){
                if(indel.SV_overlap(result_BGS_INV.back()) && indel.get_length() >= 50){
                    fprintf(stderr, "INV is removed because it overlap with other INS/DELs\n");
                    if(result_BGS_INV.back().writeVCF_final(&vcfBuffStr, NGS_read.file._hdr, NULL))
                        fprintf(stderr, "%s", vcfBuffStr.s);
                    if(indel.writeVCF_final(&vcfBuffStr, NGS_read.file._hdr, NULL))
                        fprintf(stderr, "%s", vcfBuffStr.s);
                    if(!result_BGS_INV.empty())
                    	result_BGS_INV.pop_back();
                    break;
                }
            }
        }else{
            fprintf(stderr, "GT balance %d:%d:%d %d:%d:%d ", GT_result[0][0], GT_result[0][1], GT_result[0][2], GT_result[1][0], GT_result[1][1], GT_result[1][2]);
            fprintf(stderr, "ASS BALANCE: ASS_SIGNAL_NUM %d ; ", ASS_SIGNAL_NUM);
            fprintf(stderr, "\n");
        }
	}

	void NGS_INV_CALLING_AND_GT_PROCESS(){
    	am->setRepeatMode();
    	fprintf(stderr, "INV calling process begin\n" );
       	//signals combine and store data:
    	std::vector<SVE> inv_sve_l;
    	INV_signal_combine(inv_sve_l);
        for(SVE &c_sve : inv_sve_l){
        	int ASS_SIGNAL_NUM = 0;
        	if(INV_ASSEMBLY(c_sve, ASS_SIGNAL_NUM) == false)
        		continue;
        	INV_GT(c_sve, ASS_SIGNAL_NUM);
        }
        fprintf(stderr, "INV processing Done\n");
	}

	void INV_WRITE(){
        //POST_PROCESS_INV();
        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_BGS_INV){
            fprintf(stderr, "OUTPUT BND::INV\n: \n ");
            if(sv.writeVCF_final(&vcfBuffStr, NGS_read.file._hdr, &global_SV_ID))
                fprintf(vcf_w, "%s", vcfBuffStr.s);
        }
	}

	void CCS_SIG_COLLECTION(bool print_log, BAM_handler &cur_read, std::vector<CCS_SIG> &store_sig, std::vector<CCS_SIG> &store_clip){
		Bam_file *c_b = &(cur_read.file);
		cur_read.clear();
		R_region region;
		ref->get_cur_region()->toR_region(region);
		resetRegion_ID(c_b, &region);	//reset region
		read_counter = 0;

		uint32_t total_read_len = 0;
		while (bam_next(c_b)) {
			//new sig for this read
			std::queue<int> small_sig;
			bam1_t *br = &(c_b->_brec);
			if(print_log) fprintf(stderr, "Read Begin NAME %s \n", bam_get_qname(br));
			//filter:
			//basic filter
			if(bam_is_secondary(br))
				continue;
			if(bam_is_supplementary(br))
				continue;
			if (br->core.qual < 4)
				continue;

			get_bam_seq_bin(0, br->core.l_qseq, cur_read.storeReadBuff, br);
			total_read_len += br->core.l_qseq;
			cur_read.storeCCS_ReadCore(br->core.l_qseq, cur_read.storeReadBuff, bam_get_qual(br), br->core.n_cigar, bam_get_cigar(br), br->core.pos);
			int cigar_idx = cur_read.read_list.back().cigar_index;
			int cigar_end = cigar_idx + cur_read.read_list.back().cigar_l;
			//get the reference sequence:
			//uint8_t * tseq = ref->getRefStr(br->core.pos);
			int seq_i = 0;
			int ref_i = 0;
			//uint8_t *qseq = cur_read.storeReadBuff;
			for(int i = cigar_idx; i < cigar_end; i++){
				path_segment *p = cur_read.cigar.a + i;
				if(p->type == align_t::CIGAR_SOFT_CLIP)//soft clip:
				{
					if(p->length > 1000){//at lease 1000 bp
						if(print_log) fprintf(stderr, "CLIP_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						store_clip.emplace_back();
						store_clip.back().store(CCS_SIG_CLIP_RIGHT, br->core.tid, br->core.pos + ref_i, br->core.pos + ref_i);
						if(i == cigar_idx)// in the begin:
							store_clip.back().type = CCS_SIG_CLIP_LEFT;
					}
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_INSERT)
				{
					if(p->length >= (unsigned)MIN_sv_len){
						if(print_log) fprintf(stderr, "BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						store_sig.emplace_back();
						store_sig.back().store(CCS_SIG_INS, br->core.tid, br->core.pos + ref_i, br->core.pos + ref_i);
					}else{
						//store small sigs
						//small_sig_queue_handler(print_log, br->core.tid, small_sig, br->core.pos + ref_i, p->length, store_sig);
					}
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_DELETE)
				{
					if(p->length >= (unsigned)MIN_sv_len){
						if(print_log) fprintf(stderr, "BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						store_sig.emplace_back();
						store_sig.back().store(CCS_SIG_DEL, br->core.tid, br->core.pos + ref_i, br->core.pos + ref_i + p->length);
					}else{
						//store small sigs
						//small_sig_queue_handler(print_log, br->core.tid, small_sig, br->core.pos + ref_i, p->length, store_sig);
					}
					ref_i += p->length;
				}
				else if(p->type == align_t::CIGAR_MATCH
						|| p->type == align_t::CIGAR_SEQ_MATCH
						|| p->type == align_t::CIGAR_SEQ_MISMATCH){
					//analysis the mismatch
					//small_sig_queue_handler(print_log, small_sig, POS, 0);
					ref_i += p->length;
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_HARD_CLIP){
					//DO nothing
				}
				else{
					//DO nothing
				}
			}
		}
	}

	void CCS_CALLING_MIDDLE(bool print_log, std::vector<CCS_SIG> &store_sig, bool isONT){
		fprintf(stderr, "CCS variations calling(NORMAL) process begin\n" );
		//P1: simple merging
		std::sort(store_sig.begin(), store_sig.end(), CCS_SIG::cmp_by_ref_pos);
		std::vector<CCS_SIG> &l= store_sig;
		std::vector<CCS_SIG> cmb_store_tmp;
		//simple combine:
		for(uint i = 0; i < l.size();i++){
			if(l[i].isUsed)
				continue;
			l[i].support_read_Num = 1;
			for(uint j = i+1; j < l.size();j++){
				if(l[i].similar_SIG_Normal(l[j])){
					l[i].POS = MIN(l[i].POS, l[j].POS);
					l[i].END = MAX(l[i].END, l[j].END);
					l[i].support_read_Num ++;
					l[j].isUsed = true;
				}
			}
			cmb_store_tmp.emplace_back(l[i]);
		}
		//simple combine:
		l.swap(cmb_store_tmp);
		//P2: calling (normal)
		am->setRepeatMode();
		int signal_list_size;
		if(true){
			for(CCS_SIG ss:store_sig)
				ss.show_simple();
		}
		signal_list_size = store_sig.size();
		for(int i = 0; i < signal_list_size; i++){
			set_region_addition_load_short();
			if(isONT)
				assembly_variations_ONT(print_log, store_sig[i], false, result_TGS_middle);
			else
				assembly_variations_CCS(print_log, store_sig[i], false, result_TGS_middle);
		}
	}

	void CCS_CALLING_CLIPING(bool print_log, std::vector<CCS_SIG> &store_clip, bool isONT){
		fprintf(stderr, "CCS variations calling(CLIP) process begin\n" );
		//P1： simple combine:
		std::sort(store_clip.begin(), store_clip.end(), CCS_SIG::cmp_by_ref_pos);
		std::vector<CCS_SIG> &l= store_clip;
		std::vector<CCS_SIG> cmb_store_tmp;
		//simple combine:
		for(uint i = 0; i < l.size();i++){
			if(l[i].isUsed)
				continue;
			l[i].support_read_Num = 1;
			for(uint j = i+1; j < l.size();j++){
				if(l[i].similar_SIG_CLIP(l[j])){
					l[i].POS = MIN(l[i].POS, l[j].POS);
					l[i].END = MAX(l[i].END, l[j].END);
					l[i].support_read_Num ++;
					l[j].isUsed = true;
				}
			}
			cmb_store_tmp.emplace_back(l[i]);
		}
		l.swap(cmb_store_tmp);
		//P2：pairing
		if(print_log){
			for(CCS_SIG ss:store_clip)
				ss.show_simple();
		}
		std::sort(store_clip.begin(), store_clip.end(), CCS_SIG::cmp_by_ref_pos_ignore_type);
		//clip pairing
		{//simple pairing:
			std::vector<CCS_SIG> &l= store_clip;
			std::vector<CCS_SIG> cmb_store_tmp;
			//simple combine:
			for(uint i = 0; i < l.size();i++){
				if(l[i].isUsed)
					continue;
				l[i].support_read_Num = 1;
				for(uint j = i+1; j < l.size();j++){
					if(l[i].SIG_CLIP_pairing(l[j])){
						l[i].POS = MIN(l[i].POS, l[j].POS);
						l[i].END = MAX(l[i].END, l[j].END);
						l[i].support_read_Num += l[j].support_read_Num;
						break;
					}
				}
				cmb_store_tmp.emplace_back(l[i]);
			}
			//simple combine:
			l.swap(cmb_store_tmp);
		}
		//P3： calling
		int signal_list_size = store_clip.size();
		for(int i = 0; i < signal_list_size; i++){
			set_region_addition_load_super_super_long();
			if(isONT)
				assembly_variations_ONT(print_log, store_clip[i], true, result_TGS_clip);
			else
				assembly_variations_CCS(print_log, store_clip[i], true, result_TGS_clip);
		}
	}
	void combine_CCS_clip_and_middle(std::vector<NOVA_SV_FINAL_RST_item> &result_TGS_middle, std::vector<NOVA_SV_FINAL_RST_item> &result_TGS_clip){
		//result_TGS_clip.insert(result_TGS_clip.end(), result_TGS_middle.begin(), result_TGS_middle.end());
		//return;
		std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
		//simple combine:
		for(uint i = 0; i < result_TGS_middle.size();i++){
			if(result_TGS_middle[i].SV_length < 1000){
				sv_tmp.emplace_back();
				std::swap(sv_tmp.back(), result_TGS_middle[i]);
			}
			else{
				bool SV_dup = false;
				for(uint j = 0; j < result_TGS_clip.size();j++){
					if(result_TGS_clip[j].SV_length < 1000)
						continue;
					int pos_dis = result_TGS_clip[j].st_pos - result_TGS_middle[i].st_pos;
					pos_dis = ABS(pos_dis);
					float size_dis = result_TGS_clip[j].SV_length - result_TGS_middle[i].SV_length;
					float size_dis_f = size_dis/(result_TGS_middle[i].SV_length);
					if(size_dis_f < 0.05 && size_dis_f > -0.05 && pos_dis < 200){
						SV_dup = true;
						break;
					}
				}
				if(!SV_dup){
					sv_tmp.emplace_back();
					std::swap(sv_tmp.back(), result_TGS_middle[i]);
				}else{
					fprintf(stderr, "SV is removed: ");
					result_TGS_middle[i].printSimple(stderr);
				}
			}
		}
		//simple combine:
		result_TGS_middle.clear();
		result_TGS_clip.insert(result_TGS_clip.end(), sv_tmp.begin(), sv_tmp.end());
	}

	void remove_CCS_duplications(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> &TGS_r_l){
		std::sort(TGS_r_l.begin(), TGS_r_l.end(), NOVA_SV_FINAL_RST_item::cmp_by_position_TGS);
		if(print_log){
			//set basic GT:
	        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: TGS_r_l)
	            sv.calGT_TGS();
			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: TGS_r_l){
				if(sv.writeVCF_final(&vcfBuffStr, cur_header, NULL))
					fprintf(stderr, "\nBefore remove_TGS_duplications for: %s ", vcfBuffStr.s);
			}
		}
		std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
		for(uint i = 0; i < TGS_r_l.size();i++){
			if(TGS_r_l[i].chr_ID == -1)
				continue;
			for(uint j = i+1; j < TGS_r_l.size();j++){
				if(TGS_r_l[i].st_pos != TGS_r_l[j].st_pos)
					break;
				if(TGS_r_l[i].is_same_var(TGS_r_l[j]) && TGS_r_l[i].is_same_supp_read_TGS(TGS_r_l[j])){
					TGS_r_l[j].chr_ID = -1;
				}
			}
			sv_tmp.emplace_back();
			std::swap(sv_tmp.back(), TGS_r_l[i]);
		}
		//simple combine:
		TGS_r_l.swap(sv_tmp);
	}

	void combine_Homozygous_SVs(std::vector<NOVA_SV_FINAL_RST_item> &TGS_r_l){
		std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
		//simple combine:
		for(uint i = 0; i < TGS_r_l.size();i++){
			if(TGS_r_l[i].chr_ID == -1)
				continue;
			for(uint j = i+1; j < TGS_r_l.size();j++){
				if(TGS_r_l[i].st_pos != TGS_r_l[j].st_pos)
					break;
				if(TGS_r_l[i].is_same_var(TGS_r_l[j]) && TGS_r_l[i].is_same_global_region_ID_TGS(TGS_r_l[j])){
					TGS_r_l[i].combine_supp_read_TGS(TGS_r_l[j]);
					TGS_r_l[j].chr_ID = -1;
				}
			}
			sv_tmp.emplace_back();
			std::swap(sv_tmp.back(), TGS_r_l[i]);
		}
		//simple combine:
		TGS_r_l.swap(sv_tmp);
	}

	void debug_code_load_SVs_from_vcf_f(std::vector<NOVA_SV_FINAL_RST_item> & result_sv_l, const char * vcf_fn, RefHandler *ref){
		result_sv_l.clear();
		BCF_FILE vcf_r;//vcf for read
		VCF_open_read(&vcf_r, vcf_fn);//open for read

		char *c_sv_type = (char *)malloc(1000);
		bcf_hdr_t *header = vcf_r.header;
		int *support_read_number = (int*)xcalloc(3,4);
		int32_t SV_LEN = 0;
		RefRegion * r = ref->get_cur_region();
		char *GT[1]; char GT_1[3]; GT[0] = GT_1;
		do{//read one
			bcf1_t *c_r = &( vcf_r.r);
			if(r->chr_ID != c_r->rid || c_r->pos < r->st_pos || c_r->pos > r->ed_pos)
				continue;
			//unpack the vcf data to get the alt string
			bcf_unpack(c_r, BCF_UN_STR);
			vcf_get_sv_GT(vcf_r.header, c_r, GT);
			vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
			vcf_get_sv_LENGTH(vcf_r.header, c_r, &SV_LEN);
			vcf_get_sv_SR(vcf_r.header, c_r, support_read_number);
			bool is_vntr = vcf_get_sv_flag(header, c_r, "VNTR");
			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t\n", GT_1, is_vntr, c_r->rid , c_r->pos, SV_LEN, c_sv_type, c_r->d.allele[0], c_r->d.allele[1]);// chrID+st+len
			if(c_r->d.allele[1][0] == '<' && c_r->d.allele[1][1] == 'D' && c_r->d.allele[1][2] == 'E'){
				std::string s; s += c_r->d.allele[0][0];
				for(int i = 0;i < -SV_LEN; i++)
					s += 'N';
				NOVA_SV_FINAL_RST_item::add_to_vector(result_sv_l, c_r->rid, c_r->pos + 1, c_sv_type, s.c_str(), c_r->d.allele[0], 0,  0, 0, NULL, 0, 0, 0, 0, r->st_pos);
			}else{
				NOVA_SV_FINAL_RST_item::add_to_vector(result_sv_l, c_r->rid, c_r->pos + 1, c_sv_type, c_r->d.allele[0], c_r->d.allele[1],
									0,  0, 0, NULL, 0, 0, 0, 0, r->st_pos);
			}

			result_sv_l.back().isVNTRcallerRst = is_vntr;
			if(*c_r->d.flt == 0)//filter: PASS
				result_sv_l.back().QUAL_presetVNTR = 30;
			else
				result_sv_l.back().QUAL_presetVNTR = 0;
			result_sv_l.back().set_supp_read_NGS(support_read_number[0], support_read_number[1], support_read_number[2]);
			int GT_final = 0;
			if(GT[0][0] == 4) GT_final++;
			if(GT[0][1] == 4) GT_final++;
			result_sv_l.back().setGenotype_directly(GT_final);
			//store data into:
		}while(VCF_next(&vcf_r));
		//close
		bcf_close(vcf_r.file);
	}

	void remove_similar_SVs_in_Target(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> & query, std::vector<NOVA_SV_FINAL_RST_item> & target){
		std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
		//simple combine:
		for(uint i = 0; i < query.size();i++){
			NOVA_SV_FINAL_RST_item * NGS_r = &(query[i]);
			for(uint j = 0; j < target.size();j++){
				NOVA_SV_FINAL_RST_item * TGS_r = &(target[j]);
				//position is nearby:
				int ABS_POS = ABS_U(NGS_r->st_pos, TGS_r->st_pos);
				if(NGS_r->chr_ID != TGS_r->chr_ID || ABS_POS > 300)
					continue;
				//SV type is same and SV length is similar
//						float length_rate = (float)NGS_r->SV_length/TGS_r->SV_length;
//						if(length_rate <= 0.7 || length_rate >= 1.4)
//							continue;
				if(print_log){
					fprintf(stderr, "NGS results is removed:\n");
					NGS_r->printSimple_M2(stderr);
					fprintf(stderr, "When comparing with CCS SV:\n");
					TGS_r->printSimple_M2(stderr);
				}
				NGS_r->chr_ID = -1;
				break;
			}
			if(NGS_r->chr_ID != -1){
				sv_tmp.emplace_back();
				std::swap(sv_tmp.back(), query[i]);
			}
		}
		//simple combine:
		query.swap(sv_tmp);
	}

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

	void re_genotype_NGS_USING_TGS_data(bool print_log){
		//show all the results
		if(print_log){
			fprintf(stderr, "Combine with CCS begin\n\n");
			showSVList(result_TGS_clip, "Show CCS result");
			showSVList(result_NGS_INDEL, "Show NGS result");
		}
		//merging and remove duplication NGS results that similar with CCS results
		remove_similar_SVs_in_Target(print_log, result_NGS_INDEL, result_TGS_clip);
		//
		//genotyping the NGS result using CCS data set
		for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_INDEL ){
			if(print_log) showSV(sv, "GT NGS SVs using CCS reads[Before]:\n"); //show all the SVs
			int TGS_read_number = 0;
			//get the depth of CCS reads
			for(int mode = 0; mode < 2; mode ++){
				int check_position = sv.st_pos;
				if(mode == 1){
					if(sv.SV_length > 0)	continue;
					else					check_position = sv.st_pos - sv.SV_length;
				}
				int check_edge = 500;
				RefRegion r(sv.chr_ID, check_position - check_edge, check_position + check_edge);
				if(print_log) r.show();
				//load all reads iin this regions:
				for(READ_record &rr :CCS_read.read_list){
					std::string s;
					CCS_read.load_read_CCS(s, rr, r.st_pos, r.ed_pos, false);
					if(!s.empty() && s.size() > check_edge * 1.9)//only full cover CCS reads is used
						TGS_read_number++;
				}
				if(mode == 1)
					TGS_read_number = TGS_read_number/2;
			}
			sv.setTGS_INFO(-1, 0, TGS_read_number, 0, -1, -1, -1, -1);
			int max_TGS_read_num = 3;
			if(TGS_read_number > max_TGS_read_num)	sv.setGenotype_directly(0);
			else if(TGS_read_number > 0)			sv.setGenotype_directly(1);
			if(print_log) showSV(sv,  "GT NGS SVs using CCS reads[After]:\n");
		}
		if(print_log) fprintf(stderr, "Show NGS result after merging\n\n");
	}

    void SVE_handle_region(){
        /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART1: Signal Process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
        fprintf(stderr, "Signal process begin\n" );
        clear_sve_signals_in_a_region();//S1
        if(ONT_bamFileName != NULL){//handle ONT dataset
        	return;
        	bool print_log = true;
        	result_TGS_middle.clear();
        	result_TGS_clip.clear();
    		//SIG list:
    		std::vector<CCS_SIG> store_sig;
    		std::vector<CCS_SIG> store_clip;
    		CCS_SIG_COLLECTION(print_log, ONT_read, store_sig, store_clip);
    		if(print_log) {
				fprintf(stderr, "NORMAL_SIG show:\n");	for(CCS_SIG &s : store_sig)  s.show_simple();
				fprintf(stderr, "CLIP_SIG show:\n");    for(CCS_SIG &s : store_clip) s.show_simple();
			}
    		{
    			CCS_CALLING_MIDDLE(print_log, store_sig, true);
    			CCS_CALLING_CLIPING(print_log, store_clip, true);
    		}
        }

        //handle CCS reads
        if(CCS_bamFileName != NULL){
        	bool print_log = false;
        	result_TGS_middle.clear();
        	result_TGS_clip.clear();
    		//SIG list:
    		std::vector<CCS_SIG> store_sig;
    		std::vector<CCS_SIG> store_clip;

        	CCS_SIG_COLLECTION(print_log, CCS_read, store_sig, store_clip);

			if(print_log) {
				fprintf(stderr, "NORMAL_SIG show:\n");	for(CCS_SIG &s : store_sig)  s.show_simple();
				fprintf(stderr, "CLIP_SIG show:\n");    for(CCS_SIG &s : store_clip) s.show_simple();
			}

			//test: collection signals using ONT reads

			CCS_CALLING_CLIPING(print_log, store_clip, false);
			CCS_CALLING_MIDDLE(print_log, store_sig, false);
			//combine CLIPING and MIDDLE, store the final results in the result_TGS_clip
			combine_CCS_clip_and_middle(result_TGS_middle, result_TGS_clip);
			//remove duplication SVs
			remove_CCS_duplications(true, result_TGS_clip);

			if(print_log) showSVList(result_TGS_clip, "\n\n After remove_CCS_duplications:\n");

			//combine same SVs in different haplotype-s
			//0/1 + 0/1 ---> 1/1
			combine_Homozygous_SVs(result_TGS_clip);

			if(print_log) showSVList(result_TGS_clip, "\n\n After combine_Homozygous_SVs:\n");

			//set basic GT:
	        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_TGS_clip)
	            sv.calGT_TGS();

			//re-genotyping CCS result using NGS reads when data available
	        if(NGS_bamFileName != NULL && !result_TGS_clip.empty() ){
				fprintf(stderr, "genotypingTGS_Using_NGS begin:\n");
				for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_TGS_clip){
					if(ABS(sv.SV_length) < MIN_sv_len)//only re-genotype SVs NOT shorter than MIN_sv_len
						continue;
					showSV(sv, "\nTGS re-genotype using NGS, Genotyping for:");
					int new_GT = sv.genotypingTGS_Using_NGS(false, sig_para.MaxReadLen, &NGS_read.file, &ga, ref);
					if(sv.TGS_need_filter_by_NGS() && sv.contig_not_support_by_NGS())	sv.setGenotype_directly(0);//set to 0/0
					else if(sv.TGS_need_reGT_by_NGS() && new_GT > 0 && new_GT <= 2)		sv.setGenotype_directly(new_GT);
				}
			}else fprintf(stderr, "genotypingTGS_Using_NGS SKIP\n\n");

        	if(random_phasing){
    	        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_TGS_clip)
    	        	sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
        	}
			//output the final results
	        for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_TGS_clip){
	            if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID))
	                fprintf(vcf_w, "%s", vcfBuffStr.s);
	        }
        }

        if(NGS_bamFileName != NULL){//debug code: loading NGS vcfs, mergin and output:
            //handle NGS reads: clear result list:
			result_NGS_INDEL.clear();//store the result for insertion, deletion, large tra-insertion
			result_NGS_BND.clear(); //store the result for the INV and the BNDs
			result_BGS_INV.clear();
			//calling SVs using NGS reads
//        	if(CCS_bamFileName != NULL){//DEBUG mode: when with CCS dataset, run in simple mode ,30X, simple merge:
//        		const char * vcf_fn = NULL;
//        		if(false) vcf_fn = "/home/user/zhanganqi/gyli/GC_SV2_1/CALL/HG002_30X.bench.vcf";
//        		else     vcf_fn = "/home/user/zhanganqi/gyli/GC_SV2_1/CALL/HG002_60X.bench.vcf";
//        		debug_code_load_SVs_from_vcf_f(result_NGS_INDEL, vcf_fn, ref);
//        	}
//        	else
        	{  //pure NGS mode: truly handle NGS reads
                //get signals from NGS CRAM
                GET_NGS_SIGNALS_FROM_READs();//S2:load reads and get signals
                //not used right now
                if(false)load_realignment_read_from_file();//S3:try to load realignment reads into UNMAPPED read list
                CLUSTERING_AND_COMBINE_NGS_signals();//S4: cluster signals and get SV signals
                if(true)
                	showNGS_signals();            //log: print signals
                /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART1: INDEL calling(assembly and genotyping) process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
                //results is stored in "result_INDEL"
                //call small INDELs
                NGS_SMALL_VAR_CALLING_PROCESS();/* SMALL INDEL calling process*/
                NGS_VNTR_VAR_CALLING_PROCESS();/* VNTR calling process*/
                NGS_LARGE_DELETION_CALLING_PROCESS();/* LARGE DELETION calling process*/

                REMOVE_DUPLICATED_INDEL();/*sort and removing duplication process(for INS and DEL)*/
                NGS_GENOTYPING_INDELS();/*GENOTYPING and output for INS/DEL*/
                /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART3: INV calling process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
    			NGS_INV_CALLING_AND_GT_PROCESS();
                /*～～～～～～～～～～～～～～～～～～～～～～～～～～PART2: BND calling(assembly and genotyping) process～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～～*/
                NGS_BND_CALLING_AND_GT_PROCESS();
        	}

            //remove SVs called by NGS data-set that is similar with CCS results
            if(!result_NGS_INDEL.empty() && CCS_bamFileName != NULL){
            	bool print_log = false;
            	re_genotype_NGS_USING_TGS_data(print_log);
            }

        	if(random_phasing){
    			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_INDEL ){
    				sv.randomly_phasing_Genotype(((float)global_SV_ID*1.897239));
    			}
        	}

            //Output all the NGS results
            {
            	bool print_log = true;
    			for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_INDEL ){
    				if(sv.writeVCF_final(&vcfBuffStr, cur_header, &global_SV_ID)){
    					if(print_log) fprintf(stderr, "%s", vcfBuffStr.s);
    					fprintf(vcf_w, "%s", vcfBuffStr.s);
    				}
    			}
    			INV_WRITE();
                BND_WRITE();//BND write
            }
        }
    }


private:
	/*************************************BASIC FUNCTIONs**********************************************/
	void clear_sve_signals_in_a_region(){//used only within function "process"
		for(int i = 0; i < SIG::LEN; i++) for(int j = 0; j < SV::LEN; j++) sve[i][j].clear();
		sveVNTR.clear();
		NGS_read.clear();
	}

	SveHandler (SveHandler &B);//SveHandler can`t be copied
	/*************************************SIGNALS FUNCTIONs**********************************************/
	//S1: get signals from all reads, in this step, a read pair may be used as DR or SA or UM signals, for each read signals, a SVE will be stored //signals functions
	bool is_mate_fwd(uint16_t flag) { return (!((flag & BAM_MATE_STRAND) != 0));}
	void GET_NGS_SIGNALS_FROM_READs();
	void load_realignment_read_from_file();
	void handleDRSignal(bam1_core_t *core, int middle_size);
	void storeMismatchSignals(bam1_t *br, READ_record &c_r);
	void storeClipSignals(bool isClipAtRight, uint32_t ori_pos, uint8_t read_mapq);
	void handleSASignal(READ_record& c_r, bam1_t *br);//handle a SA signal and store results
	void handleUMSignal(bam1_t *br);

	/*************************************CLUSTERING FUNCTIONs**********************************************/
	//S2: combine signal from multiple reads. In this steps, signals from multiple reads will be clustered and joint together. At last, signals combined from many reads(SVEs) will form.
	void CLUSTERING_AND_COMBINE_NGS_signals();
	int  getTopPossibilityIdx(int r_min, int r_max, SVE_L & l, bool isR1, bool is_forward, float &max_poss, std::vector<float> &dis_bp_percent);//used in sve_combine_STEP1_DR
	void single_type_sve_combine(bool print_log, SVE_L & l, int min_score, SIG::T sigT, SV::T svt);
	void combine_duplication(SVE_L & l);
	void DR_SH_signal_combining(SVE_L &SH, SVE_L &DR);//S3: Joint Calling for single sv, multi-sample and multi-signal type //r1 and r2 of SH must be within r1 and r2 of DR

	/*************************************ASSEMBLY FUNCTIONs**********************************************/
	//S4: read assembly
	bool SV_region_depth_filter(SVE & c_sve); //return whether pass assembly filter
	bool assembly_load_read(bool print_log, SVE &sv, RefRegion main, RefRegion supp);

	void assembly_load_CCS_read(bool print_log, RefRegion &main, Ass_Block & a_b, bool using_clip);
	bool assembly_load_read_ALL1(bool print_log, RefRegion r1, RefRegion r2, RefRegion main,
			RefRegion supp, int minQUAL, uint MAX_load_reads, bool usingTGS, bool usingNGS, bool using_clip);
	bool assembly_variations_ALU_AND_LONG_SV(SVE &sv, bool call_ALU_TL, bool &withRepeat);
	bool assembly_variations_INS_DEL_VNTR(SVE &sv, bool usingTGS, bool usingNGS);
	bool assembly_variations_ONT(bool print_log, CCS_SIG &sig, bool is_clip, std::vector<NOVA_SV_FINAL_RST_item> & SV_result_final);
	bool assembly_variations_CCS(bool print_log, CCS_SIG &sig, bool is_clip, std::vector<NOVA_SV_FINAL_RST_item> & SV_result_final);
	void force_calling_ALU_ins(bool print_log, SVE &sv, std::vector<AssemblyContig> &contigs, int main_st_pos);
	void combine_repeat_tail_of_contigs(bool print_log, int ori_read_number, std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, std::vector<AssemblyContig> &contigs);\
	void get_BP_from_tran_based_ins(bool print_log, std::vector<Tran_ins_suggest_pos> &magic_pos_set, std::vector<AssemblyContig> &contigs, std::vector<TRANS_INS_info> &bp_l);
	void handle_tran_based_ins(bool print_log, int ori_read_number, std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, std::vector<AssemblyContig> &contigs);
	bool read_cover_repeat_filter(bool print_log, AssemblyContig &contig, uint32_t SV_size_ori,
			std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, int ori_read_number);

	bool assembly_and_get_breakpoints_TRA(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am,
			int main_tid, int &main_pos, bool main_read_before_BP, bool main_read_forward,
			int supp_tid, int &supp_pos, bool supp_read_before_BP, bool supp_read_forward);
	int assembly_and_get_breakpoints_INV(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am, RefRegion &main_region,
			RefRegion &supp_region, int &main_BP, int &supp_BP);
	bool assembly_and_genetyping_BND(SVE &sv,
			bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward, bool main_ref_is_forward,
			bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward, bool supp_ref_is_forward,
			std::vector<NOVA_SV_FINAL_RST_item> &region_SVs, bam_hdr_t * header);//buff to store contig
	bool INV_genotyping_and_store(SVE &sv,
			bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward,
			bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward,
			bam_hdr_t * header, int *GT_result);//buff to store contig
	//S4.1: assembling contigs re-alignment
	void get_suggention_alignment_position_list(AssemblyContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f );
	void getSuggestSVlength(AssemblyContig &contig);
	void alignment_and_get_var(bool print_log, int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len, int ref_region_length);

	void set_region_addition_load_short(){
		region_addition_load = 500;
		ca.setZdrop(1200, 1200);
	}

	void set_region_addition_load_long(){
		region_addition_load = 2000;
		ca.setZdrop(2000, 2000);
	}

	void set_region_addition_load_extramly_long(){
		region_addition_load = 5000;
		ca.setZdrop(5000, 5000);
	}

	void set_region_addition_load_super_super_long(){
		region_addition_load = 50000;
		ca.setZdrop(50000, 50000);
	}

	//abandoned functions
	void jointAssemblingUmReads();	//abandoned
	void umAssembly(int beginIdx);	//abandoned
	int getTopPossibilityIdx_UM(int r_min, int r_max, float &max_poss);	//abandoned
	/*************************************VCF writing FUNCTIONs**********************************************/
};


struct REF_COMBINE{
	std::vector<uint8_t> s;
	int BP_region_st[2];
	int BP_region_ed[2];

	int s1_ref_st;
	int s2_ref_st;

	bool s1_forward;
	bool s2_forward;

	bool is_reverse;
	bool with_contig_supp;

	int main_pos;
	int supp_pos;

	int get_len(){ if(with_contig_supp) return supp_pos- main_pos; else return 0;}
	void print_ref(){
		//fprintf(stderr, "region_overlap %d combined_ref is: \n", region_overlap);
		fprintf(stderr, "combined_ref is: \n");
		for(uint i = 0; i < s.size(); i++)
			fprintf(stderr, "%c",  "ACGTNN"[s[i]]);
		fprintf(stderr, "\n");
	}

	void print_results(){
		if(with_contig_supp)
			fprintf(stderr, "Final main_pos %d supp_pos %d, length %d \n", main_pos, supp_pos, supp_pos- main_pos);
		else
			fprintf(stderr, "NO data \n");
	}

	void store_ref(
			std::vector<uint8_t> & s1, bool s1_forward, int s1_ref_st,
			std::vector<uint8_t> & s2, bool s2_forward, int s2_ref_st,
			bool is_reverse
			){
		s.clear();
		BP_region_st[0] = s.size();
		s = s1;
		BP_region_ed[0] = s.size();
		s.emplace_back(4);
		s.emplace_back(4);
		s.emplace_back(4);
		BP_region_st[1] = s.size();
		s.insert(s.end(), s2.begin(), s2.end());
		BP_region_ed[1] = s.size();

		this->s1_ref_st = s1_ref_st;
		this->s2_ref_st = s2_ref_st;
		this->s1_forward = s1_forward;
		this->s2_forward = s2_forward;

		this->is_reverse = is_reverse;
		this->with_contig_supp = false;
		if(is_reverse){
			std::vector<uint8_t> s_rev;
			for(int i = 0; i < (int)s.size(); i++){
				uint8_t c = s[s.size() - i - 1];
				if(c > 3)
					s_rev.emplace_back(4);
				else
					s_rev.emplace_back(3- s[s.size() - i - 1]);
			}

			std::swap(s,s_rev);
			this->s1_forward = !this->s2_forward;
			this->s2_forward = !this->s1_forward;
			std::swap(this->s1_ref_st,this->s2_ref_st);
		}
	}

	void check_and_store(int BP1, int BP2){
		if(BP1 >= BP_region_st[0] && BP1 <= BP_region_ed[0] && BP2 >= BP_region_st[1] && BP2 <= BP_region_ed[1] ){
			BP2 -= BP_region_st[1]; //align the BP2 to base-0
			main_pos = (s1_forward)?(s1_ref_st + BP1):(s1_ref_st + BP_region_ed[0] - BP_region_st[0] - BP1);
			supp_pos = (s2_forward)?(s2_ref_st + BP2):(s2_ref_st + BP_region_ed[1] - BP_region_st[1] - BP2);
			if(main_pos > supp_pos) std::swap(main_pos, supp_pos);
			fprintf(stderr, "\n After is main_pos(M1) %d supp_pos %d, length %d s1_ref_st %d , s2_ref_st %d , s1_forward %d , s2_forward %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, supp_pos- main_pos, s1_ref_st, s2_ref_st, s1_forward, s2_forward, BP1, BP2);
			with_contig_supp = true;
			fprintf(stderr, "check_and_store Success\n");
		}
	}

	bool check_and_store_TRAN_INS(int BP1, int BP2, bool is_right_part, int &source_BP_pos, int &target_BP_pos){
		if(BP1 >= BP_region_st[0] && BP1 <= BP_region_ed[0] && BP2 >= BP_region_st[1] && BP2 <= BP_region_ed[1] ){
			BP2 -= BP_region_st[1]; //align the BP2 to base-0
			main_pos = (s1_forward)?(s1_ref_st + BP1):(s1_ref_st + BP_region_ed[0] - BP_region_st[0] - BP1);
			supp_pos = (s2_forward)?(s2_ref_st + BP2):(s2_ref_st + BP_region_ed[1] - BP_region_st[1] - BP2);
			if(is_right_part){ source_BP_pos = main_pos; target_BP_pos = supp_pos;}
			else			 { target_BP_pos = main_pos; source_BP_pos = supp_pos;}

			if(false){
				if(is_right_part)	fprintf(stderr, "\n R(S+T): source_pos %d target_pos %d, s1_ref_st %d , s2_ref_st %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, s1_ref_st, s2_ref_st, BP1, BP2);
				else				fprintf(stderr, "\n L(T+S): target_pos %d source_pos  %d, s1_ref_st %d , s2_ref_st %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, s1_ref_st, s2_ref_st, BP1, BP2);
				fprintf(stderr, "check_and_store Success\n");
			}
			return true;
		}
		return false;
	}

	void chece_and_store_M2(int BP1, int BP2){
		int region_id = 0;
		if(BP1 <= BP_region_ed[0] && BP2 <= BP_region_ed[0]){
			region_id = 0;
		}
		else if(BP1 >= BP_region_st[1] && BP2 >= BP_region_st[1]){
			region_id = 1;
			BP1 -= BP_region_st[1]; //align the BP2 to base-0
			BP2 -= BP_region_st[1]; //align the BP2 to base-0
		}
		else
			return;
		bool is_forward = (region_id == 0)?s1_forward:s2_forward;
		int ref_st = (region_id == 0)?s1_ref_st:s2_ref_st;
		int region_len = (region_id == 0)?(BP_region_ed[0] - BP_region_st[0]):(BP_region_ed[1] - BP_region_st[1]);
		main_pos = (is_forward)?(ref_st + BP1):(ref_st + region_len - BP1);
		supp_pos = (is_forward)?(ref_st + BP2):(ref_st + region_len - BP2);
		if(main_pos > supp_pos) std::swap(main_pos, supp_pos);
		fprintf(stderr, "After is main_pos(M2) %d supp_pos %d, length %d s1_ref_st %d , s2_ref_st %d , s1_forward %d , s2_forward %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, supp_pos- main_pos, s1_ref_st, s2_ref_st, s1_forward, s2_forward, BP1, BP2);
		with_contig_supp = true;
		fprintf(stderr, "chece_and_store Success\n");
	}
};

bam_hdr_t* checkBamHeadersAreSAME(std::vector<SveHandler>& h);
void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
void ref_load_char_2_bin(bool cur_ref_is_forward, int load_len, std::vector<uint8_t> & to_bin_string, char* cur_load_ref);
void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig);
void store_bin_from_char(const char * contig_seq, int contig_seq_len, std::vector<uint8_t> &bin_contig);

#endif /* SRC_SIGNAL_SVEHANDLER_HPP_ */
