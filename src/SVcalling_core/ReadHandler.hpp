/*
 * BamHandler.hpp
 *
 *  Created on: 2020年4月28日
 *      Author: fenghe
 */

#ifndef BAMHANDLER_HPP_
#define BAMHANDLER_HPP_

extern "C"
{
	#include "../clib/bam_file.h"
	#include "../clib/vcf_lib.h"
extern uint64_t kmerMask[33];
}

#include"../SVcalling_core/RefHandler.hpp"
#include "../cpp_lib/Assembler/GreedyAssembler.hpp"
#include"../SVcalling_core/sve.hpp"
#include<vector>
#include<algorithm>

struct READ_record{
	READ_record(uint16_t flag_, uint8_t direction_, uint64_t position_, int read_l_, int cigar_l_,
			int read_index_, int quality_index_, int cigar_index_,	int soft_left_, int soft_right_, int NM_NUM_) :
			flag(flag_), direction(direction_), position(position_), read_l(read_l_), cigar_l(cigar_l_),
			read_index(read_index_), quality_index(quality_index_), cigar_index(cigar_index_),
			soft_left(soft_left_), soft_right(soft_right_), NM_NUM(NM_NUM_) {
	}
	//basic
	uint16_t 	flag;
	uint8_t 	direction; // the direction of the read, and will be direction of mate when reads are unmapped
	int		 	position; // the position of the read, and will be position of mate when reads are unmapped
	int 		read_l;
	int 		cigar_l;
	int 	 	read_index;
	int 		quality_index;
	int  		cigar_index;
	int soft_left;
	int soft_right;
	int NM_NUM;

	//std::sort(l.begin(), l.end(), SVE::cmp_by_position);
	static inline int cmp_position(const READ_record &a, const READ_record &b){
		return a.position < b.position;
	}

	void show(){
		fprintf(stderr, "Read is: [%d %d %d]\t", flag, direction, position);
		fprintf(stderr, "[len: %d CIGAR_L: %d %d]\t", read_l, cigar_l, read_index);
		fprintf(stderr, "[%d %d]\t", quality_index, cigar_index);
		fprintf(stderr, "[SL %d SR %d NM %d]\n", soft_left, soft_right, NM_NUM);
	}
};

namespace Read_type{
	enum T {SR, DR, UM, TL,TGS, unknown};
}

struct ASS_reads_info{
	ASS_reads_info(int read_in_ref_offset_, int read_st_pos_, int read_ed_pos_, int read_list_index_, int soft_left_, int soft_right_, int number_NM_, Read_type::T signal_type_){
		signal_type = signal_type_;
		read_in_ref_offset = read_in_ref_offset_;
		read_list_index = read_list_index_;
		read_st_pos = read_st_pos_;
		read_ed_pos = read_ed_pos_;
		soft_left = soft_left_;
		soft_right = soft_right_;
		number_NM = number_NM_;
	}
	int read_list_index;
	int read_in_ref_offset;
	int read_st_pos;//read_in_ref_offset - left_clip
	int read_ed_pos;//read_in_ref_offset + middle_size
	int soft_left;
	int soft_right;
	int number_NM;
	Read_type::T signal_type;
	void print(FILE* log_file, int first_read_ID){
		fprintf(log_file, "[type: ");
		switch (signal_type) {
			case Read_type::SR: fprintf(log_file, "SR"); break;
			case Read_type::DR: fprintf(log_file, "DR"); break;
			case Read_type::UM: fprintf(log_file, "UM"); break;
			case Read_type::TL: fprintf(log_file, "TL"); break;
			case Read_type::TGS: fprintf(log_file, "TGS"); break;
			default:  fprintf(log_file, "unknown"); break;
		}
		fprintf(log_file, " id: %d(%d) {%d~%d}, ref pos %d, {%d, %d, %d}]\t",
				read_list_index - first_read_ID, read_list_index, read_st_pos, read_ed_pos, read_in_ref_offset, soft_left, soft_right, number_NM);
	}
};

struct Ass_Block{
	//input read info
	std::vector<ASS_reads_info> ass_read_list;
	//input read data
	std::vector<std::string> reads;
	//output: read result
	std::vector<AssemblyContig> contigs;

	void clear(){
		ass_read_list.clear();
		reads.clear();
	}

#define MIN_BASE_QUAL 6
	void add_read_bin(Read_type::T t, uint8_t * read_str, uint8_t * qual_str, int read_len, int read_in_ref_offset_, int left_clip, int right_clip, int NM_num, int read_list_index_){
		int read_st_pos = read_in_ref_offset_ + left_clip;
		int read_ed_pos =  read_in_ref_offset_ + read_len - right_clip;
		ass_read_list.emplace_back(read_in_ref_offset_, read_st_pos, read_ed_pos, read_list_index_, left_clip, right_clip, NM_num, t);
		reads.emplace_back(); std::string &s_store = reads.back(); //string to store: string pointer type
		s_store.resize(read_len);
		int last_N_base = -1;
		for(int i = 0; i < read_len; i++){
			s_store[i] = ("ACGTN"[read_str[i]]);//store binary to acgt
			if(qual_str[i] < MIN_BASE_QUAL){
				s_store[i] = 'N';
				if(i - last_N_base < 5)
					for(int j = last_N_base+1;j<i;j++) s_store[j] = 'N';
				last_N_base = i;
			}
		}
	}

	void add_read_bin_compact(Read_type::T t, uint8_t * bam_seq, uint8_t * qual_str,
			int read_len, int read_in_ref_offset_, int left_clip, int right_clip, int NM_num, int read_list_index_, int minQUAL){
	    int read_st_pos = read_in_ref_offset_ + left_clip;
	    int read_ed_pos =  read_in_ref_offset_ + read_len - right_clip;
	    ass_read_list.emplace_back(read_in_ref_offset_, read_st_pos, read_ed_pos, read_list_index_, left_clip, right_clip, NM_num, t);
	    reads.emplace_back(); std::string &s_store = reads.back(); //string to store: string pointer type
	    s_store.resize(read_len);
	    int last_N_base = -1;
	    for(int i = 0; i < read_len; i++){
	      uint8_t b_c = bam_seqi(bam_seq, i);
	      char c_c;
	      switch(b_c)
	      {
	      case 1: c_c ='A'; break;
	      case 2: c_c ='C'; break;
	      case 4: c_c ='G'; break;
	      case 8: c_c ='T'; break;
	      default: c_c ='N'; break;
	      }
	      s_store[i] = (c_c);//store binary to acgt
	      if(qual_str[i] < minQUAL){
	        s_store[i] = 'N';
//	        if(i - last_N_base < 5)
//	          for(int j = last_N_base+1;j<i;j++) s_store[j] = 'N';
//	        last_N_base = i;
	      }
	    }
	}

	//adding trans located reads(or trans-located reads) into assembly read list, compared with normal reads, those data has no position or cigar information
	//besides, those reads are stored in CHAR
	void add_read_TL(char * read_str, uint8_t * qual_str, int read_len, bool store_in_reverse, int read_list_index_, int tid, int pos){
		for(int i = 0; i < read_len; i++)
			if(qual_str[i] < MIN_BASE_QUAL)
				read_str[i] = 'N';
		if(store_in_reverse)
			getReverseStr_char(read_str, read_len);

		ass_read_list.emplace_back(-1, tid, pos, read_list_index_, -1, -1, store_in_reverse, Read_type::TL);
		reads.emplace_back(read_str);
	}

	//adding TGS(CCS) reads into assembly read list
	void add_read_CCS(std::string & read_str, int read_list_index_, int tid, int pos){
		ass_read_list.emplace_back(-1, tid, pos, read_list_index_, -1, -1, 0, Read_type::TGS);
		reads.emplace_back(read_str);
	}

	void run_assembly(MainAssemblyHandler *am){
		//store read string
		std::swap(am->reads, reads);
		//store read info:
		uint read_size = ass_read_list.size();
		am->readInfo.clear();
		am->readInfo.resize(read_size);
		for(uint i = 0; i < read_size; i++){
			if(ass_read_list[i].signal_type == Read_type::TGS){
				am->readInfo[i].base_depth = 3;//todo::
			}
			else{
				am->readInfo[i].base_depth = 1;//todo::
			}
		}
		//run assembly
		am->assembley();
		//reset read strings
		std::swap(am->contigs, contigs);
		std::swap(am->reads, reads);//swap back the read list
	}
};

//used to store trans-location read pairs:
struct trans_Read_item{
	int32_t tid;
	int32_t mtid;
	int32_t pos;
	int32_t mpos;
	uint8_t flag;

	trans_Read_item(	int32_t tid_, int32_t mtid_, int32_t pos_, int32_t mpos_, uint8_t flag_){ tid = tid_, mtid = mtid_, pos = pos_, mpos = mpos_, flag = flag_; }
	trans_Read_item(	bam1_t *br){ tid = br->core.tid, mtid = br->core.mtid, pos = br->core.pos, mpos = br->core.mpos, flag = br->core.flag; }
	//copy:
	trans_Read_item(const trans_Read_item &b){ memcpy(this, &b, sizeof(trans_Read_item)); }
	//sort mtid:
	static inline int cmp_by_mtid(const trans_Read_item &a, const trans_Read_item &b){
		if(a.mtid == b.mtid)	return a.mpos < b.mpos;
		else							return a.mtid < b.mtid;
	}

	void prinf(FILE * log){
		fprintf(stderr, "[cur:%d:%d ; mate: %d:%d]\n", tid, pos, mtid, mpos );
	}
};

struct TransReadLoader{

	void search_read_bg_idx(int st_pos){
		int begin_idx = last_read_end_index;
		cur_bg_idx =  MAX_int32t;
//fprintf(stderr, "S411 %ld st_pos %d cur_bg_idx %d\n", trans_list.size(), st_pos, cur_bg_idx);
		if(trans_list.size() == 0) {
//			fprintf(stderr, "S4111\n");
			return;
		}
//fprintf(stderr, "S412\n");
		int trans_list_size =  (int)trans_list.size();
		if(begin_idx >= trans_list_size){ begin_idx = trans_list_size - 1; }
		bool search_forward = (trans_list[begin_idx].pos < st_pos);
//fprintf(stderr, "S413\n");
		if(search_forward){
			for(;begin_idx < trans_list_size && trans_list[begin_idx].pos < st_pos; begin_idx++);
		}else
			for(;begin_idx >= 0 && trans_list[begin_idx].pos > st_pos; begin_idx--);
//fprintf(stderr, "S414\n");
		begin_idx = MAX(begin_idx, 0);
		cur_bg_idx = begin_idx;
	}

	void search_read_ed_idx(int ed_pos){
//fprintf(stderr, "S41_1_1, cur_bg_idx\n", cur_bg_idx);
		int end_idx = cur_bg_idx + 1;
		if(end_idx < 0) end_idx = 0;
		int trans_list_size =  trans_list.size();
//fprintf(stderr, "S41_1_2, trans_list_size %d end_idx %d \n", trans_list_size, end_idx);
		for(;end_idx < trans_list_size && trans_list[end_idx].pos < ed_pos; end_idx++);
//fprintf(stderr, "S41_1_3\n");
		cur_ed_idx = end_idx;
		last_read_end_index = cur_ed_idx;
	}

	void bam_aligned_clip_analysis(bam1_t *b, int *clip_left, int *clip_right){
		*clip_left = 0;
		*clip_right = 0;
		uint32_t* bam_cigar = bam_get_cigar(b);
		//NM = MIS + INS + DEL
		int begin_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
		int end_type   = (int)(1 + (bam_cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK));
		*clip_left = 0; *clip_right = 0;
		if(	begin_type == CIGAR_SOFT_CLIP || begin_type == CIGAR_HARD_CLIP)
			*clip_left = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
		if(	end_type == CIGAR_SOFT_CLIP ||	end_type == CIGAR_HARD_CLIP)
			*clip_right = (bam_cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
	}

	int load_read(bool print_log, int tid, int st_pos, int ed_pos, Ass_Block & ab, RefHandler *ref, int MAX_load_read_number){
		if(false){//debug code
			load_read_M1(print_log, tid, st_pos, ed_pos, ab, ref, MAX_load_read_number);
			load_read_M2(print_log, tid, st_pos, ed_pos, ab, ref, MAX_load_read_number);
			ab.clear();
			load_read_M1(print_log, tid, st_pos, ed_pos, ab, ref, MAX_load_read_number);
			ab.clear();
			return 0;
		}

		if(TL_bam != NULL) return load_read_M1(print_log, tid, st_pos, ed_pos, ab, ref, MAX_load_read_number);
		else   			   return load_read_M2(print_log, tid, st_pos, ed_pos, ab, ref, MAX_load_read_number);
	}

	int load_read_M1(bool print_log, int tid, int st_pos, int ed_pos, Ass_Block & ab, RefHandler *ref, int MAX_load_read_number){
		//search reads in reads list:
		search_read_bg_idx(st_pos);
		search_read_ed_idx(ed_pos);
		fprintf(stderr, "\nload_read(M1) BEGIN\n");
		fprintf(stderr, "[TL_LOADING: st_pos  %d ed_pos %d]\t", st_pos, ed_pos);
		//copy and sort
		region_trans_list.clear();
		fprintf(stderr, "[TL_LOADING: cur_ed_idx  %d cur_bg_idx %d]\t", cur_ed_idx, cur_bg_idx);
		//if(cur_ed_idx - cur_bg_idx <= 3){return 0;}
		if(cur_ed_idx - cur_bg_idx > MAX_load_read_number)
			return 0;

		std::set<uint64_t> read_pos_list;
		for(int i = cur_bg_idx; i < cur_ed_idx; i++){
			uint64_t magic = (trans_list[i].pos); magic <<= 32; magic += trans_list[i].mpos;
			read_pos_list.emplace(magic);
		}

		R_region region;
		region.chr_ID = tid;
		region.st_pos = st_pos;
		region.ed_pos = ed_pos;
		int total_load_read = 0;
		resetRegion_ID(TL_bam, &region);	//reset region
		int UM_read_id = 0;
		fprintf(stderr, "load_read(M1) BEGIN\n");
		while (bam_next(TL_bam)) {
			bam1_t *br = &(TL_bam->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))       	continue;
			{
				uint64_t magic = (br->core.pos); magic <<= 32; magic += br->core.mpos;
				if(read_pos_list.find(magic) == read_pos_list.end()){
					continue;
				}
			}

			//reverse the TID and POS
			std::swap(br->core.pos, br->core.mpos);
			std::swap(br->core.tid, br->core.mtid);

			//load reads
			int read_len = br->core.l_qseq;
			get_bam_seq(0, read_len, storeReadBuff, br);
			//if dir == mdir: reverse the string:
			bool b_forward = bam_is_fwd_strand(br);
			bool m_forward = bam_is_mate_fwd_strand(br);
			bool store_in_reverse = (b_forward == m_forward);
			int clip_left, clip_right;
			bam_aligned_clip_analysis(br, &clip_left, &clip_right);
			int32_t r_pos = br->core.pos - clip_left;
			total_load_read++;
			ab.add_read_TL(storeReadBuff, bam_get_qual(br), read_len, store_in_reverse, UM_read_id++, br->core.tid, r_pos);
			if(false){
				//debug Code:
				char * bam_name = (char *)bam_qname(br);//string values(;, ref)
				fprintf(stderr, "[TL_LOADING(M1): :%s %d:%d @QUAL %d ; mate: %d:%d dir:(%c, %c) %s] ", bam_name, br->core.tid, br->core.pos, br->core.qual, br->core.mtid, br->core.mpos, b_forward?'F':'R', m_forward?'F':'R', storeReadBuff);
				uint32_t* bam_cigar = bam_get_cigar(br);
				uint32_t n_cigar = br->core.n_cigar;
				for (uint i = 0; i < n_cigar; ++i)
				{
					int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
					int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
					fprintf(stderr, "%d%c", c_size, "NMIDKSHPMX"[c_type]);
				}
				fprintf(stderr, "\n");
			}
		}
		fprintf(stderr, "total_load_read %d\n", total_load_read);
		return total_load_read;
	}

	int load_read_M2(bool print_log, int tid, int st_pos, int ed_pos, Ass_Block & ab, RefHandler *ref, int MAX_load_read_number){
		//search reads in reads list:
		search_read_bg_idx(st_pos);
		search_read_ed_idx(ed_pos);
		fprintf(stderr, "\nload_read(M2) BEGIN\n");
		fprintf(stderr, "[TL_LOADING: st_pos  %d ed_pos %d]\t", st_pos, ed_pos);
		//copy and sort
		region_trans_list.clear();
		fprintf(stderr, "[TL_LOADING: cur_ed_idx  %d cur_bg_idx %d]\t", cur_ed_idx, cur_bg_idx);
		if(cur_ed_idx - cur_bg_idx <= 3){
			return 0;
		}
		if(cur_ed_idx - cur_bg_idx > MAX_load_read_number)
			return 0;
		for(int i = cur_bg_idx; i < cur_ed_idx; i++)
			region_trans_list.emplace_back(trans_list[i]);
		std::sort(region_trans_list.begin(), region_trans_list.end(), trans_Read_item::cmp_by_mtid);
		//cluster:
		int region_size = region_trans_list.size();
		int total_load_read = 0;
		int check_region_size = 0;
		int check_read_num = 0;
		for(int bg_idx = 0; bg_idx < region_size;){
			trans_Read_item &ctr = region_trans_list[bg_idx];
			int try_idx = bg_idx;
			for(;try_idx < region_size; try_idx++){
				if(ctr.mtid != region_trans_list[try_idx].mtid || ctr.mpos + 20000 <  region_trans_list[try_idx].mpos)
					break;
			}
			int cluster_read_number = try_idx - bg_idx;
			fprintf(stderr, "[New region: cluster_read_number:%d]\n", cluster_read_number);
			int UM_read_id = 0;
			if(cluster_read_number >= 1){
				R_region region;
				region.chr_ID = ctr.mtid;
				region.st_pos = ctr.mpos + 1-100;
				region.ed_pos = region_trans_list[try_idx - 1].mpos + 1+100;
				int load_ref_length = 0;
				bool low_complex = false;
				if(true){
					char * c_reference = ref->load_ref_by_region(region.chr_ID, region.st_pos, region.ed_pos + 100, &load_ref_length);
					if(false)fprintf(stderr, "%s\n", c_reference);
					std::map<char, int> char_count;
					for (int i = 0; i < load_ref_length; i++)
						char_count[c_reference[i]]++;
					for(std::map<char, int>::value_type & c: char_count){
						if(c.second > (load_ref_length * 0.8 + 1))
						low_complex = true;
					}
					free(c_reference);
				}else{
					fprintf(stderr, "[SKIP low_complex check]\n");
				}

				int region_check_read = 0;
				if(!low_complex){
					fprintf(stderr, "\t[NOT low_complex ]\t");
					resetRegion_ID(bam, &region);	//reset region
					check_region_size++;
					//reference check:
					while (bam_next(bam)) {
						check_read_num++;
						region_check_read++;
						if(region_check_read > 3000)//skip region with too many reads
							break;
						bam1_t *br = &(bam->_brec);
						if(bam_is_secondary(br))		continue;
						if(bam_is_supplementary(br))	continue;
						if(bam_is_duplicate(br))       	continue;
						if(br->core.mtid == tid && br->core.mpos > st_pos && br->core.mpos < ed_pos){
							//load reads
							int read_len = br->core.l_qseq;
							get_bam_seq(0, read_len, storeReadBuff, br);
							//if dir == mdir: reverse the string:
							bool b_forward = bam_is_fwd_strand(br);
							bool m_forward = bam_is_mate_fwd_strand(br);
							bool store_in_reverse = (b_forward == m_forward);
							int clip_left, clip_right;
							bam_aligned_clip_analysis(br, &clip_left, &clip_right);
							int32_t r_pos = br->core.pos - clip_left;
							ab.add_read_TL(storeReadBuff, bam_get_qual(br), read_len, store_in_reverse, UM_read_id++, br->core.tid, r_pos);
							if(total_load_read++ > MAX_load_read_number){
								fprintf(stderr, "[TL_LOADING: total load read is over max load %d, skip others and return]\n", total_load_read);
								return 0;
							}
							if(false){
								//debug Code:
								char * bam_name = (char *)bam_qname(br);//string values(;, ref)
								fprintf(stderr, "[TL_LOADING(M2): :%s %d:%d @QUAL %d ; mate: %d:%d dir:(%c, %c) %s] ", bam_name, br->core.tid, br->core.pos, br->core.qual, br->core.mtid, br->core.mpos, b_forward?'F':'R', m_forward?'F':'R', storeReadBuff);
								uint32_t* bam_cigar = bam_get_cigar(br);
								uint32_t n_cigar = br->core.n_cigar;
								for (uint i = 0; i < n_cigar; ++i)
								{
									int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
									int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
									fprintf(stderr, "%d%c", c_size, "NMIDKSHPMX"[c_type]);
								}
								fprintf(stderr, "\n");
							}
						}
					}
					fprintf(stderr, "check_read_num: %d total_load_read %d\n", check_read_num, total_load_read);
				}else{
					fprintf(stderr, "check_read_num: %d total_load_read %d\n", check_read_num, total_load_read);
				}
			}
			bg_idx += cluster_read_number;
		}
		fprintf(stderr, "Region low_complex\n");
		return total_load_read;
	}

	void init(Bam_file *bam_, Bam_file *TL_bam_){
		storeReadBuff = (char *)xcalloc(2000, 1);
		bam = bam_;
		TL_bam = TL_bam_;
	}

	void clear(){
		last_read_end_index = 0;
		trans_list.clear();
	}
	int cur_bg_idx;
	int cur_ed_idx;
	int last_read_end_index;
	std::vector<trans_Read_item> trans_list;
	std::vector<trans_Read_item> region_trans_list;
	//pointer
	Bam_file *bam;
	Bam_file *TL_bam;

	char * storeReadBuff;
};

struct BAM_handler{
public:
	typedef std::vector<READ_record> READ_LIST;
	void init(char * bamFileName, char *realignment_bamFileName, char *TL_read_Filename, char * ref_file_name){
		const int MAX_READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
		sr_read_index = new int [MAX_READ_INDEX_SIZE];
		um_read_index = new int [MAX_READ_INDEX_SIZE];
		storeReadBuff = new uint8_t[10000000];//10MBP at MAX
		kv_init(base);
		kv_init(quality);
		kv_init(cigar);
		bam_file_open(bamFileName, ref_file_name, NULL, &file);

		if(realignment_bamFileName != NULL)
			bam_file_open(realignment_bamFileName, ref_file_name, NULL, &re_alignment_file);
		else{
			memset(&re_alignment_file, 0, sizeof(Bam_file));
		}

		if(TL_read_Filename != NULL){
			bam_file_open(TL_read_Filename, ref_file_name, NULL, &TL_file);
			tr_loader.init(&file, &TL_file);
		}else{
			memset(&TL_file, 0, sizeof(Bam_file));
			tr_loader.init(&file, NULL);
		}

	}

	void destroy(){
		delete[]sr_read_index;
		delete[]um_read_index;
		delete[]storeReadBuff;
		bam_file_close(&file);
		if(re_alignment_file._hdr != NULL)
			bam_file_close(&re_alignment_file);
	}
	bool clip_low_quality_Filter(bam1_t *br, int read_len, int soft_left, int soft_right);
	bool storeReadSR(bam1_t *br, int soft_left, int soft_right, int gap_mismatch_inside);// store normal SR signals
	void storeReadUM(bam1_t *br, uint8_t *query);
	void storeCCS_ReadCore(int readLen, uint8_t *seq, uint8_t * qual, int cigarLen, uint32_t* bam_cigar, int32_t pos);
	void storeReadCore(	Read_type::T t, int readLen, uint8_t *seq, uint8_t * qual,	int cigarLen,
			uint32_t* bam_cigar, uint16_t flag, int32_t pos,	int soft_left, int soft_right, int gap_mismatch_inside);

	void getReadStr(uint8_t * to, READ_record& c_r, int bg, int len){
		memcpy(to, base.a + c_r.read_index + bg, len);
	}
	inline uint8_t* getReadStr(READ_record& c_r){ return base.a + c_r.read_index;}
	inline uint8_t* getQualStr(READ_record& c_r){ return quality.a + c_r.quality_index;}
	inline path_segment* getCIGARStr(READ_record& c_r){ return cigar.a + c_r.cigar_index;}

	void load_read_CCS(std::string &s_store, READ_record& c_r,
			 int ref_pos_st, int ref_pos_ed, bool using_clip){
		s_store.clear();
		if(c_r.position > ref_pos_ed || c_r.position + c_r.read_l + 10000 < ref_pos_st){
			return;
		}
		int MIN_QUAL_TGS = 0;
		//store ref position to read position:
		//
		path_segment* p = getCIGARStr(c_r);
		uint8_t * read_str = getReadStr(c_r);
		uint8_t * qual_str = getQualStr(c_r);

		uint32_t n_cigar = c_r.cigar_l;
		int seq_i = 0;
		int ref_i = c_r.position;
		for (uint i = 0; i < n_cigar; ++i)
		{
			int c_type = p[i].type;
			int c_size = p[i].length;
			switch (c_type){
			case CIGAR_MATCH: case CIGAR_SEQ_MATCH: case CIGAR_SEQ_MISMATCH:
				for(int i = 0; i < c_size; i++, seq_i++, ref_i++){
					if(ref_i >= ref_pos_st && ref_i < ref_pos_ed){
						uint8_t c = "ACGTNN"[read_str[seq_i]];
						if(qual_str[seq_i] < MIN_QUAL_TGS)
							c = 'N';
						if(c == 'N'){
							fprintf(stderr, " ");
						}
						s_store.push_back(c);
					}
				}
				break;
			case CIGAR_INSERT:
				for(int i = 0; i < c_size; i++, seq_i++){
					if(ref_i >= ref_pos_st && ref_i < ref_pos_ed){
						uint8_t c = "ACGTNN"[read_str[seq_i]];
						if(read_str[seq_i] < MIN_QUAL_TGS)
							c = 'N';
						s_store.push_back(c);
					}
				}
				break; //do nothing
			case CIGAR_DELETE:
				for(int i = 0; i < c_size; i++, ref_i++);
				break;
			case CIGAR_SOFT_CLIP:
				for(int i = 0; i < c_size; i++, seq_i++)
				{
					if(using_clip){
						uint8_t c = "ACGTNN"[read_str[seq_i]];
						if(read_str[seq_i] < MIN_QUAL_TGS)
							c = 'N';
						s_store.push_back(c);
					}
				}
				break;
			case CIGAR_HARD_CLIP:
				break;
			default:	break;
			}
		}
	}
	bool isSameHeader(const BAM_handler &B) const;
	static bool pass_compact_filter(uint8_t *s, const int len);
	//bool clip_AAA_TailFilter(uint8_t *query, const int read_len);//when a clip string clip to an "AAAAAAAA..." tail, return true
	bool clip_AAA_TailFilter(uint8_t *read_bin, const int read_len, int soft_left, int soft_right);

	void clear(){
		read_list.clear();
		um_read_list.clear();
		base.n = 0;
		quality.n = 0;
		cigar.n = 0;
		index_already_built = false;
	}

	void build_read_index_SINGLE(READ_LIST &l, int *index, int region_st_index){
		for(int i = 0; i < READ_INDEX_SIZE; i++)
			index[i] = -1;
		int read_ID = 0;
		for(auto r: l){
			int c_index = (r.position - region_st_index)/ 100;
			if(c_index >= 0 && c_index < READ_INDEX_SIZE && index[c_index] == -1)
				index[c_index] = read_ID;
			read_ID++;
		}
	}

	inline void build_read_index(int region_st_index){
		if(index_already_built) return;
		build_read_index_SINGLE(read_list, sr_read_index, region_st_index);
		build_read_index_SINGLE(um_read_list, um_read_index, region_st_index);
		index_already_built = true;
	}

	int search_reads(Read_type::T t, Ass_Block & ab, int st_pos, int ed_pos, int max_load, int region_st_index){
		//init:
		if(t != Read_type::SR && t != Read_type::UM) return 0;
		if(index_already_built == false) build_read_index(region_st_index);
		if(st_pos > ed_pos) return 0;


		//select read and index list
		READ_LIST *c_read_l = NULL; int *index = NULL;
		if(t == Read_type::SR){ c_read_l = &read_list; index = sr_read_index; }
		else				  {	c_read_l = &um_read_list; index = um_read_index; }
		READ_LIST &l = *c_read_l;

		//loading read signals
		int c_index = (st_pos - region_st_index) / 100;

		if(c_index < 0) c_index = 0;
		if(c_index >= READ_INDEX_SIZE) c_index = READ_INDEX_SIZE - 1;

		while(c_index < READ_INDEX_SIZE && index[c_index] == -1)
			c_index++;
		if(c_index == READ_INDEX_SIZE) return 0;//no data, skip: C1
		int search_st = index[c_index];


		int list_size = l.size();
		if(list_size == 0) return 0;//no data, skip: C2

		if(search_st >= list_size) return 0;//no data, skip: C3
		//search the beginning
		//fprintf(stderr, "c_index %d search_st %d st_pos %d\n", c_index, search_st, st_pos);

		for(;search_st < list_size; search_st ++){
			//fprintf(stderr, "XX l[search_st].position %d\n", l[search_st].position);
			if(l[search_st].position >= st_pos)
				break;
		}

		//search for the end
		int total_load_num = 0;
		for(auto r = l.begin() + search_st; r->position <= ed_pos && search_st < list_size; r++, search_st ++){
			ab.add_read_bin(t, getReadStr(*r), getQualStr(*r), r->read_l, r->position - r->soft_left,  r->soft_left,  r->soft_right ,r->NM_NUM,r - l.begin());
			if(++total_load_num >= max_load) break;
		}


		return total_load_num;
	}

	READ_LIST read_list;
	READ_LIST um_read_list;

	TransReadLoader tr_loader;
	//std::vector<trans_Read_item> trans_read_list;//used to store trans-location read pairs:

	Bam_file file;
	Bam_file TL_file;
	Bam_file re_alignment_file;

	uint8_t * storeReadBuff;
private:
	BAM_handler(const BAM_handler &b);//can`t be copy

	//read lists and index
	bool index_already_built = false;
	const static int READ_INDEX_SIZE = SEGMENT_LEN/100 + 1;
	int * sr_read_index;
	int * um_read_index;

	//buffs
	kvec_T(uint8_t, BASE_STR);//store read string
	kvec_T(uint8_t, QUALITY_STR);

public:
	BASE_STR base;
	QUALITY_STR quality;
	path_t cigar;
};

#endif /* BAMHANDLER_HPP_ */
