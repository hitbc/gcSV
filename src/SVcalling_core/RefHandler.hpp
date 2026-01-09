/*
 * RefHandler.hpp
 *
 *  Created on: 2020-4-28
 *      Author: fenghe
 */

#ifndef SIGNAL_REFHANDLER_HPP_
#define SIGNAL_REFHANDLER_HPP_

extern "C"
{
	#include "../clib/utils.h"
	//#include "clib/bam_file.h"
	extern uint64_t kmerMask[33];
}
extern unsigned char charToDna5n[256];

#include <algorithm>
#include<array>
#include<unordered_map>

#include "../cpp_lib/cpp_utils.hpp"
#include "../cpp_lib/RefRegion.hpp"

#define REF_BLOCK_LEN 4000 // at most 8 k data
#define REF_BLOCK_NUM 500 // 4K * 500 = 2M
#define HASH_LEN 14 // 7mer
#define KMER_LEN 7 // 7mer = 14/2
#define SEGMENT_LEN 5000000 //5M

//------------------------index for reference----------------------------//
struct REF_HASH{
	REF_HASH(uint16_t key_, uint16_t start_pos_):key(key_), start_pos(start_pos_){}
	uint16_t key;//use 7 mer as key
	uint16_t next = 0;//store 16K data at most
	uint16_t start_pos;//position from "ref_st", store 16K data at most
};

struct REF_index{
		hash_T(REF_HASH, REF_HASH_T);

		void init(){
			hash_t_init(REF_HASH, &hash, HASH_LEN);
		}

		void destroy(){
			 hash_t_destroy(&hash);
		}

		REF_HASH *getHashByKey(uint16_t kmer) const {
			REF_HASH *rst;
			hash_t_search(&hash, kmer, rst);
			return rst;
		}
		REF_HASH *getNextHash(REF_HASH * old_hash)const {
			REF_HASH *rst;
			hash_t_search_next(&hash, old_hash, rst);//search in ref index:
			return rst;
		}

		int32_t ref_st = 0; //local in a segmentation
		int32_t ref_ed = 0; //local in a segmentation
		REF_HASH_T hash;

private:
		REF_index(const REF_index & i);//can`t be copy
	};

struct RefLocalRepeatItem{
	int repeat_level;
};

struct RefLocalRepeatItemR{
	int chrID;
	int position;
	int repeat_level;
	RefLocalRepeatItemR(int chrID, int position, int repeat_level){
		this->chrID = chrID;
		this->position = position;
		this->repeat_level = repeat_level;
	}
	void show(){
		fprintf(stderr, "chrID %d , position %d , repeat_level  %d \n", chrID, position, repeat_level);
	}
};

struct RefLocalRepeat{
	std::unordered_map<uint64_t, RefLocalRepeatItem> *repeatInfoMap;
	int region_length = 0;
	int region_step_size = 0;

	void get_middle(RefLocalRepeatItemR& r, int & middle_pos, int &range){
		middle_pos = r.position + region_length/2;
		range = region_length / 2;
	}

	void load(char * ref_i_fn){
		if(repeatInfoMap == NULL)
			repeatInfoMap = new std::unordered_map<uint64_t, RefLocalRepeatItem>();
		std::vector<std::string> tmp;
		std::vector<std::string> item_value;
		load_string_list_from_file(ref_i_fn, tmp);
		bool isHeader = true;
		for(std::string &s :tmp){
			split_string(item_value, s.c_str(), "\t");
			xassert(item_value.size() == 3, " ");
			if(isHeader){
				region_length = atoi(item_value[0].c_str());
				region_step_size = atoi(item_value[1].c_str());
				isHeader = false;
			}
			else{
				int chrID = atoi(item_value[0].c_str());
				int pos = atoi(item_value[1].c_str());
				uint64_t pos_global = pos + ((uint64_t)chrID << 32);
				int repeat_level = atoi(item_value[2].c_str());
				RefLocalRepeatItem n = {repeat_level};
				if(repeat_level >= 40){
					repeatInfoMap->insert({pos_global, n});
				}
			}
		}
	}

	void search(int chrID, int region_bg, int region_ed, std::vector<RefLocalRepeatItemR> &r){
		r.clear();
		int bgn_pos = region_bg  - region_bg % region_step_size;
		int end_pos = region_ed  - region_ed % region_step_size + region_step_size;
		for(int i = bgn_pos; i <= end_pos; i++){
			uint64_t pos_global = i + ((uint64_t)chrID << 32);
			std::unordered_map<uint64_t, RefLocalRepeatItem>::iterator it = repeatInfoMap->find(pos_global);
			if(it !=repeatInfoMap->end()){
				r.emplace_back(chrID, i, it->second.repeat_level);
			}
		}
	}
};

struct RefIndexItem{
	int position;
	int count;

	void add_count(){
		count ++;
	}

	void set_new(int position){
		count = 1;
		this->position = position;
	}
};
struct RefHandler{

private:
	//basic parameters
	int st_chr_ID; int st_pos;
	int ed_chr_ID; int ed_pos;
	char*	referenceFilename;

	bam_hdr_t* bam_header;
	faidx_t * c_ref_idx;

	//reference regions
	RefRegion curRefRegion;//current reference region
	RefRegion cur_r;//region for the index, it must be within the current reference region

	//reference string
	struct REF_STR{ size_t n; size_t m; uint8_t *a; } ref;	//a buff to store the current reference region

	std::unordered_map<int32_t, RefIndexItem> ref_kmer_index;


public:
	int REF_IDX_KMER_LEN;
	uint64_t REF_IDX_KMER_MASK;
	void init(
			int st_chr_ID_, int st_pos_,
			int ed_chr_ID_, int ed_pos_,
			char* referenceFilename_){
		st_chr_ID = st_chr_ID_;
		st_pos =  st_pos_;
		ed_chr_ID = ed_chr_ID_;
		ed_pos = ed_pos_;
		referenceFilename = referenceFilename_;

		cur_r.set(st_chr_ID, st_pos - SEGMENT_LEN, st_pos);
		curRefRegion.set(-1, -1, -1);
		c_ref_idx = NULL;
		if(c_ref_idx == NULL)//load index (build new index when without index file), open reference file
			c_ref_idx = reference_index_load(referenceFilename);
		ref.n = ref.m = 0;
		ref.a = NULL;
		REF_IDX_KMER_LEN = 16;
		extern uint64_t kmerMask[33];
		REF_IDX_KMER_MASK = kmerMask[REF_IDX_KMER_LEN];
	}

	void destory(){
		if(ref.a != NULL)free(ref.a);
		reference_index_free(c_ref_idx);//Free fai index for fna file, and close reference file
	}

	faidx_t * getIdx(){
		return c_ref_idx;
	}

	int get_chr_ID(){ return cur_r.chr_ID; }

	void registHeader(bam_hdr_t* bam_header_){bam_header = bam_header_;}

	char * load_ref_by_region(int tid, int st_pos, int ed_pos, int * max_load){
		char ref_region_char[1024];
		sprintf(ref_region_char, "%s:%d-%d", bam_header->target_name[tid], st_pos, ed_pos);
		return fai_fetch(c_ref_idx, ref_region_char, max_load);
	}

	bool load_seg_index(){
		//get new current handle region:
		if(cur_r.st_pos + SEGMENT_LEN > (int)bam_header->target_len[cur_r.chr_ID] - 1)
			cur_r.set(cur_r.chr_ID + 1, 0, SEGMENT_LEN);
		else
			cur_r.set(cur_r.chr_ID, cur_r.st_pos + SEGMENT_LEN, cur_r.st_pos + SEGMENT_LEN + SEGMENT_LEN - 1);
		if(cur_r.chr_ID < bam_header->n_targets)
			cur_r.ed_pos = MIN(cur_r.ed_pos, (int)bam_header->target_len[cur_r.chr_ID] - 1);
		if(cur_r.chr_ID == ed_chr_ID)
			cur_r.ed_pos = MIN(cur_r.ed_pos, ed_pos);

		//when reaching end of
		if(cur_r.chr_ID >= bam_header->n_targets || //reach end of all target
				cur_r.chr_ID > ed_chr_ID || //reach end of parameter
				(cur_r.chr_ID == ed_chr_ID && cur_r.st_pos > ed_pos))
			return false;

		if(cur_r.chr_ID != curRefRegion.chr_ID){
			if(!load_reference(cur_r.chr_ID))
				return false;
		}

		xassert(cur_r.Within(curRefRegion), "The region for building index must be within the current reference region!");
		return true;
	}

	uint8_t * getRefStr(int bg_pos){	return ref.a + bg_pos;}
	RefRegion* get_cur_region(){ return &cur_r;}
	char * get_refFileName(){return referenceFilename;}

	void build_ref_index(){
		//cur_r.show();
		uint8_t * buff_bin = getRefStr(cur_r.st_pos);
		std::unordered_map<int32_t, RefIndexItem>::iterator it;
		int kmer_number = cur_r.getLen() - REF_IDX_KMER_LEN + 1;
		uint64_t kmer = bit2_nextKmer_init(buff_bin, REF_IDX_KMER_LEN);
		uint64_t MASK = REF_IDX_KMER_MASK;
		ref_kmer_index.clear();
		for(int i = 0; i < kmer_number; i++){
			kmer = bit2_nextKmerMASK( buff_bin + i, kmer, REF_IDX_KMER_LEN);
			if(i % 3 == 0){
				it = ref_kmer_index.find(kmer);
				if(it!=ref_kmer_index.end())
					it->second.add_count();
				else
					ref_kmer_index[kmer].set_new(i + cur_r.st_pos);
			}
		}
	}

	//search the index and get the position
	int search_ref_idx(int32_t kmer_query){
		std::unordered_map<int32_t, RefIndexItem>::iterator it;
		it = ref_kmer_index.find(kmer_query);
		if(it!=ref_kmer_index.end() && (it->second.count == 1)){
			return it->second.position;
		}
		return -1;
	}

private:

	//load reference string and store data into "ref"
	bool load_reference(int chr_ID){

		if(chr_ID >= bam_header->n_targets) return false;

		int load_ref_length = 0;
		xassert(bam_header != NULL, "BAM header should be registered for a REF_BUFF before loading reference!\n");
		if(c_ref_idx == NULL)//load index (build new index when without index file), open reference file
			c_ref_idx = reference_index_load(referenceFilename);
		char * c_reference = get_reference_region(c_ref_idx, bam_header->target_name[chr_ID], &load_ref_length);
		if(ref.m < (size_t)load_ref_length)
			kv_resize(uint8_t, ref, (size_t)load_ref_length);
		ref.n = load_ref_length;
		//change char to bin
		uint8_t *char_ref = ref.a;
		for(int i = 0; i < load_ref_length; i++)
			char_ref[i] = charToDna5n[(int8_t)c_reference[i]];
		free(c_reference);

		curRefRegion.set(chr_ID, 0, load_ref_length - 1);

		return true;
	}

};
struct KMER_ERROR_PROFILE_HANDLER{
	//load:
	//
	float* kmer_error_profile = NULL;
	int kmer_error_profile_kmer_length = 0;
	int kmer_error_profile_total_kmer_number = 0;//change to struct.....

	uint64_t KMER_ALL_A_HOMO;
	uint64_t KMER_ALL_T_HOMO;
	int indexing_edge_len;
	void load_data(const char * KMER_ERROR_PROFILE_File){
		std::vector<std::string> load_line;
		std::vector<std::string> item_value;
		load_string_list_from_file(KMER_ERROR_PROFILE_File, load_line);
		size_t line_num = load_line.size();//skip the header line
		bool is_the_first_line = true;
		//skip the head line
		extern uint8_t charToDna5n[];
		for(uint64_t i = 1; i < line_num; i++){
			split_string(item_value, load_line[i].c_str(), "\t");//skip the header line
			std::string & kmer = item_value[0];
			float prob = atof(item_value[1].c_str());
			if(is_the_first_line){
				kmer_error_profile_kmer_length = kmer.size();
				kmer_error_profile_total_kmer_number = 0x1 << (kmer_error_profile_kmer_length*2);
				kmer_error_profile = (float*)xcalloc(kmer_error_profile_total_kmer_number, sizeof(float));
				is_the_first_line = false;
			}
			int kmer_bin = char2Kmer(kmer.c_str(), kmer_error_profile_kmer_length, charToDna5n);
			xassert(kmer_error_profile_total_kmer_number > kmer_bin, "");
			kmer_error_profile[kmer_bin] = prob;
		}
		KMER_ALL_A_HOMO = 0;
		KMER_ALL_T_HOMO = kmer_error_profile_total_kmer_number - 1;
		indexing_edge_len = 25;
	}

	//input: the contig and the CONTIG_IDX_KMER_LEN
	//output: show witch kmer need to be indexing
	//kmer_need_indexing: 'P' polishing; 'E' edge; 'N' do nothing
	void search_profile_and_indexing_region(std::vector<uint8_t> & bin_string, std::string & kmer_need_indexing, int CONTIG_IDX_KMER_LEN){
		//todo::
		//S1: get all kmers
		int kmer_number_indexing = bin_string.size() - CONTIG_IDX_KMER_LEN + 1;
		kmer_need_indexing.resize(kmer_number_indexing);
		//for each kmer
		extern uint64_t kmerMask[33];
		uint64_t CONTIG_IDX_KMER_MASK;
		CONTIG_IDX_KMER_MASK = kmerMask[kmer_error_profile_kmer_length];
		uint64_t MASK = CONTIG_IDX_KMER_MASK;
		//set the kmer need polishing
		//long AAA strings add, additional add kmer < 0.8
		uint8_t * buff_bin = &(bin_string[0]);
		uint64_t kmer = bit2_nextKmer_init(buff_bin, kmer_error_profile_kmer_length);
		for(int i = 0; i < kmer_number_indexing; i++){
			kmer = bit2_nextKmerMASK( buff_bin + i, kmer, kmer_error_profile_kmer_length);
			if(kmer == KMER_ALL_A_HOMO || KMER_ALL_T_HOMO == kmer || kmer_error_profile[kmer] < 0.5)
				kmer_need_indexing[i] = 'P';
			else
				kmer_need_indexing[i] = 'N';
		}
		//set the kmer need indexing
		for(int i = 0; i < kmer_number_indexing; i++){
			if(kmer_need_indexing[i] == 'P'){
				int min_reverse = MAX(0,i - indexing_edge_len); //reverse
				for(int j = i - 1; j >= min_reverse && kmer_need_indexing[j] == 'N'; j--)
					kmer_need_indexing[j] = 'E';
				int max_forward = MIN(kmer_number_indexing - 1, i + indexing_edge_len); 				//forward
				for(int j = i + 1; j <= max_forward && kmer_need_indexing[j] == 'N'; j++)
					kmer_need_indexing[j] = 'E';
			}
		}
	}
};



#endif /* SIGNAL_REFHANDLER_HPP_ */
