/*
 * RefHandler.hpp
 *
 *  Created on: 2020年4月28日
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

#include<array>
#include<unordered_map>

#include "../cpp_lib/cpp_utils.hpp"
#include "../cpp_lib/RefRegion.hpp"

#define REF_BLOCK_LEN 4000 // at most 8 k data
#define REF_BLOCK_NUM 500 // 4K * 500 = 2M
#define HASH_LEN 14 // 7mer
#define KMER_LEN 7 // 7mer = 14/2
#define SEGMENT_LEN 2000000 //2M

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

public:

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

#endif /* SIGNAL_REFHANDLER_HPP_ */
