/*
 * Contig_polishing.hpp
 *
 *  Created on: 2025-3-29
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_CONTIG_POLISHING_HPP_
#define SVCALLING_CORE_CONTIG_POLISHING_HPP_

#include "../cpp_lib/cpp_utils.hpp"
//#include "../SVcalling_core/SveHandler.hpp"
extern "C"
{
#include "../clib/utils.h"
#include "../clib/vcf_lib.h"
}

struct ContigIndexItem{
	int position; int count;
	void add_count(){count ++;}
	void set_new(int position){ count = 1; this->position = position;}
};

struct READ_CONTIG_MATCH{
	int count = 0;
	int min_contig_position = MAX_int32t;
	int max_contig_position = 0;
	int min_read_position = MAX_int32t;
	int max_read_position = 0;

	void set(int contig_position, int NGS_read_position){
		count++;
		min_contig_position = MIN(min_contig_position, contig_position);
		max_contig_position = MAX(max_contig_position, contig_position);
		min_read_position = MIN(min_read_position, NGS_read_position);
		max_read_position = MAX(max_read_position, NGS_read_position);
	}

	void show(){
		fprintf(stderr, "Anchored KMER number(MAX MEM) %5d, min_contig_position %5d, max_contig_position %5d, min_read_position %5d, max_read_position %5d\t",
			count, min_contig_position, max_contig_position, min_read_position, max_read_position);
	}
};


struct NGS_PLOSHING_STRING{
	struct read_info{
		int read_ID;
		int read_position;
		void set(int read_ID, int read_position){
			this->read_ID = read_ID;
			this->read_position = read_position;
		}
		void show(){
			fprintf(stderr, "%5d@%5d:\t",
					read_ID, read_position);
		}
	};

	std::map<std::vector<uint8_t>, std::vector<read_info>> rl;
	std::map<std::vector<uint8_t>, std::vector<read_info>>::iterator con_it;
	int contig_positon = -1;
	std::vector<uint8_t> s;
	//std::string s;
	float total_read_number = -1;

	void add_signal(bool print_log, uint8_t * buff_bin, int read_ID, int read_position, int contig_positon, int kmer_number, int kmer_len){
		this->contig_positon = contig_positon;
		s.clear();
		s.resize(kmer_number + kmer_len);
		memcpy(&(s[0]), buff_bin + read_position, kmer_number + kmer_len);
		//for(int i = 0; i < kmer_number + kmer_len; i++)	s[i] = buff_bin[read_position + i];

		if(print_log && false){
			fprintf(stderr, "[add_signal:]read_ID %d, read_position %d, string :\n", read_ID, read_position);
			for(uint i = 0; i < s.size(); i++)
				fprintf(stderr, "%c", "ACGTNNNN"[s[i]]);
			fprintf(stderr, "\n");
		}

		rl[s].emplace_back();
		rl[s].back().set(read_ID, read_position);
	}
//
	void cal_read_number(){
		if(total_read_number >= 0)
			return;
		total_read_number = 0;
		std::map<std::vector<uint8_t>, std::vector<read_info>>::iterator it = rl.begin();
		for(;it != rl.end(); it++)
			total_read_number += it->second.size();
	}

	void show_one_consensus_string(std::map<std::vector<uint8_t>, std::vector<read_info>>::iterator it){
		std::vector<read_info> &ri_l = it->second;
		const std::vector<uint8_t> & s = it->first;
		fprintf(stderr, "update string: ");
		for(uint i = 0; i < s.size(); i++)
			fprintf(stderr, "%c", "ACGTNNNN"[s[i]]);
		fprintf(stderr, "\tread num %ld: ", ri_l.size());
		for(read_info & i: ri_l)
			i.show();
		fprintf(stderr, "\n");
	}

	void show(){
		fprintf(stderr, "contig_positon %d: \n", contig_positon);
		std::map<std::vector<uint8_t>, std::vector<read_info>>::iterator it = rl.begin();
		for(;it != rl.end(); it++){
			show_one_consensus_string(it);
		}
		fprintf(stderr, "\n");
	}

	bool has_consensus_string(){
		cal_read_number();
		std::map<std::vector<uint8_t>, std::vector<read_info>>::iterator it = rl.begin();
		for(;it != rl.end(); it++){
			std::vector<read_info> &ri_l = it->second;
			if(ri_l.size() >= 3 && ri_l.size() >= total_read_number *0.7){
				con_it = it;
				return true;
			}
		}
		return false;
	}
};

struct HOMO_COUNT{
	int read_pos;
	int read_id;
	int homo_len;
	int contig_position;
	char c;
	void set(int read_pos, int read_id, int homo_len, int contig_position, char c){
		this->read_pos = read_pos;
		this->read_id = read_id;
		this->homo_len = homo_len;
		this->contig_position = contig_position;
		this->c = c;
	}

	void show(){
		fprintf(stderr, "contig_position %d, len %d, \"%c\",read_id %d, read pos %d\n", contig_position, homo_len, c, read_id, read_pos);
	}

	//sort by the contig position
	static inline int cmp_by_contig_pos(const HOMO_COUNT &a, const HOMO_COUNT &b){
		if(a.contig_position != b.contig_position)
			return a.contig_position < b.contig_position;
		return a.homo_len < b.homo_len;
	}
};

struct HOMO_ANA_ITEM{
	//store the results
	int suggest_homo_size;
	int min_contig_position;
	int max_contig_position;
	int support_read_count;
	void set(int suggest_homo_size, int min_contig_position, int max_contig_position, int support_read_count){
		this->suggest_homo_size = suggest_homo_size;
		this->min_contig_position = min_contig_position;
		this->max_contig_position = max_contig_position;
		this->support_read_count = support_read_count;
	}
	void show(){
		fprintf(stderr, "Contig_position %d~%d, count %d \t", min_contig_position, max_contig_position, support_read_count);
		fprintf(stderr, "suggest_homo_size %d \n", suggest_homo_size);
	}
};

struct CONTIG_POLISHING_HANDLER{
	//indexing the bin-contig:
	int CONTIG_IDX_KMER_LEN;
	uint64_t CONTIG_IDX_KMER_MASK;

	std::string kmer_need_indexing;
	std::map<int32_t, ContigIndexItem> contig_kmer_index;

	std::map<int, READ_CONTIG_MATCH> read_contig_match_l;
	std::map<int, NGS_PLOSHING_STRING> ngs_polishing_string_l;

	std::vector<HOMO_COUNT> homo_count_list;

	uint8_t storeReadBuff[1024];

	KMER_ERROR_PROFILE_HANDLER *kmer_profile_handler;

	Bam_file *c_b = NULL;

	void init(KMER_ERROR_PROFILE_HANDLER *kmer_profile_handler, Bam_file *c_b){
		extern uint64_t kmerMask[33];
		CONTIG_IDX_KMER_LEN = 13;//4^13=2^26=64M
		CONTIG_IDX_KMER_MASK = kmerMask[CONTIG_IDX_KMER_LEN];
		this->kmer_profile_handler = kmer_profile_handler;
		this->c_b = c_b;
	}

	//kmer_need_indexing: 'P' polishing; 'E' edge; 'N' do nothing
	void contig_index_building(bool print_log, std::vector<uint8_t> & bin_contig){
		int kmer_number = bin_contig.size() - CONTIG_IDX_KMER_LEN + 1;
		uint64_t MASK = CONTIG_IDX_KMER_MASK;
		//long AAA strings add, additional add kmer < 0.8
		kmer_profile_handler->search_profile_and_indexing_region(bin_contig, kmer_need_indexing, CONTIG_IDX_KMER_LEN);

		if(print_log) fprintf(stderr, "kmer_need_indexing %s \n ", kmer_need_indexing.c_str());

		uint8_t * buff_bin = &(bin_contig[0]);
		std::map<int32_t, ContigIndexItem>::iterator it;
		uint64_t kmer = bit2_nextKmer_init(buff_bin, CONTIG_IDX_KMER_LEN);
		contig_kmer_index.clear();
		for(int i = 0; i < kmer_number; i++){
			kmer = bit2_nextKmerMASK( buff_bin + i, kmer, CONTIG_IDX_KMER_LEN);
			if(kmer_need_indexing[i] != 'N'){
				it = contig_kmer_index.find(kmer);
				if(it!=contig_kmer_index.end())
					it->second.add_count();
				else
					contig_kmer_index[kmer].set_new(i);
			}
		}
	}

	void read_contig_kmer_match(uint8_t * buff_bin, int read_len, uint64_t MASK){
		std::map<int32_t, ContigIndexItem>::iterator it_idx;
		uint64_t kmer = bit2_nextKmer_init(buff_bin, CONTIG_IDX_KMER_LEN);
		int read_kmer_number = read_len - CONTIG_IDX_KMER_LEN + 1;
		//the result list:
		read_contig_match_l.clear();
		for(int i = 0; i < read_kmer_number; i++){
			kmer = bit2_nextKmerMASK( buff_bin + i, kmer, CONTIG_IDX_KMER_LEN);
			it_idx = contig_kmer_index.find(kmer);
			//only UNIQ match will be used
			if(it_idx!=contig_kmer_index.end() && it_idx->second.count == 1){
				int position_diff = it_idx->second.position - i;
				read_contig_match_l[position_diff].set(it_idx->second.position, i);
			}
		}
	}

	READ_CONTIG_MATCH *get_max_count_match(bool print_log, int read_ID, int &anchor_mem_number, int & total_anchor_kmer_number, int &read_2_contig_align_position){
		//analysis the mapping results:
		std::map<int, READ_CONTIG_MATCH>::iterator read_match_it = read_contig_match_l.begin();
		anchor_mem_number = 0;

		READ_CONTIG_MATCH * max_count_match = NULL;
		for(;read_match_it != read_contig_match_l.end(); read_match_it++){
			total_anchor_kmer_number += read_match_it->second.count;
			if(read_match_it->second.count > 25){//at least anchor 25 kmer
				if(print_log && false){
					fprintf(stderr, "read_ID %d read-contig alignment position is %d\t", read_ID, read_match_it->first);
					read_match_it->second.show();
				}
				anchor_mem_number++;
				if(anchor_mem_number < read_match_it->second.count){
					anchor_mem_number = read_match_it->second.count;
					max_count_match = &(read_match_it->second);
					read_2_contig_align_position = read_match_it->first;
				}
			}
		}
		return max_count_match;
	}

	void collect_NGS_polishing_string(bool print_log, R_region &NGS_load_region){
		resetRegion_ID(c_b, &NGS_load_region);   //reset region
		int read_ID = -1;
		uint64_t MASK = CONTIG_IDX_KMER_MASK;
		ngs_polishing_string_l.clear();
		homo_count_list.clear();
		while (bam_next(c_b)) {
			bam1_t *br = &(c_b->_brec);
			if(bam_is_duplicate(br)) 	        		continue;
			if(bam_is_secondary(br))	                continue;
			if(bam_is_supplementary(br))                continue;
			//if(bam_is_unmapped(br))	                 	continue;
			if(br->core.qual < 20)	                	continue;
			//load NGS read string
			const int read_len = br->core.l_qseq;
			xassert(read_len < 1024, "");
			if(get_bam_seq_bin(0, read_len, storeReadBuff, br) == false)
				continue;//get string failed, return
			 //search the read in the index
			uint8_t * buff_bin = storeReadBuff;
			read_ID++;
			//read alignment
			read_contig_kmer_match(buff_bin, read_len, MASK);
			int anchor_mem_number = 0;
			int total_anchor_kmer_number = 0;
			int read_2_contig_align_position = 0;
			READ_CONTIG_MATCH * max_count_match = get_max_count_match(print_log, read_ID, anchor_mem_number, total_anchor_kmer_number, read_2_contig_align_position);
			//define the read alignment position and read alignment coverage
			//get the string of read in the polishing regions
			if(anchor_mem_number == 0)//do nothing when no alignment results
				continue;
			xassert(max_count_match != NULL, "");
			//special : when no-linear alignment, todo::
			if(anchor_mem_number == 1){//linear alignment

			}
			else{//no-linear alignment

			}
			if(print_log && false){
				fprintf(stderr, "Read is anchored to contig: \t total_anchor_kmer_number %d, ", total_anchor_kmer_number);
				max_count_match->show();
				fprintf(stderr, " the read string is: ");
				for(int i = 0; i < read_len; i++){
					fprintf(stderr, "%c", "ACGTNNNNN"[buff_bin[i]]);
				}
				fprintf(stderr, "\n");
			}

			//store the read alignment result to list:
			//store the polishing regions
			//kmer_need_indexing: 'P' polishing; 'E' edge; 'N' do nothing
			for(int i = max_count_match->min_contig_position; i < max_count_match->max_contig_position; i++){
				if(kmer_need_indexing[i] == 'P'){
					int j = 1;
					for(;j + i < max_count_match->max_contig_position && kmer_need_indexing[j + i] == 'P'; j++);
					int read_position = i - read_2_contig_align_position;
					//show the kmer in the read:
					//store the search result:
					ngs_polishing_string_l[i].add_signal(print_log, buff_bin, read_ID, read_position, i, j, CONTIG_IDX_KMER_LEN);
					i += j;
				}
			}

			//store the homopolymer-count
			int MAX_SAME = 1; //AAAACCCT --> ACT
			for(int i = 1; i < read_len; i++){
				if(buff_bin[i] == buff_bin[i - 1])
					MAX_SAME ++;
				else{
					if(MAX_SAME >= 8){
						//store the data:
						homo_count_list.emplace_back();
						homo_count_list.back().set(i, read_ID, MAX_SAME, read_2_contig_align_position + i - MAX_SAME, "ACGTNNN"[buff_bin[i - 1]]);
					}
					MAX_SAME = 1;
				}
			}

		}
		//end of each read
	}
	//end of the function:
	//SRS_CONSENSUS_p, CONTIG_p
	int get_change_base_number( const uint8_t * s1_p, const uint8_t * s2_p, int size_update){
		//same check:
		int change_base_n = 0;
		for(int i = 0; i < size_update; i++){
			if(s1_p[i] != s2_p[i]){
				change_base_n++;
				fprintf(stderr, "New change_base: %d %d %d \n", i, s1_p[i], s2_p[i]);
			}
		}
		return change_base_n;
	}

	void show(std::map<int, NGS_PLOSHING_STRING>::iterator it, const uint8_t * s1_p, const uint8_t * s2_p, int size_update, int change_base_n){
		fprintf(stderr, "Original contig string: \t");
		for(int i = 0; i < size_update; i++)
			fprintf(stderr, "%c", "ACGTNNNN"[s2_p[i]]);
		fprintf(stderr, " \t");
		fprintf(stderr, "change_base_n %d\t", change_base_n);
		it->second.show_one_consensus_string(it->second.con_it);
	}

	void update_wrong_base(const uint8_t * s1_p, uint8_t * s2_p, int size_update){
		//re-place the wrong strings, update the wrong bases
		for(int i = 0; i < size_update; i++)
			s2_p[i] = s1_p[i];
	}

	void HOMO_signal_classification(std::vector<HOMO_ANA_ITEM> &homo_ana_list){
		//classification
		int min_contig_position = -100;
		int read_count = 1;
		std::map<int,int> len_count; len_count.clear();
		homo_ana_list.clear();
		for(uint i = 0; i < homo_count_list.size() + 1; i++){
			if(i == homo_count_list.size() || min_contig_position + 10 < homo_count_list[i].contig_position){
				if(min_contig_position >= 0){
					xassert(i >= 1, "");
					int suggest_homo_size = -1;
					std::map<int,int>::iterator it = len_count.begin();
					for(;it != len_count.end(); it++){
						float NGS_support_rate = (float)it->second / read_count;
						if(false) fprintf(stderr, "len %d, count %d(%f)\t", it->first, it->second, NGS_support_rate);
						if(NGS_support_rate > 0.7 && it->second >= 3){
							suggest_homo_size = it->first;
						}
					}
					homo_ana_list.emplace_back();
					homo_ana_list.back().set(suggest_homo_size, min_contig_position, homo_count_list[i - 1].contig_position, read_count);
				}
				//reset the counter
				if(i < homo_count_list.size()){
					min_contig_position = homo_count_list[i].contig_position;
					read_count = 1;
					len_count.clear();
					len_count[homo_count_list[i].homo_len] = 1;
				}
			}
			else{
				read_count++;
				xassert(i >= 0 && i < homo_count_list.size(), "");
				if(len_count.find(homo_count_list[i].homo_len) == len_count.end()){
					len_count[homo_count_list[i].homo_len] = 1;
				}else{
					len_count[homo_count_list[i].homo_len] ++;
				}
			}
		}
	}

	void HOMO_ERROR_update(bool print_log, std::vector<HOMO_ANA_ITEM> &homo_ana_list, std::vector<uint8_t> & bin_contig){
		int add_offset = 0;//todo::
		for(HOMO_ANA_ITEM & ha: homo_ana_list){
			ha.show();
			if(ha.suggest_homo_size > 0){
				int middle_pos_ref = (ha.min_contig_position + ha.max_contig_position) /2;
				//search forward
				uint8_t c = bin_contig[middle_pos_ref];
				int contig_size = bin_contig.size();
				int contig_home_st_pos = middle_pos_ref;
				int home_len_in_contig = 0;
				for(int i = middle_pos_ref; i < middle_pos_ref + ha.suggest_homo_size && i < contig_size; i++){
					if(c == bin_contig[i]){
						home_len_in_contig++;
						contig_home_st_pos--;
					}else
						break;
				}

				for(int i = middle_pos_ref - 1 ; i > middle_pos_ref - ha.suggest_homo_size && i >= 0; i--){
					if(c == bin_contig[i]){
						home_len_in_contig++;
					}else
						break;
				}
				fprintf(stderr, "[HOMO:]Original contig string: \t");
				for(int i = -10; i < ha.suggest_homo_size + 10; i++)
					fprintf(stderr, "%c", "ACGTNNNN"[bin_contig[i + ha.min_contig_position]]);
				fprintf(stderr, " home_len_in_contig %d\n", home_len_in_contig);
				//todo:: homo-polymer update the original string:
				{
					int change_base_number = ha.suggest_homo_size - home_len_in_contig;
					add_offset += change_base_number;
					//remove the original
					//remove the string:
					if(change_base_number == 0){
						/*DO NOTHING*/
					}
					if(change_base_number > 0){
						std::vector<uint8_t> add_s;
						for(int i = 0; i < change_base_number; i++)	add_s.emplace_back(c);
						bin_contig.insert(bin_contig.begin() + middle_pos_ref, add_s.begin(), add_s.end());
					}else if(change_base_number < 0){
						xassert(contig_home_st_pos >= 0 && contig_home_st_pos + change_base_number < (int)bin_contig.size(), "");
						bin_contig.erase(bin_contig.begin() + contig_home_st_pos, bin_contig.begin() + contig_home_st_pos + change_base_number);
					}
				}
			}
		}
	}

	void contig_correction(std::vector<uint8_t> & bin_contig, bool print_log){
		//part 1: read pile up:
		if(print_log) fprintf(stderr, "NGS correction (M1): read pile up:\n");
		for(std::map<int, NGS_PLOSHING_STRING>::iterator it = ngs_polishing_string_l.begin(); it != ngs_polishing_string_l.end(); it++){
			//show all the results:
			if(print_log && false) it->second.show();
			//analysis the results
			fprintf(stderr, "contig_positon %d \t", it->second.contig_positon);
			if(false == it->second.has_consensus_string())
				fprintf(stderr, "NO consensus_string \n");
			else {
				const uint8_t * s1_p = &(it->second.con_it->first[0]);//SRS_CONSENSUS_p
				uint8_t * s2_p = &(bin_contig[it->second.contig_positon]);//CONTIG_p
				int len = it->second.s.size();//region_size
				int change_base_n = get_change_base_number(s1_p, s2_p, len);
				show(it, s1_p, s2_p, len, change_base_n);
				if(change_base_n > 0 && change_base_n <= 2)//update SNPs
					update_wrong_base(s1_p, s2_p, len);
			}
		}

		//part 2: homo-polymer  correction
		if(print_log) fprintf(stderr, "NGS correction (M2): homo-polymer correction:\n");
		std::sort(homo_count_list.begin(), homo_count_list.end(), HOMO_COUNT::cmp_by_contig_pos);
		if(print_log && false) for(uint i = 0; i < homo_count_list.size(); i++) homo_count_list[i].show();

		std::vector<HOMO_ANA_ITEM> homo_ana_list;
		HOMO_signal_classification(homo_ana_list);
		HOMO_ERROR_update(print_log, homo_ana_list, bin_contig);
	}

	void run(std::vector<uint8_t> & bin_contig, R_region &NGS_load_region){
		bool print_log = true;
		//region:
		//get the profile
		contig_index_building(print_log, bin_contig);
		collect_NGS_polishing_string(print_log, NGS_load_region);
		contig_correction(bin_contig, print_log);
	}
};

#endif /* SVCALLING_CORE_CONTIG_POLISHING_HPP_ */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
