/*
 * assembler.hpp
 *
 *  Created on: 2021-8-9
 *      Author: fenghe
 */

#ifndef CPP_LIB_ASSEMBLER_ASSEMBLER_HPP_
#define CPP_LIB_ASSEMBLER_ASSEMBLER_HPP_

#include <iosfwd>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include "../../clib/utils.h"

#include "map_define.hpp"

//read action record where the read was add to the contig
struct AssemblyReadAction{
	AssemblyReadAction(	int position_, int read_ID_, bool isAdd_){
		setBasic(position_, read_ID_, isAdd_);
	}
	void setBasic(int position_, int read_ID_, bool isAdd_){
		position_in_contig = position_;
		read_ID = read_ID_;
		 isAdd = isAdd_;
	}
	int position_in_contig;
	int read_ID;
	bool isAdd;//true for adding read action; false for delete read action

	//used in contig alignment
	int position_read;
	int suggest_contig_offset_in_ref;
	int search_len;//debugcode
	int wrong_base;

	void set_read_pos(std::string &read_seq, std::string &contig_seq, int ass_begin_offset_in_contig, int word_len, std::set<int> &remove_read_set, int read_in_ref_offset){

		if(isAdd == false)
			{ remove_read_set.emplace(read_ID); return;}
		else
			remove_read_set.erase(read_ID);

		position_read = -1;
		int read_kmer_num = read_seq.size() - word_len + 1;
		search_len = 0;
		if(position_in_contig < 0){
			for(int i = read_kmer_num - 1; i >= 0; i--){
				int cmp = read_seq.compare(i, word_len, contig_seq, position_in_contig - ass_begin_offset_in_contig, word_len);
				if(cmp == 0){
					position_read = i;
					break;
				}
				search_len++;
			}
		}else{
			for(int i = 0; i < read_kmer_num; i++){
				int cmp = read_seq.compare(i, word_len, contig_seq, position_in_contig - ass_begin_offset_in_contig, word_len);
				if(cmp == 0){
					position_read = i;
					break;
				}
				search_len++;
			}
		}

		suggest_contig_offset_in_ref = read_in_ref_offset - (position_in_contig - ass_begin_offset_in_contig - position_read);

		xassert(position_read != -1, "");
		//position_ref -=	ref_begin_pos;
	}
	void print(FILE *output, bool is_from_main_SV) const {
		fprintf(output,
			"[CP: %d RID:%d %c from_MAIN: %s ",
			position_in_contig, read_ID, isAdd?'+':'-', is_from_main_SV?"YES":"NO");

		if(isAdd)
			fprintf(output,
				"RP: %d SL: %d SUG: %d WB %d ",
				position_read, search_len, suggest_contig_offset_in_ref, wrong_base);
		fprintf(output,	"]\n");
	}
	static inline int cmp_by_position(const AssemblyReadAction &a, const AssemblyReadAction &b){
		int contig_pos_A = a.position_in_contig - a.position_read;
		int contig_pos_B = b.position_in_contig - b.position_read;
		return contig_pos_A < contig_pos_B;
	}

};

//Suggest
struct SUGGEST_POS_LIST_ITEM{
	SUGGEST_POS_LIST_ITEM(int suggest_pos_){
		suggest_pos = suggest_pos_;
		read_count = 0;
		low_wrong_base_read_number = 0;//
		high_wrong_base_read_number = 0;
		read_start_pos_sum = 0;
		ave_read_start = 0;
	}

	void add_read_start_pos(int read_start_pos, int wrong_base){
		if(wrong_base < 3)
			low_wrong_base_read_number ++;
		else if(wrong_base < 200)
			high_wrong_base_read_number ++;
		else return;

		read_count++;
		read_start_pos_sum += read_start_pos;
	}

	void count_ave_read_start(){ ave_read_start = (float)read_start_pos_sum/((float)(read_count) + 0.1); }

	void printf(FILE * log){
		if(ave_read_start < -5000 || ave_read_start > 5000)
			fprintf(log, "Fatal ERROR, ave_read_start wrong");
		fprintf(log, "[SG: %d, read count: %d, ave_read_start: %f]\t", suggest_pos, read_count, ave_read_start);
	}

	static inline int cmp_by_read_count(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.read_count > b.read_count; }
	static inline int cmp_by_read_position(const SUGGEST_POS_LIST_ITEM &a, const SUGGEST_POS_LIST_ITEM &b){ return a.ave_read_start < b.ave_read_start; }

	int suggest_pos;
	int read_count;
	int low_wrong_base_read_number;
	int high_wrong_base_read_number;
	int read_start_pos_sum;
	float ave_read_start;

};

struct AssemblyContig {
	bool contig_is_filtered = false;
	std::string seq;  ///< contigsequence
	unsigned seedReadCount = 0;  ///< no of reads containing the seeding kmer
	std::set<unsigned> supportReads;
	std::set<unsigned> rejectReads;//read were removed in the building of contig
	std::set<int> remove_read_set;//read were removed in the final read set
	pos_t conservativeRange_bgn = 0; ///< subsection of the contig with conservative coverage
	pos_t conservativeRange_end = 0;
	int ending_reason[2]; //index: 0 for left and 1 for right; data : 0 for reason 0[maxBaseCount < o.minCoverage] and 1 for reason 1[ after seeing one repeat word]
	int new_support_read = 0;
	unsigned kmerLength;
	std::vector<AssemblyReadAction> actions;
	std::vector<SUGGEST_POS_LIST_ITEM> SUGGEST_pos_list;
	int ass_begin_offset_in_contig;// the offset of begin of contig compare with the start assembly kmer

	void debug_print(FILE *output) const {

		fprintf(output,
				"CONTIG size: [%ld] seedCount: [%d] supportReads: [%ld] "
				"ending_reason: [%d %d]"
				"seq:\n",
				seq.size(), seedReadCount, supportReads.size(), ending_reason[0], ending_reason[1]);

		const char *seq_ = seq.c_str();
		static const unsigned rowSize(200);
		static const unsigned sectionSize(100);
		assert(nullptr != seq_);
		const unsigned seqLen(strlen(seq_));
		for (unsigned i(0); i < seqLen; ++i) {
			if (i) {
				if (0 == (i % rowSize))
					fprintf(output, "\n");
				else if (0 == (i % sectionSize))
					fprintf(output, " ");
			}
			fprintf(output, "%c", seq_[i]);
		}

		fprintf(output, "\nsupportReads: ");
		for (auto &r : supportReads)
			fprintf(output, "%d ", r);
		fprintf(output, "\n");
		if(!rejectReads.empty())
		{
			fprintf(output, "\nrejectReads: ");
			for (auto &r : rejectReads)
				fprintf(output, "%d ", r);
			fprintf(output, "\n");
		}
	}
};

/**************************************************************/
/// Information added to each read in the process of assembly
struct AssemblyReadInformation {
	AssemblyReadInformation(bool isPseudo_ = false) : isPseudo(isPseudo_) {}
	//bool isUsed = false;
	/// If true, the read was an assembled contig
	bool isPseudo = false;
	int base_depth = -1;
	/// Index of the contigs that this read is used in
	//std::vector<unsigned> contigIds;
};

/// Input parameters for IterativeAssembler
///
struct AssemblerOptions {
	AssemblerOptions() {}
	/// the symbol set used during assembly
	std::string alphabet = "ACGT";
	/// initial word (kmer) length
	unsigned minWordLength = 26;
	unsigned maxWordLength = 121;//todo::must be less than read length
	unsigned word_length_list[128] = { 26, 31, 36, 41, 51, 61, 71, 81, 101, 121, 141, 999};

	unsigned minCoverage = 2;/// min. coverage required for contig extension
	unsigned minConservativeCoverage = 2;	/// coverage required for conservative contig sub-range
	unsigned minUnusedReads = 3;/// min. number of unused reads to enable search for more contigs
	unsigned minSupportReads = 3;	/// min. number of reads required to start assembly
	unsigned maxAssemblyCount = 5;	/// Max. number of assembly returned for a given set of reads
	//bool reject_read_reused = false;/// reuse rejected reads
};

struct MainAssemblyHandler {
	///-----------------------------input--------------------------------------------//
	std::vector<std::string> reads;
	///-------------------------------output------------------------------------------//
	std::vector<AssemblyContig> contigs;
	bool reachMaxKmerLength;

	void clear() {
		reads.clear();
		readInfo.clear();
	}
	void assembley();
	void setRepeatMode(){
		//o.reject_read_reused = true;
		o.maxAssemblyCount = 5;
	}
	void setNormalMode(){
		//o.reject_read_reused = false;
		o.maxAssemblyCount = 10;
	}

	void run_ass_test(){
		//std::swap(am->reads, reads);
		reads.clear();

reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTG");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCG");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTG");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTG");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGAT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAA");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTT");
reads.emplace_back("ATGGGACTGGCTGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTT");
reads.emplace_back("TGACAGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("AGCTGACCTCCCTCCTGTGCCCAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("TGACCTCCCTCCTGTGTCCACAGAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("CCTCCCTCCTGTGCCTAGAGAGCCACAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("GAGCCATAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");
reads.emplace_back("CAGCCTGCTCGGTGGCCCGCATGGGATGAAACCTTT");


		assembley();
		for (AssemblyContig &contig : contigs) {
			if(contig.ending_reason[0] != 0 ||  contig.ending_reason[1] != 0) {continue;}
			fprintf(stderr,	"CONTIG seq:  [%s]\n ", contig.seq.c_str());
		}
	}

	std::vector<AssemblyReadInformation> readInfo;
private:

	AssemblerOptions o;

	//************************get kmer count***************************************/
	//buffs for getKmerCounts
	std::set<std::string> getKmerCounts_readWords_BUFF;
	void getKmerCounts(const unsigned wordLength, str_uint_map_t &wordCount,
			str_set_uint_map_t &wordSupportReads);

	//****************************getRepeatKmers*******************************************/
	str_pair_uint_map_t getRepeatKmers_wordIndices;
	std::vector<std::string> getRepeatKmers_wordStack;
	void getRepeatKmers(const str_uint_map_t &wordCount,
			std::set<std::string> &repeatWords);

	//**************************buildContigs********************************/
	str_uint_map_t buildContigs_wordCount;
	// records the supporting reads for each kmer
	str_set_uint_map_t buildContigs_wordSupportReads;
	// identify repeat kmers (i.e. circles from the de bruijn graph)
	std::set<std::string> buildContigs_repeatWords;
	std::set<std::string> buildContigs_unusedWords;
	bool buildContigs(const unsigned wordLength);

	/************************selectContigs*******************************/
	std::vector<AssemblyContig> tmpContigs;
	std::set<unsigned> selectContigs_usedReads;
	std::set<unsigned> selectContigs_usedPseudoReads;
	std::set<unsigned> selectContigs_contigs2Remove;
	std::set<unsigned> selectContigs_newSupportReads;
	void selectContigs(const unsigned normalReadCount);

	/***********************************Walk**************************************/
	//buffs for walk
	std::set<std::string> walk_wordsInContig;
	std::set<unsigned> walk_maxWordReads;
	std::set<unsigned> walk_maxContigWordReads;
	std::set<unsigned> walk_previousWordReads;
	std::set<unsigned> walk_supportReads2Remove;
	std::set<unsigned> walk_rejectReads2Add;
	std::set<unsigned> walk_contigWordReads;
	std::set<unsigned> walk_sharedReads;
	std::set<unsigned> walk_toRemove;
	std::set<unsigned> walk_toAdd;
	std::set<unsigned> walk_sharedReads_alleles;
	std::set<unsigned> walk_toUpdate;
	bool walk(const std::string &seed, const unsigned wordLength,
			const str_uint_map_t &wordCount,
			const str_set_uint_map_t &wordReads,
			const std::set<std::string> &repeatWords,
			std::set<std::string> &unusedWords, AssemblyContig &contig);

};

#endif /* CPP_LIB_ASSEMBLER_ASSEMBLER_HPP_ */
