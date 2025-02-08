/*
 * fa_stat.hpp
 *
 *  Created on: 2023年12月29日
 *      Author: fenghe
 */

#ifndef SV_SPA_CORE_FA_STAT_HPP_
#define SV_SPA_CORE_FA_STAT_HPP_

extern "C"
{
#include "../clib/bam_file.h"
}

struct FA_STATUS{

	int N_count(const std::string &word){
		int N = 0;
		const char *s = word.c_str();
		int n = word.size();
		for(int i = 0; i < n; i++){
			if(s[i] == 'N')
				N++;
		}
		return N;
	}

	int run_ana(int argc, char *argv[]){
		//bs->load(nova_para->statsFileName);
		//	//PART5: for each choromosome
		if(argc != 1){
			fprintf(stderr, "WARNING: the reference file is needed\n\n");
			fprintf(stderr, "Function: Indexing the reference for SV calling\n");
			fprintf(stderr, "Usage: fa_stat [reference] > ref.idx\n");
			return -1;
		}
		char * fa_fn_in = argv[0];//separate by ','
		//char * bed_fn_out = argv[2];
		//FILE * bed_file = xopen(bed_fn_out, "w");
		//default parameters
		int REGION_LEN = 600;
		int REGION_STEP_LEN = 300;
		int MIN_K_size = 20;
		int MAX_K_size = 180;
		int K_step_LEN = 20;

		//test file exist
		FILE* try_open = xopen(fa_fn_in, "r");
		fclose(try_open);
		faidx_t * c_ref_idx = reference_index_load(fa_fn_in);
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(26, N_seq);

		//header line
		fprintf(stdout, "%d\t%d\t%d\n", REGION_LEN, REGION_STEP_LEN, 0);
		std::set<std::string> kmerCountingSet;
		for(int i = 0; i < N_seq; i++){
			int len = 0;
			const char * chrName = faidx_iseq(c_ref_idx, i);
			char * c_reference = fai_fetch(c_ref_idx, chrName, &len);
			int beginPos = 0; int endPos = 0; int LENGTH = REGION_LEN;
			int chrLength = faidx_seq_len(c_ref_idx, chrName);
			std::string regionString;
			for(;beginPos < chrLength; beginPos += REGION_STEP_LEN){
				endPos = beginPos + LENGTH;
				endPos = MIN(endPos, chrLength);
				regionString.clear();
				for(int j = beginPos; j < endPos; j++){
					regionString += c_reference[j];
				}
				bool withloop = false;
				int kmerSizeF = 0;
				for(int kmerSize = MIN_K_size; kmerSize <= MAX_K_size; kmerSize += K_step_LEN){
					kmerCountingSet.clear();
					const std::string &seq = regionString;
					const unsigned seqLen = seq.size();
					// track all words from the read, including repetitive k
					unsigned kmer_number = seqLen - kmerSize + 1;
					for (unsigned j = 0; j < kmer_number; ++j) {
						const std::string word(seq.substr(j, kmerSize));
						if (N_count(word) > 0 )
							continue;
						std::set<std::string>::iterator it = kmerCountingSet.find(word);
						if(it == kmerCountingSet.end()){
							kmerCountingSet.emplace(word);
						}else{
							withloop = true;
							kmerSizeF = kmerSize;
							break;
						}
					}
					if(!withloop){//am.DBGTest(kmerSize)){
						break;
					}
				}
				//debug show results
				//if(REPEAT_DEBUG) {for( DBG_INFO_item_ * kmerIdx_i : id2info) kmerIdx_i->showData();}
				if(withloop){
					fprintf(stdout, "%d\t%d\t%d\n",i, beginPos, kmerSizeF);
				}
			}
			free(c_reference);
		}
		//fclose(bed_file);
		return 0;
	}

};

#endif /* SV_SPA_CORE_FA_STAT_HPP_ */
