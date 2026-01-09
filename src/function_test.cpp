/*
 * function_test.cpp
 *
 *  Created on: 2022-6-16
 *      Author: fenghe
 */

#ifndef FUNCTION_TEST_CPP_
#define FUNCTION_TEST_CPP_

#include "index_building/haplotype.hpp"
#include "function_test.hpp"
#include "mm_aln/deBGA_index.hpp"

int get_var_ALT_REF_same_length(){
	std::string varREF = "CAAAAAAAAA";
	std::string varALT = "C";
	HAP_VAR_ITEM var_v(46, varREF, varALT);
	std::string REF;
	REF = "CTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGA";

	fprintf(stderr, "%s \n", REF.c_str());
	fprintf(stderr, "%s ", var_v.REF.c_str());
	fprintf(stderr, "%s ", var_v.ALT.c_str());
	fprintf(stderr, "\n");

	//M1
	std::string ALT = REF;
	ALT.erase(ALT.begin() + var_v.ref_pos, ALT.begin() + var_v.ref_pos + var_v.REF.size());
	ALT.insert(ALT.begin() + var_v.ref_pos, var_v.ALT.begin(), var_v.ALT.end());

	fprintf(stderr, "%s\n" ,REF.c_str());
	fprintf(stderr, "%s\n" ,ALT.c_str());

	int SAM_len1 = 0;
	for(int i = var_v.ref_pos + 1; i < 300; i++){
		if(REF[i] == ALT[i])
			SAM_len1++;
		else
			break;
	}
	//M2
	int SAM_len2 = 0;
	int var_len = var_v.get_len();
	if(var_len < 0)//del
	{
		int max_search = REF.size() + var_len;
		for(int i = var_v.ref_pos + 1; i < max_search; i++){
			if(REF[i] == REF[i - var_len])
				SAM_len2++;
			else
				break;
		}
	}else if( var_len > 0){ //ins
		bool is_end = false;
		for(int i = 0; i < var_len; i++){
			if(REF[var_v.ref_pos + i + 1] == var_v.ALT[i + 1])
				SAM_len2++;
			else{
				is_end = true;
				break;
			}
		}
		if(!is_end){
			int max_search = REF.size() - var_len;
			for(int i = var_v.ref_pos + 1; i < max_search; i++){
				if(REF[i] == REF[i + var_len])
					SAM_len2++;
				else
					break;
			}
		}
	}else
		SAM_len2 = 0;

	xassert(SAM_len2 == SAM_len1, "FATAL_error\n");
	return SAM_len2;

}


struct Repeat_region{ int st; int ed; Repeat_region(int st_, int ed_):st(st_),ed(ed_){} };

//void clear_repeat_region(uint repeat_region_size){
//	//memset(&(repeat_region[0]), 0, repeat_region_size);
//	//already_clear = true;
//}
//
//	void simple_repeat_check(){
//		std::string  ref = "TAAATTATATAAATATAATATATATTTTATTATATAATATAATATATATTATATAAATATAATATATAAATTATATAATATAATATATATTATATAATATAATATATTTTATTATATAAATATATATTATATTATATAATATATATTTTATTATAT";
//		uint8_t *repeat_region = (uint8_t *)xcalloc(500,1);
//		bool with_data = false;
//		std::vector<Repeat_region> r;
//		r.clear();
//		const char * s = ref.c_str();
//		uint max_check_size = 25;
//		uint min_region_size = 3;
//		uint max_repeat_len = 0;
//		for(uint i = 1; i < max_check_size + 1; i++){
//			for(uint j = i; j < ref.size(); j++){
//				if(s[j] == s[j - i] && s[j] != 'N'){
//					max_repeat_len++;
//				}else{
//					if(max_repeat_len >= min_region_size){
//						uint region_ed = j;
//						uint region_bg = (region_ed-1)-((max_repeat_len) - 1);
//						with_data = true;
//						for(uint k = region_bg; k < region_ed; k++)
//							repeat_region[k] = 1;
//						if(true){
//							for(int i = region_bg; i <= region_ed; i++)
//								fprintf(stderr, "%c", ref[i]);
//							fprintf(stderr, "max_repeat_len %d, i %d , [%d, %d) \n", max_repeat_len, i, region_bg, region_ed);
//						}
//					}
//					max_repeat_len = 0;
//				}
//			}
//			if(max_repeat_len >= min_region_size){
//				uint region_ed = ref.size();
//				uint region_bg = (region_ed-1)-((max_repeat_len) - 1);
//				with_data = true;
//				for(uint k = region_bg; k < region_ed; k++)
//					repeat_region[k] = 1;
//				if(true){
//					for(int i = region_bg; i <= region_ed; i++)
//						fprintf(stderr, "%c", ref[i]);
//					fprintf(stderr, "max_repeat_len %d, i %d , [%d, %d) \n", max_repeat_len, i, region_bg, region_ed);
//				}
//			}
//			max_repeat_len = 0;
//		}
//		//finally:
//		if(with_data){
//			//get data
//			if(false)
//				fprintf(stderr, "REF：%s\n", ref.c_str());
//			uint region_bg = 0;
//			uint region_ed = 0;
//			for(uint j = 0; j < ref.size(); j++){
//				if(repeat_region[j] == 1){
//					region_ed = j+1;
//				}else{
//					if(region_ed > region_bg){
//						r.emplace_back(region_bg, region_ed);
//						if(true)
//							fprintf(stderr, "load region [%d, %d) \n", region_bg, region_ed);
//					}
//					region_bg = region_ed = j+1;
//				}
//			}
//			//final one
//			if(region_ed > region_bg){
//				r.emplace_back(region_bg, region_ed);
//				if(true)
//					fprintf(stderr, "load region [%d, %d) \n", region_bg, region_ed);
//			}
//			clear_repeat_region(ref.size());
//			with_data = false;
//		}
//	}

int function_text1(int argc, char *argv[]){
	//part5: reference
	Simple_ref_handler ref;
	ref.load_reference_from_file("/media/fenghe/data2/4_20/test_data/GRCh38_chr1_2.fa");
	//dump bin reference
	FILE * bin_ref_f = xopen("/media/fenghe/data2/4_20/index_dir/ref.bin", "wb");
	FILE * bin_ref_info_f = xopen("/media/fenghe/data2/4_20/index_dir/ref_info.bin", "wb");
	ref.dump_bin_ref(bin_ref_f, bin_ref_info_f);
	fclose(bin_ref_f);
	fclose(bin_ref_info_f);
	std::string rst;
	ref.load_ref_from_buff(0, 200000, 1000,  rst);
	fprintf(stderr, "%s\n", rst.c_str());

	Simple_ref_handler ref1;
	ref1.load_bin_ref("/media/fenghe/data2/4_20/index_dir/");
	ref1.load_ref_from_buff(0, 200000, 1000,  rst);
	fprintf(stderr, "%s\n", rst.c_str());

	return 0;

}

int function_text2(int argc, char *argv[]){
	std::string INDEX_PATH = argv[optind];
	MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
    //idx->load_all_index(INDEX_PATH.c_str(), NULL, true);
	idx->load_part_index(INDEX_PATH.c_str(), 11445375, 11445379);

    Simple_ref_handler ref1;
	ref1.load_bin_ref(INDEX_PATH.c_str());
	hap_string_loader_single_thread hl_r1;
	hap_string_loader_single_thread hl_DEBUG;
	std::string ref;
	FILE * outf = stdout;
	fprintf(stdout, "Total_wb_size %ld\n", idx->wb_info_size);
	char *var_string = (char *)xcalloc(1000000, 1);
	std::map<std::string, uint32_t> var_map;

	//Simple_repeat_checker src;
	//src.init(10000);

	for(uint wb_id = 0; wb_id < idx->wb_info_size; wb_id++){
		fprintf(stderr, " %d\n", wb_id);

		window_block_info &c_wb_info = idx->wb_info[wb_id];
		if(c_wb_info.is_SV()){ break; }

		ref1.load_ref_from_buff(c_wb_info.chrID, c_wb_info.region_st, c_wb_info.region_length, ref);
		String_list_and_var_list & A = hl_r1.get_string_list_and_var_list(wb_id, ref, (idx));
		A.print(stderr);
			continue;
	}
	return 0;
}


int function_text(int argc, char *argv[]){
	//get_var_ALT_REF_same_length();
	//simple_repeat_check();
	//return 0;
	std::string SV_S;
	fprintf(stderr, "hap_N %d\n", 1);
	Seed_handler seq_1_seed;
	seq_1_seed.init(10, 21, false, 500, 500);

	SV_S = "ATCACCTTCACAATGGAGTAC";
	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
	seq_1_seed.show_read_mm();

	//part 3: search the mm index of SV
//	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
//	seq_1_seed.show_read_mm();

	MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
	deBGA_INDEX *unitig_idx = (deBGA_INDEX *)new (deBGA_INDEX);
	//idx->load_all_index("/media/fenghe/data2/4_20/index_dir", "/media/fenghe/data3/panSVR_data_3-29/deBG_index/", false);
    //idx->load_all_index("/home/user/zhanganqi/wbenchmark/panGenomeKmer/test/index/2.5W/indexw5m22", "/home/user/zhanganqi/wbenchmark/debga/index/GRCh38_full_analysis_set_plus_decoy_hla_index_v1/", true);
    idx->load_all_index("/home/user/zhanganqi/wbenchmark/panGenomeKmer/cost/1KG/window10kmer21_0729/",
    		"/home/user/hitbio/data/commondata/index_GRCh38_plus_decoy_hla/", false);
    //load deBGA index
    //todo::

    //1
	SV_S = "GTTCAATTACCCCAACACAGGCAAAATTTCTTACTAACTCAAAAGGAACTGAATTCACATATAAAGACAAACGTGGGATTCTGTTGAGTTCTGTTGGACACAGAACACAATTGCTCAAGCACTGGTTGGGCACCTGTATTCTAATAGCTC";
	fprintf(stderr, "%s\n", SV_S.c_str() );
	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
	seq_1_seed.show_read_mm();

    seq_1_seed.search_seed(idx, true, false);

	fprintf(stderr, "read 1 SAME dir \n\n");
	for(auto & a : seq_1_seed.same_dir_seed_ref){
		a.print(stderr);
		uint32_t unitig_ID = (a).unitig_ID;
		uint64_t ref_pos_n = unitig_idx->get_ref_N(unitig_ID);
		uint32_t unitig_offset = (a).unitig_offset - (21 - 1);//0-base k_t=21

		fprintf(stderr, "A: %u %u\n", unitig_offset, unitig_ID);

		uint64_t global_offset = unitig_idx->get_global_offset(unitig_ID, 0, 0) + 1;//0-base

		fprintf(stderr, "A: g: %u\n", global_offset);
		uint32_t unitig_length = unitig_idx->get_UNITIG_length(unitig_ID);
		//idx->ref.load_ref_from_buff_global_offset(global_offset, unitig_length, buffer_seq, &(idx->unitig_idx));//1-base

		uint32_t uni_offset_s_l = unitig_offset;
		uint32_t uni_offset_s_r = unitig_length - uni_offset_s_l - 21;//我不能理解具体含义
		fprintf(stderr, "A: l: %u r: %u length: %d\n", uni_offset_s_l, uni_offset_s_r, unitig_length);
	}

	fprintf(stderr, "read 1 NOT SAME dir \n\n");
	for(auto & a : seq_1_seed.not_same_dir_seed_ref){
		uint32_t unitig_ID = (a).unitig_ID;
		uint64_t ref_pos_n = unitig_idx->get_ref_N(unitig_ID);
		uint32_t unitig_offset = (a).unitig_offset - (21 - 1);//0-base k_t=21

		fprintf(stderr, "A: %u %u\n", unitig_offset, unitig_ID);

		uint64_t global_offset = unitig_idx->get_global_offset(unitig_ID, 0, 0) + 1;//0-base

		fprintf(stderr, "A: g: %u\n", global_offset);
		uint32_t unitig_length = unitig_idx->get_UNITIG_length(unitig_ID);
		//idx->ref.load_ref_from_buff_global_offset(global_offset, unitig_length, buffer_seq, &(idx->unitig_idx));//1-base

		uint32_t uni_offset_s_l = unitig_offset;
		uint32_t uni_offset_s_r = unitig_length - uni_offset_s_l - 21;//我不能理解具体含义
		fprintf(stderr, "A: l: %u r: %u length: %d\n", uni_offset_s_l, uni_offset_s_r, unitig_length);
	}
	//return 0;
	//2
	SV_S = "CACAGCTGAAGATCACTTATAGTCAATAACACTTGTCAAACTTTAATTTTAAATTCATTAAGTCTCTCCCTAATGTATCGTACTTTTTAATAAAGGAAACTGTCAGATACGCTATAATACAGAGGAAAACGTTTATATACTCTTTCATTA";
	fprintf(stderr, "%s\n", SV_S.c_str() );
	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
	seq_1_seed.show_read_mm();

    seq_1_seed.search_seed(idx, true, false);

	fprintf(stderr, "read 1 SAME dir \n\n");
	for(auto & a : seq_1_seed.same_dir_seed_ref){
		a.print(stderr);
	}

	fprintf(stderr, "read 1 NOT SAME dir \n\n");
	for(auto & a : seq_1_seed.not_same_dir_seed_ref){
		a.print(stderr);
	}
	//3
	SV_S = "TCCCAGGAGAGAAGGGATATCTTCTTCTATGCCCCTACTATGAACTAGCATATCCTTTGAGTCCTTTCTTTGTATGAAGCACTGGGGAAATTTAATCAACTTCATAATCACTTATAATTCTCACAGTAATTTTATGAGATAAGTTCTATT";
	fprintf(stderr, "%s\n", SV_S.c_str() );
	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
	seq_1_seed.show_read_mm();

    seq_1_seed.search_seed(idx, true, false);

	fprintf(stderr, "read 1 SAME dir \n\n");
	for(auto & a : seq_1_seed.same_dir_seed_ref){
		a.print(stderr);
	}

	fprintf(stderr, "read 1 NOT SAME dir \n\n");
	for(auto & a : seq_1_seed.not_same_dir_seed_ref){
		a.print(stderr);
	}

	//4
	SV_S = "ATCTCTAATATTATCTCATTAAGTGCCCTCTAACTACAGAATGACTCCTCCTTCTGGGTTTCCTCTGCCTGAATGGCACCACCATTCAGTGGCCCACGACAAATGTGGACCTCATTTGGAACTTCTCCCTTTCCCTCATGTTTCATTTTC";
	fprintf(stderr, "%s\n", SV_S.c_str() );
	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
	seq_1_seed.show_read_mm();

    seq_1_seed.search_seed(idx, true, false);

	fprintf(stderr, "read 1 SAME dir \n\n");
	for(auto & a : seq_1_seed.same_dir_seed_ref){
		a.print(stderr);
	}

	fprintf(stderr, "read 1 NOT SAME dir \n\n");
	for(auto & a : seq_1_seed.not_same_dir_seed_ref){
		a.print(stderr);
	}

    unitig_idx->show_all_unitig_global_offset(stderr);

	return 0;
//
	hap_string_loader_single_thread hl_r1;
//	uint32_t * hap_n_stat = (uint32_t *)xcalloc(500,  4);
//
//	Simple_ref_handler ref1;
//	ref1.load_bin_ref("/home/user/zhanganqi/wbenchmark/panGenomeKmer/test/index/2.5W/indexw5m22/");
//	//ref1.load_bin_ref("/media/fenghe/data2/4_20/index_dir");
//	std::string ref_str;
	int window_ID = hl_r1.get_windows_ID(21, 20005776, idx);
//	int hap_N = hl_r1.get_hap_N(window_ID, idx);
//	ref1.load_ref_from_buff(idx->wb_info[window_ID].chrID, idx->wb_info[window_ID].region_st, idx->wb_info[window_ID].region_length,  ref_str);
////	for(int i = 0; i < hap_N;i++){
////    	hl_r1.get_string(window_ID, i, ref_str, idx);
////    	hl_r1.printf_modify_info_after_get_string();
////	}
//
//	String_list_and_var_list & A = hl_r1.get_string_list_and_var_list(window_ID, ref_str, (idx));
//	A.print();
//
//	return 0;
//	double cpu_time = cputime();
//    //count the wb data
//    for(int i = 0; i < idx->wb_info_size; i++){
//    	if(idx->wb_info[i].flag == 1){
//    		break;
//    	}
//    	ref1.load_ref_from_buff(idx->wb_info[i].chrID, idx->wb_info[i].region_st, idx->wb_info[i].region_length,  ref_str);
//    	std::vector<std::string> & A = hl_r1.get_string_list(i, ref_str, idx);
//
////		for(auto & S: A){
////			fprintf(stderr, "S: %s\n", S.c_str());
////		}
//
//    	//int hap_N = hl_r1.get_hap_N(i, idx);
//    	//if(hap_N < 500){
//    	//	hap_n_stat[hap_N]++;
//    	//}else{
//    	//	hap_n_stat[499]++;
//    	//}
//    }
//
//    fprintf(stderr, "Classify CPU: %.3f sec\n", cputime() - cpu_time);
//
//    for(int i = 0; i < 500; i++){
//    	fprintf(stderr, "hap_n_stat i %d N %d \n", i, hap_n_stat[i]);
//    }
//	std::string SV_S = "TAACAAACCTGCACATGTACCTCCTGAACCTAAAATAGAAGTTGAGAAAAGAAAAAGAAGATATCAGTTTCCTATAGGGAGCTCCAGGTTTGTTCCAAACAGTAGAGAGACAAACCTATCAAGGATGAAATATATACAAGAGTAGGCACA";
//
//	seq_1_seed.init(10, 22, false, 500, 500);
//	//part 3: search the mm index of SV
//	//std::string SV_S = "CTTTCTTTTTTTTTTTTTTTTTTTTTTGAGGAGTTCCTTGTCGCCGCTGGGTGGCGGCGCGATTGCTCCTGCAGCTCCGCCCCCGTCCCCATTCCTGCCTCGCCTCCCAAGTACTGGACTCAGCGCCCCCTCGCCCGGCTAATTTTTGTATTTTTAGTAAGACGTTTCCGTTTAGCGGGGTTCGATCTCTGACTTCGTGTCCTCCGCCTCGCTCCCAGTGTGATTACAGCTGACCACCCCCCCAG";
//	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
//	seq_1_seed.search_seed(idx, true, false);
//
//	fprintf(stderr, "read 1 SAME dir \n\n");
//
//	count_hit_number_ref1(seq_1_seed.same_dir_seed_ref, idx, stderr, 500);
//	fprintf(stderr, "read 1 NOT SAME dir \n\n");
//
//	count_hit_number_ref1(seq_1_seed.not_same_dir_seed_ref, idx, stderr, 500);
//
//	SV_S = "CAACATGTGATAACTATGTGAAGACATGAATATGTTAATTAGCTTGATAGTGAAAATCATTTTAATATATACTTGATGGGCGCAGGGTGAGAGGAGAGTGAGGGTCAAAACACTACCTGCTGGGGACTATGCTCACTAGCTGGGTGATGA";
//	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
//	seq_1_seed.search_seed(idx, true, false);
//
//	fprintf(stderr, "read 2 SAME dir \n\n");
//
//	count_hit_number_ref1(seq_1_seed.same_dir_seed_ref, idx, stderr, 500);
//
//	fprintf(stderr, "read 2 NOT SAME dir \n\n");
//
//	count_hit_number_ref1(seq_1_seed.not_same_dir_seed_ref, idx, stderr, 500);
//

//
//	for(auto & a: seq_1_seed.same_dir_seed_ref){
//		uint64_t global_offset = unitig_idx->get_global_offset(a.unitig_ID, a.unitig_offset - 21, 0) + 1;//0-base
//		int unitig_length = unitig_idx->get_UNITIG_length(a.unitig_ID);
//		std::string B;
//		idx->ref.load_ref_from_buff_global_offset(global_offset, unitig_length, B, &(idx->unitig_idx));//1-base
//		fprintf(stderr, "B %s\n", B.c_str());
//	}
//
//
////	return 1;
//
//
//	//seed index pointer
//
//	int chr_ID, POS;
//	unitig_idx->get_chr_ID_POS(308924372, chr_ID, POS);
//
//	//part 1: search the SNP/INDEL wb
//	window_ID = hl_r1.get_windows_ID(1, 76538473, idx);
//	window_block_info &c_wb_info = idx->wb_info[window_ID];
//	c_wb_info.show_wb(stderr);
//	hap_N = hl_r1.get_hap_N(window_ID, idx);
//	fprintf(stderr, "hap_N %d\n", hap_N);
//
	std::string ref;

	for(int i = 0; i < idx->wb_data_size; i++){
		window_block_info &c_wb_info = idx->wb_info[i];
		fprintf(stderr, "ID %d\n", i);
		c_wb_info.show_wb(stdout);
		idx->ref.load_ref_from_buff(c_wb_info.chrID, c_wb_info.region_st, c_wb_info.region_length, ref);
		String_list_and_var_list & A = hl_r1.get_string_list_and_var_list(i, ref, (idx));
		A.print(stdout);
	}

//	//part 2: search the SV wb
//	std::vector<uint32_t> & SV_wb_list = hl_r1.search_SV_IN_region(1, 7650000, 7960000, (idx));
//	for(uint32_t SV_wb: SV_wb_list){
//		window_block_info &c_wb_info = idx->wb_info[SV_wb];
//		c_wb_info.show_wb(stderr);
//		fprintf(stderr, "c_wb_info.is_SV()? %d\n", c_wb_info.is_SV());
//		ref.clear();
//		std::vector<char> & S = hl_r1.get_string(SV_wb, 0, ref,  &(idx));
//		fprintf(stderr, "SV: S: %s\n", &(S[0]));
//	}
//
//	seq_1_seed.init(5, 28, false, 500, 500);
//	//part 3: search the mm index of SV
//	//std::string SV_S = "CTTTCTTTTTTTTTTTTTTTTTTTTTTGAGGAGTTCCTTGTCGCCGCTGGGTGGCGGCGCGATTGCTCCTGCAGCTCCGCCCCCGTCCCCATTCCTGCCTCGCCTCCCAAGTACTGGACTCAGCGCCCCCTCGCCCGGCTAATTTTTGTATTTTTAGTAAGACGTTTCCGTTTAGCGGGGTTCGATCTCTGACTTCGTGTCCTCCGCCTCGCTCCCAGTGTGATTACAGCTGACCACCCCCCCAG";
//	seq_1_seed.read_mm(SV_S.c_str(), SV_S.size());
//	seq_1_seed.search_seed(idx, false, false);
//	//part 4: search the mm index of the SNP
//	for(auto & seed : seq_1_seed.same_dir_seed_alt){
//		seed.print(stderr);
//		idx->wb_info[seed.window_ID].show_wb(stderr);
//	}
//	for(auto & seed : seq_1_seed.not_same_dir_seed_alt){
//		seed.print(stderr);
//		idx->wb_info[seed.window_ID].show_wb(stderr);
//	}
//	//part 4: search the mm index of UNITIG
//	std::string SV_S_1 = "TGCCTGTTTCTCCACAAAGTGTTTACTTTTGGATTTTTGCCAGTCTAACAGGTGAAGCCCTGGAGATTCTTATTAGTGATTTGGGCTGGGGCCTGGCCATGTGTATTTTTTTAAATTTCCACTGATGATTTTGCTGCATG";
//	seq_1_seed.read_mm(SV_S_1.c_str(), SV_S_1.size());
//	seq_1_seed.search_seed(idx, false, false);
//	//part 4: search the mm index of the SNP
//	for(auto & seed : seq_1_seed.same_dir_seed_ref){
//		seed.print(stderr);
//		uint64_t ref_N = unitig_idx->get_ref_N(seed.unitig_ID);
//		for(uint64_t i = 0; i <ref_N; i++){
//			uint64_t gloabl_offset = unitig_idx->get_global_offset(seed.unitig_ID, seed.unitig_offset, i);
//			int ref_chr_ID; int ref_POS;
//				unitig_idx->get_chr_ID_POS(gloabl_offset, ref_chr_ID, ref_POS);
//				fprintf(stderr, "UNITIG: gloabl_offset %ld, ref_chr_ID %d, ref_POS %d, s.mm.read_pos %d unitig_offset %d\n",  gloabl_offset, ref_chr_ID, ref_POS, seed.mm.read_pos, seed.unitig_offset);
//		}
//		idx->wb_info[seed.unitig_ID].show_wb(stderr);
//	}
//	for(auto & seed : seq_1_seed.not_same_dir_seed_ref){
//		seed.print(stderr);
//		uint64_t ref_N = unitig_idx->get_ref_N(seed.unitig_ID);
//		for(uint64_t i = 0; i <ref_N; i++){
//			uint64_t gloabl_offset = unitig_idx->get_global_offset(seed.unitig_ID, seed.unitig_offset, i);
//			int ref_chr_ID; int ref_POS;
//				unitig_idx->get_chr_ID_POS(gloabl_offset, ref_chr_ID, ref_POS);
//				fprintf(stderr, "UNITIG: gloabl_offset %ld, ref_chr_ID %d, ref_POS %d, s.mm.read_pos %d unitig_offset %d\n",  gloabl_offset, ref_chr_ID, ref_POS, seed.mm.read_pos, seed.unitig_offset);
//		}
//		idx->wb_info[seed.unitig_ID].show_wb(stderr);
//	}

	return 0;
}

int minimizer_generate(int argc, char *argv[]) {
	//loading strings
	char * alt_str_fn = argv[1];
	//minimizer settings
	int minimizer_len = atol(argv[2]);
	int window_size = atol(argv[3]);

	gzFile alt_str_file = gzopen(alt_str_fn, "rb");

	char *temp = new char[MAX_LINE_LENGTH];//10M
	char *analysis_line = new char[MAX_LINE_LENGTH];//10M
	std::vector<std::string> item_value;
	std::vector<mm128_t> mm_v;
	Minimizer_generater mg;
	mg.init(window_size, minimizer_len, false);
	int total_minimizer_num = 0;
	int total_alt_string_len = 0;
	int total_alt_string_number = 0;

	while(true){
		bool is_the_final_block = false;
		if(NULL == gzgets(alt_str_file, analysis_line, MAX_LINE_LENGTH))
			is_the_final_block = true;
		else{
			//remove the final \n
			analysis_line[strlen(analysis_line) - 1] = 0;//todo:: need check whether to remove
			if(analysis_line[0] == '>')//skip contig NAME line in .fa format
				analysis_line[0] = 0;
			split_string(item_value, temp, analysis_line, "N");
		}

		if(is_the_final_block)
			break;
		else{
			for(auto & alt_str : item_value){
				mm_v.clear();
				mg.mm_sketch(alt_str.c_str(), alt_str.size(), mm_v);
				int mini_number = mm_v.size();
				total_minimizer_num += mini_number;
				total_alt_string_len += alt_str.size();
				total_alt_string_number ++;
				//if(total_alt_string_number % 10 == 0){ fprintf(stderr, "[total_alt_string_number: %d] [total_alt_string_len:  %d][total_minimizer_num: %d]\n", total_alt_string_number, total_alt_string_len, total_minimizer_num);	}
				if(false) {
					fprintf(stderr, "string [%s] has %d minimizer\n", alt_str.c_str(), mini_number);
					mg.print_minimizer_list(mm_v);
				}
			}
		}
	}
	fprintf(stderr, "[total_minimizer_num]: %d\n", total_minimizer_num);
	gzclose(alt_str_file);
	return 0;
}
#endif /* FUNCTION_TEST_CPP_ */
