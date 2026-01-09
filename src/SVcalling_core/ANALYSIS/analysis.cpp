/*
 * analysis.cpp
 *
 *  Created on: 2025/7/9
 *      Author: fenghe
 */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <math.h>
#include "cpp_lib/get_option_cpp.hpp"
#include "cpp_lib/RefRegion.hpp"
#include "analysis.hpp"
#include "VNTR_analysis.hpp"
#include "Genome_context_analysis.hpp"
#include "ASM_analysis_handler.hpp"

#include "FC_joint_calling.hpp"

//VNTR---------------------------------------------------------
int vntr_analysis_run(int argc, char *argv[]){
	vntr_analysis_handler r;
	return r.vntr_analysis(argc, argv);
}

int bed_in_repeat_masker_region_run(int argc, char *argv[]){
	vntr_analysis_handler r;
	return r.bed_in_repeat_masker_region(argc, argv);
}

//CPX
int complex_sv_region_detection_run(int argc, char *argv[]){
	COMPLEX_SV_region r;
	return r.complex_sv_region_detection(argc, argv);
}

//CPX
int complex_sv_region_detection_500W_run(int argc, char *argv[]){
	COMPLEX_SV_region r;
	return r.complex_sv_region_detection_500W(argc, argv);
}

int vcf_2_contig_run(int argc, char *argv[]){
	COMPLEX_SV_region r;
	return r.vcf_2_contig(argc, argv);
}

int contig_file_split_run(int argc, char *argv[]){
	COMPLEX_SV_region r;
	return r.contig_file_split(argc, argv);
}

//REF---------------------------------------------------------
int getReverseStr_run(int argc, char *argv[]){
	Ref_region_handler r;
	return r.getReverseStr_simple(argc, argv);
}

int ref_split_run(int argc, char *argv[]){
	Ref_region_handler r;
	return r.ref_split(argc, argv);
}

int get_ref_region_run(int argc, char *argv[]){
	Ref_region_handler r;
	return r.get_ref_region(argc, argv);
}

int dump_ref_by_region_run(int argc, char *argv[]){
	Ref_region_handler r;
	return r.dump_ref_by_region(argc, argv);
}

//GC analysis---------------------------------------------------------
int genomic_contest_AnalysisNGS_run(int argc, char *argv[]){
	Genomic_context_Analysis_handler r;
	return r.genomic_contest_AnalysisNGS(argc, argv);
}

int genomic_contest_AnalysisCCS_run(int argc, char *argv[]){
	Genomic_context_Analysis_handler r;
	return r.genomic_contest_AnalysisCCS(argc, argv);
}

int heatMapData_run(int argc, char *argv[]){
	Genomic_context_Analysis_handler r;
	return r.heatMapData(argc, argv);
}

int heatMapDataSum_run(int argc, char *argv[]){
	Genomic_context_Analysis_handler r;
	return r.heatMapDataSum(argc, argv);
}

//REFERENCE STAT analysis---------------------------------------------------------
int fa_stat_full_run(int argc, char *argv[]){
	FA_stat_handler r;
	return r.fa_stat_full(argc, argv);
}

int string_local_repeat_run(int argc, char *argv[]){
	FA_stat_handler r;
	return r.string_local_repeat(argc, argv);
}

//BAM  analysis---------------------------------------------------------
int alignment_Entropy_Analysis_run(int argc, char *argv[]){
	Alignment_Entropy_Analysis r;
	return r.ana_run_Entropy(argc, argv);
}

int alignment_SMALL_VAR_DENSITY_Analysis_run(int argc, char *argv[]){
	SMALL_VAR_DENSITY r;
	return r.ana_run_SMALL_VAR_DENSITY(argc, argv);
}

//FC---------------------------------------------------------
int bed_merging_run(int argc, char *argv[]){
	BED_MERGING_handler r;
	return r.bed_merging(argc, argv);
}

int FC_joint_ana_run(int argc, char *argv[]){
	FC_joint_calling_handler r;
	return r.ana_main(argc, argv);
}

int align_main_run(int argc, char *argv[]){
	FC_joint_calling_handler r;
	return r.align_main(argc, argv);
}

//T2T calling---------------------------------------------------------

int T2T_calling_ana_run(int argc, char *argv[]){
	ASM_ANALYSIS_HANDLER r;
	return r.run(argc, argv);
}


int t2t_split_alignment_analysis_run(int argc, char *argv[]){
	ASM_ANALYSIS_HANDLER r;
	return r.t2t_split_alignment_analysis(argc, argv);
}

//------------------------------MAIN--------------------------------------//
int analysis_main(int argc, char *argv[])
{
	COMMAND_HANDLER ch;
	ch.set_main_command("analysis");

	ch.add_function("vntr_analysis", "", vntr_analysis_run);
	ch.add_function("bed_in_repeat_masker_region", "", bed_in_repeat_masker_region_run);

	ch.add_function("csv_region", "usage: tools csv_region [VCF_f] [Distance] [min_SV_len]", complex_sv_region_detection_run);
	ch.add_function("csv_region_500w", "usage: tools csv_region [VCF_f] [Distance] [min_SV_len]", complex_sv_region_detection_500W_run);
	ch.add_function("contig_file_split", "usage: ?????", contig_file_split_run);

	ch.add_function("dump_ref_by_region", "", dump_ref_by_region_run);
	ch.add_function("get_ref_region", "get_ref_region hg38.fa chr1:100000-2000000", get_ref_region_run);
	ch.add_function("ref_split", "Split reference file by chr_ID [input.fa]", ref_split_run);
	ch.add_function("getReverseStr", "Get Reverse Str of a DNA string [string (ACGT.....)]", getReverseStr_run);

	ch.add_function("genomic_contest_AnalysisNGS", "", genomic_contest_AnalysisNGS_run);
	ch.add_function("genomic_contest_AnalysisCCS", "", genomic_contest_AnalysisCCS_run);
	ch.add_function("genomic_contest_heatMapData", "", heatMapData_run);
	ch.add_function("genomic_contest_heatMapDataSum", "", heatMapDataSum_run);

	ch.add_function("fa_stat_full", "", fa_stat_full_run);
	ch.add_function("string_local_repeat", "", string_local_repeat_run);

	ch.add_function("vcf_2_contig", "", vcf_2_contig_run);
	ch.add_function("entropy_analysis", "", alignment_Entropy_Analysis_run);
	ch.add_function("small_var_den", "", alignment_SMALL_VAR_DENSITY_Analysis_run);

	ch.add_function("bed_merging", "tools bed_merging [ref.fa] [input.bed]", bed_merging_run);
	ch.add_function("FC_joint_ana", "tools FC_joint_ana [contig.fa] [RepeatMasker.out] [TRF.anno] [TRUR/FALSE]", FC_joint_ana_run);
	ch.add_function("align_seq", "tools align_seq", align_main_run);

	ch.add_function("T2T_calling_ana", "tools T2T_calling_ana", T2T_calling_ana_run);
	ch.add_function("t2t_split_alignment_analysis", "", t2t_split_alignment_analysis_run);

	return ch.run(argc, argv);

//	//not used functions
//
//
//


//

//	ch.add_function("vcf_add_alt_string", "", vcf_add_alt_string);
//	//ch.add_function("randomGenerateSV", "Generate 20000 random SV [sam_header_fn.sam, ref.fa, int_rand_seed]", randomGenerateSV);
//	//ch.add_function("ROC_PR", "", analysis_ROC_PR);
//	//ch.add_function("signalPattenAnalysis", "", signalPattenAnalysis);
//	//ch.add_function("UCSC_mappibility_binary_dump", "", UCSC_mappibility_binary_dump);

//	ch.add_function("combine_sort_vcf", "Combine vcf records from multiple files", combine_sort_vcf);
//	ch.add_function("isize_count", "Count ISIZE for a bam file [input.bam]", isize_count);
//	ch.add_function("bam2Fastq", "Convert bam into fastq [input.bam, output.fq]", bam2Fastq);
//	ch.add_function("bamDump", "Dump first N record in a bam file [input.bam, ref.fa, output.bam, N]", bamDump);
//	ch.add_function("gz_head", "Get the N characters from offset P of a XXX.gz file.[input.gz, N, P]", gz_head);
//	ch.add_help_msg_back("zlib only support SEEK_SET, and read from begin P is at most 200M, not support from-where = 2: SEEK_END");
//	ch.add_function("read_ACGT_analysis", "Analysis the acgt distribution of cram read", read_ACGT_analysis);
//	ch.add_function("vcf_dump", "Filter and dump the SV items in VCF file using 'sample_ID'， ‘SV_TYPE‘ or ’chrID’", vcf_dump);
//	ch.add_help_msg_back("[in_fn, out_fn, sample_ID, SV_TYPE, chrID], set 'ALL' for options the get all sample or SV type");

	//not used functions
	//else if (strcmp(argv[1],"ref_dump") == 0) 		ref_dump(argv[2], atoi(argv[3]));
	//else if (strcmp(argv[1],"SH_bam") == 0) 		get_all_record_with_SH(argv[2], argv[3]);
	//else if (strcmp(argv[1],"DR_bam") == 0) 		get_all_record_with_DR(argv[2], argv[3]);
	//else if (strcmp(argv[1],"SA_bam") == 0) 		get_all_record_with_SA(argv[2], argv[3]);
	//else if (strcmp(argv[1],"vcf_nstd_dump") == 0) 	vcf_nstd_dump(argv[2], argv[3]);
	//else if (strcmp(argv[1],"vcf_manta_dump") == 0) vcf_manta_dump(argv[2], argv[3]);
	//else if (strcmp(argv[1],"vcf_jlra_dump") == 0)  vcf_jlra_dump(argv[2], argv[3]);
	//else if (strcmp(argv[1],"vcf_GIAB_getSV") == 0) vcf_GIAB_getSV(argv[2], argv[3]);
	//else if (strcmp(argv[1],"vcf_SURVIVOR_getSV") == 0) vcf_SURVIVOR_getSV(argv[2], argv[3]);
	//else if (strcmp(argv[1],"SVRegionCombine") == 0) SVRegionCombine(argv[2], atoi(argv[3]));
	//else if (strcmp(argv[1],"GIAB_SV_region_full_test") == 0) GIAB_SV_region_full_test(argv[2], argv[3], argv[4], atoi(argv[5]));
	//GIAB test phase 2

	//else if (strcmp(argv[1],"bamDump_discard_both_unmapped") == 0) bamDump_discard_both_unmapped(argv[2], argv[3]);
	//else if (strcmp(argv[1],"bamDump_discard_map_len_100") == 0) bamDump_discard_map_len_100(argv[2], argv[3], atoi(argv[4]));
	//else if (strcmp(argv[1],"SP_region_result_full_test") == 0) SP_region_result_full_test(argv[2]);
	//else if (strcmp(argv[1],"vcf_sample") == 0)		vcf_sample(argv[2], argv[3], argv[4], argv[5], argv[6]);
	//else if (strcmp(argv[1],"vcf_compare") == 0) 	vcf_compare(argv[2], argv[3], argv[4], argv[5], argv[6]);
		 //SV test phase 3
	//else if (strcmp(argv[1],"liftoverBuildingIndex") == 0)	liftoverBuildingIndex(argv[2], argv[3]);
	//else if (strcmp(argv[1],"liftoverSearch") == 0)	liftoverSearchRegion(argv[2]);
	//else if (strcmp(argv[1],"liftoverRead") == 0)	liftoverRead(argv[2], argv[3]);
		 //SV calling: phase 4
	//else if (strcmp(argv[1],"simple_kmer_counter") == 0)	simple_kmer_counter(argv[2]);

	return 0;

}

