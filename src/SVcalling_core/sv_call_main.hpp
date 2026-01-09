/*
 * sv_call_main.hpp
 *
 *  Created on: 2020-4-27
 *      Author: fenghe
 */

#include <getopt.h>
#include <string.h>
#include <SVcalling_core/SV_core.hpp>
#include <string>
#include <iostream>
#include <vector>

#include "../cpp_lib/get_option_cpp.hpp"
#include "../SVcalling_core/bam_stat.hpp"

struct NOVA_SV_PARA{
	//Necessary
	char*	referenceFilename;
	char*   SRS_BAM_F;
	char*   LRS_BAM_F;

	char*   LRS_Tumor_BAM_F;
	char*   LRS_Normal_BAM_F;

	char * FC_BED_F;
	char * OUT_MODE_char;
	int OUT_MODE_int;
	//char * FC_Contig_L;

	char* 	outputFile;

	char*	realigned_Filename;
	char*	TL_read_Filename;
	char* 	statsFileName;
	char* 	ref_idx_FileName;

	char*   KMER_ERROR_PROFILE_File;

	bool 	is_compression;
	int 	MIN_SV_len;

	char*   preset_platform_str;
	int   preset_platform;

	char*   INNER_TEST_FLAG;

	bool not_output_vcf_header;
	int random_phasing;
	int output_small_var;
	bool force_calling_ALU;

	int st_chr_ID; int st_pos;
	int ed_chr_ID; int ed_pos;

	int get_option(int argc, char *argv[]){
		char default_outputFile[] = "stdout";
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  call  [Options]\n");
	    l.add_title_string("\n\n");
	    l.add_title_string("     Necessary parameters:\n");
	    l.add_title_string("     \t\tAt least one BAM/CRAM is needed;The reference file is necessary.");
	    l.add_option("reference",  		'r', "the reference (.fa)", false, ""); l.set_arg_pointer_back((void *)&referenceFilename);
	    l.add_help_msg_back("\tGermline SVs calling:");
	    l.add_option("LRS_BAM", 'l', "The long read sequencing (LRS: pacBIO-HiFi, ONT-hac or ONT-sup)(At least Q26) BAM file", false, "");
	    l.set_arg_pointer_back((void *)&LRS_BAM_F);
		l.add_option("NGS_BAM", 'n', "The NGS BAM file", false, ""); l.set_arg_pointer_back((void *)&SRS_BAM_F);
		l.add_help_msg_back("\tSomtic SVs calling:(only support LRS (ONT-Q20 or HiFi) right now)");
		l.add_option("LRS_Tumor_BAM", 'X', "The tumor BAM file", false, ""); l.set_arg_pointer_back((void *)&LRS_Tumor_BAM_F);
		l.add_option("LRS_Normal_BAM", 'x', "The normal BAM file", false, ""); l.set_arg_pointer_back((void *)&LRS_Normal_BAM_F);
		l.add_help_msg_back("");
		l.add_help_msg_back("Optional parameters:");
		l.add_help_msg_back("\tRegion parameters:(SV calling in a specific REGION)");
		l.add_help_msg_back("\t\tExample 1: using: '-S 0 -E 0 -s 100000 -F 200000' to call SVs in region: chr1:100000-200000");
		l.add_help_msg_back("\t\tExample 2: using: '-S 0 -E 23'(default) to call SVs in the whole-genome: from chr1 to chrY");
		l.add_option("st_chr_ID", 		'S', "The start chr_ID(0-base) for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_chr_ID);
		l.add_option("ed_chr_ID", 		'E', "The end chr_ID(0-base) for SV calling", true, 23); l.set_arg_pointer_back((void *)&ed_chr_ID);
		l.add_option("st_pos", 			's', "The start position for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_pos);
		l.add_option("ed_pos", 			'F', "The end position for SV calling", true, 500000000); l.set_arg_pointer_back((void *)&ed_pos);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tOutput parameters:");
		l.add_option("output_format",   'B', "The result format. One of VCF or BED", true, "VCF"); l.set_arg_pointer_back((void *)&OUT_MODE_char);
		l.add_option("random_phasing", 	'f', "Randomly phasing all SVs (0/1-->0|1 or 1|0; 1/1-->1|1)", true, 1); l.set_arg_pointer_back((void *)&random_phasing);
		l.add_option("MIN_SV_len", 		'm', "MIN length of output SVs", true, 45); l.set_arg_pointer_back((void *)&MIN_SV_len);
		l.add_option("Small_var", 		'v', "Output small VARs around SVs(useless in MODE: ERR_PRONE)", true, 1); l.set_arg_pointer_back((void *)&output_small_var);
		l.add_option("output",  		'o', "the output file (.vcf)", true, default_outputFile); l.set_arg_pointer_back((void *)&outputFile);
		l.add_option("no-header", 		'H', "Not output VCF header"); l.set_arg_pointer_back((void *)&not_output_vcf_header);
		l.add_help_msg_back("\tLRS calling optional parameters:");
		l.add_help_msg_back("\t\t\"ONT_Q20\" and \"HIFI\" can be used when error rate is less than 1%.");
		l.add_help_msg_back("\t\t\"ERR_PRONE\" can be used when error rate is 1%~6%.");
		l.add_option("LRS_preset",   'p', "one of \"ONT_Q20\", \"ERR_PRONE\" or \"HIFI\".", true, "ONT_Q20"); l.set_arg_pointer_back((void *)&preset_platform_str);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tSRS calling optional parameters:");
		l.add_option("ref_idx",  		'I', "The reference index used for NGS-base calling(generate this file using 'ngs_fa_stat' command)", false, ""); l.set_arg_pointer_back((void *)&ref_idx_FileName);
		l.add_option("TL_BAM",    		'L', "The NGS trans-location BAM file(.bam)(generate this file using 'ngs_trans_reads' command)", false, ""); l.set_arg_pointer_back((void *)&TL_read_Filename);
		l.add_option("status",  		'T', "The NGS status file name [.json](generate this file using 'ngs_bam_stat' command)", false, ""); l.set_arg_pointer_back((void *)&statsFileName);
		l.add_option("r_BAM",   'R', "(Optional) The re-aligned NGS reads.", false, ""); l.set_arg_pointer_back((void *)&realigned_Filename);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tLRS+SRS hybrid optional parameters:");
		l.add_option("KMER_PROFILE", 'k', "The KMER ERROR rate profile file of LRS [used for SRS polishing LRS]", false, ""); l.set_arg_pointer_back((void *)&KMER_ERROR_PROFILE_File);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tContigs Force generating and Joint-Calling in [BED] regions");
		l.add_help_msg_back("\t\tOnly support LRS_BAM (-l) right now");
		l.add_help_msg_back("\t\tEach regions in BED file should be shorter then 50,000 bp");
		l.add_option("FC_BED", 'b', "[BED]The Contigs Force generation regions(In BED format)", false, ""); l.set_arg_pointer_back((void *)&FC_BED_F);
		//l.add_option("FC_Contig_L", 'J', "The Contigs File list to be joint-calling.", false, ""); l.set_arg_pointer_back((void *)&FC_Contig_L);

		l.add_help_msg_back("\tOthers");
		l.add_option("INNER_TEST", 'Y', "INNER_TEST", false, ""); l.set_arg_pointer_back((void *)&INNER_TEST_FLAG);

		force_calling_ALU = true;
		is_compression = false;
		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);
		if (argc - optind < 0)
			return l.output_usage();
		else if(referenceFilename == NULL){
			fprintf(stderr, "FATAL ERROR: Reference file (.fa/.fna/.fasta) is needed as input.\n");
			return l.output_usage();
		}
		else if(SRS_BAM_F == NULL && LRS_BAM_F == NULL && LRS_Tumor_BAM_F == NULL){
			fprintf(stderr, "FATAL ERROR: At least one BAM/CRAM (-l -n or -X) is needed as input.\n");
			return l.output_usage();
		}

		if(		strcmp(preset_platform_str, "ONT_Q20") == 0
				|| strcmp(preset_platform_str, "ont_Q20") == 0
		)
			preset_platform = PRESET_ONT_Q20;
		else if(strcmp(preset_platform_str, "HIFI") == 0
				|| strcmp(preset_platform_str, "CCS") == 0
				|| strcmp(preset_platform_str, "HiFi") == 0
				|| strcmp(preset_platform_str, "ccs") == 0
		)
			preset_platform = PRESET_CCS;
		else if(strcmp(preset_platform_str, "ASM") == 0
				|| strcmp(preset_platform_str, "asm") == 0
		)
			preset_platform = PRESET_ASM;
		else if(strcmp(preset_platform_str, "ERR_PRONE") == 0
				|| strcmp(preset_platform_str, "err_prone") == 0
		)
			preset_platform = PRESET_ERR;
		else	/*default*/
			preset_platform = PRESET_ONT_Q20;

		//NOT output small var when ERROR_PRONE
		if(preset_platform == PRESET_ERR)
			output_small_var = false;

		OUT_MODE_int = OUT_MODE_UNSET;
		if(strcmp(OUT_MODE_char, "PURE_STR") == 0){
			OUT_MODE_int = OUT_MODE_PURE_STR;
		}else if(strcmp(OUT_MODE_char, "VCF") == 0){
			OUT_MODE_int  = OUT_MODE_VCF;
		}
		return 0;
	}
};

class NOVA_SV_handler{

public:
	void run(NOVA_SV_PARA *nova_para){
		RefHandler refHandler;// = (RefHandler *)xcalloc(1,sizeof(RefHandler));
		refHandler.init(nova_para->st_chr_ID, nova_para->st_pos, nova_para->ed_chr_ID, nova_para->ed_pos,
				nova_para->referenceFilename);
		RefLocalRepeat *rlr = NULL;
		BAM_STATUS *bs = NULL;
		if(nova_para->SRS_BAM_F != NULL){
			rlr = (RefLocalRepeat *)xcalloc(1,sizeof(RefLocalRepeat));
			rlr->load(nova_para->ref_idx_FileName);
			bs = (BAM_STATUS *)xcalloc(1, sizeof(BAM_STATUS));
			if(nova_para->statsFileName != NULL){
				fexist_check(nova_para->statsFileName);
				bs->json_load(nova_para->statsFileName);
			}else
				bs->generate_state(nova_para->referenceFilename, nova_para->SRS_BAM_F, false);
		}
		KMER_ERROR_PROFILE_HANDLER kmer_profile_handler;
		if(nova_para->KMER_ERROR_PROFILE_File != NULL){
			kmer_profile_handler.load_data(nova_para->KMER_ERROR_PROFILE_File);
		}
		//SveHandler *sve_h = (SveHandler *)xcalloc(1, sizeof(SveHandler));
		SV_CALLING_Handler *sve_h = new SV_CALLING_Handler;
		{
			if(nova_para->FC_BED_F != NULL){
				sve_h->fc_lrs_handler.load_region_list(nova_para->FC_BED_F, refHandler.getIdx());
			}
		}
		sve_h->init(&refHandler, &kmer_profile_handler, rlr, nova_para->SRS_BAM_F, nova_para->realigned_Filename, nova_para->TL_read_Filename,
				nova_para->LRS_BAM_F, nova_para->LRS_Tumor_BAM_F, nova_para->LRS_Normal_BAM_F, nova_para->INNER_TEST_FLAG,
				bs, nova_para->outputFile, nova_para->is_compression, nova_para->MIN_SV_len,
				!nova_para->not_output_vcf_header, nova_para->force_calling_ALU,
				nova_para->random_phasing, nova_para->output_small_var, nova_para->preset_platform,
				nova_para->FC_BED_F, nova_para->OUT_MODE_int);
		//	//PART5: for each choromosome
		while(refHandler.load_seg_index())//get a new ref block in the reference, 2M bases per block
		{
			//build index for long SV detection
			//RefRegion * r = refHandler->get_cur_region();
			char log_title[1024]; sprintf(log_title, "[refHandler.build_ref_index]");
			double __cpu_time = cputime(); double __real_time = realtime();
			if(nova_para->OUT_MODE_int == OUT_MODE_VCF)//build index only when generating vcfs
				refHandler.build_ref_index();
			fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
			//get fc_sig for this region:
			sve_h->SV_calling_in_region_main_pipeline();//S1: in step1: we clear all signals in the SVE list
		}

		sve_h->distory();
	}
};
