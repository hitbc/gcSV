/*
 * sv_main.cpp
 *
 *  Created on: 2020年4月27日
 *      Author: fenghe
 */

#include <getopt.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>

#include "../cpp_lib/get_option_cpp.hpp"
#include "../SVcalling_core/bam_stat.hpp"
#include "../SVcalling_core/SveHandler.hpp"

struct NOVA_SV_PARA{
	char*	referenceFilename;
	char*	realigned_Filename;
	char*	TL_read_Filename;
	char* 	outputFile;
	char* 	statsFileName;
	char* 	ref_idx_FileName;
	char*   NGS_BAM_File;
	bool 	is_compression;
	int 	MIN_SV_len;

	bool not_output_vcf_header;
	bool random_phasing;
	bool output_small_var;
	bool force_calling_ALU;

	int st_chr_ID; int st_pos;
	int ed_chr_ID; int ed_pos;

	//cuteSV-FC
	char*   CCS_BAM_F;
	char*   ONT_BAM_F;

	int get_option(int argc, char *argv[]){
		char default_outputFile[] = "stdout";
		options_list l;
		l.show_command(stderr, argc + 1, argv - 1);
	    l.add_title_string("\n");
	    l.add_title_string("  Usage:     ");  l.add_title_string(PACKAGE_NAME);  l.add_title_string("  call  [Options]\n");
	    l.add_title_string("\n\n");
	    l.add_title_string("     Necessary parameters:\n");
	    l.add_title_string("     \tAt least one BAM/CRAM is needed;\n");
	    l.add_title_string("     \tThe reference file is necessary.");
	    l.add_option("CCS_BAM", 'c', "The CCS BAM file", false, ""); l.set_arg_pointer_back((void *)&CCS_BAM_F);
	    l.add_option("ONT_BAM", 'O', "The ONT(At least Q26) BAM file", false, ""); l.set_arg_pointer_back((void *)&ONT_BAM_F);
		l.add_option("NGS_BAM", 'n', "The NGS BAM file", false, ""); l.set_arg_pointer_back((void *)&NGS_BAM_File);
		l.add_option("reference",  		'r', "the reference (.fa)", false, ""); l.set_arg_pointer_back((void *)&referenceFilename);
		l.add_help_msg_back("");
		l.add_help_msg_back("Optional parameters:");
		l.add_help_msg_back("\tRegion parameters:(SV calling only in a specific REGION)");
		l.add_help_msg_back("\t\tExample 1: using: '-S 0 -E 0 -s 100000 -F 200000' to call SVs in region: chr1:100000-200000");
		l.add_help_msg_back("\t\tExample 2: using: '-S 0 -E 23'(default) to call SVs in the whole-genome: from chr1 to chrY");
		l.add_option("st_chr_ID", 		'S', "The start chr_ID(0-base) for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_chr_ID);
		l.add_option("ed_chr_ID", 		'E', "The end chr_ID(0-base) for SV calling", true, 23); l.set_arg_pointer_back((void *)&ed_chr_ID);
		l.add_option("st_pos", 			's', "The start position for SV calling", true, 0); l.set_arg_pointer_back((void *)&st_pos);
		l.add_option("ed_pos", 			'F', "The end position for SV calling", true, 500000000); l.set_arg_pointer_back((void *)&ed_pos);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tOutput parameters:");
		l.add_option("random_phasing", 	'f', "Randomly phasing all SVs randomly (0/1-->0|1 or 1|0; 1/1-->1|1)"); l.set_arg_pointer_back((void *)&random_phasing);
		l.add_option("MIN_SV_len", 		'm', "MIN length of output SVs", true, 30); l.set_arg_pointer_back((void *)&MIN_SV_len);
		l.add_option("Small_var", 		'v', "Output small VARs around SVs [useful when applying haplotype-based benchmarking like vcfdist]", true, false); l.set_arg_pointer_back((void *)&output_small_var);
		l.add_option("output",  		'o', "the output file (.vcf)", true, default_outputFile); l.set_arg_pointer_back((void *)&outputFile);
		l.add_option("no-header", 		'H', "Not output VCF header"); l.set_arg_pointer_back((void *)&not_output_vcf_header);
		l.add_help_msg_back("");
		l.add_help_msg_back("\tNGS calling optional parameters:");
		l.add_option("ref_idx",  		'I', "The reference index used for NGS-base calling(generate this file using 'ngs_fa_stat' command)", false, ""); l.set_arg_pointer_back((void *)&ref_idx_FileName);
		l.add_option("TL_BAM",    		'L', "The NGS trans-location BAM file(.bam)(generate this file using 'ngs_trans_reads' command)", false, ""); l.set_arg_pointer_back((void *)&TL_read_Filename);
		l.add_option("status",  		'T', "The NGS status file name [.json](generate this file using 'ngs_bam_stat' command)", false, ""); l.set_arg_pointer_back((void *)&statsFileName);
		l.add_option("r_BAM",   'R', "(Optional) The re-aligned NGS reads.", false, ""); l.set_arg_pointer_back((void *)&realigned_Filename);
		l.add_help_msg_back("");

		force_calling_ALU = true;
		is_compression = false;
		if(l.default_option_handler(argc, argv)) {
			l.show_c_value(stderr);
			return 1;
		}
		l.show_c_value(stderr);
		if (argc - optind < 0)
			return l.output_usage();
		else if(CCS_BAM_F == NULL && NGS_BAM_File == NULL && ONT_BAM_F == NULL){
			fprintf(stderr, "At least one BAM/CRAM is needed as input.\n");
			return l.output_usage();
		}else if(CCS_BAM_F != NULL && ONT_BAM_F != NULL){
			fprintf(stderr, "At most one TGS (CCS or ONT) BAM/CRAM can be used as input.\n");
				return l.output_usage();
		}
		return 0;
	}
};

class NOVA_SV_handler{

public:
	void run(NOVA_SV_PARA *nova_para){
		RefHandler *refHandler = (RefHandler *)xcalloc(1,sizeof(RefHandler));
		refHandler->init(nova_para->st_chr_ID, nova_para->st_pos, nova_para->ed_chr_ID, nova_para->ed_pos,
				nova_para->referenceFilename);
		RefLocalRepeat *rlr = NULL;
		BAM_STATUS *bs = NULL;
		if(nova_para->NGS_BAM_File != NULL){
			rlr = (RefLocalRepeat *)xcalloc(1,sizeof(RefLocalRepeat));
			rlr->load(nova_para->ref_idx_FileName);
			bs = (BAM_STATUS *)xcalloc(1, sizeof(BAM_STATUS));
			if(nova_para->statsFileName != NULL){
				fexist_check(nova_para->statsFileName);
				bs->json_load(nova_para->statsFileName);
			}else
				bs->generate_state(nova_para->referenceFilename, nova_para->NGS_BAM_File, false);
		}

		SveHandler *sve_h = (SveHandler *)xcalloc(1, sizeof(SveHandler));
		sve_h->init(refHandler, rlr, nova_para->NGS_BAM_File, nova_para->realigned_Filename, nova_para->TL_read_Filename,
				nova_para->CCS_BAM_F, nova_para->ONT_BAM_F,
				bs, nova_para->outputFile, nova_para->is_compression, nova_para->MIN_SV_len,
				!nova_para->not_output_vcf_header, nova_para->force_calling_ALU, nova_para->random_phasing, nova_para->output_small_var);
		//	//PART5: for each choromosome
		while(refHandler->load_seg_index())//get a new ref block in the reference, 2M bases per block
		{
			//get fc_sig for this region:
			sve_h->SVE_handle_region();//S1: in step1: we clear all signals in the SVE list
		}

		sve_h->distory();
	}
};
