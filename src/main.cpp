
#include<stdio.h>
#include<string.h>
#include <SVcalling_core/sv_call_main.hpp>
#include <time.h>
#include <vector>
#include <map>
#include "cpp_lib/Assembler/GreedyAssembler.hpp"
#include "cpp_lib/get_option_cpp.hpp"
#include "SVcalling_core/bam_stat.hpp"
#include "SVcalling_core/GT_cal.hpp"
#include "SVcalling_core/fa_stat.hpp"
#include "SVcalling_core/forceCalling.hpp"
#include "SVcalling_core/NGS_merging.hpp"

extern "C"
{
#include "clib/utils.h"
#include "clib/desc.h"
}

int analysis_main(int argc, char *argv[]);

int forceCalling_main(int argc, char *argv[]){
	FC_handler FC;
	FC.forceCalling_run(argc, argv);
	return 0;
}

int NGS_merging_main(int argc, char *argv[]){
	NGS_SV_MERGING_HANDLER SM;
	SM.SV_merging(argc, argv);
	return 0;
}

int SVLoci_main(int argc, char *argv[]){

	NOVA_SV_PARA * nova_para = (NOVA_SV_PARA *)xcalloc(1, sizeof(NOVA_SV_PARA));
	if(nova_para->get_option(argc, argv) != 0)
		return 0;
	NOVA_SV_handler h;
	h.run(nova_para);
	return 0;
}

int transLocation_read_main(int argc, char *argv[]){

	char *reference_fn =  argv[1];
	char *input_bam_fn =  argv[2];
	char *output_bam_fn =  argv[3];

	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file

	if (NULL != reference_fn)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, reference_fn);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(input_file, referenceFilenameIndex);
		xassert(ret == 0, "Failed to use reference for BAM/CRAM file");
	}

	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br_i = {0};
	while(sam_read1(input_file, header, &br_i) >=0)//read record{
	{
		bam1_t *br = &br_i;
		if(bam_is_secondary(br))		continue;
		if(bam_is_supplementary(br))	continue;
		if(bam_is_duplicate(br))       continue;
		if(br->core.tid == -1 || br->core.mtid == -1)  continue;
		if(br->core.mpos == br->core.pos)  continue;//skip unmapped reads

		if(br->core.isize == 0 || br->core.isize > 5000 || br->core.isize < -5000){//condition for TL reads
			std::swap(br->core.pos, br->core.mpos);
			std::swap(br->core.tid, br->core.mtid);
			xassert(sam_write1(output_file, header, br) >= 0, "");//write record
		}
	} 
	//close file
	hts_close(input_file);
	hts_close(output_file);
	return 0;
}

int get_fa_statistics(int argc, char *argv[]){
	FA_STATUS *fs = (FA_STATUS *)xcalloc(1, sizeof(FA_STATUS));
	return fs->run_ana(argc - 1, argv + 1);
}

int get_bam_statistics(int argc, char *argv[]){
	BAM_STATUS *bs = (BAM_STATUS *)xcalloc(1, sizeof(BAM_STATUS));
	bs->GEN_STAT_and_json_dump(argv[optind], argv[optind + 1]);

	return 0;
}

int main(int argc, char *argv[])
{

	//fprintf(stderr, "\n\n main version V2.01\n\n");
	COMMAND_HANDLER ch;
	ch.add_function("call", "GERMLINE or SOMAITC SV calling using (LRS or SRS or HYBRID(LRS+SRS))BAM file", SVLoci_main);
	ch.add_function("ngs_bam_stat", "Get SRS-BAM ISIZE distribution.", get_bam_statistics);
	ch.add_help_msg_back("Usage: srs_bam_stat reference.fa SRS.bam/cram > STAT.json");
	ch.add_function("srs_fa_stat", "Indexing the reference for SRS SV calling.", get_fa_statistics);
	ch.add_help_msg_back("Usage: srs_fa_stat reference.fa > ref.idx");
	ch.add_function("srs_trans_reads", "Extract trans_read from SRS BAM", transLocation_read_main);
	ch.add_help_msg_back("Usage: srs_trans_reads reference.fa input.bam output.bam");
	ch.add_function("ngs_merging", "Force calling tools for NGS dataset", NGS_merging_main);
	ch.add_function("ngs_fc", "Force calling tools for NGS dataset", forceCalling_main);
	//ch.add_function("fc", "Small tools", analysis_main);
	ch.add_function("tools", "Small tools", analysis_main);
	return ch.run(argc, argv);
}

