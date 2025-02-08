#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <math.h>
#include "cpp_lib/RefRegion.hpp"
#include "cpp_lib/get_option_cpp.hpp"
#include "cpp_lib/cpp_utils.hpp"
#include "NGS_merging.hpp"

extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/vcf_lib.h"
}
//#----------------------------------------------------------
typedef struct
{
	char ID[64];
	char type[32];
	int rid;
	int sample;
	int st;
	int ed;
}VCF_ANALYSIS;
kvec_T(VCF_ANALYSIS, VCF_A_L);

#define INSERT_MAX 100000

int isize_count(int argc, char *argv[])
{
	char* bam_file_name = argv[1];
	//open bam file
	FILE* try_open = xopen(bam_file_name, "r");
	fclose(try_open);
	Bam_file bf;
	bam_file_open(bam_file_name, NULL, NULL, &bf);

	int insert_size[INSERT_MAX];
	for(int i = 0; i < INSERT_MAX; i++)
		insert_size[i] = 0;
	//count
	int read_number = 0;
	int high_mapq_read = 0;
	int overflow = 0;
	while(bam_next(&bf))
	{
		read_number++;
		if(read_number % 100000 == 0)
			fprintf(stderr, "loading read_number:%d\r", read_number);//-----------------------
		bam1_t *br = &(bf._brec);
		if(br->core.qual < 15)
			continue;
		high_mapq_read++;
		int i_size = ABS(br->core.isize);
		if(i_size < INSERT_MAX)
			insert_size[i_size] ++;
		else
			overflow++;
	}

	printf(
			"Read number: %d\n"
			"high_mapq_read: %d\n"
			"insert size 0: %d\n"
			"overflow: %d\n"
			"\n\n",
			read_number,
			high_mapq_read,
			insert_size[0],
			overflow
			);

	int sum = 0;
	for(int i = 0; i < INSERT_MAX; i++)
	{
		sum += insert_size[i];
		if(insert_size[i] != 0)
			printf(
				"insert size: %d\t"
				"Number: %d\t"
				"Percent: %f\t"
				"Sum percent %f\n",
				i,
				insert_size[i],
				(float)insert_size[i]/high_mapq_read,
				(float)sum/high_mapq_read
				);
	}
	return 0;
}

static void analysis_usage()
{
	fprintf(stderr, "analysis_usage\n\n");
}

static int ref_split(int argc, char *argv[])
{
	char *fasta_file = argv[1];
	gzFile fp = xzopen(fasta_file, "r");
	kstream_t *_fp = ks_init(fp);

	kseq_t temp = {0};
	temp.f = _fp;
	while( kseq_read(&temp) >= 0)
	{
		char out_file_name[1024];
		strcpy(out_file_name, temp.name.s);
		strcat(out_file_name, ".fa");
		FILE *out_file = xopen(out_file_name, "w");
		printf(">%s %s\n", temp.name.s, temp.comment.s);
		fprintf(out_file, ">%s %s\n", temp.name.s, temp.comment.s);
		//print content
		for(unsigned int i = 0; i < temp.seq.l; i++)
		{
			fprintf(out_file, "%c", temp.seq.s[i]);
			if(i %70 == 69)
				fprintf(out_file, "\n");
		}
		fclose(out_file);
	}
	return 0;
}

void ref_dump(char *fasta_file, int maxChr)
{
	gzFile fp = xzopen(fasta_file, "r");
	kstream_t *_fp = ks_init(fp);

	kseq_t temp = {0};
	temp.f = _fp;
	int chr = 0;
	while( kseq_read(&temp) >= 0 && chr++ < maxChr)
	{
		char out_file_name[1024];
		strcpy(out_file_name, temp.name.s);
		strcat(out_file_name, ".fna");
		FILE *out_file = xopen(out_file_name, "w");
		printf(">%s %s\n", temp.name.s, temp.comment.s);
		fprintf(out_file, ">%s %s\n", temp.name.s, temp.comment.s);
		//print content
		for(unsigned int i = 0; i < temp.seq.l; i++)
		{
			fprintf(out_file, "%c", temp.seq.s[i]);
			if(i %70 == 69)
				fprintf(out_file, "\n");
		}
		fclose(out_file);
	}
}

//todo
#define MAX_LINE_LENGTH 1000
void getSVRegion(char * region_fn, std::vector<R_region>&regionList, std::vector<std::string> &id2Name){
	//get sv region:
	char *temp = new char[MAX_LINE_LENGTH];//1M
	std::map<std::string, uint32_t> name2ID;

	std::ifstream regionListFile(region_fn);
	while(true){
		regionListFile.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)
			break;
		//get chrID
		char *token = strtok(temp, "\t");
		std::string strID(token); int intID = 0;
		std::map<std::string, uint32_t>::iterator it = name2ID.find(strID);
		if(it!=name2ID.end()){
			intID = it->second;
		}
		else{
			intID = name2ID.size();
			name2ID[strID] = name2ID.size();
			id2Name.emplace_back(strID);
		}

		//get refSt
		token = strtok(NULL, "\t");
		uint32_t regionSt = strtoul(token, NULL, 10);
		//get refEd
		token = strtok(NULL, "\t");
		uint32_t regionEnd = strtoul(token, NULL, 10);
		//
		regionList.emplace_back();
		R_region & r = regionList.back();
		r.chr_ID = intID; r.st_pos = regionSt; r.ed_pos = regionEnd;
	}
	regionListFile.close();
	delete [] temp;
}

void dump_ref_by_region(char *fasta_file, char * region_fn)//todo:::::
{
	//load regions
	std::vector<R_region> regionList;
	std::vector<std::string> id2Name;
	getSVRegion(region_fn, regionList, id2Name);

	//open fasta file
//	gzFile fp = xzopen(fasta_file, "r");
//	kstream_t *_fp = ks_init(fp);

	//open fai file
	faidx_t * fai = fai_load(fasta_file);

	//load regions
	for(R_region &r :regionList){
		char reg[1024]; int load_len = 0;
		sprintf(reg, "%s:%d-%d", id2Name[r.chr_ID].c_str(), r.st_pos, r.ed_pos);
		//char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);
		char *ref = fai_fetch(fai, reg, &load_len);

		printf(">%s_%d_%d len_%d_%d_%d\n",  id2Name[r.chr_ID].c_str(), r.st_pos, r.ed_pos, load_len, load_len/70, load_len % 70);
		for(int i = 0; i < load_len; i++){
			printf("%c", ref[i]);
			if(i %70 == 69)
				printf("\n");
		}
		if((load_len - 1) % 70 != 69)
			printf("\n");
		free(ref);
	}

	fai_destroy(fai);
}

void get_all_record_with_SH(char *input_bam_fn, char * output_bam_fn)
{
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br;
	while(sam_read1(input_file, header, &br) >=0)//read record
	{
		int soft_left, soft_right;
		if(bam_has_SH_cigar(&br, &soft_left, &soft_right))
			xassert(sam_write1(output_file, header, &br) >= 0, "");//write record
	}
	//close file
	hts_close(input_file);
	hts_close(output_file);
}

void get_all_record_with_DR(char *input_bam_fn, char * output_bam_fn)
{
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br;
	while(sam_read1(input_file, header, &br) >=0)//read record
	{
		int min = 100, max = 1200;
		if(bam_is_DR_signal(&br, min, max))
			xassert(sam_write1(output_file, header, &br) >= 0, "");//write record
	}
	//close file
	hts_close(input_file);
	hts_close(output_file);
}

void get_all_record_with_SA(char *input_bam_fn, char * output_bam_fn)
{
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br;
	while(sam_read1(input_file, header, &br) >=0)//read record
	{
		if(bam_get_string_tag(&br, "SA") != NULL)
			xassert(sam_write1(output_file, header, &br) >= 0, "");//write record
	}
	//close file
	hts_close(input_file);
	hts_close(output_file);
}

void get_all_record_with_SIGNAL(char *input_bam_fn, char * output_bam_fn)
{
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br;
	while(sam_read1(input_file, header, &br) >=0)//read record
	{
		int min = 100, max = 1200;
		if(bam_is_DR_signal(&br, min, max))
			xassert(sam_write1(output_file, header, &br) >= 0, "");//write record
	}
	//close file
	hts_close(input_file);
	hts_close(output_file);
}

static int bam2Fastq(int argc, char *argv[]){
	char *input_bam_fn = argv[1];
	char * output_fastq_fn = argv[2];
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	FILE *output_file = xopen(output_fastq_fn, "w");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	bam1_t br = {0};
	char * seq = new char [10000];
	uint8_t * qual = new uint8_t [10000];
	while(sam_read1(input_file, header, &br) >=0){//read record
		const int read_len = br.core.l_qseq;
		get_bam_seq(0, read_len, seq, &br);//store in binary format
		get_bam_quality_str(0, read_len, qual, &br);
		fprintf(output_file, ""
				"@%s\n"
				"%s\n"
				"+%s\n"
				"%s\n",
				bam_qname(&br), seq, bam_qname(&br), qual);
	}
	fprintf(output_file, "\n");
	//close file
	hts_close(input_file);
	fclose(output_file);
	delete [] seq;
	delete [] qual;
	return 0;
}

static int bamDump(int argc, char *argv[]){

	char *input_bam_fn =  argv[1];
	char * reference_fn =  argv[2];
	char * output_bam_fn =  argv[3];
	int readN = atoi(argv[4]);

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
	bam1_t br = {0}; int readNum = 0;
	while(sam_read1(input_file, header, &br) >=0 && readNum++ < readN)//read record{
	{
		xassert(sam_write1(output_file, header, &br) >= 0, "");//write record
	}
	//close file
	hts_close(input_file);
	hts_close(output_file);
	return 0;
}

void bamDump_discard_both_unmapped(char *input_bam_fn, char * output_bam_fn){

	fprintf(stderr, "version 1.07 \n");

	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br1 = {0}; bam1_t br2 = {0};
	int readNum = 0; int load_read_number = 0;
	char readName[1024]; int sam_rst1 = 0; int sam_rst2 = 0;
	while(1){//read record
		//load SAM 1 & 2
		do{	sam_rst1 = sam_read1(input_file, header, &br1);	} while( (sam_rst1 >= 0) && (bam_is_secondary(&br1) || bam_is_supplementary(&br1) || bam_is_duplicate(&br1)));
		do{	sam_rst2 = sam_read1(input_file, header, &br2);	} while( (sam_rst2 >= 0) && (bam_is_secondary(&br2) || bam_is_supplementary(&br2) || bam_is_duplicate(&br2)));
		if(sam_rst1 < 0 || sam_rst2 < 0)		break;
		//SAME read name
		bool namesAreSame = (0 == strcmp((char *)bam_qname(&br1), (char *)bam_qname(&br2))); xassert(namesAreSame, "");
		//analysis read name
		strcpy(readName, (char *)bam_qname(&br1));

		// get basic information for that read:
		//HISEQ1:93:H2YHMBCXX:1:1203:3107:29946__15_46406331__15_46406262_MQ1_56_MQ2_56_UM1M_UM2M_ISIZE_296_FR_IM1_47_IM2_2_CL1_0_CL2_17_
		//get read name:
		char * token = strtok(readName, "_");
		token = strtok(NULL, "_"); //int chr_ID_1 = atoi(token);
		token = strtok(NULL, "_"); //int chr_pos1 = atoi(token);
		token = strtok(NULL, "_"); //int chr_ID_2 = atoi(token);
		token = strtok(NULL, "_"); //int chr_pos2 = atoi(token);
		token = strtok(NULL, "_"); //nothing
		token = strtok(NULL, "_"); //int mapq1 = atoi(token);
		token = strtok(NULL, "_"); //nothing
		token = strtok(NULL, "_"); //int mapq2 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_");
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); //int isize = atoi(token);
		token = strtok(NULL, "_"); //char direction1 = token[0]; char direction2 = token[1];
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); //int SNP1 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); //int SNP2 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); //int CL1 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); //int CL2 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); int score1 = atoi(token);
		token = strtok(NULL, "_");
		token = strtok(NULL, "_"); int score2 = atoi(token);

		int c_score1 = getScoreByCigar_BWA_MEM_LIKE(&br1);
		int c_score2 = getScoreByCigar_BWA_MEM_LIKE(&br2);

		bool score_is_bigger = false;
		if(score1 + score2 < c_score1 + c_score2 + 3)
			score_is_bigger = true;

		load_read_number += 2;
		if(load_read_number % 10000 == 0){
			fprintf(stderr, "signal read number : %d\t", load_read_number);
			fprintf(stderr, "signal read number : %d\r", readNum);
		}

		if(score_is_bigger && !(bam_is_unmapped(&br1) && bam_is_unmapped(&br2))){
			readNum +=2;
			if(readNum % 10000 == 0)
				fprintf(stderr, "signal read number : %d\r", readNum);

			xassert(sam_write1(output_file, header, &br1) >= 0, "");//write record
			xassert(sam_write1(output_file, header, &br2) >= 0, "");//write record
		}
	}
	fprintf(stderr, "total signal read number : %d\n", readNum);
	//close file
	hts_close(input_file);
	hts_close(output_file);
}

void bamDump_discard_map_len_100(char *input_bam_fn, char * output_bam_fn, int MIN_len){

	fprintf(stderr, "version 1.01 \n");

	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	htsFile *output_file = hts_open(output_bam_fn, "wb");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	xassert(sam_hdr_write(output_file, header) >=0, "write header wrong!");//write header
	bam1_t br1 = {0};
	bam1_t br2 = {0};

	int totalLenCount[600] = {0};

	int readNum = 0;
	while(sam_read1(input_file, header, &br1) >=0 && sam_read1(input_file, header, &br2) >=0){//read record
		bool namesAreSame = (0 == strcmp((char *)bam_qname(&br1), (char *)bam_qname(&br2)));
		xassert(namesAreSame, "");
		int soft_left1; int soft_right1, read_len1 = 0;
		if(!bam_is_unmapped(&br1)){
			bam_has_SH_cigar(&br1, &soft_left1, &soft_right1);
			read_len1 = br1.core.l_qseq - soft_left1 - soft_right1;
		}

		int soft_left2; int soft_right2, read_len2 = 0;
		if(!bam_is_unmapped(&br2)){
			bam_has_SH_cigar(&br2, &soft_left2, &soft_right2);
			read_len2 =  br2.core.l_qseq - soft_left2 - soft_right2;
		}

		int total_len = read_len1 + read_len2;
		totalLenCount[total_len] ++;
		if(total_len < MIN_len) continue;

		 readNum ++;
		 xassert(sam_write1(output_file, header, &br1) >= 0, "");//write record
		 xassert(sam_write1(output_file, header, &br2) >= 0, "");//write record
	}
	fprintf(stderr, "total signal read number : %d\n", readNum);
	for(int i = 0; i < 600; i ++){
		fprintf(stderr, "len : %d count: %d \n", i, totalLenCount[i]);
	}

	//close file
	hts_close(input_file);
	hts_close(output_file);
}

void SP_region_result_full_test(char *input_bam_fn){

	fprintf(stderr, "version 1.03 \n");

	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file

	bam_hdr_t *header = sam_hdr_read(input_file);//read header

	bam1_t br = {0};

	int readNum = 0;
	int c_tid = -1;
	uint32_t regionSt = 0, regionEnd = 0;
	uint32_t chrID_ori = 0;
	int region_length = 1;
	int c_read_number = 0;
	uint64_t region_total_length = 0;
	bool record_finish = false;
	bool regionFinish = false;
	int newMapped = 0;
	int regionID = 1;

	int totalLenCount[600] = {0};

	while(1){//read record

		if(sam_read1(input_file, header, &br) < 0)
			record_finish = true;
		regionFinish = (c_tid != br.core.tid)?true:false;
		if((record_finish || regionFinish) && (c_tid != -1) ){
			//read depth:
			float readdpth = (float)region_total_length/ region_length;
			fprintf(stderr, "(%d) %d:%d-%d\t readDepth\t%f\tnewMapped\t%f\n", regionID++, chrID_ori, regionSt, regionEnd, readdpth, (float)newMapped/c_read_number);
		}

		if(!record_finish && regionFinish){
			c_tid = br.core.tid;
			//2_91600041_91601285
			//get the region length:
			char * region_name = header->target_name[c_tid];
			char * token = strtok(region_name, "_");
			chrID_ori = strtoul(token, NULL, 10);
			//get refSt
			token = strtok(NULL, "_");
			regionSt = strtoul(token, NULL, 10);
			//get refEd
			token = strtok(NULL, "\0");
			regionEnd = strtoul(token, NULL, 10);
			region_length = regionEnd - regionSt;
			c_read_number = 0;
			region_total_length = 0;
			newMapped = 0;
		}

		if(record_finish) break;

		readNum++;

		// get basic information for that read:
		//HISEQ1:93:H2YHMBCXX:1:1203:3107:29946__15_46406331__15_46406262_MQ1_56_MQ2_56_UM1M_UM2M_ISIZE_296_FR_IM1_47_IM2_2_CL1_0_CL2_17_
		//get read name:

		char * readName = (char *)bam_qname(&br);
		char * token = strtok(readName, "_");
		token = strtok(NULL, "_");// int chr_ID_1 = atoi(token);
		token = strtok(NULL, "_");// int chr_pos1 = atoi(token);
		token = strtok(NULL, "_");// int chr_ID_2 = atoi(token);
		token = strtok(NULL, "_");// int chr_pos2 = atoi(token);
		token = strtok(NULL, "_"); //nothing
		token = strtok(NULL, "_");// int mapq1 = atoi(token);
		token = strtok(NULL, "_"); //nothing
		token = strtok(NULL, "_");// int mapq2 = atoi(token);
		token = strtok(NULL, "_"); bool unmapped1 = (token[3] == 'U');
		token = strtok(NULL, "_"); bool unmapped2 = (token[3] == 'U');
//		token = strtok(NULL, "_");
//		token = strtok(NULL, "_"); int isize = atoi(token);
//		token = strtok(NULL, "_"); char direction1 = token[0]; char direction2 = token[1];
//		token = strtok(NULL, "_");
//		token = strtok(NULL, "_"); int SNP1 = atoi(token);
//		token = strtok(NULL, "_");
//		token = strtok(NULL, "_"); int SNP2 = atoi(token);
//		token = strtok(NULL, "_");
//		token = strtok(NULL, "_"); int CL1 = atoi(token);
//		token = strtok(NULL, "_");
//		token = strtok(NULL, "_"); int CL2 = atoi(token);

		//current stats:
		bool isFirstRead = bam_is_first(&br);
		bool oriUM = (isFirstRead)?unmapped1:unmapped2;

		bool UM = bam_is_unmapped(&br);
		//clip:
		int soft_left; int soft_right;
		bam_has_SH_cigar(&br, &soft_left, &soft_right);
		int	total_clip_len = (soft_left + soft_right);

		//map length
		int readMapLen = br.core.l_qseq - total_clip_len;
		readMapLen = MAX(0, readMapLen);
		region_total_length += readMapLen;
		xassert(readMapLen >= 0 && readMapLen < 300, "");
		totalLenCount[readMapLen] ++;

		//analysis
		c_read_number++;
		if(!UM && (oriUM)){
			newMapped++;
		}
	}

	for(int i = 0; i < 600; i ++){
		fprintf(stderr, "len : %d count: %d \n", i, totalLenCount[i]);
	}
	//close file
	hts_close(input_file);
}


void bam2Info(char *input_bam_fn, char * output_fn){
	htsFile *input_file = hts_open(input_bam_fn, "rb");//open input file
	FILE *output_file = xopen(output_fn, "w");//open output file
	bam_hdr_t *header = sam_hdr_read(input_file);//read header
	bam1_t br = {0};
	char * seq = new char [10000];
	uint8_t * qual = new uint8_t [10000];
	while(sam_read1(input_file, header, &br) >=0){//read record
		const int read_len = br.core.l_qseq;
		get_bam_seq(0, read_len, seq, &br);//store in binary format
		get_bam_quality_str(0, read_len, qual, &br);
		fprintf(output_file, ""
				"@%s\n"
				"%s\n"
				"+%s\n"
				"%s\n",
				bam_qname(&br), seq, bam_qname(&br), qual);
	}
	fprintf(output_file, "\n");
	//close file
	hts_close(input_file);
	fclose(output_file);
	delete [] seq;
	delete [] qual;
}

void vcf_sample(char *fn_in, char * fn_out, char * sample_name, char * sv_type, char * CHROM_ID)
{
	bool is_compression = false;//BCF or VCF
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read
	BCF_FILE vcf_w;//vcf for write
	VCF_open_write(&vcf_w, fn_out, is_compression);
	bcf_hdr_write(vcf_w.file, vcf_r.header);
	char *c_sample = (char *)malloc(1000);
	char *c_sv_type = (char *)malloc(1000);
	int all_sample = false;
	int all_type = false;
	int all_chrom = false;
	if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)
		all_sample = true;
	if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)
		all_type = true;
	if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)
		all_chrom = true;
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
		if(all_sample == false)
		{
			vcf_get_sample(vcf_r.header, c_r, c_sample);
			if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
				continue;
		}
		if(all_type == false)
		{
			vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
			if(strcmp(c_sv_type, sv_type) != 0)
				continue;
		}
		if(all_chrom == false)
		{
			if(c_r->rid != strtol(CHROM_ID,0,10))
				continue;
		}
		if(c_r->rid > 24)//discard decoy sequence
			continue;
		bcf_write(vcf_w.file, vcf_r.header, c_r);
		//vcf_get_genotype(vcf_r.header, c_r);
	}
	//close
	bcf_close(vcf_r.file);
	bcf_close(vcf_w.file);
}

int vcf_dump(int argc, char *argv[]){
	char *fn_in = argv[1];
	char * fn_out = argv[2];
	char * sample_name = argv[3];
	char * sv_type = argv[4];
	char * CHROM_ID = argv[5];
	vcf_sample(fn_in, fn_out, sample_name, sv_type, CHROM_ID);
	return 0;
}

struct SV_INFO_vcf_add_alt_string{
	std::vector<std::string> info;
	int sample_ID;
	int SV_len;
	int chrID;
	int pos;
	int sample_numer;//todo::
	int max_same_sv_len_sample_numer;//todo::

	std::string ALT_str;
	std::string REF_str;
	int true_sample_pos;

	void store_basic(std::vector<std::string> &item_value, faidx_t *c_ref_idx){
		for(int j = 0; j < 9; j ++)
			info.emplace_back(item_value[j]);
		chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
		pos = atoi(item_value[1].c_str());
	}
};

struct SAMPLE_INFO_vcf_add_alt_string{
	std::string sample_name;
	std::vector<std::vector<int>> SV_IDs;
};

int vcf_add_alt_string(int argc, char *argv[]){
	//char *joint_vcf_fn = argv[1];
	faidx_t *c_ref_idx = reference_index_load("/media/fenghe/Data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa");
	const char *joint_vcf_fn = "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/filltag_S5.vcf";
	std::vector<SV_INFO_vcf_add_alt_string> SV_info;
	std::vector<SAMPLE_INFO_vcf_add_alt_string> sample_info;
	//load joint data
	{
		std::vector<std::string> joint_vcf_data;
		load_string_list_from_file_MAX_line(joint_vcf_fn, joint_vcf_data, 5000);
		size_t line_num = joint_vcf_data.size();//skip the header line
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_2;
		std::map <int, int> SVLEN_counter;
		for(uint64_t i = 0; i < line_num; i++){
			//show header
			const char * line_data = joint_vcf_data[i].c_str();
			if(line_data[0] == '#' && line_data[1] == '#'){
				printf("%s\n", joint_vcf_data[i].c_str());
			}
			else if(line_data[0] == '#'){
				split_string(item_value, joint_vcf_data[i].c_str(), "\t");//skip the header line
				for(int j = 0; j < 9; j ++)
					printf("%s\t", item_value[j].c_str());
				printf("SAMPLE\n");
				//store sample name
				for(uint j = 9; j < item_value.size(); j ++){
					sample_info.emplace_back();
					sample_info.back().sample_name = item_value[j];
					sample_info.back().SV_IDs.resize(24);
				}
			}
			else{
				split_string(item_value, joint_vcf_data[i].c_str(), "\t");//skip the header line
				SV_info.emplace_back();
				SV_info.back().store_basic(item_value,c_ref_idx);
				//get the sample with max same SV length number
				SVLEN_counter.clear();
				for(uint j = 9; j < item_value.size(); j++){
					split_string(item_value_2, item_value[j].c_str(), ":");//skip the header line
					if(strcmp(item_value_2[0].c_str(), "./.") != 0){
						int SVLEN = atoi(item_value_2[1].c_str());
						std::map<int, int >::iterator it = SVLEN_counter.find(SVLEN);
						if(it != SVLEN_counter.end()){
							it->second++;
						}else
							SVLEN_counter[SVLEN] = 1;
					}
				}
				int max_SV_len = -1;
				int max_count = -1;
				for(std::map<int, int >::iterator it = SVLEN_counter.begin(); it != SVLEN_counter.end(); it++){
					if(max_count < it->second){
						max_count = it->second;
						max_SV_len = it->first;
					}
				}
				SV_info.back().SV_len = max_SV_len;
				for(uint j = 9; j < item_value.size(); j ++){
					split_string(item_value_2, item_value[j].c_str(), ":");//skip the header line
					if(strcmp(item_value_2[0].c_str(), "./.") != 0){
						int SVLEN = atoi(item_value_2[1].c_str());
						if(SVLEN == max_SV_len){
							SV_info.back().sample_ID = j - 9;
							sample_info[SV_info.back().sample_ID].SV_IDs[SV_info.back().chrID].emplace_back(SV_info.size() - 1);
							break;
						}
					}
				}
			}
		}
	}

	//recover from single sample vcf data:
	{
		//load single sample vcf
		for(uint i = 0; i < sample_info.size(); i++){
			//if(i > 30) break;
			char vcf_fn[1024];
			sprintf(vcf_fn, "/media/fenghe/MyPassport_4T/bamdata/1KGP_gcSV/sorted_vcf/%s_sv_spa.sort.vcf.gz", sample_info[i].sample_name.c_str());
			fprintf(stderr, "Current handle %s\n", vcf_fn);
			BCF_FILE input_vcf;

			VCF_open_read(&input_vcf, vcf_fn);
			bcf_hdr_t *vcf_header = input_vcf.header;
			int recover_SV_num = 0;
			//load each SVs
			do
			{
				bcf1_t *c_r = &( input_vcf.r);
				bcf_unpack(c_r, BCF_UN_INFO);
				int sv_chr_ID = c_r->rid;
				//vcf_read.header->
				//faidx_iseq(c_ref_idx, SV_CHR_ID);
				int chrID = faidx_get_chrID(c_ref_idx, bcf_hdr_id2name(vcf_header, sv_chr_ID), NULL, 0);
				if(chrID >= 24)
					continue;
				int SV_POS = c_r->pos;
				int SV_length = 0;
				vcf_get_sv_LENGTH(input_vcf.header, c_r, &SV_length);
				//search
				std::vector<int> &SV_list = sample_info[i].SV_IDs[chrID];
				//search if this SV is needed
				for(int SV_ID: SV_list){
					if(ABS_U(SV_info[SV_ID].pos,SV_POS) < 1000){
						if(ABS(SV_length) == ABS(SV_info[SV_ID].SV_len)){
							//store the ALT string; //store the REF string
							SV_info[SV_ID].REF_str.append(c_r->d.allele[0]);
							SV_info[SV_ID].ALT_str.append(c_r->d.allele[1]);
							SV_info[SV_ID].true_sample_pos = SV_POS;
							recover_SV_num++;
						}
					}
				}
			}while(VCF_next(&input_vcf));//read one
			fprintf(stderr, "\n recover_SV_num %d\n", recover_SV_num);
			bcf_close(input_vcf.file);
		}

		//output the final result
		{
			//std::vector<SV_INFO_vcf_add_alt_string> SV_info;
			int SV_number = SV_info.size();
			for(int i = 0; i < SV_number; i++){
				SV_INFO_vcf_add_alt_string & sv = SV_info[i];
				std::vector<std::string> & info = sv.info;
				bool STR_reset = (!sv.ALT_str.empty());
				int store_pos = (STR_reset)?(sv.true_sample_pos + 1):(sv.pos);
				printf("%s\t%d\t%s\t", info[0].c_str(), store_pos, info[2].c_str());//chrID, pos, name
				if(STR_reset){
					printf("%s\t%s\t", sv.REF_str.c_str(), sv.ALT_str.c_str());//chrID
				}else
					printf("%s\t%s\t", info[3].c_str(), info[4].c_str());//chrID
				// sample_info[sv.sample_ID].sample_name.c_str()
				printf("%s\t%s\t%s\tGT\t0/1\n", info[5].c_str(), info[6].c_str(), info[7].c_str() );//chrID, pos, name
			}
		}

	}
	fprintf(stderr, "\nEND\n");
	return 0;
}

struct BED_ITEM{
	int chrID;
	int pos;
	int end;

	std::string info1;
	std::string info2;

	int repeat_len;
	std::string repeat_str;
	//for VNTR count:
	int min_sv_len;
	int max_sv_len;

	int repeat_range;

	int bed_original_index = 0;
	int bed_repeat_range_index = 0;

	void load_data(std::vector<std::string> &item_value, faidx_t *c_ref_idx, std::vector<std::string> &item_value_buff){
		chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
		pos = atoi(item_value[1].c_str());
		end = atoi(item_value[2].c_str());
		info1 = item_value[3];
		info2 = item_value[4];
		split_string(item_value_buff, info2.c_str(), ")");
		repeat_len = item_value_buff[0].size() - 1;
		repeat_str.clear();
		repeat_str.append(item_value_buff[0].c_str() + 1);
	}

	void init_SV_min_max(){
		min_sv_len = MAX_int32t;
		max_sv_len = -MAX_int32t;
		repeat_range = -1;
	}

	void store_SV_min_max(int SV_len){
		min_sv_len = MIN(min_sv_len, SV_len);
		max_sv_len = MAX(max_sv_len, SV_len);
	}

	void store_repeat_range(){
		if(min_sv_len == MAX_int32t)
			repeat_range = -1;
		else
			repeat_range = (max_sv_len - min_sv_len) / repeat_len;
	}
	void show(){
		fprintf(stderr, "BED1: %d\t%d\t%d\t%s\t%s\t%d\t", chrID, pos, end, info1.c_str(), info2.c_str(), repeat_len);
		fprintf(stderr, "BED2: %d\t%d\t%d\t%d\t%d\n", min_sv_len, max_sv_len, repeat_range, bed_original_index, bed_repeat_range_index);
	}

	static inline int cmp_by_repeat_range(const BED_ITEM &a, const BED_ITEM &b){
		if(a.repeat_range != b.repeat_range)
			return a.repeat_range < b.repeat_range;
		else if(a.chrID != b.chrID)
			return a.chrID < b.chrID;
		else
			return a.pos < b.pos;
	}

	static inline int cmp_by_bed_original_index(const BED_ITEM &a, const BED_ITEM &b){
		if(a.bed_original_index != b.bed_original_index)
			return a.bed_original_index < b.bed_original_index;
		else xassert(0, "");
	}
};

#define BED_IDX_BIN_SIZE 1000

struct BED_IDX{
	struct BLOCK_IDX{
		std::vector<int> overlap_bed_id;
	};

	struct CHR_IDX{
		int chrID;
		int chr_len;
		std::vector<BLOCK_IDX> pos_idx;
		void initBlock(int block_size){
			int total_block_length = chr_len/block_size + 1;
			pos_idx.resize(total_block_length);
		}
	};

	std::vector<CHR_IDX> chr_idx;
	void init(faidx_t *c_ref_idx ){
		//load reference index
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int i = 0; i < N_seq; i++){
			chr_idx.emplace_back();
			chr_idx.back().chrID = i;
			chr_idx.back().chr_len = faidx_seq_len(c_ref_idx, faidx_iseq(c_ref_idx, i));
			chr_idx.back().initBlock(BED_IDX_BIN_SIZE);
		}
	}
	void store_bed_item(int bed_id, BED_ITEM & b){
		int bg = b.pos/BED_IDX_BIN_SIZE;
		int ed = b.end/BED_IDX_BIN_SIZE;
		for(int i = bg; i <= ed; i++){
			chr_idx[b.chrID].pos_idx[i].overlap_bed_id.emplace_back(bed_id);
		}
	}

	int search_bed_full_cover_SV(int sv_chr, int sv_pos, int sv_end, std::vector<BED_ITEM> & bed_list){
		int bg = sv_pos/BED_IDX_BIN_SIZE;
		for(int bed_i: chr_idx[sv_chr].pos_idx[bg].overlap_bed_id){
			if(sv_pos >= bed_list[bed_i].pos && sv_end <= bed_list[bed_i].end)
				return bed_i;
		}
		return -1;
	}
};

struct SV_in_repeat{
	int SV_id;
	int repeatID;
	int SV_len;
	static inline int cmp_by_repeatID(const SV_in_repeat &a, const SV_in_repeat &b){
		if(a.repeatID != b.repeatID)
			return a.repeatID < b.repeatID;
		else if(a.SV_len != b.SV_len)
			return a.SV_len < b.SV_len;
		else
			return a.SV_id < b.SV_id;
	}
};

int bed_in_repeat_masker_region(int argc, char *argv[]){
	bool print_log = true;
	//parameters
	char * ref_fn = argv[1];
	char * analysis_bed_fn_in = argv[2];//separate by ','
	char * repeat_masker_region_bed_fn = argv[3];

	faidx_t *c_ref_idx = NULL;
	BED_IDX bed_idx;
	std::vector<BED_ITEM> bed_l;
	//load reference index
	{
		c_ref_idx = reference_index_load(ref_fn);
		bed_idx.init(c_ref_idx);
	}
	//load bed file
	{
		std::vector<std::string> bed_record_l;
		load_string_list_from_file_MAX_line(repeat_masker_region_bed_fn, bed_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &bed_r :bed_record_l){
			split_string(item_value, bed_r.c_str(), "\t");
			bed_l.emplace_back();
			bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
			if(bed_l.back().chrID >= 24){
				bed_l.pop_back();
			}else
				bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
		}
	}
	//load input bed file
	{
		std::vector<std::string> bed_record_l;
		load_string_list_from_file_MAX_line(analysis_bed_fn_in, bed_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &bed_r :bed_record_l){
			split_string(item_value, bed_r.c_str(), "\t");
			if(item_value.size() < 3){
				fprintf(stderr, "BED_in_range: #\n");
			}
			int chrID = faidx_get_chrID(c_ref_idx, item_value[0].c_str(), NULL, 0);
			int pos = atoi(item_value[1].c_str());
			int end = atoi(item_value[2].c_str());
			if(chrID < 0){
				fprintf(stderr, "BED_in_range: %s\n", item_value[0].c_str());
			}else{
				int bed_search_idx = bed_idx.search_bed_full_cover_SV(chrID, pos, end, bed_l);
				//output:
				if(bed_search_idx == -1)
					fprintf(stderr, "BED_in_range: %d\t%d\t%d\tNA\n", chrID, pos, end);
				else
					fprintf(stderr, "BED_in_range: %d\t%d\t%d\t%s\t%s\n", chrID, pos, end, bed_l[bed_search_idx].info1.c_str(), bed_l[bed_search_idx].info2.c_str());
			}
		}
	}
	return 0;
}

int SV_in_simple_repeat(int argc, char *argv[]){
	bool print_log = true;
	//parameters
	char * ref_fn = argv[1];
	char * vcf_fn_in = argv[2];//separate by ','
	char * vntr_region_bed_fn = argv[3];
	char * sample_region_fn = argv[4];

	std::vector<Merging_SV_info> SV_info;
	std::vector<BED_ITEM> bed_l;
	std::vector<SV_in_repeat> sv_repeat;

	std::vector<int> sample_region;

	faidx_t *c_ref_idx = NULL;
	BED_IDX bed_idx;

	//load sample region file:
	{
		std::vector<std::string> sample_record_l;
		load_string_list_from_file_MAX_line(sample_region_fn, sample_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &s_r :sample_record_l){
			split_string(item_value, s_r.c_str(), "\t");
			int main_region_id = atoi(item_value[6].c_str()) - 1;
			sample_region.emplace_back(main_region_id);
		}
	}

	//load reference index
	{
		c_ref_idx = reference_index_load(ref_fn);
		bed_idx.init(c_ref_idx);
	}

	//load vcf binary
	{
		std::vector<std::string> vcf_record_l;
		load_string_list_from_file_MAX_line(vcf_fn_in, vcf_record_l, -1);//todo::
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &vcf_r :vcf_record_l){
			split_string(item_value, vcf_r.c_str(), "\t");
			SV_info.emplace_back();
			SV_info.back().load_data(item_value, item_value_buff);
		}
	}

	//load bed file
	{
		std::vector<std::string> bed_record_l;
		load_string_list_from_file_MAX_line(vntr_region_bed_fn, bed_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &bed_r :bed_record_l){
			split_string(item_value, bed_r.c_str(), "\t");
			bed_l.emplace_back();
			bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
			if(bed_l.back().chrID >= 24){
				bed_l.pop_back();
			}else
				bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
		}
	}

	//search vcf in bed
	{
		for(int SV_id = 0; SV_id < SV_info.size(); SV_id++){
			Merging_SV_info &cur_sv_info = SV_info[SV_id];
			int sv_chr = cur_sv_info.chrID;
			int sv_pos = cur_sv_info.pos;
			int sv_end = (cur_sv_info.SV_len > 0)?sv_pos:(sv_pos - cur_sv_info.SV_len);
			int bed_search_idx = bed_idx.search_bed_full_cover_SV(sv_chr, sv_pos, sv_end, bed_l);
			if(bed_search_idx != -1){
				//analysis if the SV is a true VNTR;
				std::string repeat_string = bed_l[bed_search_idx].repeat_str;
				int max_dis = bed_l[bed_search_idx].repeat_len * 5;
				std::string check_str;
				if(cur_sv_info.ALT_str.size() > cur_sv_info.REF_str.size())
					check_str = cur_sv_info.ALT_str;
				else
					check_str = cur_sv_info.REF_str;

				bool pass_filter = false;
				{
					int total_gap_length = 0;
					int found_pos = 0;
					int repeat_number = 0;
					while(found_pos != -1){
						int found = check_str.find(repeat_string, found_pos + 1);
						//fprintf(stderr, "%ld\t", found);
						int gap_len = 0;
						if(found != -1){
							repeat_number++;
							gap_len = found-found_pos;
						}
						else
							gap_len = check_str.size()-found_pos;
						if(gap_len > max_dis){
							//fprintf(stderr, "GAP:%d, %d, %d\n", found_pos, found, gap_len);
							total_gap_length += gap_len;
						}
						found_pos = found;
					}
					int final_sv_len = check_str.size() - 1 - total_gap_length;
					fprintf(stderr, "Final repeat_number: %d total_gap_length %d final_sv_len %d %s %s %s\n", repeat_number, total_gap_length, final_sv_len, cur_sv_info.REF_str.c_str(), cur_sv_info.ALT_str.c_str(), repeat_string.c_str());
					if(repeat_number != 0 && final_sv_len >= 30 && final_sv_len*2 > check_str.size()){
						pass_filter = true;
					}else
						pass_filter = false;
				}

				if(pass_filter){
					sv_repeat.emplace_back();
					sv_repeat.back().SV_id = SV_id;
					sv_repeat.back().repeatID = bed_search_idx;
					sv_repeat.back().SV_len = cur_sv_info.SV_len;
				}
//				if(12422321 == cur_sv_info.pos){
//					fprintf(stderr, " ");
//				}
			}else{
				//cur_sv_info.show();
				//fprintf(stderr, "NO data\n");
			}
		}

		for(BED_ITEM & b : bed_l)
			b.init_SV_min_max();

		//analysis:
		std::sort(sv_repeat.begin(), sv_repeat.end(), SV_in_repeat::cmp_by_repeatID);
		int old_repeatID = -1;
		for(SV_in_repeat & sr : sv_repeat){
			if(old_repeatID != sr.repeatID){
				bed_l[sr.repeatID].show();
			}
			old_repeatID = sr.repeatID;
			SV_info[sr.SV_id].show();
			//set the max and min SV size:
			bed_l[sr.repeatID].store_SV_min_max(SV_info[sr.SV_id].SV_len);
		}

		//set the repeat_range and the original_index for each BED
		int bed_original_index = 0;
		for(BED_ITEM & b : bed_l){
			b.store_repeat_range();
			b.bed_original_index = bed_original_index;
			bed_original_index ++;
		}

		//sort by repeat range
		std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_repeat_range);
		//set the repeat range index
		int bed_repeat_range_index = 0;
		for(BED_ITEM & b : bed_l){
			if(b.repeat_range != -1){
				b.bed_repeat_range_index = bed_repeat_range_index;
				bed_repeat_range_index ++;
			}
		}

		//sort by original index
		std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_bed_original_index);

		std::vector<int> sample_region_counter;
		sample_region_counter.resize(5);
		//set the repeat range index for all SVs
		old_repeatID = -1;
		for(SV_in_repeat & sr : sv_repeat){
			if(old_repeatID != sr.repeatID){
				fprintf(stderr, "bed_l show final\t");
				bed_l[sr.repeatID].show();
			}
			old_repeatID = sr.repeatID;
			int repeat_unit_n = SV_info[sr.SV_id].SV_len/bed_l[sr.repeatID].repeat_len;
			fprintf(stderr, "SV_in_range: %d\t%d\t%d\t", bed_l[sr.repeatID].bed_original_index, bed_l[sr.repeatID].bed_repeat_range_index, repeat_unit_n);
			SV_info[sr.SV_id].count_sample_region(sample_region, sample_region_counter);
			SV_info[sr.SV_id].show();
		}
	}
	return 0;
}

//gcSV tools vntr_analysis ref.fa single_sample_vcf_list.txt VNTR_region.bed
//1kgp_sample_info.txt

int vntr_analysis(int argc, char *argv[]){
	bool print_log = true;
	//parameters
	char * ref_fn = argv[1];
	char * single_sample_vcf_list = argv[2];//separate by ','
	char * vntr_region_bed_fn = argv[3];
	char * sample_region_fn = argv[4];

	faidx_t *c_ref_idx = NULL;
	std::vector<Merging_SV_info> final_SV_l;
	std::vector<BED_ITEM> bed_l;
	std::vector<SV_in_repeat> sv_repeat;

	std::vector<std::string> sample_file_name;//vcf file name for each samples
	std::vector<int> sample_region;

	BED_IDX bed_idx;//simple repeat BED file
	//char *joint_vcf_fn = argv[1];
	std::vector<SV_store_index> SV_idx;
		//index
		//load sample list
	//S1: load reference index
	{
		c_ref_idx = reference_index_load(ref_fn);
		bed_idx.init(c_ref_idx);
		NGS_SV_MERGING_HANDLER::init_SV_index(SV_idx, c_ref_idx);
	}

	//load SVs from all single sample vcf data:
	{
		//S2: load sample file names
		load_string_list_from_file_MAX_line(single_sample_vcf_list, sample_file_name, 1000);
		//load single sample vcf
		int total_SV_num = 0;
		char *c_sv_type = (char *)malloc(1000);
		uint32_t * GT_output_buff = (uint32_t *) xcalloc(10,sizeof(uint32_t));
		for(uint sample_ID = 0; sample_ID < sample_file_name.size(); sample_ID++){
			//if(sample_ID > 10) break;
			fprintf(stderr, "Current handle %d: %s\n",sample_ID, sample_file_name[sample_ID].c_str());
			BCF_FILE input_vcf;

			VCF_open_read(&input_vcf, sample_file_name[sample_ID].c_str());
			bcf_hdr_t *vcf_header = input_vcf.header;
			int new_SV_num = 0;
			//load each SVs
			do
			{
				bcf1_t *c_r = &( input_vcf.r);
				bcf_unpack(c_r, BCF_UN_INFO);
				int sv_chr_ID = c_r->rid;
				//vcf_read.header->
				//faidx_iseq(c_ref_idx, SV_CHR_ID);
				int chrID = faidx_get_chrID(c_ref_idx, bcf_hdr_id2name(vcf_header, sv_chr_ID), NULL, 0);
				if(chrID >= 24)
					continue;
				vcf_get_sv_type(vcf_header, c_r, c_sv_type);
				if(strcmp(c_sv_type, "BND") == 0)
					continue;
				if(*c_r->d.flt != 0 && *c_r->d.flt != 13)
					continue;
				int SV_POS = c_r->pos;
				int SV_length = 0;
				vcf_get_sv_LENGTH(input_vcf.header, c_r, &SV_length);
				std::string REF_str; std::string ALT_str;
				REF_str.append(c_r->d.allele[0]);
				ALT_str.append(c_r->d.allele[1]);
				int mdat;
				//dat = 0/0	0/2	2/2	2/4	4/4	4/6	5/6	5/7	6/6	7/7
				//GT  = ./.	./0	0/0	0/1	1/1	1/2	1/2	1|2	2/2	2|2
				bcf_get_genotypes(vcf_header,c_r,(void**)&GT_output_buff, &mdat);
				//store the SVs to SV list
				bool is_new = true;
				bool high_complex_region = (SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.size() > 100);
				if(high_complex_region){
					fprintf(stderr, " ");
				}
				for(Merging_SV_info & sv: SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs){
					if(sv.sameWith(SV_POS, ALT_str, REF_str,high_complex_region, SV_length)){
						sv.supportSAMPLE_ID.emplace_back();
						sv.supportSAMPLE_ID.back().store(sample_ID, GT_output_buff[0], GT_output_buff[1], c_r->qual);
						is_new = false;
						break;
					}
				}
				if(is_new){
					new_SV_num ++;
					SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.emplace_back();
					Merging_SV_info & sv = SV_idx[chrID].b[SV_POS/SV_MERGE_IDX_B_SIZE].SVs.back();
					sv.store_basic(SV_POS, ALT_str, REF_str, SV_length);
					sv.supportSAMPLE_ID.emplace_back();
					sv.supportSAMPLE_ID.back().store(sample_ID, GT_output_buff[0], GT_output_buff[1], c_r->qual);
				}
			}while(VCF_next(&input_vcf));//read one
			total_SV_num += new_SV_num;
			fprintf(stderr, "\n new_SV_num %d total_SV_num %d\n", new_SV_num, total_SV_num);
			bcf_close(input_vcf.file);
		}

		//output the final result
		{
			//std::vector<SV_INFO_vcf_add_alt_string> SV_info;
			for(uint chrID = 0; chrID < SV_idx.size(); chrID++){
				for(uint b_id = 0; b_id < SV_idx[chrID].b.size(); b_id++){
					//for each block
					std::vector<Merging_SV_info> SV_l = SV_idx[chrID].b[b_id].SVs;
					// sort:
					std::sort(SV_l.begin(), SV_l.end(), Merging_SV_info::cmp_by_pos);
					final_SV_l.emplace_back();
					for(uint SV_idx_in_B = 0; SV_idx_in_B < SV_l.size(); SV_idx_in_B++){
						final_SV_l.emplace_back();
						std::swap(final_SV_l.back(), SV_l[SV_idx_in_B]);
					}
				}
			}
		}
	}

	//load sample region file:
	{
		std::vector<std::string> sample_record_l;
		load_string_list_from_file_MAX_line(sample_region_fn, sample_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &s_r :sample_record_l){
			split_string(item_value, s_r.c_str(), "\t");
			int main_region_id = atoi(item_value[6].c_str()) - 1;
			sample_region.emplace_back(main_region_id);
		}
	}

	//load vntr bed file
	{
		std::vector<std::string> bed_record_l;
		load_string_list_from_file_MAX_line(vntr_region_bed_fn, bed_record_l, -1);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_buff;
		for(std::string &bed_r :bed_record_l){
			split_string(item_value, bed_r.c_str(), "\t");
			bed_l.emplace_back();
			bed_l.back().load_data(item_value, c_ref_idx, item_value_buff);
			if(bed_l.back().chrID >= 24){
				bed_l.pop_back();
			}else
				bed_idx.store_bed_item(bed_l.size() - 1, bed_l.back());
		}
	}

	//search vcf in bed
	{
		for(int SV_id = 0; SV_id < final_SV_l.size(); SV_id++){
			Merging_SV_info &cur_sv_info = final_SV_l[SV_id];
			int sv_chr = cur_sv_info.chrID;
			int sv_pos = cur_sv_info.pos;
			int sv_end = (cur_sv_info.SV_len > 0)?sv_pos:(sv_pos - cur_sv_info.SV_len);
			int bed_search_idx = bed_idx.search_bed_full_cover_SV(sv_chr, sv_pos, sv_end, bed_l);
			if(bed_search_idx != -1){
				//analysis if the SV is a true VNTR;
				std::string repeat_string = bed_l[bed_search_idx].repeat_str;
				int max_dis = bed_l[bed_search_idx].repeat_len * 5;
				std::string check_str;
				if(cur_sv_info.ALT_str.size() > cur_sv_info.REF_str.size())
					check_str = cur_sv_info.ALT_str;
				else
					check_str = cur_sv_info.REF_str;

				bool pass_filter = false;
				{
					int total_gap_length = 0;
					int found_pos = 0;
					int repeat_number = 0;
					while(found_pos != -1){
						int found = check_str.find(repeat_string, found_pos + 1);
						//fprintf(stderr, "%ld\t", found);
						int gap_len = 0;
						if(found != -1){
							repeat_number++;
							gap_len = found-found_pos;
						}
						else
							gap_len = check_str.size()-found_pos;
						if(gap_len > max_dis){
							//fprintf(stderr, "GAP:%d, %d, %d\n", found_pos, found, gap_len);
							total_gap_length += gap_len;
						}
						found_pos = found;
					}
					int final_sv_len = check_str.size() - 1 - total_gap_length;
					fprintf(stderr, "Final repeat_number: %d total_gap_length %d final_sv_len %d %s %s %s\n", repeat_number, total_gap_length, final_sv_len, cur_sv_info.REF_str.c_str(), cur_sv_info.ALT_str.c_str(), repeat_string.c_str());
					if(repeat_number != 0 && final_sv_len >= 30 && final_sv_len*2 > check_str.size()){
						pass_filter = true;
					}else
						pass_filter = false;
				}

				if(pass_filter){
					sv_repeat.emplace_back();
					sv_repeat.back().SV_id = SV_id;
					sv_repeat.back().repeatID = bed_search_idx;
					sv_repeat.back().SV_len = cur_sv_info.SV_len;
				}
			}else{/*DO NOTHING*/}
		}

		for(BED_ITEM & b : bed_l)
			b.init_SV_min_max();

		//analysis:
		std::sort(sv_repeat.begin(), sv_repeat.end(), SV_in_repeat::cmp_by_repeatID);
		int old_repeatID = -1;
		for(SV_in_repeat & sr : sv_repeat){
			if(old_repeatID != sr.repeatID){
				bed_l[sr.repeatID].show();
			}
			old_repeatID = sr.repeatID;
			final_SV_l[sr.SV_id].show();
			//set the max and min SV size:
			bed_l[sr.repeatID].store_SV_min_max(final_SV_l[sr.SV_id].SV_len);
		}

		//set the repeat_range and the original_index for each BED
		int bed_original_index = 0;
		for(BED_ITEM & b : bed_l){
			b.store_repeat_range();
			b.bed_original_index = bed_original_index;
			bed_original_index ++;
		}

		//sort by repeat range
		std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_repeat_range);
		//set the repeat range index
		int bed_repeat_range_index = 0;
		for(BED_ITEM & b : bed_l){
			if(b.repeat_range != -1){
				b.bed_repeat_range_index = bed_repeat_range_index;
				bed_repeat_range_index ++;
			}
		}

		//sort by original index
		std::sort(bed_l.begin(), bed_l.end(), BED_ITEM::cmp_by_bed_original_index);

		std::vector<int> sample_region_counter;
		sample_region_counter.resize(5);
		//set the repeat range index for all SVs
		old_repeatID = -1;
		for(SV_in_repeat & sr : sv_repeat){
			if(old_repeatID != sr.repeatID){
				fprintf(stderr, "bed_l show final\t");
				bed_l[sr.repeatID].show();
			}
			old_repeatID = sr.repeatID;
			int repeat_unit_n = final_SV_l[sr.SV_id].SV_len/bed_l[sr.repeatID].repeat_len;
			fprintf(stderr, "SV_in_range: %d\t%d\t%d\t", bed_l[sr.repeatID].bed_original_index, bed_l[sr.repeatID].bed_repeat_range_index, repeat_unit_n);
			final_SV_l[sr.SV_id].count_sample_region(sample_region, sample_region_counter);
			final_SV_l[sr.SV_id].show();
		}
	}

	fprintf(stderr, "\nEND\n");
	return 0;


}
#define MIN_SV_LEN_GIAB 50
//get the SV in the original GIAB vcf file(mixed with all types of variation) and store it into output file
//file type: GIAB CCS DATA: deepvariant_.GRCh38_15kb_37X_SequelII.vcf
//Only get SV when : the ref_len >= 50 base pair or the Alt_len >= 50 base pair
void vcf_GIAB_getSV(char *fn_in, char * fn_out)
{
	bool is_compression = false;//BCF or VCF
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read
	BCF_FILE vcf_w;//vcf for write
	VCF_open_write(&vcf_w, fn_out, is_compression);
	bcf_hdr_write(vcf_w.file, vcf_r.header);
	char *c_sample = (char *)malloc(1000);
	char *c_sv_type = (char *)malloc(1000);
	int all_sample = false;
	int all_type = false;
	int all_chrom = false;
	char sample_name[10] = "ALL";
	char sv_type[10] = "ALL";
	char CHROM_ID[10] = "ALL";
	if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)	all_sample = true;
	if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)			all_type = true;
	if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)		all_chrom = true;
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		//unpack the vcf data to get the alt string
		bcf_unpack(c_r, BCF_UN_STR);
		//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
		if(all_sample == false)
		{
			vcf_get_sample(vcf_r.header, c_r, c_sample);
			if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
				continue;
		}
		if(all_type == false)
		{
			vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
			if(strcmp(c_sv_type, sv_type) != 0)
				continue;
		}
		if(all_chrom == false)
		{
			if(c_r->rid != strtol(CHROM_ID,0,10))
				continue;
		}
		//check the alleles
		uint32_t refLen = strlen(c_r->d.allele[0]);
		uint32_t maxAlleleLen = 0;
		for(uint32_t i = 1; i < c_r->n_allele; i++){
			uint32_t allele_len = strlen(c_r->d.allele[i]);
			maxAlleleLen = MAX(maxAlleleLen, allele_len);
		}
		if(refLen >= MIN_SV_LEN_GIAB || maxAlleleLen >= MIN_SV_LEN_GIAB){
			fprintf(stderr, "%d\t%d\t%d\t%d\t\n", c_r->rid, c_r->pos, refLen, maxAlleleLen);
			bcf_write(vcf_w.file, vcf_r.header, c_r);
		}

		//vcf_get_genotype(vcf_r.header, c_r);
	}
	//close
	bcf_close(vcf_r.file);
	bcf_close(vcf_w.file);
}

//get the SV in the original GIAB vcf file(mixed with all types of variation) and store it into output file
//file type: GIAB CCS DATA: deepvariant_.GRCh38_15kb_37X_SequelII.vcf
//Only get SV when : the ref_len >= 50 base pair or the Alt_len >= 50 base pair
void vcf_SURVIVOR_getSV(char *fn_in, char *sv_type)
{
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read

	char *c_sample = (char *)malloc(1000);
	char *c_sv_type = (char *)malloc(1000);
	int all_sample = false;
	int all_type = false;
	int all_chrom = false;
	char sample_name[10] = "ALL";
	char CHROM_ID[10] = "ALL";
	bcf_hdr_t *header = vcf_r.header;

	std::set<std::string> SV_TYPE_SET;

	if(strcmp("all", sample_name) == 0 || strcmp("ALL", sample_name) == 0)	all_sample = true;
	if(strcmp("all", sv_type) == 0 || strcmp("ALL", sv_type) == 0)			all_type = true;
	if(strcmp("all", CHROM_ID) == 0 || strcmp("ALL", CHROM_ID) == 0)		all_chrom = true;

	do//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		//unpack the vcf data to get the alt string
		bcf_unpack(c_r, BCF_UN_STR);
		//int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
		if(all_sample == false)
		{
			vcf_get_sample(vcf_r.header, c_r, c_sample);
			if(*c_sample != 0 && strcmp(c_sample, sample_name) != 0)
				continue;
		}
		if(all_type == false)
		{
			vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
			std::string str_SV_TYPE(c_sv_type);
			SV_TYPE_SET.emplace(str_SV_TYPE);
			if(strcmp(c_sv_type, sv_type) != 0)
				continue;
		}
		if(all_chrom == false)
		{
			if(c_r->rid != strtol(CHROM_ID,0,10))
				continue;
		}
		//get the start position:
		int32_t stPos = c_r->pos;
		//get the end position:
		int32_t edPos = 0;
		vcf_get_sv_END(vcf_r.header, c_r, &edPos);
		const char *chrID = bcf_hdr_id2name(header, c_r->rid);
		int32_t SV_LEN = edPos - stPos;
		if(SV_LEN < 0 || SV_LEN > 10000)
			SV_LEN = 0;
		printf("%s\t%d\t%d\t\n",chrID , stPos, SV_LEN);// chrID+st+len

	}while(VCF_next(&vcf_r));

	for(auto & t :SV_TYPE_SET)
		fprintf(stderr, "%s\t",t.c_str());

	//close
	bcf_close(vcf_r.file);
}

void SVRegionCombine(char* regionListFN, int SV_EDGE_LENGTH){
	//step 1: load region, 其中每一行的c_r->pos - 1000 作为开始 position； c_r->pos + refLen + 1000 作为结束 position
	char *temp = new char[MAX_LINE_LENGTH];//1M
	std::vector<RefRegion>regionList;
	std::ifstream regionListFile(regionListFN);
	std::map<std::string, uint32_t> name2ID;
	std::vector<std::string> id2Name;

	while(true){
		regionListFile.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0)
			break;
		//get chrID
		char *token = strtok(temp, "\t");
		std::string strID(token); int intID = 0;
		std::map<std::string, uint32_t>::iterator it = name2ID.find(strID);
		if(it!=name2ID.end()){
			intID = it->second;
		}
		else{
			intID = name2ID.size();
			name2ID[strID] = name2ID.size();
			id2Name.emplace_back(strID);
		}

		//get refSt
		token = strtok(NULL, "\t");
		uint32_t refSt = strtoul(token, NULL, 10);
		//get refEd
		token = strtok(NULL, "\t");
		uint32_t refLen = strtoul(token, NULL, 10);
		regionList.emplace_back(intID, refSt - SV_EDGE_LENGTH, refSt + refLen + SV_EDGE_LENGTH);
	}
	regionListFile.close();
	delete [] temp;
	std::sort(regionList.begin(), regionList.end(), RefRegion::cmp_by_pos);
	//merge region
	auto r_ed = regionList.end();//for each sve
	int total_output_ITEM = 0;
	for(auto r = regionList.begin(); r < r_ed;)	{
		auto r_try = r + 1;
		for(; r_try < r_ed && r->region_overlap(*r_try); r_try++)//use sve as main, than try to combine
			r->Combine(*r_try, true);
		printf("%s\t%d\t%d\t\n", id2Name[r->chr_ID].c_str(), r->st_pos, r->ed_pos);
		total_output_ITEM++;
		r = r_try;
	}
	fprintf(stderr, "total_output_ITEM: %d\n", total_output_ITEM);
}

#define MAX_read_LEN 256
void GIAB_SV_region_full_test(char *input_bam_fn, char *ref_fn, char * region_fn, int readLen)
{
	xassert(readLen < MAX_read_LEN, "MAX read length: 250 bp");
	fprintf(stderr, "\n\n V1.02\n\n");
	//read bam/cram file:
	Bam_file c_b;
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
	bam_hdr_t* hdr = c_b._hdr;

	//get sv region:
	char *temp = new char[MAX_LINE_LENGTH];//1M
	std::vector<R_region>regionList;
	std::ifstream regionListFile(region_fn);
	while(true){
		regionListFile.getline(temp, MAX_LINE_LENGTH);
		if(*temp == 0 || *temp == 'N')
			break;
		//get chrID
		char *token = strtok(temp, "\t");
		const char *strID = token;
		int chrID = bam_name2id(hdr, strID) ;

		//get refSt
		token = strtok(NULL, "\t");
		uint32_t regionSt = strtoul(token, NULL, 10);
		//get refEd
		token = strtok(NULL, "\t");
		uint32_t regionEnd = strtoul(token, NULL, 10);
		//
		regionList.emplace_back();
		R_region & r = regionList.back();
		r.chr_ID = chrID; r.st_pos = regionSt; r.ed_pos = regionEnd;
	}
	regionListFile.close();
	delete [] temp;

	char SA_TAG[3] = "SA";
	for(R_region & r: regionList){
		//reset region
		resetRegion_ID(&c_b, &r);
		//reset calculators
		uint32_t readNum = 0;
		//mapQ, region 0~60, if mapQ > 60; MAPQ[61]++
		uint32_t MAPQCount[62] = {0};
		uint32_t MAPQ_less_than30 = 0;
		//insert size: when it is [0 ~1023]: 0~15 store in insertSize[0], 16~31 in insertSize[1].....
		uint32_t insertSizeCount[64] = {0};
		uint32_t insertSizeOver1K = 0;//the read of insert size > 1023
		uint64_t totalInsert_SIZE = 0;
		uint32_t insertSizeOver100K = 0;//the read of insert size > 10000
		//mate unmapped:
		uint32_t mate_unmapped_count = 0;
		// mate_different_chromsome
		uint32_t mate_different_chromsome_count = 0;

		uint32_t insertSizeFF = 0;//both forward:
		uint32_t insertSizeRR = 0;//both reverse:
		uint32_t insertSizeRF = 0;//the read of smaller position is reverse, the other is forward
		uint32_t insertSizeFR = 0;//the read of smaller position is forward, the other is reverse

		//SA signal
		uint32_t SA_read_len[MAX_read_LEN] = {0};
		//SNP and INDEL
		uint32_t SNP_INDEL_len[MAX_read_LEN] = {0};
		uint32_t SNP_INDEL_total = 0;
		//soft/hard clip
		uint32_t clip_len_left[MAX_read_LEN] = {0};
		uint32_t clip_len_right[MAX_read_LEN] = {0};
		uint32_t total_clip_len = 0;
		//flags
		uint32_t flag_count[MAX_uint16_t] = {0};

		// analysis the reads
		while (bam_next(&c_b)) {
			bam1_t *b = &(c_b._brec);
			readNum++;
			fprintf(stderr, "%d\n", readNum);
			//mapq:
			uint8_t mapq = b->core.qual;
			if(mapq > 60)	MAPQCount[61]++;
			else{			MAPQCount[mapq]++;	if(mapq < 30)	MAPQ_less_than30 ++;	}
			bool abnormal_insert_size = false;
			//mate read:
			if(bam_is_mate_unmapped(b)) 	{abnormal_insert_size = true; mate_unmapped_count ++;}
			if(b->core.tid != b->core.mtid)	{abnormal_insert_size = true; mate_different_chromsome_count ++;}
			//insert size
			if(!abnormal_insert_size){
				int32_t insertSizeOri = b->core.isize;
				uint32_t insertSize = ABS(insertSizeOri);
				if(insertSize <= 1023) insertSizeCount[insertSize/16]++;
				else insertSizeOver1K++;

				if(insertSize <= 100000) { totalInsert_SIZE += insertSize;}
				else insertSizeOver100K++;

			}
			//read pair orientation
			bool direction = bam_is_fwd_strand(b);
			bool mateDirection = bam_is_mate_fwd_strand(b);
			if(b->core.pos > b->core.mpos) std::swap(direction, mateDirection);

			if		(direction == FORWARD && mateDirection == FORWARD) insertSizeFF++;
			else if (direction == FORWARD && mateDirection == REVERSE) insertSizeFR++;
			else if (direction == REVERSE && mateDirection == FORWARD) insertSizeRF++;
			else if (direction == REVERSE && mateDirection == REVERSE) insertSizeRR++;

			//SA signal:
			const char * SA_String = bam_get_string_tag(b, SA_TAG);
			if(SA_String != NULL){
				uint32_t SA_len = strlen(SA_String);
				SA_read_len[SA_len]++;
			}

			//SNP and INDEL
			int INDEL_len = 0;
			int NM_len = 0;
			bam_get_INDEL_NM(b,&INDEL_len, &NM_len);

			SNP_INDEL_len[NM_len]++;
			SNP_INDEL_total+= NM_len;

			//clip:
			int soft_left; int soft_right;
			bam_has_SH_cigar(b, &soft_left, &soft_right);
			clip_len_left[soft_left]++;
			clip_len_right[soft_right]++;
			total_clip_len += (soft_left + soft_right);

			//flag
			flag_count[b->core.flag]++;
		}

		// out put information for that region
		//title:
		printf("[%d:%d~%d]", r.chr_ID, r.st_pos, r.ed_pos);
		printf("\t[RN:%d]", readNum);
		float readDepth = (float)readLen * readNum / (r.ed_pos - r.st_pos);
		printf("\t[DP:%f]", readDepth);
		if(readNum == 0) readNum = 1;
		printf("\t[MUM:%d, %f%%]", mate_unmapped_count,(float)mate_unmapped_count*100/readNum );
		printf("\t[MAPQ<30:%f%%]", (float)MAPQ_less_than30*100/readNum);
		printf("\t[AVG.ISIZE(<100000):%f]", (float)totalInsert_SIZE/(readNum - insertSizeOver100K));
		printf("\t[>1023:%f%%]", (float)insertSizeOver1K*100/readNum);

		printf("\t[AVG.SNP_INDEL:%f]", (float)SNP_INDEL_total/(readNum));
		printf("\t[AVG.CLIP:%f]", (float)total_clip_len/(readNum));
		printf("\t[AVG.FR:%f%%]", (float)insertSizeFR*100/(readNum));

		//mapq:
		printf("\tMAPQ: ");
		for(int i = 0; i < 61; i++)
			if(MAPQCount[i] > 0)
				printf("[%d:%d] ", i, MAPQCount[i]);
		if(MAPQCount[61] > 0)
			printf("[> 60:%d] ", MAPQCount[61]);
		//mate_unmapped_count:

		//mate_different_chromsome_count:
		printf("\t[MDC:%d]", mate_different_chromsome_count);
		//insert size:
		printf("\tIS: ");
		for(int i = 0; i < 64; i++)
			if(insertSizeCount[i] > 0)
				printf("[%d:%d] ", i, insertSizeCount[i]);
		if(insertSizeOver1K > 0)
			printf("[> 1023:%d] ", insertSizeOver1K);
		//orientation
		printf("\t[FF:%d]", insertSizeFF);
		printf("\t[FR:%d]", insertSizeFR);
		printf("\t[RF:%d]", insertSizeRF);
		printf("\t[RR:%d]", insertSizeRR);
		// SA len
		printf("\tSA_L: ");
		for(int i = 0; i < MAX_read_LEN; i++)
			if(SA_read_len[i] > 0)
				printf("[%d:%d] ", i, SA_read_len[i]);
		//SNP INDEL
		printf("\tSNP_INDEL: ");
		for(int i = 0; i < MAX_read_LEN; i++)
			if(SNP_INDEL_len[i] > 0)
				printf("[%d:%d] ", i, SNP_INDEL_len[i]);
		//clip_len_left
		printf("\tCL: ");
		for(int i = 0; i < MAX_read_LEN; i++)
			if(clip_len_left[i] > 0)
				printf("[%d:%d] ", i, clip_len_left[i]);
		//clip_len_right
		printf("\tCR: ");
		for(int i = 0; i < MAX_read_LEN; i++)
			if(clip_len_right[i] > 0)
				printf("[%d:%d] ", i, clip_len_right[i]);
		//flags
		printf("\tFLAG: ");
		for(int i = 0; i < MAX_uint16_t; i++)
			if(flag_count[i] > 0)
				printf("[%d:%d] ", i, flag_count[i]);
		//END

		printf("\n");
	}
	//close file
	bam_file_close(&c_b);
}

struct ASS_INFO{
	ASS_INFO(){}
	ASS_INFO(const char * name_){strcpy(name, name_); }
	uint32_t length = 0;
	char name[100];
	void print(FILE * o, char endl){
		fprintf(o, "[%s %d]\t%c", name, length, endl);
	}
};

struct LIFTOVER{
	LIFTOVER(
			uint32_t tid_ass_,  uint32_t pos_ass_st_,  uint32_t pos_ass_ed_,
			uint32_t tid_hs37_, uint32_t pos_hs37_st_, uint32_t pos_hs37_ed_,
			uint32_t length_, bool direction_):
				tid_ass(tid_ass_),   pos_ass_st(pos_ass_st_),   pos_ass_ed(pos_ass_ed_),
				tid_hs37(tid_hs37_), pos_hs37_st(pos_hs37_st_), pos_hs37_ed(pos_hs37_ed_),
				length(length_), direction(direction_){}
	LIFTOVER():
		tid_ass(0),  pos_ass_st(0),  pos_ass_ed(0),
		tid_hs37(0), pos_hs37_st(0), pos_hs37_ed(0),
		length(0), direction(0){}

	uint32_t tid_ass;
	uint32_t pos_ass_st;
	uint32_t pos_ass_ed;

	uint32_t tid_hs37;
	uint32_t pos_hs37_st;
	uint32_t pos_hs37_ed;

	uint32_t length;

	bool direction;

	void print(FILE * o, std::vector<ASS_INFO> & assemble_name_str_list){
		fprintf(o,
				"\t\t[ASS:INFO %s %d ]"
				" [ %d %d %d ]"
				"[ %d %d %d ] "
				"%d %c\n",
				assemble_name_str_list[tid_ass].name, assemble_name_str_list[tid_ass].length,
				tid_ass,  pos_ass_st,  pos_ass_ed,
				tid_hs37, pos_hs37_st, pos_hs37_ed,
				length, (direction==FORWARD)?'+':'-');
	}

	static inline int cmp_by_ref_pos(const LIFTOVER &a, const LIFTOVER &b){
		if(a.tid_hs37 != b.tid_hs37)
			return a.tid_hs37 < b.tid_hs37;
		else
			return a.pos_hs37_st < b.pos_hs37_st;
	}

	static inline int cmp_by_ass_pos(const LIFTOVER &a, const LIFTOVER &b){//todo::
		if(a.tid_ass != b.tid_ass)
			return a.tid_ass < b.tid_ass;
		else
			return a.pos_ass_st < b.pos_ass_st;
	}

};

struct liftoverRegion_direction{

	liftoverRegion_direction(
			uint32_t tid_ass_,  uint32_t pos_ass_st_,  uint32_t pos_ass_ed_,
			bool direction_):
				tid_ass(tid_ass_),   pos_ass_st(pos_ass_st_),   pos_ass_ed(pos_ass_ed_),
				direction(direction_){}
	liftoverRegion_direction(){}

	uint32_t tid_ass = 0;
	uint32_t pos_ass_st  = 0;
	uint32_t pos_ass_ed = 0;

	bool direction = 0;

	static inline int cmp_by_ass_pos(const liftoverRegion_direction &a, const liftoverRegion_direction &b){//todo::
		if(a.tid_ass != b.tid_ass)
			return a.tid_ass < b.tid_ass;
		else
			return a.pos_ass_st < b.pos_ass_st;
	}
};

struct LIFTOVER_index{
	std::map<std::string, uint32_t> assembly_name_2_ID;
	std::vector<ASS_INFO> id2Name;
	std::vector<LIFTOVER> liftIndex;

	std::vector<liftoverRegion_direction> LD;

	void dump(char * fn){
		FILE * f = xopen(fn, "wb");
		uint64_t id2Name_size = id2Name.size();
		fwrite(&id2Name_size, sizeof(uint64_t), 1, f);
		fwrite(&(id2Name[0]), sizeof(ASS_INFO), id2Name_size, f);

		uint64_t liftIndex_size = liftIndex.size();
		fwrite(&liftIndex_size, sizeof(uint64_t), 1, f);
		fwrite(&(liftIndex[0]), sizeof(LIFTOVER), liftIndex_size, f);

		//LD
		uint64_t LD_size = LD.size();
		fwrite(&LD_size, sizeof(uint64_t), 1, f);
		fwrite(&(LD[0]), sizeof(liftoverRegion_direction), LD_size, f);

		fclose(f);
	}

	void load(char * fn){
		FILE * f = xopen(fn, "rb");

		uint64_t id2Name_size = 0;
		xread(&id2Name_size, sizeof(uint64_t), 1, f);
		id2Name.clear();
		id2Name.resize(id2Name_size);
		xread(&(id2Name[0]), sizeof(ASS_INFO), id2Name_size, f);

		uint64_t liftIndex_size = 0;
		xread(&liftIndex_size, sizeof(uint64_t), 1, f);
		liftIndex.clear();
		liftIndex.resize(liftIndex_size);
		xread(&(liftIndex[0]), sizeof(LIFTOVER), liftIndex_size, f);

		//LD
		uint64_t LD_size = 0;
		xread(&LD_size, sizeof(uint64_t), 1, f);
		LD.clear();
		LD.resize(liftIndex_size);
		xread(&(LD[0]), sizeof(liftoverRegion_direction), LD_size, f);

		//debug code:
		//for(auto l: liftIndex)	{l.print();}
		fclose(f);
	}

};

void liftoverBuildingIndex(char *input_bam_fn, char *liftOverData_fn){
	fprintf(stderr, "\n\n V1.09\n\n");
		//read bam/cram file:
	Bam_file c_b;
	bam_file_open(input_bam_fn, NULL, NULL, &c_b);
	bam_hdr_t* hdr = c_b._hdr;
	uint64_t readNum = 0;

	//print header
	for(int i = 0; i < hdr->n_targets; i++)
		fprintf(stderr, "%d %s %d \n", i, hdr->target_name[i], hdr->target_len[i]);

	bam1_t b = {0};//BAM record for the first read in a pair
	// analysis the reads'
	//reset region
	//resetRegion_ID(&c_b, &r);
	int sam_rst1 = 0;

	LIFTOVER_index idx;
	std::map<std::string, uint32_t> & assembly_name_2_ID = idx.assembly_name_2_ID;
	std::vector<ASS_INFO> &id2Name = idx.id2Name;
	std::vector<LIFTOVER> &liftIndex = idx.liftIndex;
	std::vector<liftoverRegion_direction> &LD =  idx.LD;

	while (1){
		//load SAM 1 & 2
		readNum++; // debug code //if(readNum == 10) break;
		sam_rst1 = sam_read1(c_b._hfp, hdr, &b);
		if(sam_rst1 < 0) break;
		//load the read name
		std::string strID((char *)bam_qname(&b)); int assemblyID = 0;
		std::map<std::string, uint32_t>::iterator it = assembly_name_2_ID.find(strID);
		if(it!=assembly_name_2_ID.end()){assemblyID = it->second;}
		else{
			assemblyID = assembly_name_2_ID.size();
			assembly_name_2_ID[strID] = assembly_name_2_ID.size();
			id2Name.emplace_back(strID.c_str());
		}
		if(!bam_is_secondary(&b) && !bam_is_supplementary(&b) && !bam_is_duplicate(&b))
			id2Name[assemblyID].length = b.core.l_qseq;
		//get cigar
		bool direction1 = bam_is_fwd_strand(&b);

		uint32_t tid = b.core.tid;
		uint64_t st_pos_in_ref = b.core.pos;
		uint64_t st_pos_in_assembly = 0;
		uint32_t* bam_cigar = bam_get_cigar(&b);

		for (uint32_t i = 0; i < b.core.n_cigar; ++i)
		{
			int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			switch (type)
			{
			case CIGAR_MATCH:
			case CIGAR_SEQ_MATCH:
				liftIndex.emplace_back(
					assemblyID, st_pos_in_assembly, st_pos_in_assembly + length,
					tid, st_pos_in_ref, st_pos_in_ref + length,
					length, direction1);
				//liftIndex.back().print(stderr);
				st_pos_in_ref += length; st_pos_in_assembly += length;
				break;
			case CIGAR_INSERT:
				st_pos_in_assembly += length;
				break;
			case CIGAR_DELETE:
				st_pos_in_ref += length;
				break;
			case CIGAR_SKIP:
				break;
			case CIGAR_SOFT_CLIP:
			case CIGAR_HARD_CLIP:
				st_pos_in_assembly += length;
				break;
			case CIGAR_PAD:
				break;
			case CIGAR_SEQ_MISMATCH:
				break;
			default:
				break;
			}
		}
		int cigar_st_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
		int cigar_st_length = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
		uint32_t region_st_in_ass = (cigar_st_type == CIGAR_SOFT_CLIP || cigar_st_type == CIGAR_HARD_CLIP)?cigar_st_length:0;
		LD.emplace_back(assemblyID, region_st_in_ass, liftIndex.back().pos_ass_ed,direction1);
	}
	//close file
	bam_file_close(&c_b);

	//modify LD
	for(auto &l : LD){
		if(l.direction == REVERSE){
			uint32_t ass_len = id2Name[l.tid_ass].length;
			l.pos_ass_st = ass_len - l.pos_ass_st;
			l.pos_ass_ed = ass_len - l.pos_ass_ed;
			std::swap(l.pos_ass_st, l.pos_ass_ed);
		}
	}
	//sort by pos
	std::sort(LD.begin(), LD.end(), liftoverRegion_direction::cmp_by_ass_pos);

	//modify the position for REVERSE records
//	for(auto &l : liftIndex){
//		if(l.direction == REVERSE){
//			uint32_t ass_len = id2Name[l.tid_ass].length;
//			l.pos_ass_st = ass_len - l.pos_ass_st;
//			l.pos_ass_ed = ass_len - l.pos_ass_ed;
//			std::swap(l.pos_ass_st, l.pos_ass_ed);
//		}
//	}

	//store lift over data
	idx.dump(liftOverData_fn);
}

void liftoverSearchRegion(char *liftOverData_fn){
		//read bam/cram file:
	LIFTOVER_index idx;
	idx.load(liftOverData_fn);

	std::vector<ASS_INFO> &id2Name = idx.id2Name;
	std::vector<LIFTOVER> &liftIndex = idx.liftIndex;

	uint64_t ass_total_len = 0;
	for(auto & a: id2Name){	a.print(stderr, '\n');	ass_total_len += a.length;}
	fprintf(stderr, "[total lenghth: %ld]\n", ass_total_len);
	std::sort(liftIndex.begin(), liftIndex.end(), LIFTOVER::cmp_by_ref_pos);
	//for(auto & a: liftIndex){	a.print(stderr, id2Name);}
	std::vector<RefRegion> r_list;

	r_list.emplace_back(1, 1013668, 1014819);
	r_list.emplace_back(7, 699616, 700715);
	r_list.emplace_back(12, 129674126, 129675212);
	r_list.emplace_back(20, 62757546, 62759452);
	r_list.emplace_back(2, 91600041, 91601285);
	r_list.emplace_back(7, 103549585, 103550702);
	r_list.emplace_back(19, 8772442, 8773733);
	r_list.emplace_back(7, 104895073, 104896074);
	r_list.emplace_back(21, 11141346, 11142347);
	r_list.emplace_back(10, 71752046, 71753047);
	r_list.emplace_back(18, 77712216, 77713220);
	r_list.emplace_back(10, 42673500, 42674501);
	for(auto & r : r_list){
		int chrID_ref = r.chr_ID - 1;
		int st_pos_ref = r.st_pos;
		int ed_pos_ref = r.ed_pos;

		auto l_bg = liftIndex.begin();
		auto l_ed = liftIndex.end();
		auto l_c_bg = l_bg;
		for(; l_c_bg < l_ed; l_c_bg++){
			if(l_c_bg->tid_hs37 != chrID_ref)			continue;
			if(l_c_bg->pos_hs37_st < st_pos_ref)			continue;
			else										break;
		}
		auto l_c_ed = l_c_bg;
		for(; l_c_ed < l_ed; l_c_ed++){
			if(l_c_ed->tid_hs37 != chrID_ref)				break;
			if(l_c_ed->pos_hs37_st < ed_pos_ref)				continue;
			else											break;
		}

		LIFTOVER new_l_bg = l_c_bg[-1];

		fprintf(stderr, "Region in hs37d5:\t");
		r.print(stderr);

		int offset_bg = st_pos_ref - new_l_bg.pos_hs37_st;
		new_l_bg.pos_hs37_st += offset_bg;
		new_l_bg.pos_ass_st += offset_bg;

		LIFTOVER new_l_ed = l_c_ed[-1];
		int offset_ed = ed_pos_ref - new_l_ed.pos_hs37_st;
		new_l_ed.pos_hs37_st += offset_ed;
		new_l_ed.pos_ass_st += offset_ed;

		if(new_l_bg.direction == REVERSE){
			//fprintf(stderr, "\t\t before reverse BG:\t");
			//new_l_bg.print(stderr, id2Name);
			uint32_t length = id2Name[new_l_bg.tid_ass].length;
			new_l_bg.pos_ass_st = length - new_l_bg.pos_ass_st;
		}

		if(new_l_ed.direction == REVERSE){
			//fprintf(stderr, "\t\t before reverse ED:\t");
			//new_l_bg.print(stderr, id2Name);
			//new_l_ed.print(stderr, id2Name);
			std::swap(new_l_bg, new_l_ed);
			uint32_t length = id2Name[new_l_ed.tid_ass].length;
			new_l_ed.pos_ass_st = length - new_l_ed.pos_ass_st;
		}

		if(new_l_bg.direction == REVERSE && new_l_ed.direction == REVERSE){		std::swap(new_l_bg, new_l_ed);	}

		fprintf(stderr, "\t\t region liftover begin:\t");
		new_l_bg.print(stderr, id2Name);
		fprintf(stderr, "\t\t region liftover end:\t");
		new_l_ed.print(stderr, id2Name);
		fprintf(stderr, "\t\t final region: [%s:%d-%d]\n", id2Name[new_l_bg.tid_ass].name, new_l_bg.pos_ass_st, new_l_ed.pos_ass_st);
	}
}

void liftoverRead(char *liftOverData_fn, char * liftoverBam_fn){
	//load lift over file

	LIFTOVER_index idx;
	idx.load(liftOverData_fn);
	std::vector<ASS_INFO> &id2Name = idx.id2Name;
	std::vector<LIFTOVER> &liftIndex = idx.liftIndex;
	std::vector<liftoverRegion_direction> &LD = idx.LD;

	//lift over file check/sort and building index
	uint64_t ass_total_len = 0;
	for(auto & a: id2Name){	a.print(stderr, '\n');	ass_total_len += a.length;}
	fprintf(stderr, "[total lenghth: %ld]\n", ass_total_len);
	std::sort(liftIndex.begin(), liftIndex.end(), LIFTOVER::cmp_by_ass_pos); //sort by assembly position
	//for(auto & a: liftIndex){	a.print(stderr, id2Name);}
	//build simple index for the lift index
	std::vector<uint32_t> ass_ID_2_st_pos_in_liftover;
	ass_ID_2_st_pos_in_liftover.emplace_back(0);
	uint32_t old_ass_ID = MAX_uint32_t;
	//may be duplication in some position
	for(uint32_t i = 0; i < liftIndex.size(); i++){
		auto & a = liftIndex[i];
		if(old_ass_ID != a.tid_ass){
			old_ass_ID = a.tid_ass;
			uint32_t old_index = ass_ID_2_st_pos_in_liftover.back();
			ass_ID_2_st_pos_in_liftover.resize(old_ass_ID + 1, old_index);
			ass_ID_2_st_pos_in_liftover[ass_ID_2_st_pos_in_liftover.size() - 1] = i;
		}
	}
	ass_ID_2_st_pos_in_liftover.emplace_back(liftIndex.size());

	//pos in read ----> pos in assmbly ----- > pos in hs37d5
	//read bam/cram file:
	//read bam/cram file:
	Bam_file c_b;
	bam_file_open(liftoverBam_fn, NULL, NULL, &c_b);
	bam_hdr_t* hdr = c_b._hdr;
	uint64_t readNum = 0;
	bam1_t b; int sam_read_rst = 0;

	uint32_t cigar_buff[1000];
	uint32_t cigar_len = 0;

	while (1){
		readNum++; // debug code //if(readNum == 10) break;
		sam_read_rst = sam_read1(c_b._hfp, hdr, &b);
		if(sam_read_rst < 0) break;

		int st_pos_read_in_ass = b.core.pos;
		int ass_tid_in_read = b.core.tid;
		//get end pos:
		uint32_t* bam_cigar = bam_get_cigar(&b);
		int ed_pos_read_in_ass = st_pos_read_in_ass + bam_cigar2rlen( b.core.n_cigar, bam_cigar);

		//search direction in LD list
		bool direction_in_ass = true;
		for(liftoverRegion_direction &l: LD ){
			if(l.tid_ass != ass_tid_in_read)				continue;
			if(l.pos_ass_st < st_pos_read_in_ass)			continue;
			else{direction_in_ass = l.direction; break;}
		}

		//load cigar
		//reverse CIGAR when needed
		cigar_len = b.core.n_cigar;
		bam_cigar = bam_get_cigar(&b);

		//copy the cigar
		if(direction_in_ass == FORWARD)
			for (uint32_t i = 0; i < cigar_len; ++i)
				cigar_buff[i] = bam_cigar[i];
		else{
			for (uint32_t i = 0; i < cigar_len; ++i)
				cigar_buff[i] = bam_cigar[cigar_len - i - 1];
			uint32_t ass_len = id2Name[ass_tid_in_read].length;
			st_pos_read_in_ass = ass_len - st_pos_read_in_ass;
			ed_pos_read_in_ass = ass_len - ed_pos_read_in_ass;
			std::swap(st_pos_read_in_ass, ed_pos_read_in_ass);
		}

		//get beginning lift over block
		//search the lift over index
		auto l_bg = liftIndex.begin() + ass_ID_2_st_pos_in_liftover[ass_tid_in_read];
		auto l_ed = liftIndex.end() + ass_ID_2_st_pos_in_liftover[ass_tid_in_read + 1];
		//for start position
		auto l_c_bg = l_bg;
		for(; l_c_bg < l_ed; l_c_bg++)
			if(l_c_bg->direction != direction_in_ass || l_c_bg->pos_ass_st < st_pos_read_in_ass)		continue;
			else										break;
		//modify the CIGAR

		//todo::::



	}



//	for(auto & r : r_list){
//		int chrID_ref = r.chr_ID - 1;
//		int st_pos_ref = r.st_pos;
//		int ed_pos_ref = r.ed_pos;
//
//		auto l_bg = liftIndex.begin();
//		auto l_ed = liftIndex.end();
//		auto l_c_bg = l_bg;
//		for(; l_c_bg < l_ed; l_c_bg++){
//			if(l_c_bg->tid_hs37 != chrID_ref)			continue;
//			if(l_c_bg->pos_hs37 < st_pos_ref)			continue;
//			else										break;
//		}
//		auto l_c_ed = l_c_bg;
//		for(; l_c_ed < l_ed; l_c_ed++){
//			if(l_c_ed->tid_hs37 != chrID_ref)				break;
//			if(l_c_ed->pos_hs37 < ed_pos_ref)				continue;
//			else											break;
//		}
//
//		LIFTOVER new_l_bg = l_c_bg[-1];
//
//		fprintf(stderr, "Region in hs37d5:\t");
//		r.print(stderr);
//
//		int offset_bg = st_pos_ref - new_l_bg.pos_hs37;
//		new_l_bg.pos_hs37 += offset_bg;
//		new_l_bg.pos_ass += offset_bg;
//
//		LIFTOVER new_l_ed = l_c_ed[-1];
//		int offset_ed = ed_pos_ref - new_l_ed.pos_hs37;
//		new_l_ed.pos_hs37 += offset_ed;
//		new_l_ed.pos_ass += offset_ed;
//
//		if(new_l_bg.direction == REVERSE){
//			fprintf(stderr, "\t\t before reverse BG:\t");
//			new_l_bg.print(stderr, id2Name);
//			uint32_t length = id2Name[new_l_bg.tid_ass].length;
//			new_l_bg.pos_ass = length - new_l_bg.pos_ass;
//		}
//
//		if(new_l_ed.direction == REVERSE){
//			fprintf(stderr, "\t\t before reverse ED:\t");
//			new_l_bg.print(stderr, id2Name);
//			new_l_ed.print(stderr, id2Name);
//			std::swap(new_l_bg, new_l_ed);
//			uint32_t length = id2Name[new_l_ed.tid_ass].length;
//			new_l_ed.pos_ass = length - new_l_ed.pos_ass;
//		}
//
//		if(new_l_bg.direction == REVERSE && new_l_ed.direction == REVERSE){		std::swap(new_l_bg, new_l_ed);	}
//
//		fprintf(stderr, "\t\t region liftover begin:\t");
//		new_l_bg.print(stderr, id2Name);
//		fprintf(stderr, "\t\t region liftover end:\t");
//		new_l_ed.print(stderr, id2Name);
//		fprintf(stderr, "\t\t final region: [%s:%d-%d]\n", id2Name[new_l_bg.tid_ass].name, new_l_bg.pos_ass, new_l_ed.pos_ass);
//	}
}


void print_one_vcf_analysis(FILE *stream, VCF_ANALYSIS *c_va)
{
	fprintf(stream,
			"%d\t"
			"%s\t"
			"%d\t"
			"%d\t"
			"%d\t"
			"%s\n",
			c_va->sample,
			c_va->type,
			c_va->rid,
			c_va->st,
			c_va->ed,
			c_va->ID
			);
}
int VCF_ANALYSIS_cmp_sample_type_pos(const void*a_,const void*b_)
{
	VCF_ANALYSIS *a = (VCF_ANALYSIS *)a_;
	VCF_ANALYSIS *b = (VCF_ANALYSIS *)b_;

	if(a->sample != b->sample)
		return a->sample > b->sample;
	int type_cmp = strcmp(a->type, b->type);
	if(type_cmp != 0)
		return (type_cmp > 0);
	if(a->rid != b->rid)
		return a->rid > b->rid;
	return a->st > b->st;

}

//------------------------------NSTD VCF--------------------------------------//
int get_one_vcf_nstd(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
{
	char sample[128];
	vcf_get_sample(header, vcf, sample);
	int sample_int = 100;
	if(strcmp(sample, "HG00512") == 0)
		sample_int = 0;
	else if(strcmp(sample, "HG00513") == 0)
		sample_int = 1;
	else if(strcmp(sample, "HG00514") == 0)
		sample_int = 2;
	if(sample_int > 2)
		return false;
	vcf_a->sample = sample_int;
	vcf_get_sv_type(header, vcf, vcf_a->type);
	vcf_get_sv_END(header, vcf, &(vcf_a->ed));
	strcpy(vcf_a->ID, vcf->d.id);
	vcf_a->rid = vcf->rid;
	vcf_a->st = vcf->pos;
	return true;
}

void vcf_nstd_load(char *fn_in, VCF_A_L *vcf_l)
{
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read
	bcf_hdr_t *header = vcf_r.header;
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		VCF_ANALYSIS *c_va = NULL;
		kv_pushp_2(VCF_ANALYSIS, vcf_l, c_va);
		if(! get_one_vcf_nstd(c_r, header, c_va))
			vcf_l->n--;
	}
	qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
	bcf_close(vcf_r.file);	//close
}

void vcf_nstd_dump(char *fn_in, char *fn_out)
{
	VCF_A_L vcf_l = {0};
	FILE * out = xopen(fn_out, "w");
	vcf_nstd_load(fn_in, &vcf_l);
	for(unsigned int i = 0; i < vcf_l.n; i++)
		print_one_vcf_analysis(out, &(vcf_l.a[i]));
	fclose(out);
}

//------------------------------MANTA VCF--------------------------------------//

int get_one_vcf_manta(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
{
	vcf_get_sv_type(header, vcf, vcf_a->type);
	vcf_get_sv_END(header, vcf, &(vcf_a->ed));
	strcpy(vcf_a->ID, vcf->d.id);
	vcf_a->rid = vcf->rid;
	vcf_a->st = vcf->pos;
	return true;
}

void vcf_manta_load(char *fn_in, VCF_A_L *vcf_l)
{
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read
	bcf_hdr_t *header = vcf_r.header;
	char *GT[3]; char GT_1[3];char GT_2[3];char GT_3[3]; GT[0] = GT_1; GT[1] = GT_2; GT[2] = GT_3;
	int nsmpl = bcf_hdr_nsamples(header);
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		VCF_ANALYSIS c_va;
		get_one_vcf_manta(c_r, header, &c_va);

		int N_GT = vcf_get_sv_GT(header, c_r, GT)/nsmpl;
		for(int i = 0; i < N_GT; i++)
		{
			if(GT[i][0] == 4 || GT[i][1] == 4)
			{
				c_va.sample = i;
				kv_push_2(VCF_ANALYSIS, vcf_l, c_va);
			}
		}
	}
	qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
	//close
	bcf_close(vcf_r.file);
}

void vcf_manta_dump(char *fn_in, char *fn_out)
{
	VCF_A_L vcf_l = {0};
	FILE * out = xopen(fn_out, "w");
	vcf_manta_load(fn_in, &vcf_l);
	for(unsigned int i = 0; i < vcf_l.n; i++)
		print_one_vcf_analysis(out, &(vcf_l.a[i]));
	fclose(out);
}

//------------------------------JRLA VCF--------------------------------------//

int get_one_vcf_jlra(bcf1_t *vcf, bcf_hdr_t *header, VCF_ANALYSIS *vcf_a)
{
	char sample[128];
	vcf_get_sample(header, vcf, sample);
	int sample_int = 100;
	if(strcmp(sample, "0") == 0)
		sample_int = 0;
	else if(strcmp(sample, "1") == 0)
		sample_int = 1;
	else if(strcmp(sample, "2") == 0)
		sample_int = 2;
	if(sample_int > 2)
		return false;
	vcf_a->sample = sample_int;
	vcf_get_sv_type(header, vcf, vcf_a->type);
	vcf_get_sv_END(header, vcf, &(vcf_a->ed));
	strcpy(vcf_a->ID, vcf->d.id);
	vcf_a->rid = vcf->rid;
	vcf_a->st = vcf->pos;
	return true;
}

void vcf_jlra_load(char *fn_in, VCF_A_L *vcf_l)
{
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_in);//open for read
	bcf_hdr_t *header = vcf_r.header;
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		VCF_ANALYSIS *c_va = NULL;
		kv_pushp_2(VCF_ANALYSIS, vcf_l, c_va);
		if(! get_one_vcf_jlra(c_r, header, c_va))
			vcf_l->n--;
	}
	qsort(vcf_l->a, vcf_l->n, sizeof(VCF_ANALYSIS), VCF_ANALYSIS_cmp_sample_type_pos);
	//close
	bcf_close(vcf_r.file);
}

void vcf_jlra_dump(char *fn_in, char *fn_out)
{
	VCF_A_L vcf_l = {0};
	FILE * out = xopen(fn_out, "w");
	vcf_jlra_load(fn_in, &vcf_l);
	for(unsigned int i = 0; i < vcf_l.n; i++)
		print_one_vcf_analysis(out, &(vcf_l.a[i]));
	fclose(out);
}

//------------------------------COMPARE DEL--------------------------------------//

//compare two vcf record, return overlap rate
float vcf_compare_one(VCF_ANALYSIS *v1, VCF_ANALYSIS *v2)
{
	int MIN_ED = MIN(v1->ed, v2->ed);
	int MAX_ST = MAX(v1->st, v2->st);
	int overlap = MIN_ED - MAX_ST;
	if(overlap < 0)
		return 0;
	int length1 = v1->ed - v1->st;
	int length2 = v2->ed - v2->st;
	int MAX_LENGTH = MAX(length1, length2);
	return ((float)overlap)/MAX_LENGTH;
}

//to compare data from one sample, one a sv type and one chromosome
void vcf_compare_sample_type_chromosome(VCF_ANALYSIS *vl1, int n1, VCF_ANALYSIS *vl2, int n2,
		FILE* flog, int * overlap_90, int *overlap_50)
{
	int overlap_90_percent = 0;
	int overlap_50_percent = 0;

	//for v1
	for(int i = 0; i < n1; i++)
	{
		VCF_ANALYSIS *v1 = vl1 + i;
		float max_overlap = 0;
		VCF_ANALYSIS c_MAX;
		for(int j = 0; j < n2; j++)
		{
			float overlap = vcf_compare_one(v1, vl2 + j);
			if(overlap > max_overlap)
			{
				max_overlap = overlap;
				memcpy(&c_MAX, vl2 + j, sizeof(VCF_ANALYSIS));
			}
		}
		print_one_vcf_analysis(flog, v1);
		if(max_overlap > 0.01)
		{
			fprintf(flog, "\t\t%f\t", max_overlap);
			print_one_vcf_analysis(flog, &c_MAX);
			if(max_overlap > 0.9)
				overlap_90_percent++;
			if(max_overlap > 0.5)
				overlap_50_percent++;
		}
	}
	*overlap_90 = overlap_90_percent;
	*overlap_50 = overlap_50_percent;
}

void vcf_compare(char *fn_1, char * file_type1, char * fn_2, char *file_type2, char * fn_log)
{
	//open log file
	FILE *flog = xopen(fn_log, "w");
	//open files
	VCF_A_L vcf_l1 = {0};
	VCF_A_L vcf_l2 = {0};
	for(int loop = 0; loop < 2; loop++)
	{
		char *fn = (loop == 0)?fn_1:fn_2;
		char *ft = (loop == 0)?file_type1:file_type2;
		VCF_A_L *vcf_l = (loop == 0)?(&vcf_l1):(&vcf_l2);
		if(strcmp(ft, "jlra") == 0)
			vcf_jlra_load(fn, vcf_l);
		if(strcmp(ft, "manta") == 0)
			vcf_manta_load(fn, vcf_l);
		if(strcmp(ft, "nstd") == 0)
			vcf_nstd_load(fn, vcf_l);
	}
	//for one sample
	int sample_ID = 0;
	unsigned int sample_st = 0, sample_ed = 0;
	unsigned int sample_st2 = 0, sample_ed2 = 0;
	while(1)
	{
		//get st/ed and sample ID for method1
		if(sample_st >= vcf_l1.n)//end of data
			break;
		sample_ID = vcf_l1.a[sample_st].sample;
		for(sample_ed = sample_st; sample_ed < vcf_l1.n && vcf_l1.a[sample_ed].sample == sample_ID; sample_ed++);
		//get st/ed for method2 for sample ID
		int old_st = sample_st2, old_ed = sample_ed2;
		for(; sample_st2 < vcf_l2.n && vcf_l2.a[sample_st2].sample != sample_ID; sample_st2++);
		for(sample_ed2 = sample_st2; sample_ed2 < vcf_l2.n && vcf_l2.a[sample_ed2].sample == sample_ID; sample_ed2++);
		if(sample_st2 == sample_ed2) {sample_st2 = old_st;sample_ed2 = old_ed;}//reset when no results

		//for one type
		char type_ID[32];
		int type_st = sample_st, type_ed = type_st;
		int type_st2 = sample_st2, type_ed2 = type_st2;
		while(1)
		{
			//get st/ed and type ID for method1
			if(type_st >= sample_ed)//end of data
				break;
			strcpy(type_ID, vcf_l1.a[type_st].type);
			for(type_ed = type_st; type_ed < sample_ed && (strcmp(vcf_l1.a[type_ed].type, type_ID) == 0); type_ed++);
			//get st/ed for method2 for type ID
			int old_st = type_st2, old_ed = type_ed2;
			for(; type_st2 < sample_ed2 && (strcmp(vcf_l2.a[type_st2].type, type_ID) != 0); type_st2++);
			for(type_ed2 = type_st2; type_ed2 < sample_ed2 &&  (strcmp(vcf_l2.a[type_ed2].type, type_ID) == 0); type_ed2++);
			if(type_st2 == type_ed2) {type_st2 = old_st;type_ed2 = old_ed;}//reset when no results
			//for one chromosome
			//for one type
			int chr_ID;
			int chr_st = type_st, chr_ed = chr_st;
			int chr_st2 = type_st2, chr_ed2 = chr_st2;
			while(1)
			{
				//get st/ed and type ID for method1
				if(chr_st >= type_ed)//end of data
					break;
				chr_ID = vcf_l1.a[chr_st].rid;
				for(chr_ed = chr_st; chr_ed < type_ed && vcf_l1.a[chr_ed].rid == chr_ID; chr_ed++);
				//get st/ed for method2 for type ID
				int old_st = chr_st2, old_ed = chr_ed2;
				for(; chr_st2 < type_ed2 && vcf_l2.a[chr_st2].rid != chr_ID; chr_st2++);
				for(chr_ed2 = chr_st2; chr_ed2 < type_ed2 && vcf_l2.a[chr_ed2].rid == chr_ID; chr_ed2++);
				if(chr_st2 == chr_ed2) {chr_st2 = old_st;chr_ed2 = old_ed;}//reset when no results
				//compare
				int overlap_90_percent, overlap_50_percent;
				int number_method1 = chr_ed - chr_st;
				int number_method2 = chr_ed2 - chr_st2;
				vcf_compare_sample_type_chromosome(vcf_l1.a + chr_st, number_method1,
						vcf_l2.a + chr_st2, number_method2, flog, &overlap_90_percent, &overlap_50_percent);
				//print results
				fprintf(stderr,
						"sample_ID: %d\t"
						"type_ID:%s\t"
						"chr_ID:%d\t"
						"number_method1:%d\t"
						"number_method2:%d\t"
						"overlap_90_percent:%d\t"
						"overlap_50_percent:%d\t"
						"SEN_90:%f\t"
						"SEN_90:%f\t"
						"\n",
						sample_ID,
						type_ID,
						chr_ID,
						number_method1,
						number_method2,
						overlap_90_percent,
						overlap_50_percent,
						((float)overlap_90_percent)/number_method1,
						((float)overlap_90_percent)/number_method1);
				//end
				chr_st = chr_ed;
				chr_st2 = chr_ed2;
			}
			//end
			type_st = type_ed;
			type_st2 = type_ed2;
		}
		//end
		sample_st = sample_ed;
		sample_st2 = sample_ed2;
	}
	fclose(flog);
}

int getReverseStr(int argc, char *argv[]){
	char * str = argv[1];
	int len = strlen(str);
	getReverseStr_char(str, len);
	fprintf(stderr, "%s\n", str);
	return 0;
}

#define KMER_COUNT_LEN 20
extern uint64_t kmerMask[33];
void simple_kmer_counter(char *fn){
	fprintf(stderr, "version 1.00\n");

	FILE *f_read = xopen(fn, "r");
	char buff_[1000];
	uint8_t buff_bin[1000];
	uint8_t buff_bin_rev[1000];
	std::map<uint64_t, int> kmer_set;
	std::map<uint64_t, int>::iterator kmer_set_it;
	std::map<uint64_t, int> global_kmer_set;

	uint64_t MASK = kmerMask[KMER_COUNT_LEN];
	int read_number = 0;
	while(fgets(buff_, 1000, f_read)){
		read_number++;
		if(read_number % 100 == 0)
			fprintf(stderr, "%d\r", read_number);
		//find \t
		int start_idx = 0;
		while(buff_[start_idx++] != '\t');
		char * str_p = buff_ + start_idx;
		 int string_n = 0;
		for(;str_p[string_n] != '\n';string_n++)
		{
			switch(str_p[string_n]){
			case 'A': buff_bin[string_n] = 0; break;
			case 'C': buff_bin[string_n] = 1; break;
			case 'G': buff_bin[string_n] = 2; break;
			case 'T': buff_bin[string_n] = 3; break;
			default : xassert(0, "");
			}
		}
		//reverse string:
		for(int i = 0; i < string_n; i++){
			buff_bin_rev[string_n - i - 1] = 3 - buff_bin[i];
		}
		uint64_t kmer     = bit2_nextKmer_init(buff_bin, KMER_COUNT_LEN);
		uint64_t kmer_rev = bit2_nextKmer_init(buff_bin_rev, KMER_COUNT_LEN);
		int kmer_number = string_n - KMER_COUNT_LEN + 1;
		kmer_set.clear();
		for(int i = 0; i < kmer_number; i++){
			kmer     = bit2_nextKmerMASK( buff_bin     + i, kmer, KMER_COUNT_LEN);
			kmer_set_it = kmer_set.find(kmer);		if(kmer_set_it!=kmer_set.end()){kmer_set_it->second ++;	}	else{kmer_set[kmer] = 1;	}
		}
		for(int i = 0; i < kmer_number; i++){
			kmer_rev = bit2_nextKmerMASK( buff_bin_rev + i, kmer_rev, KMER_COUNT_LEN);
			kmer_set_it = kmer_set.find(kmer_rev);	if(kmer_set_it!=kmer_set.end()){kmer_set_it->second ++;	}	else{kmer_set[kmer_rev] = 1;}
		}
		kmer     = bit2_nextKmer_init(buff_bin, KMER_COUNT_LEN);
		kmer_rev = bit2_nextKmer_init(buff_bin_rev, KMER_COUNT_LEN);

		for(kmer_set_it = kmer_set.begin(); kmer_set_it != kmer_set.end(); kmer_set_it++){
			if(kmer_set_it->second >= 3){
				auto g_it = global_kmer_set.find(kmer_set_it->first);
				if(g_it!=global_kmer_set.end()){g_it->second += kmer_set_it->second;}	else{global_kmer_set[kmer_set_it->first] = kmer_set_it->second;	}
			}
		}
	}

	for(auto it = global_kmer_set.begin(); it != global_kmer_set.end(); it++){
		char kmerStr[100];
		uint64_t kmer = it->first;
		for(int i = 0; i < KMER_COUNT_LEN; i++){
			kmerStr[i] = "ACGT"[(kmer >> ((KMER_COUNT_LEN - 1 - i) * 2)) & 0x3];
		}
		kmerStr[KMER_COUNT_LEN] = 0;
		fprintf(stdout, "%d\t%s\n", it->second, kmerStr);
	}

	fclose(f_read);
}

static int gz_head(int argc, char *argv[]){
	const char * fn = argv[1];
	int total_read_len = atoi(argv[2]);
	int begin_offset = atoi(argv[3]);

	gzFile g = xzopen(fn, "rb");
	char *buff = (char *)xmalloc(total_read_len + 1);

	//文件头0(SEEK_SET)，当前位置1(SEEK_CUR)，文件尾2(SEEK_END)）
	gzseek(g, begin_offset, 0);

	err_gzread(g, buff, total_read_len);
	fprintf(stdout, "%s", buff);
	free(buff);
	gzclose(g);
	return 0;
}

int read_ACGT_analysis(int argc, char *argv[]){
	const char * input_bam_fn = argv[1]; const char * ref_fn = argv[2];
	Bam_file c_b;
	memset(&c_b, 0, sizeof(Bam_file));
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
	bam_hdr_t* hdr = c_b._hdr;

	bam1_t b = {0};//BAM record for the first read in a pair

	char * seq_buff = (char *)xmalloc(10000);
	int * analysis_cnt = (int *)xcalloc(10000, sizeof(int));
	int read_ID = 0;
	while (sam_read1(c_b._hfp, hdr, &b) >= 0){

		const int read_len = b.core.l_qseq;
		get_bam_seq(0, read_len, seq_buff, &b);//store in binary format
		int A_count = 0;
		int C_count = 0;
		int G_count = 0;
		int T_count = 0;
		int N_count = 0;
		for(int i = 0; i < read_len; i++){
			switch(seq_buff[i])
			{
			case 'A': case 'a': A_count++; break;
			case 'C': case 'c': C_count++; break;
			case 'G': case 'g': G_count++; break;
			case 'T': case 't': T_count++; break;
			default: N_count++; break;
			}
		}

		if(N_count > 25){
			fprintf(stderr, "[N_count %d @ %d]%s\n", N_count, read_ID ,seq_buff );
		}

//		if(A_count >= 75 || C_count >= 75 || G_count >= 75 || T_count >= 75){
//			fprintf(stdout, "[%d %d %d %d]@ %d %s\n", A_count, C_count, G_count, T_count, read_ID, seq_buff );
//		}

		if((G_count >= 75 && C_count < 20) || (C_count >= 75 && G_count < 20)){
			fprintf(stdout, "[%d %d %d %d]@ ID %d @ pos %d:%d %s\n", A_count, C_count, G_count, T_count, read_ID,b.core.tid, b.core.pos, seq_buff );
		}

		A_count = A_count / 16;
		C_count = C_count / 16;
		G_count = G_count / 16;
		T_count = T_count / 16;

		int count_number = 0;
		xassert(A_count < 10, ""); count_number*= 10; count_number += A_count;
		xassert(C_count < 10, ""); count_number*= 10; count_number += C_count;
		xassert(G_count < 10, ""); count_number*= 10; count_number += G_count;
		xassert(T_count < 10, ""); count_number*= 10; count_number += T_count;

		analysis_cnt[count_number] ++;
		if(read_ID % 1000000 == 0){
			for(int i = 0; i < 10000; i++){
				if(analysis_cnt[i] > 0){
					fprintf(stderr, "[%d %d]\n", i, analysis_cnt[i] );
				}
			}
		}
		read_ID++;
	}

	bam_file_close(&c_b);

	return 0;
}

//random generate 10000 deletion + 10000 insertion + 500 duplication
int randomGenerateSV(int argc, char *argv[]){
	char * header_fn = argv[1];
	const char * ref_fn = argv[2];
	int seed_random = atoi(argv[3]);
	srand(seed_random);

	gzFile fp = xzopen(ref_fn, "r");
	kstream_t *_fp = ks_init(fp);

	kseq_t ref_seq = {0};
	ref_seq.f = _fp;

	//load header file
	htsFile *header_file = hts_open(header_fn, "r");//open output file
	bam_hdr_t *header = sam_hdr_read(header_file);
	hts_close(header_file);

	int64_t total_length = 0;
	int32_t n_targets = header->n_targets;
	n_targets = MIN(23, n_targets);


	fprintf(stdout, "##fileformat=VCFv4.2\n");
	fprintf(stdout, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	fprintf(stdout, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n");
	fprintf(stdout, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
	fprintf(stdout, "##INFO=<ID=SIM_INS_ST,Number=.,Type=Integer,Description=\"Random insertion string original position.(compared with POS)\">\n");
	fprintf(stdout, "##INFO=<ID=SVANN,Number=.,Type=String,Description=\"Repeat annotation of structural variant\">\n");
	fprintf(stdout, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n");
	fprintf(stdout, "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n");
	fprintf(stdout, "##INFO=<ID=MATEDIST,Number=1,Type=Integer,Description=\"Distance to the mate breakend for mates on the same contig\">\n");
	fprintf(stdout, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(stdout, "##INFO=<ID=SHADOWED,Number=0,Type=Flag,Description=\"CNV overlaps with or is encapsulated by deletion\">\n");
	fprintf(stdout, "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(stdout, "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(stdout, "##ALT=<ID=CNV,Description=\"Copy number variable region\">\n");
	fprintf(stdout, "##FILTER=<ID=Decoy,Description=\"Variant involves a decoy sequence\">\n");
	fprintf(stdout, "##FILTER=<ID=NearReferenceGap,Description=\"Variant is near (< 1000 bp) from a gap (run of >= 50 Ns) in the reference assembly\">\n");
	fprintf(stdout, "##FILTER=<ID=NearContigEnd,Description=\"Variant is near (< 1000 bp) from the end of a contig\">\n");
	fprintf(stdout, "##FILTER=<ID=InsufficientStrandEvidence,Description=\"Variant has insufficient number of reads per strand (< 0).\">\n");
	fprintf(stdout, "##FILTER=<ID=NotFullySpanned,Description=\"Duplication variant does not have any fully spanning reads.\">\n");
	fprintf(stdout, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(stdout, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth per allele\">\n");
	fprintf(stdout, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">\n");
	fprintf(stdout, "##FORMAT=<ID=SAC,Number=.,Type=Integer,Description=\"Number of reads on the forward and reverse strand supporting each allele including reference\">\n");
	fprintf(stdout, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n");
	fprintf(stdout, "##reference=file:%s\n", ref_fn);
	for(int i = 0; i < n_targets; i++){
		fprintf(stdout, "##contig=<ID=%s,length=%d>\n", header->target_name[i], header->target_len[i]);
	}
	fprintf(stdout, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	demo_SAMPLE\n");

	for(int i = 0; i < n_targets; i++){
		total_length += header->target_len[i];
	}

	int del_ID = 0;
	int ins_ID = 0;
	for(int chr_ID = 0; chr_ID < n_targets; chr_ID++){
		//load reference
		xassert( kseq_read(&ref_seq) >= 0, "");

		int64_t chr_len = header->target_len[chr_ID];
		int del_num = ((10000 * chr_len) / total_length) + 1;
		int ins_num = ((10000 * chr_len) / total_length) + 1;
		//int dup_num = (500 * chr_len / total_length) + 1;

		//generate deletion
		for(int sv_ID = 0; sv_ID < del_num; sv_ID++){
			int st_pos = rand() % chr_len;
			int length = 0;
			while(rand() % 3 > 0){ //66% X10
				length += (rand() % 9) + 1; length *= 10;
			}
			if(st_pos + length + 1 > chr_len || length < 50 || length > 3000){ sv_ID--; continue; }
			if(ref_seq.seq.s[st_pos] == 'N') { sv_ID--; continue; }
			fprintf(stdout, "%s\t%d\trandom.DEL.%d\t", header->target_name[chr_ID], st_pos, del_ID++);
			for(int i = 0; i < length + 1; i++){
				fprintf(stdout, "%c", ref_seq.seq.s[st_pos + i - 1]);
			}
			fprintf(stdout, "\t%c\t.\tPASS\tSVTYPE=DEL;END=%d;SVLEN=%d\tGT:AD:DP\t1/1:0,31:31\n", ref_seq.seq.s[st_pos - 1], st_pos + length, length);
		}
		//INS:
		for(int sv_ID = 0; sv_ID < ins_num; sv_ID++){
			int st_pos = rand() % chr_len;
			int length = 0;
			while(rand() % 3 > 0){ 	length += (rand() % 9) + 1; length *= 10; }//66% X10
			if(st_pos + length > chr_len || length < 50 || length > 3000){ sv_ID--; continue; }
			xassert(length < 3000, "");
			if(ref_seq.seq.s[st_pos] == 'N') { sv_ID--; continue; }
			fprintf(stdout, "%s\t%d\trandom.INS.%d\t", header->target_name[chr_ID], st_pos , ins_ID++);
			fprintf(stdout, "%c\t%c", ref_seq.seq.s[st_pos - 1], ref_seq.seq.s[st_pos - 1]);
			int insert_seq_st = 0;
			do{
				int levle = 0;
				while(rand() % 3 > 0 && levle++ < 6){ 	insert_seq_st += (rand() % 9) + 1; insert_seq_st *= 10; }//66% X10
				if(rand()%2 == 0)
					insert_seq_st = -insert_seq_st;
			}while(insert_seq_st == 0 || st_pos + insert_seq_st < 0 || st_pos + insert_seq_st + length > chr_len || ref_seq.seq.s[st_pos + insert_seq_st] == 'N');

			for(int i = 0; i < length; i++){
				fprintf(stdout, "%c", ref_seq.seq.s[st_pos + insert_seq_st + i - 1]);
			}
			fprintf(stdout, "\t.\tPASS\tSVTYPE=INS;END=%d;SIM_INS_ST=%d;SVLEN=%d\tGT:AD:DP\t1/1:0,31:31\n", st_pos, insert_seq_st, length);
		}
	}

	return 0;
}

struct VCF_COM_Record{
	std::string sample;
	bcf1_t 		r;
	VCF_COM_Record(){
		memset(&r, 0, sizeof(bcf1_t));
	}
	static inline int cmp_by_pos(const VCF_COM_Record &a, const VCF_COM_Record &b){
		//var basic
		if(a.r.rid != b.r.rid) 	return a.r.rid < b.r.rid;
		if(a.r.pos != b.r.pos) 	return a.r.pos < b.r.pos;
		return a.sample.compare(b.sample);
	}
};

int combine_sort_vcf(int argc, char *argv[]){
	char * vcf_fn_in = argv[1];//separate by ','
	char * vcf_fn_out = argv[2];

	//get bam file list
	std::vector<std::string> vcf_files_names;
	split_string(vcf_files_names, vcf_fn_in, ",");
	std::vector<BCF_FILE> bam_list;
	for(std::string &bam_fn:  vcf_files_names){
		bam_list.emplace_back();
		memset(&bam_list.back(), 0, sizeof(BCF_FILE));
		VCF_open_read(&bam_list.back(),bam_fn.c_str());
		//bcf_hdr_append(bam_list.back().header, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
	}

	std::vector<VCF_COM_Record> vcf_list;
	uint32_t load_size = 0;

	for(BCF_FILE & bcf_f: bam_list){
		std::string sample_name(bcf_f.header->samples[0]);
		do{
			if(vcf_list.size() < load_size + 1)
				vcf_list.emplace_back();
			 vcf_list[load_size].sample = sample_name;
		}while(VCF_next_dump(&bcf_f, &(vcf_list[load_size].r)) && vcf_list[load_size].r.rid == 0 && ++load_size);
	}

	std::sort(vcf_list.begin(), vcf_list.end(), VCF_COM_Record::cmp_by_pos);

	//write::
	//open write file
	BCF_FILE vcf_out;
	VCF_open_write(&vcf_out, vcf_fn_out, false);

	bcf_hdr_t *write_header = bam_list[0].header;
	//bcf_hdr_append(write_header, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
	//fprintf(stderr, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample Strings\">\n");
    vcf_hdr_write(vcf_out.file, write_header);

    //int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);

	//append sample info for each line:
	for(VCF_COM_Record & vcf_r: vcf_list){
	    //int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);
		if ( !(vcf_r.r.unpacked & BCF_UN_ALL) ) bcf_unpack(&vcf_r.r, BCF_UN_ALL);
		bcf_add_id(write_header, &vcf_r.r, vcf_r.sample.c_str());
		//bcf_update_info(write_header, &vcf_r.r, "SAMPLE",vcf_r.sample.c_str() , vcf_r.sample.size(), BCF_HT_STR);
	    vcf_write(vcf_out.file, write_header, &vcf_r.r);
	}
	//close::
	for(BCF_FILE & bcf_f: bam_list){
		bcf_close(bcf_f.file);
	}
	return 0;
}

int analysis_ROC_PR(int argc, char *argv[]){
	struct CVS_item{
		int label;
		float score;
		CVS_item(int label, float score){
			this->label = label;
			this->score = score;
		}
		static inline int cmp_by_score(const CVS_item &a, const CVS_item &b){
			//var basic
			return a.score > b.score;
		}
	};

	char * csv_fn_in = argv[1];//separate by ','
	//char * vcf_fn_out = argv[2];
	//load data:
	std::vector<std::string> csv_data_str;
	load_string_list_from_file(csv_fn_in, csv_data_str);
	size_t line_num = csv_data_str.size() - 1;//skip the header line
	std::vector<CVS_item> cvs_l;
	for(uint64_t i = 0; i < line_num; i++){
		std::vector<std::string> item_value;
		split_string(item_value, csv_data_str[i + 1].c_str(), ",");//skip the header line
		cvs_l.emplace_back(atoi(item_value[1].c_str()), atof(item_value[2].c_str()));
	}

	//sort the cvs_l
	std::sort(cvs_l.begin(), cvs_l.end(), CVS_item::cmp_by_score);

	if(false){
		for(uint64_t i = 0; i < line_num; i++)
			fprintf(stderr, "%d ",cvs_l[i].label);
		fprintf(stderr, "\n");

		for(uint64_t i = 0; i < line_num; i++)
			fprintf(stderr, "%f ",cvs_l[i].score);
		fprintf(stderr, "\n");
	}

//	struct RP_item{
//		float precision;
//		float recall;
//		float thresholds;
//	};
//	std::vector<RP_item> rp_l;
	fprintf(stdout, "SEN,PRE,THRED\n");
	{
		int total_base_number = 0;
		for(uint64_t i = 0; i < line_num; i++){
			if(cvs_l[i].label == 1){
				total_base_number++;
			}
		}
		//fprintf(stdout, "0,1,10000000\n");
		for(uint64_t i = 0; i < line_num; i++){
			if(i < line_num - 1  && cvs_l[i].score == cvs_l[i+1].score){//skip duplication score
				continue;
			}
			int TP = 0;
			int FP = 0;
			for(uint64_t j = 0; j <= i; j++){
				if(cvs_l[j].label == 1){TP++;}
				else				   {FP++;}
			}
			if(i == line_num - 1)
				FP=1000000000;
			float SEN = (float)TP/total_base_number;
			float PRE = (float)TP/(TP+FP);
			float THRED=cvs_l[i].score;
			if(PRE > 0.0001)
				fprintf(stdout, "%f,%f,%f\n", SEN, PRE, THRED);
		}
	}
	return 0;
}

void setR_Region(std::vector<R_region>& SV_region_l, int rid, int st, int end){
	SV_region_l.emplace_back();
	SV_region_l.back().chr_ID = rid;
	SV_region_l.back().st_pos = st;
	SV_region_l.back().ed_pos = end;
}

struct SV_basic_infomation{
	int pos_begin;
	std::string SVType;
	int SV_length;
	int maxVNTR_value;
	int isHIGH;

	int TGS_CIGAR_SIG_N;
	int TGS_CLIP_SIG_N;

	void store_SV_Data(char * type, int len){
		this->SV_length = len;
		SVType.clear();
		SVType.append(type);
	}
	void storemaxVNTR_value(int maxVNTR_value){
		this->maxVNTR_value = MAX(this->maxVNTR_value, maxVNTR_value);
	}
};

//#define ANALYSIS_BLOCK_SIZE 500

struct REF_BASIC_DATA{
	uint32_t ref_length;
	uint8_t * mappbility_value;
	const char * chrName;
	int chrID;
	int centromeres_region_bg;
	int centromeres_region_ed;
	int total_block_length = 0;
	std::vector<SV_basic_infomation> region_withSV;

	REF_BASIC_DATA(int chrID, uint32_t ref_length, const char * chrName){
		this->ref_length = ref_length;
		this->mappbility_value = (uint8_t * )xcalloc(this->ref_length, 1);
		this->chrName = chrName;
		centromeres_region_bg = 0;
		centromeres_region_ed = 0;
		this->chrID = chrID;
	}

	REF_BASIC_DATA(){
		this->ref_length = 0;
		this->mappbility_value = NULL;
		this->chrName = NULL;
		centromeres_region_bg = 0;
		centromeres_region_ed = 0;
		this->chrID = 0;
	}

	void initRegionWithSV(int ana_block_step_size){
		int total_block_length = ref_length/ana_block_step_size + 1;
		region_withSV.resize(total_block_length);
		int position_bg_cur = 0;
		for(SV_basic_infomation & r : region_withSV){
			r.pos_begin = position_bg_cur;
			position_bg_cur += ana_block_step_size;
			r.maxVNTR_value = 0;
			r.isHIGH = false;
			//init TGS counter
			r.TGS_CIGAR_SIG_N = 0;
			r.TGS_CLIP_SIG_N = 0;
		}
	}
};

#define MIN_COMPACT_LEN 8
bool pass_compact_filter(uint8_t *s, const int len) {
	int compact_str_len = 1; //AAAACCCT --> ACT
	for (int i = 1; i < len; i++)
		if (s[i - 1] != s[i])
			compact_str_len++;
	if (compact_str_len < MIN_COMPACT_LEN)
		return false;
	return true;
}

//when a clip string clip to an "AAAAAAAA..." tail, return false
bool clip_AAA_TailFilter(uint8_t *read_bin, const int read_len, int soft_left, int soft_right) {
	int middle_len = read_len - soft_left - soft_right;
	if(middle_len < 20) return false;
	//left check:
	if(soft_left > 0){
		int base_ACGT_left[4] = { 0 };
		uint8_t * st_base = read_bin + soft_left;
		int MAX_SAME = 1;
		for (int i = 1; i < 20; i++){
			if(st_base[i] == st_base[i - 1])
				MAX_SAME ++;
			else
				MAX_SAME = 1;
		}
		for (int i = 0; i < 20; i++)
			base_ACGT_left[st_base[i]]++;
		if((base_ACGT_left[0] > 15 || base_ACGT_left[3] > 15) && MAX_SAME > 9)
			return false;
	}
	if(soft_right > 0){
		//right check:
		int base_ACGT_right[4] = { 0 };
		uint8_t * st_base = read_bin + read_len - soft_right - 20;
		int MAX_SAME = 1;
		for (int i = 1; i < 20; i++){
			if(st_base[i] == st_base[i - 1])
				MAX_SAME ++;
			else
				MAX_SAME = 1;
		}
		for (int i = 0; i < 20; i++)
			base_ACGT_right[st_base[i]]++;
		if((base_ACGT_right[0] > 15 || base_ACGT_right[3] > 15) && MAX_SAME > 9)
			return false;
	}
	return true;
}

float get_ave_quality_value1(int st_pos, int end_pos, bam1_t* _bp)
{
	uint32_t total_qual = 0;
	uint8_t* bam_quality = bam_get_qual(_bp);
	end_pos = MIN(end_pos, _bp->core.l_qseq);
	if(st_pos > end_pos){ return -100; }
	for(int i = st_pos; i <= end_pos; i++)
		total_qual += bam_quality[i];
	return (float)total_qual / (end_pos - st_pos + 1);
}

bool clip_low_quality_Filter(bam1_t *br, int read_len, int soft_left, int soft_right) {
	float mapped_part_ave_quality_value = get_ave_quality_value1(soft_left, read_len - soft_right, br);
	float clip_part_ave_quality_value_left;
	float clip_part_ave_quality_value_right;
	if(soft_left > 0){
		clip_part_ave_quality_value_left = get_ave_quality_value1(0, soft_left, br);
		if(clip_part_ave_quality_value_left < 0.6 * mapped_part_ave_quality_value)
			soft_left = 0;
	}
	if(soft_right > 0){
		clip_part_ave_quality_value_right = get_ave_quality_value1(read_len - soft_right - 1, read_len - 1, br);
		if(clip_part_ave_quality_value_right < 0.6 * mapped_part_ave_quality_value)
			soft_right = 0;
	}
	if(soft_left > 0 || soft_right > 0)
		return true;
	return false;
}

void get_block_ID_by_pos(int & id1, int &id2, int pos, int ana_block_length, int ana_block_step_size){
	id1 = (pos - ana_block_length + ana_block_step_size)/ana_block_step_size;
	id1 = MAX(0,id1);
	id2 = pos/ana_block_step_size;
}

void init_and_load_local_repeat_data(char * input_VNTR_region_file, std::vector<REF_BASIC_DATA> &reference_list,
		int & ana_block_length, int & ana_block_step_size, faidx_t *c_ref_idx){
	//load VNTR repeat level
	std::vector<std::string> VNTR_bed_l;
	load_string_list_from_file(input_VNTR_region_file, VNTR_bed_l);
	std::vector<std::string> item_value;

	bool isHeader = true;
	for(std::string &s :VNTR_bed_l){
		split_string(item_value, s.c_str(), "\t");
		if(isHeader){
			ana_block_length = atoi(item_value[0].c_str());
			ana_block_step_size = atoi(item_value[1].c_str());
			isHeader = false;
			{
				std::vector<std::string> item_value;
				//load reference index
				int N_seq = faidx_nseq(c_ref_idx);
				N_seq = MIN(24, N_seq);
				for(int i = 0; i < N_seq; i++){
					char fn[1024];
					sprintf(fn, "%d_dump_file.bin", i);
					reference_list.emplace_back();
					reference_list.back().chrID = i;
					reference_list.back().ref_length = faidx_seq_len(c_ref_idx, faidx_iseq(c_ref_idx, i));
					reference_list.back().centromeres_region_bg = 0;
					reference_list.back().centromeres_region_ed = 0;
					reference_list.back().initRegionWithSV(ana_block_step_size);
				}
			}
		}
		else{
			int chrID = atoi(item_value[0].c_str());
			int pos = atoi(item_value[1].c_str());
			int maxVNTR_value = atoi(item_value[2].c_str());
			if(chrID <= 23){
				reference_list[chrID].region_withSV[pos/ana_block_step_size].storemaxVNTR_value(maxVNTR_value);
			}
		}
	}
}

void load_vcf_data_file(std::vector<REF_BASIC_DATA> &reference_list, char * vcf_fn_in, int ana_block_length, int ana_block_step_size){
	//open vcf files
	BCF_FILE vcf_read;//vcf for read
	VCF_open_read(&vcf_read, vcf_fn_in);//open for read
	char *c_sv_type = (char *)malloc(1000);
	int SV_CHR_ID;
	int SV_POS;
	int SV_END;
	int SV_length;
	do{
		bcf1_t *c_r = &( vcf_read.r);
		if(c_r->d.flt != NULL)
			*c_r->d.flt = 0;
		//unpack the vcf data to get the Filter
		bcf_unpack(c_r, BCF_UN_INFO);
		//vcf filters
		if(c_r->d.flt != NULL && *c_r->d.flt != 0)//filter: PASS
			continue;

		if(c_r->rid > 23)
			continue;

		SV_CHR_ID = c_r->rid;
		SV_POS = c_r->pos;
		fprintf(stderr, "%d\t%d\n",SV_CHR_ID,SV_POS);
		vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
		vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
		vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

		if(SV_length < 50)
			continue;

		if(strcmp(c_sv_type, "INS") == 0)
			SV_END = SV_POS;
		else
			SV_END = SV_POS + SV_length;

		if(strcmp(c_sv_type, "BND") == 0)
			continue;

		std::vector<int> SV_BP_pos_l;

		if(ABS_U(SV_POS, SV_END) < 10){
			int SV_POS_middle = (SV_POS + SV_END) / 2;
			SV_BP_pos_l.emplace_back(SV_POS_middle);
		}else{
			SV_BP_pos_l.emplace_back(SV_POS);
			SV_BP_pos_l.emplace_back(SV_END);
		}
		//
		for(int SV_BP_pos:SV_BP_pos_l){
			//get SV reion for this BP
			int store_block_id1 = 0, store_block_id2 = 0;
			get_block_ID_by_pos(store_block_id1, store_block_id2, SV_BP_pos, ana_block_length, ana_block_step_size);
			//store data:
			reference_list[SV_CHR_ID].region_withSV[store_block_id1].store_SV_Data(c_sv_type, SV_length);
			reference_list[SV_CHR_ID].region_withSV[store_block_id2].store_SV_Data(c_sv_type, SV_length);

		}
	}
	while(VCF_next(&vcf_read));//read one
	bcf_close(vcf_read.file);
}

void load_HIGH_CON_bed(std::vector<REF_BASIC_DATA> &reference_list, char * HIGH_CONFIDENCE_region_file, faidx_t *c_ref_idx, int ana_block_step_size){
	std::vector<std::string> HIGH_CON_bed_l;
	load_string_list_from_file(HIGH_CONFIDENCE_region_file, HIGH_CON_bed_l);
	std::vector<std::string> item_value;
	for(std::string &s :HIGH_CON_bed_l){
		split_string(item_value, s.c_str(), "\t");
		const char *chrNAME = item_value[0].c_str();
		int pos_bg = atoi(item_value[1].c_str());
		int pos_ed = atoi(item_value[2].c_str());
		int chrID = faidx_get_chrID(c_ref_idx, chrNAME, NULL, 0);
		if(chrID <= 23){
			for(int pos = pos_bg; pos < pos_ed; pos+=100)
				reference_list[chrID].region_withSV[pos/ana_block_step_size].isHIGH = true;
		}
	}
}

int genomic_contest_AnalysisCCS(int argc, char *argv[]){
	//parameters
	char * ref_fn = argv[1];
	char * vcf_fn_in = argv[2];//separate by ','
	char * input_bam_fn = argv[3];
	char * input_VNTR_region_file = argv[4];
	char * HIGH_CONFIDENCE_region_file = argv[5];
	bool used_BLANK = (atoi(argv[6]) == 1);

	//open bam file
	Bam_file c_b;
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);

	std::vector<REF_BASIC_DATA> reference_list;
	faidx_t *c_ref_idx = reference_index_load(ref_fn);

	int ana_block_length;
	int ana_block_step_size;
	init_and_load_local_repeat_data(input_VNTR_region_file, reference_list, ana_block_length, ana_block_step_size, c_ref_idx);
	load_vcf_data_file(reference_list, vcf_fn_in,ana_block_length, ana_block_step_size);
	load_HIGH_CON_bed(reference_list, HIGH_CONFIDENCE_region_file, c_ref_idx, ana_block_step_size);

	//load read signals from reads
	std::vector<int> mismatch_position;
	//store data:
	bool print_log = false;
	int MIN_sv_len = 30;
	path_segment* path = (path_segment* )calloc(10000, sizeof(path_segment));
	{
		//set the region:
		int read_ID = 0;
		while (bam_next(&c_b)) {
			bam1_t *br = &(c_b._brec);
			read_ID ++;
			if(bam_is_secondary(br))
				continue;
			if(bam_is_supplementary(br))
				continue;
			if (br->core.qual < 4)
				continue;
			if(read_ID % 10000 == 0){
				fprintf(stderr, "%d\t%d\t%d\n", read_ID, br->core.tid, br->core.pos);
			}
			//if(br->core.tid > 0 || br->core.pos > 1000000 )				break;
			if(br->core.tid > 23)
				break;

			int seq_i = 0;
			int ref_i = 0;
			int cigar_idx = 0;
			int n_cigar = br->core.n_cigar;
			uint32_t* bam_cigar = bam_get_cigar(br);
			for (unsigned int i = 0; i < n_cigar; ++i){
				path[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
				path[i].type   = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			}

			//uint8_t *qseq = cur_read.storeReadBuff;
			for(int i = 0; i < n_cigar; i++){
				path_segment *p = path + i;
				if(p->type == align_t::CIGAR_SOFT_CLIP)//soft clip:
				{
					if(p->length > 1000){//at lease 1000 bp
						if(print_log) fprintf(stderr, "CLIP_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						int tid = br->core.tid;
						int breakpoint = br->core.pos + ref_i;
						//store CLIP signals
						int store_block_id1 = 0, store_block_id2 = 0;
						get_block_ID_by_pos(store_block_id1, store_block_id2, breakpoint, ana_block_length, ana_block_step_size);
						reference_list[tid].region_withSV[store_block_id1].TGS_CLIP_SIG_N++;
						if(store_block_id2 != store_block_id1){
							reference_list[tid].region_withSV[store_block_id2].TGS_CLIP_SIG_N++;
						}
					}
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_INSERT){
					if(p->length >= (unsigned)MIN_sv_len){
						if(print_log) fprintf(stderr, "BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						int tid = br->core.tid;
						int breakpoint = br->core.pos + ref_i;
						int store_block_id1 = 0, store_block_id2 = 0;
						get_block_ID_by_pos(store_block_id1, store_block_id2, breakpoint, ana_block_length, ana_block_step_size);
						reference_list[tid].region_withSV[store_block_id1].TGS_CIGAR_SIG_N++;
						if(store_block_id2 != store_block_id1){
							reference_list[tid].region_withSV[store_block_id2].TGS_CIGAR_SIG_N++;
						}
					}else{
						//store small sigs
						//small_sig_queue_handler(print_log, br->core.tid, small_sig, br->core.pos + ref_i, p->length, store_sig);
					}
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_DELETE)
				{
					if(p->length >= (unsigned)MIN_sv_len){
						if(print_log) fprintf(stderr, "BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid, br->core.pos + ref_i, p->length, pathSTR[p->type], bam_get_qname(br));
						int tid = br->core.tid;
						{
							int breakpoint = br->core.pos + ref_i;
							int store_block_id1 = 0, store_block_id2 = 0;
							get_block_ID_by_pos(store_block_id1, store_block_id2, breakpoint, ana_block_length, ana_block_step_size);
							reference_list[tid].region_withSV[store_block_id1].TGS_CIGAR_SIG_N++;
							if(store_block_id2 != store_block_id1){
								reference_list[tid].region_withSV[store_block_id2].TGS_CIGAR_SIG_N++;
							}
						}
						{
							int breakpoint = br->core.pos + ref_i + p->length;
							int store_block_id1 = 0, store_block_id2 = 0;
							get_block_ID_by_pos(store_block_id1, store_block_id2, breakpoint, ana_block_length, ana_block_step_size);
							reference_list[tid].region_withSV[store_block_id1].TGS_CIGAR_SIG_N++;
							if(store_block_id2 != store_block_id1){
								reference_list[tid].region_withSV[store_block_id2].TGS_CIGAR_SIG_N++;
							}
						}
					}else{
						//store small sigs
						//small_sig_queue_handler(print_log, br->core.tid, small_sig, br->core.pos + ref_i, p->length, store_sig);
					}
					ref_i += p->length;
				}
				else if(p->type == align_t::CIGAR_MATCH
						|| p->type == align_t::CIGAR_SEQ_MATCH
						|| p->type == align_t::CIGAR_SEQ_MISMATCH){
					//analysis the mismatch
					//small_sig_queue_handler(print_log, small_sig, POS, 0);
					ref_i += p->length;
					seq_i += p->length;
				}
				else if(p->type == align_t::CIGAR_HARD_CLIP){
					//DO nothing
				}
				else{
					//DO nothing
				}
			}
		}
	}

	{
		int region_ID = -1;
		std::vector<std::string> item_value;
		//load reference index
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int i = 0; i < N_seq; i++){
			REF_BASIC_DATA & u = reference_list[i];
			for(SV_basic_infomation & sr: u.region_withSV){
				region_ID++;
				//select regions:
				if(!sr.SVType.empty() || (used_BLANK)){//SV or using BLANK
					//output:
					int st_pos = sr.pos_begin - 50;
					int ed_pos = sr.pos_begin + ana_block_length + 50;
					// analysis the reads
					// out put information for that region
					//SV INFO
					bool blank_site =(sr.SVType.empty());
					if(blank_site)
						sr.SVType.append("BLANK");
					if(blank_site)
						printf("!0\tBLANK\t");
					else
						printf("!0\tSV_IN\t");
					//2
					printf("!1\t%d\t%d\t%s\t%d\t", u.chrID, sr.pos_begin, sr.SVType.c_str(), region_ID);
					//REGION:7
					printf("!2\t%d\t%d\t", st_pos, ed_pos);
					//10
					printf("!3\t%d\t%d\t", sr.maxVNTR_value, sr.isHIGH);
					//BAISC:14
					printf("!4\tCIGAR\t%d\tCLIP\t%d\t", sr.TGS_CIGAR_SIG_N, sr.TGS_CLIP_SIG_N);
					//16
					printf("\n");
				}
			}
		}
	}

	//close files
	bam_file_close(&c_b);
	return 0;
}

int genomic_contest_AnalysisNGS(int argc, char *argv[]){
	//parameters
	char * ref_fn = argv[1];
	char * vcf_fn_in = argv[2];//separate by ','
	char * input_bam_fn = argv[3];
	char * input_VNTR_region_file = argv[4];
	char * HIGH_CONFIDENCE_region_file = argv[5];
	bool used_BLANK = (atoi(argv[6]) == 1);

	std::vector<REF_BASIC_DATA> reference_list;
	faidx_t *c_ref_idx = reference_index_load(ref_fn);

	int ana_block_length;
	int ana_block_step_size;
	init_and_load_local_repeat_data(input_VNTR_region_file, reference_list, ana_block_length, ana_block_step_size, c_ref_idx);
	load_vcf_data_file(reference_list, vcf_fn_in,ana_block_length, ana_block_step_size);
	load_HIGH_CON_bed(reference_list, HIGH_CONFIDENCE_region_file, c_ref_idx,ana_block_step_size);
	//open vcf files
	Bam_file c_b;
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);

	uint8_t read_str_bin[512];
  	std::vector<int> mismatch_position;
	//store data:
	{
		int region_ID = -1;
		std::vector<std::string> item_value;
		//load reference index
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int chrID = 0; chrID < N_seq; chrID++){
			int ref_len = 0;
			char * ref_char = fai_fetch(c_ref_idx, faidx_iseq(c_ref_idx, chrID), &ref_len);
			REF_BASIC_DATA & u = reference_list[chrID];
			for(SV_basic_infomation & sr: u.region_withSV){
				region_ID++;
				//SV or using BLANK
				if(!sr.SVType.empty() || (used_BLANK)){
					//SV or using BLANK
					R_region analysis_region;
					analysis_region.chr_ID = u.chrID;

					analysis_region.st_pos = sr.pos_begin - 50;
					analysis_region.ed_pos = sr.pos_begin + ana_block_length + 50;
					analysis_region.st_pos = MAX(0, analysis_region.st_pos);
					resetRegion_ID(&c_b, &analysis_region);
					//reset calculators
					uint32_t readNum = 0;
					uint64_t total_readLen = 0;
					//clip signals
					uint32_t total_clip_len = 0;
					uint32_t readNum_clip = 0;
					uint32_t readNum_clip_0_10 = 0;
					uint32_t readNum_clip_11_50 = 0;
					uint32_t readNum_clip_over_50 = 0;
					//DRP signals
					uint32_t DRP_read_num = 0;//both forward:
					uint32_t insertSizeFF_read_num = 0;//both forward:
					uint32_t insertSizeRR_read_num = 0;//both reverse:
					uint32_t insertSizeRF_read_num = 0;//the read of smaller position is reverse, the other is forward
					uint32_t insertSizeOver1K_read_num = 0;//the read of insert size > 1023
					uint32_t insertSizeOver100K_read_num = 0;//the read of insert size > 10000
					uint32_t mate_unmapped_read_num = 0;
					uint32_t mate_different_chromsome_read_num = 0;
					//NM signals
					//SNP and INDEL
					uint32_t total_INDEL_len = 0;
					uint32_t total_NM_len = 0;
					uint32_t NM_over2_read_num = 0;
					uint32_t NM_over10_read_num = 0;
					uint32_t IXE_over2_read_num = 0;
					//MAPQ
					uint32_t total_MAPQ = 0;
					uint32_t MAPQ_less_than30_read_num = 0;
					uint32_t MAPQ_less_than10_read_num = 0;
					uint32_t MAPQ_is_0_read_num = 0;
					// analysis the reads
					while (bam_next(&c_b)) {
						bam1_t *b = &(c_b._brec);
						if(bam_is_secondary(b))
							continue;
						if(bam_is_supplementary(b))
							continue;
						if(bam_is_supplementary(b))
							continue;
						//basic
						int read_len = b->core.l_qseq;
						get_bam_seq_bin(0, b->core.l_qseq, read_str_bin, b);
						{
							readNum++;
							total_readLen += read_len;
						}
						//clip analysis:
						{
							int soft_clip_len_left; int soft_clip_len_right;
							bam_has_SH_cigar(b, &soft_clip_len_left, &soft_clip_len_right);
							if(soft_clip_len_left != 0 || soft_clip_len_right != 0){
								if(b->core.qual < 10){
									soft_clip_len_left = 0; soft_clip_len_right= 0;
								}
								//assume the quality of CLIP signals
								if (!pass_compact_filter(read_str_bin, read_len)){
									soft_clip_len_left = 0; soft_clip_len_right= 0;
								}
								else if (!clip_AAA_TailFilter(read_str_bin, read_len, soft_clip_len_left, soft_clip_len_right)){
									soft_clip_len_left = 0; soft_clip_len_right= 0;
								}
								else if (!clip_low_quality_Filter(b, read_len, soft_clip_len_left, soft_clip_len_right)){
									soft_clip_len_left = 0; soft_clip_len_right= 0;
								}
							}

							int read_total_clip_len = soft_clip_len_left + soft_clip_len_right;
							total_clip_len += read_total_clip_len;
							if(read_total_clip_len != 0)
								readNum_clip++;
							if(read_total_clip_len == 0){/*Do nothing*/}
							else if(read_total_clip_len < 11)
								readNum_clip_0_10++;
							else if(read_total_clip_len < 51)
								readNum_clip_11_50++;
							else
								readNum_clip_over_50++;
						}

						//DRP signals
						if(b->core.qual > 10){
							//read pair orientation
							bool read_direction = bam_is_fwd_strand(b);
							bool mate_direction = bam_is_mate_fwd_strand(b);
							if(b->core.pos > b->core.mpos) std::swap(read_direction, mate_direction);
							//mate read:
							bool mate_mapped = !bam_is_mate_unmapped(b);
							if(!mate_mapped) 	mate_unmapped_read_num ++;

							if(mate_mapped){
								if		(read_direction == FORWARD && mate_direction == FORWARD) insertSizeFF_read_num++;
								else if (read_direction == REVERSE && mate_direction == FORWARD) insertSizeRF_read_num++;
								else if (read_direction == REVERSE && mate_direction == REVERSE) insertSizeRR_read_num++;
							}
							//Trans
							bool mateInSame_chromsome = b->core.tid == b->core.mtid;
							if(!mateInSame_chromsome)	mate_different_chromsome_read_num ++;
							//insert size
							int32_t insertSizeOri = b->core.isize;
							uint32_t insertSize = ABS(insertSizeOri);

							if(insertSize > 1000) insertSizeOver1K_read_num++;
							if(insertSize > 100000) insertSizeOver100K_read_num++;

							//normal reads
							if(mate_mapped && read_direction == FORWARD && mate_direction == REVERSE && mateInSame_chromsome && insertSize <= 1000)
								{ /*DO NOTHING*/}
							else
								DRP_read_num++;
						}
						//SNP and INDEL
						{
							int INDEL_len = 0;
							int NM_len = 0;
							bam_get_INDEL_NM(b,&INDEL_len, &NM_len);
							int independence_event_number = 0;

							if(NM_len != 0 || INDEL_len != 0){
								char * tseq = ref_char + b->core.pos;
								uint8_t * qseq = read_str_bin;
								uint8_t * qqual = bam_get_qual(b);
								mismatch_position.clear();

								int seq_i = 0;
								int ref_i = 0;

								uint32_t* bam_cigar = bam_get_cigar(b);
								uint32_t n_cigar = b->core.n_cigar;
								for (uint i = 0; i < n_cigar; ++i)
								{
									int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
									int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
									switch (c_type){
									case CIGAR_MATCH: case CIGAR_SEQ_MATCH:
										for(int i = 0; i < c_size; i++, seq_i++, ref_i++)
											if(("ACGTNNNNNNNNN"[qseq[seq_i]] != tseq[ref_i]) && (qqual[seq_i] >= 20))
												mismatch_position.emplace_back(seq_i);
										break;
									case CIGAR_INSERT:
										for(int i = 0; i < c_size; i++, seq_i++)
											mismatch_position.emplace_back(seq_i);
										break; //do nothing
									case CIGAR_DELETE:
										for(int i = 0; i < c_size; i++, ref_i++)
											mismatch_position.emplace_back(seq_i);
										break;
									case CIGAR_SOFT_CLIP:
									case CIGAR_HARD_CLIP:
										break;
									default:	break;
									}
								}
								NM_len = mismatch_position.size();
								if(!mismatch_position.empty()){
									//check independence event number
									independence_event_number = 1;
									int mismatch_position_size = mismatch_position.size();
									for(int i = 0; i < mismatch_position_size - 1 ; i++)
										if(mismatch_position[i] + 1 != mismatch_position[i + 1] && mismatch_position[i] != mismatch_position[i + 1])
											independence_event_number ++;
								}
							}

							total_INDEL_len += INDEL_len;
							total_NM_len += NM_len;

							if(NM_len > 2)
								NM_over2_read_num++;
							if(NM_len > 10)
								NM_over10_read_num++;
							//~~~~~~~~~
							if(independence_event_number > 2)
								IXE_over2_read_num ++;
						}
						//MAPQ
						{
							uint8_t mapq = b->core.qual;
							total_MAPQ += mapq;
							if(mapq < 30)	MAPQ_less_than30_read_num++;
							if(mapq < 10)	MAPQ_less_than10_read_num++;
							if(mapq == 0)	MAPQ_is_0_read_num++;
						}
					}

					// out put information for that region
					//SV INFO
					bool blank_site =(sr.SVType.empty());
					if(blank_site)
						sr.SVType.append("BLANK");
					if(blank_site)
						printf("!0\tBLANK\t");
					else
						printf("!0\tSV_IN\t");
					//2
					printf("!1\t%d\t%d\t%s\t%d\t", u.chrID, sr.pos_begin, sr.SVType.c_str(), region_ID);
					//REGION:7
					printf("!2\t%d\t%d\t", analysis_region.st_pos, analysis_region.ed_pos);
					//10
					printf("!3\t%d\t%d\t", sr.maxVNTR_value,sr.isHIGH);
					//BAISC:14
					float depth = total_readLen/( analysis_region.ed_pos - analysis_region.st_pos + 148);
					printf("!4\t%d\t%ld\t%f\t", readNum, total_readLen, depth);
					if(readNum == 0) readNum = 1;
					//clip analysis:17
					printf("!5\t%d\t%d\t%d\t%d\t%d\t",total_clip_len ,readNum_clip ,readNum_clip_0_10, readNum_clip_11_50, readNum_clip_over_50);
					//DRP signals:23
					printf("!6\t%d\t%d\t%d\t%d\t", DRP_read_num, insertSizeFF_read_num, insertSizeRR_read_num, insertSizeRF_read_num);
					//30
					printf("!7\t%d\t%d\t%d\t%d\t", insertSizeOver1K_read_num, insertSizeOver100K_read_num,
							mate_unmapped_read_num, mate_different_chromsome_read_num);
					//SNP and INDEL:33
					printf("!8\t%d\t%d\t%d\t%d\t%d\t", total_INDEL_len, total_NM_len, NM_over2_read_num, NM_over10_read_num, IXE_over2_read_num);
					//MAPQ:39
					printf("!9\t%d\t%d\t%d\t%d\t", total_MAPQ, MAPQ_less_than30_read_num, MAPQ_less_than10_read_num, MAPQ_is_0_read_num);
					//46
					printf("\n");
				}
			}
		}
	}

	//close files
	bam_file_close(&c_b);
	return 0;
}
/*/*
int genomic_contest_Analysis(int argc, char *argv[]){
	//parameters
	char * ref_fn = argv[1];
	char * vcf_fn_in = argv[2];//separate by ','
	char * input_bam_fn = argv[3];
	char * ucsc_mappibility_binary_data_dump_path = argv[4];
	int ucsc_mappibility_key_value = atoi(argv[5]);
	char * input_centromeres_bed_fn = argv[6];
	char * input_VNTR_region_file = argv[7];
	char * HIGH_CONFIDENCE_region_file = argv[8];
	int used_chrID = atoi(argv[9]);
	bool used_BLANK = (atoi(argv[10]) == 1);

	//open vcf files
	BCF_FILE vcf_read;//vcf for read
	VCF_open_read(&vcf_read, vcf_fn_in);//open for read
	//open bam file
	Bam_file c_b;
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
	bam_hdr_t* bam_header = c_b._hdr;
	//VCF info buffs
	bcf_hdr_t * vcf_header = vcf_read.header;
	//bam info buffs

	std::vector<REF_BASIC_DATA> UCSC_mappibility_ref_l;
	faidx_t *c_ref_idx = reference_index_load(ref_fn);
	{
		std::vector<std::string> centromeres_bed_l;
		std::vector<std::string> item_value;
		//input_centromeres_bed_fn
		load_string_list_from_file(input_centromeres_bed_fn, centromeres_bed_l);
		//load reference index
		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int i = 0; i < N_seq; i++){
			char fn[1024];
			sprintf(fn, "%d_dump_file.bin", i);
			//mappibility file
			UCSC_mappibility_ref_l.emplace_back();
			UCSC_mappibility_ref_l.back().chrID = i;
			UCSC_mappibility_ref_l.back().ref_length =
					vector_load_bin(ucsc_mappibility_binary_data_dump_path, fn, (void **)&(UCSC_mappibility_ref_l.back().mappbility_value));
			//centromeres file
			split_string(item_value, centromeres_bed_l[i].c_str(), "\t");
			//BED begin and end
			UCSC_mappibility_ref_l.back().centromeres_region_bg = atoi(item_value[1].c_str());
			UCSC_mappibility_ref_l.back().centromeres_region_ed = atoi(item_value[2].c_str());
			//init region with SV
			UCSC_mappibility_ref_l.back().initRegionWithSV();
		}
	}
	//load all vcf data and store to the "RegionWithSV" list
	//
	{
		char *c_sv_type = (char *)malloc(1000);
		int SV_CHR_ID;
		int SV_POS;
		int SV_END;
		int SV_length;
		do
		{
			bcf1_t *c_r = &( vcf_read.r);
			//unpack the vcf data to get the Filter
			bcf_unpack(c_r, BCF_UN_INFO);
			//vcf filters
			if(*c_r->d.flt != 0)//filter: PASS
				continue;
			if(c_r->rid > 23)
				continue;
			//vcf get data


			SV_CHR_ID = c_r->rid;
			SV_POS = c_r->pos;
			vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
			vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
			vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

			if(strcmp(c_sv_type, "BND") == 0)
				continue;

			std::vector<int> SV_BP_pos_l;

			if(ABS_U(SV_POS, SV_END) < 10){
				int SV_POS_middle = (SV_POS + SV_END) / 2;
				SV_BP_pos_l.emplace_back(SV_POS_middle);
			}else{
				SV_BP_pos_l.emplace_back(SV_POS);
				SV_BP_pos_l.emplace_back(SV_END);
			}
			//
			for(int SV_BP_pos:SV_BP_pos_l){
				//get SV reion for this BP
				int store_block_id = SV_BP_pos/ANALYSIS_BLOCK_SIZE;
				//store data:
				UCSC_mappibility_ref_l[SV_CHR_ID].region_withSV[store_block_id].storeData(c_sv_type, SV_length);
			}
		}while(VCF_next(&vcf_read));//read one
	}

	int VNTR_region_length;
	int VNTR_region_step_size;
	//load VNTR repeat level
	{
		std::vector<std::string> VNTR_bed_l;
		load_string_list_from_file(input_VNTR_region_file, VNTR_bed_l);
		std::vector<std::string> item_value;

		bool isHeader = true;
		for(std::string &s :VNTR_bed_l){
			split_string(item_value, s.c_str(), "\t");
			if(isHeader){
				VNTR_region_length = atoi(item_value[0].c_str());
				VNTR_region_step_size = atoi(item_value[1].c_str());
				isHeader = false;
			}
			else{
				int chrID = atoi(item_value[0].c_str());
				int pos = atoi(item_value[1].c_str());
				int maxVNTR_value = atoi(item_value[2].c_str());
				if(chrID <= 23){
					for (int i = pos; i < pos + VNTR_region_length; i+=100)
					{
						UCSC_mappibility_ref_l[chrID].region_withSV[i/ANALYSIS_BLOCK_SIZE].storemaxVNTR_value(maxVNTR_value);
					}
				}
			}
		}
	}

	//store the high confidence region data
	{

		std::vector<std::string> HIGH_CON_bed_l;
		load_string_list_from_file(HIGH_CONFIDENCE_region_file, HIGH_CON_bed_l);
		std::vector<std::string> item_value;
		for(std::string &s :HIGH_CON_bed_l){
			split_string(item_value, s.c_str(), "\t");
			const char *chrNAME = item_value[0].c_str();
			int pos_bg = atoi(item_value[1].c_str());
			int pos_ed = atoi(item_value[2].c_str());
			int chrID = 0;
			//chr1
			if(chrNAME[3] == 'X'){
				chrID = 22;
			}else if(chrNAME[3] == 'Y'){
				chrID = 23;
			}else{
				chrID = atoi(chrNAME + 3) - 1;
			}
			if(chrID <= 23){
				for(int pos = pos_bg; pos < pos_ed; pos+=100)
					UCSC_mappibility_ref_l[chrID].region_withSV[pos/ANALYSIS_BLOCK_SIZE].isHIGH = true;
			}
		}
	}

	uint8_t read_str_bin[512];
  	std::vector<int> mismatch_position;
	//store data:
	int region_ID = -1;
	{
		REF_BASIC_DATA & u = UCSC_mappibility_ref_l[used_chrID];
		//load reference
		int ref_len = 0;
		char * ref_char = fai_fetch(c_ref_idx, faidx_iseq(c_ref_idx, u.chrID), &ref_len);
		for(SV_basic_infomation & sr: u.region_withSV){
			region_ID++;
			//select regions:
			if(!sr.SVType.empty() || (used_BLANK)){//SV or using BLANK
				R_region analysis_region;
				analysis_region.chr_ID = u.chrID;

				analysis_region.st_pos = sr.pos_begin - 50;
				analysis_region.ed_pos = sr.pos_begin + ANALYSIS_BLOCK_SIZE + 50;
				//skip N regions
				if(ref_char[analysis_region.st_pos] == 'N')
					continue;
				resetRegion_ID(&c_b, &analysis_region);
				uint8_t * mv = u.mappbility_value;
				float total_mappibility = 0;
				for(int i = analysis_region.st_pos; i < analysis_region.ed_pos; i++){
					total_mappibility += mv[i];
				}
				float average_mappibility = (float)total_mappibility/ucsc_mappibility_key_value/(analysis_region.ed_pos - analysis_region.st_pos);

				bool withIncentromeres_region = false;
				if(u.centromeres_region_bg < sr.pos_begin &&
						u.centromeres_region_ed > sr.pos_begin){
					withIncentromeres_region = true;
				}

				//reset calculators
				uint32_t readNum = 0;
				uint64_t total_readLen = 0;
				//clip signals
				uint32_t total_clip_len = 0;
				uint32_t readNum_clip = 0;
				uint32_t readNum_clip_0_10 = 0;
				uint32_t readNum_clip_11_50 = 0;
				uint32_t readNum_clip_over_50 = 0;
				//DRP signals
				uint32_t DRP_read_num = 0;//both forward:

				uint32_t insertSizeFF_read_num = 0;//both forward:
				uint32_t insertSizeRR_read_num = 0;//both reverse:
				uint32_t insertSizeRF_read_num = 0;//the read of smaller position is reverse, the other is forward

				uint32_t insertSizeOver1K_read_num = 0;//the read of insert size > 1023
				uint32_t insertSizeOver100K_read_num = 0;//the read of insert size > 10000
				uint32_t mate_unmapped_read_num = 0;
				uint32_t mate_different_chromsome_read_num = 0;
				//NM signals
				//SNP and INDEL
				uint32_t total_INDEL_len = 0;
				uint32_t total_NM_len = 0;
				uint32_t NM_over2_read_num = 0;
				uint32_t NM_over10_read_num = 0;

				uint32_t IXE_over2_read_num = 0;

				//MAPQ
				uint32_t total_MAPQ = 0;
				uint32_t MAPQ_less_than30_read_num = 0;
				uint32_t MAPQ_less_than10_read_num = 0;
				uint32_t MAPQ_is_0_read_num = 0;
				// analysis the reads

				while (bam_next(&c_b)) {
					bam1_t *b = &(c_b._brec);
					if(bam_is_secondary(b))
						continue;
					if(bam_is_supplementary(b))
						continue;
					if(bam_is_supplementary(b))
						continue;
					//basic
					int read_len = b->core.l_qseq;
					get_bam_seq_bin(0, b->core.l_qseq, read_str_bin, b);
					{
						readNum++;
						total_readLen += read_len;
					}
					//clip analysis:
					{
						int soft_clip_len_left; int soft_clip_len_right;
						bam_has_SH_cigar(b, &soft_clip_len_left, &soft_clip_len_right);
						if(soft_clip_len_left != 0 || soft_clip_len_right != 0){
							if(b->core.qual < 10){
								soft_clip_len_left = 0; soft_clip_len_right= 0;
							}
							//assume the quality of CLIP signals
							if (!pass_compact_filter(read_str_bin, read_len)){
								soft_clip_len_left = 0; soft_clip_len_right= 0;
							}
							else if (!clip_AAA_TailFilter(read_str_bin, read_len, soft_clip_len_left, soft_clip_len_right)){
								soft_clip_len_left = 0; soft_clip_len_right= 0;
							}
							else if (!clip_low_quality_Filter(b, read_len, soft_clip_len_left, soft_clip_len_right)){
								soft_clip_len_left = 0; soft_clip_len_right= 0;
							}
						}

						int read_total_clip_len = soft_clip_len_left + soft_clip_len_right;
						total_clip_len += read_total_clip_len;
						if(read_total_clip_len != 0)
							readNum_clip++;
						if(read_total_clip_len == 0){}//Do nothing
						else if(read_total_clip_len < 11)
							readNum_clip_0_10++;
						else if(read_total_clip_len < 51)
							readNum_clip_11_50++;
						else
							readNum_clip_over_50++;
					}

					//DRP signals
					if(b->core.qual > 10){
						//read pair orientation
						bool read_direction = bam_is_fwd_strand(b);
						bool mate_direction = bam_is_mate_fwd_strand(b);
						if(b->core.pos > b->core.mpos) std::swap(read_direction, mate_direction);
						//mate read:
						bool mate_mapped = !bam_is_mate_unmapped(b);
						if(!mate_mapped) 	mate_unmapped_read_num ++;

						if(mate_mapped){
							if		(read_direction == FORWARD && mate_direction == FORWARD) insertSizeFF_read_num++;
							else if (read_direction == REVERSE && mate_direction == FORWARD) insertSizeRF_read_num++;
							else if (read_direction == REVERSE && mate_direction == REVERSE) insertSizeRR_read_num++;
						}
						//Trans
						bool mateInSame_chromsome = b->core.tid == b->core.mtid;
						if(!mateInSame_chromsome)	mate_different_chromsome_read_num ++;
						//insert size
						int32_t insertSizeOri = b->core.isize;
						uint32_t insertSize = ABS(insertSizeOri);

						if(insertSize > 1000) insertSizeOver1K_read_num++;
						if(insertSize > 100000) insertSizeOver100K_read_num++;

						//normal reads
						if(mate_mapped && read_direction == FORWARD && mate_direction == REVERSE && mateInSame_chromsome && insertSize <= 1000)
							{ }
						else
							DRP_read_num++;
					}

					//SNP and INDEL
					{
						int INDEL_len = 0;
						int NM_len = 0;
						bam_get_INDEL_NM(b,&INDEL_len, &NM_len);
						int independence_event_number = 0;

						if(NM_len != 0 || INDEL_len != 0){
							char * tseq = ref_char + b->core.pos;
							uint8_t * qseq = read_str_bin;
							uint8_t * qqual = bam_get_qual(b);
							mismatch_position.clear();

							int seq_i = 0;
							int ref_i = 0;

							uint32_t* bam_cigar = bam_get_cigar(b);
							uint32_t n_cigar = b->core.n_cigar;
							for (uint i = 0; i < n_cigar; ++i)
							{
								int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
								int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
								switch (c_type){
								case CIGAR_MATCH: case CIGAR_SEQ_MATCH:
									for(int i = 0; i < c_size; i++, seq_i++, ref_i++)
										if(("ACGTNNNNNNNNN"[qseq[seq_i]] != tseq[ref_i]) && (qqual[seq_i] >= 20))
											mismatch_position.emplace_back(seq_i);
									break;
								case CIGAR_INSERT:
									for(int i = 0; i < c_size; i++, seq_i++)
										mismatch_position.emplace_back(seq_i);
									break; //do nothing
								case CIGAR_DELETE:
									for(int i = 0; i < c_size; i++, ref_i++)
										mismatch_position.emplace_back(seq_i);
									break;
								case CIGAR_SOFT_CLIP:
								case CIGAR_HARD_CLIP:
									break;
								default:	break;
								}
							}
							NM_len = mismatch_position.size();
							if(!mismatch_position.empty()){
								//check independence event number
								independence_event_number = 1;
								int mismatch_position_size = mismatch_position.size();
								for(int i = 0; i < mismatch_position_size - 1 ; i++)
									if(mismatch_position[i] + 1 != mismatch_position[i + 1] && mismatch_position[i] != mismatch_position[i + 1])
										independence_event_number ++;
							}
						}

						total_INDEL_len += INDEL_len;
						total_NM_len += NM_len;

						if(NM_len > 2)
							NM_over2_read_num++;
						if(NM_len > 10)
							NM_over10_read_num++;
						//~~~~~~~~~
						if(independence_event_number > 2)
							IXE_over2_read_num ++;
					}

					//MAPQ
					{
						uint8_t mapq = b->core.qual;
						total_MAPQ += mapq;
						if(mapq < 30)	MAPQ_less_than30_read_num++;
						if(mapq < 10)	MAPQ_less_than10_read_num++;
						if(mapq == 0)	MAPQ_is_0_read_num++;
					}
				}

				// out put information for that region
				//title:
					bool isLOW_SIG = (((float)(readNum_clip_11_50 +  readNum_clip_over_50) / readNum) < 0.04);
					isLOW_SIG &= (((float)(DRP_read_num) / readNum) < 0.04);
					isLOW_SIG &= (sr.maxVNTR_value < 20);
					isLOW_SIG &= (!withIncentromeres_region);
					isLOW_SIG &= (sr.isHIGH);
					//if(!isLOW_SIG)	continue;
				//SV INFO
				bool blank_site =(sr.SVType.empty());
				if(blank_site)
					sr.SVType.append("BLANK");
				if(blank_site)
					printf("!0\tBLANK\t");
				else
					printf("!0\tSV_IN\t");
				//2
				printf("!1\t%d\t%d\t%s\t%d\t", u.chrID, sr.pos_begin, sr.SVType.c_str(), region_ID);
				//REGION:7
				printf("!2\t%d\t%d\t", analysis_region.st_pos, analysis_region.ed_pos);
				//10
				printf("!3\t%f\t%d\t%d\t%d\t", average_mappibility, withIncentromeres_region,sr.maxVNTR_value,sr.isHIGH);
				//BAISC:14
				float depth = total_readLen/( analysis_region.ed_pos - analysis_region.st_pos + 148);
				printf("!4\t%d\t%ld\t%f\t", readNum, total_readLen, depth);
				if(readNum == 0) readNum = 1;
				//clip analysis:19
				printf("!5\t%d\t%d\t%d\t%d\t%d\t",total_clip_len ,readNum_clip ,readNum_clip_0_10, readNum_clip_11_50, readNum_clip_over_50);
				//DRP signals:25
				printf("!6\t%d\t%d\t%d\t%d\t", DRP_read_num, insertSizeFF_read_num, insertSizeRR_read_num, insertSizeRF_read_num);
				//30
				printf("!7\t%d\t%d\t%d\t%d\t", insertSizeOver1K_read_num, insertSizeOver100K_read_num,
						mate_unmapped_read_num, mate_different_chromsome_read_num);
				//SNP and INDEL:35
				printf("!8\t%d\t%d\t%d\t%d\t%d\t", total_INDEL_len, total_NM_len, NM_over2_read_num, NM_over10_read_num, IXE_over2_read_num);
				//MAPQ:41
				printf("!9\t%d\t%d\t%d\t%d\t", total_MAPQ, MAPQ_less_than30_read_num, MAPQ_less_than10_read_num, MAPQ_is_0_read_num);
				//46
				printf("\n");
			}
		}
	}

	//close files
	bam_file_close(&c_b);
	bcf_close(vcf_read.file);
	return 0;
}
*/

#define HEAT_ROW_NUM 100
#define HEAT_COL_NUM 300
int heatMapData(int argc, char *argv[]){

	char * GCA_fn = argv[1];
	int YValue_method = atoi(argv[2]);
	int used_row_Num = atoi(argv[3]);
	bool using_only_high_region = (atoi(argv[4]) == 1);

	//load all data
	std::vector<std::string> data;
	std::vector<std::string> item_value;
	load_string_list_from_file(GCA_fn, data);
	//init
	struct HAET_MAP_RST{
		int region_SV = 0;
		int region_BLANK = 0;
		float mappibility_bg = 0;
		float mappibility_ed = 0;

		float YValue_bg = 0;
		float YValue_ed = 0;
	};
	std::vector<std::vector<HAET_MAP_RST>> heat_map;
	heat_map.resize(HEAT_ROW_NUM);
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		hr.resize(HEAT_COL_NUM);
	}

	int TB = 0;
	int TS = 0;
	//
	long int data_index = 0;
	for(std::string & l:data){
		data_index++;
		//if(data_index % 100000 == 0){
		//	fprintf(stderr, "%ld\n", data_index);
		//}
		split_string(item_value, l.c_str(), "\t");
		bool withIncentromeres_region = (item_value[12].compare("1") == 0);

		//todo::
		{
			int chrID = atoi(item_value[5].c_str());
			int POS = atoi(item_value[6].c_str());
			if(chrID ==0 && POS > 121751000 && POS < 125184587){ withIncentromeres_region = true;}
			if(chrID ==1 && POS > 89816500 && POS < 94090557){ withIncentromeres_region = true;}
			if(chrID ==2 && POS > 90533500 && POS < 93655574){ withIncentromeres_region = true;}
			if(chrID ==3 && POS > 49091500 && POS < 51743951){ withIncentromeres_region = true;}
			if(chrID ==4 && POS > 46485901 && POS < 50059807){ withIncentromeres_region = true;}
			if(chrID ==6 && POS > 151855000 && POS < 158595000){ withIncentromeres_region = true;}
			if(chrID ==8 && POS > 43236168 && POS < 45518558){ withIncentromeres_region = true;}
			if(chrID ==9 && POS > 38485000 && POS < 41914000){ withIncentromeres_region = true;}
			if(chrID ==16 && POS > 21894500 && POS < 26885980){ withIncentromeres_region = true;}
			if(chrID ==19 && POS > 26436233 && POS < 31470500){ withIncentromeres_region = true;}
			if(chrID ==21 && POS > 11212000 && POS < 16567000){ withIncentromeres_region = true;}
			//if(chrID ==22){ withIncentromeres_region = true;}
			if(chrID ==23){ withIncentromeres_region = true;}
		}
		if(withIncentromeres_region)
			continue;
		bool withInHIGH_region = (item_value[14].compare("1") == 0);
		if(!withInHIGH_region && using_only_high_region)
			continue;
		int readNum = atoi(item_value[16].c_str());
		if(readNum < 50) continue;
		//if(readNum < 50 || readNum > 800) continue; 
		//
		bool isSV = (item_value[1].compare("SV_IN") == 0);

		float mappibility = atof(item_value[11].c_str());
		int   VNTR_CPX = atoi(item_value[13].c_str());

		//select the Y value:
		int CLIP_ALL = atoi(item_value[22].c_str()) + atoi(item_value[23].c_str()) + atoi(item_value[24].c_str());
		int CLIP_over10 = atoi(item_value[23].c_str()) + atoi(item_value[24].c_str());
		int CLIP_over50 = atoi(item_value[24].c_str());
		int DRP_ALL = atoi(item_value[26].c_str());
		int CIGAR_ALL = atoi(item_value[40].c_str());

		int YValue = 0;
		switch(YValue_method){
		case 0:	YValue = atoi(item_value[used_row_Num].c_str()); break;
		case 1:	YValue = (CLIP_over10 + DRP_ALL)/2; break;
		case 2:	YValue = MAX(CLIP_over10, DRP_ALL); break;
		case 3:	YValue = MAX(CLIP_over10, DRP_ALL); YValue = MAX(YValue, CIGAR_ALL); break;
		case 4:	YValue = CLIP_over10; break;
		case 5:	YValue = CLIP_over50; break;
		}

		//fprintf(stderr, "[%d\t%f\t%d\t]", readNum,mappibility,YValue);
		//YValue_idx: M1
		int YValue_idx = 0;
		if(false){
			float YValue_float = (float)YValue/readNum;
			YValue_idx = YValue_float*HEAT_ROW_NUM;
			if(YValue_float == 1){
				YValue_idx = HEAT_ROW_NUM - 1;
			}
		}else{
			YValue_idx = YValue;
			if(YValue_idx >= HEAT_ROW_NUM){
				YValue_idx = HEAT_ROW_NUM - 1;
			}
		}

		if(YValue_idx > HEAT_ROW_NUM){
			fprintf(stderr, "ERROR: Error data in line %s\n ", l.c_str());
			continue;
		}

		int VNTR_CPX_idx = VNTR_CPX;
		if(VNTR_CPX_idx >= HEAT_COL_NUM){
			VNTR_CPX_idx = HEAT_COL_NUM - 1;
		}
		int mappibility_idx = mappibility*HEAT_COL_NUM;
		if(mappibility == 1){
			mappibility_idx = HEAT_COL_NUM - 1;
		}
//		if(isSV){
//			heat_map[YValue_idx][mappibility_idx].region_SV++;
//		}else
//			heat_map[YValue_idx][mappibility_idx].region_BLANK++;

		if(isSV){
			heat_map[YValue_idx][VNTR_CPX_idx].region_SV++;
		}else{
			heat_map[YValue_idx][VNTR_CPX_idx].region_BLANK++;
			if(YValue_idx > 60 ) 
				fprintf(stderr, "%s\n", l.c_str());	
		}

		if(isSV)
			TS++;
		else
			TB++;
		if(false && YValue_idx ==0 && mappibility_idx == 0){
			for(std::string & i:item_value){
				fprintf(stderr, "%s\t", i.c_str());
			}
			fprintf(stderr, "\n");
		}
	}

	fprintf(stderr, "TS %d, TB %d\n\n", TS, TB);
	fprintf(stdout, "heatMapData\t");
	for(uint i = 0; i < HEAT_COL_NUM; i++){
		fprintf(stdout, "M_%d\t",i);
	}
	fprintf(stdout, "\n");
	uint i = 0;
	//P_log
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		fprintf(stdout, "Y_%d\t", i++);
		for(HAET_MAP_RST & h: hr){
			float total_region_num = (h.region_SV+h.region_BLANK);
			float P = (float)h.region_SV/(total_region_num);
			if((h.region_BLANK == 0 && h.region_SV == 0)){
				fprintf(stdout, "NA\t");	
			}else{
				float P1 = -10*log10((double)P);
				fprintf(stdout, "%.4f\t", P1);	
			}
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n\n\n");
	//P_ori
	i = 0;
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		fprintf(stdout, "Y_%d\t", i++);
		for(HAET_MAP_RST & h: hr){
			float total_region_num = (h.region_SV+h.region_BLANK);
			float P = (float)h.region_SV/(total_region_num);
			if((h.region_BLANK == 0 && h.region_SV == 0)){
				fprintf(stdout, "NA\t");
			}else{
				fprintf(stdout, "%.8f\t", P);
			}
			
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n\n\n");
	//total_region_num only
	i = 0;
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		fprintf(stdout, "Y_%d\t", i++);
		for(HAET_MAP_RST & h: hr){
			int total_region_num = (h.region_SV+h.region_BLANK);
			if((h.region_BLANK == 0 && h.region_SV == 0)){
				fprintf(stdout, "NA\t");
			}else
			fprintf(stdout, "%d\t", total_region_num);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n\n\n");

	//region_SV only
	i = 0;
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		fprintf(stdout, "Y_%d\t", i++);
		for(HAET_MAP_RST & h: hr){
			fprintf(stdout, "%d\t", h.region_SV);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n\n\n");

	//region_BLANK only
	i = 0;
	for(std::vector<HAET_MAP_RST> & hr: heat_map){
		fprintf(stdout, "Y_%d\t", i++);
		for(HAET_MAP_RST & h: hr){
			fprintf(stdout, "%d\t", h.region_BLANK);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n\n\n");

	return 0;
}

int heatMapDataSum(int argc, char *argv[]){

	char * dummary_fn_l_f = argv[1];
	uint64_t result_1[HEAT_ROW_NUM][HEAT_COL_NUM] = {0};
	uint64_t result_2[HEAT_ROW_NUM][HEAT_COL_NUM] = {0};

	std::vector<std::string> fn_l;
	std::vector<std::string> item_value;
	load_string_list_from_file(dummary_fn_l_f, fn_l);
	for(std::string & fn:fn_l){
		std::vector<std::string> data_f;
		
		load_string_list_from_file(fn.c_str(), data_f);
		fprintf(stderr, "Load data %s %ld\n", fn.c_str(), data_f.size());
		//result
		for(int i = 310; i < 410; i++){
			std::string & data_line = data_f[i];
			split_string(item_value, data_line.c_str(), "\t");
			for(int j = 1; j < HEAT_COL_NUM + 1; j++){
				int d = atoi(item_value[j].c_str());
				result_1[i-310][j-1] += d;
			}
		}
		for(int i = 413; i < 513; i++){
			std::string & data_line = data_f[i];
			split_string(item_value, data_line.c_str(), "\t");
			for(int j = 1; j < HEAT_COL_NUM + 1; j++){
				int d = atoi(item_value[j].c_str());
				result_2[i-413][j-1] += d;
			}
		}
	}

	for(int i = 0; i < 100; i++){
		fprintf(stdout, "Y_%d\t", i);
		for(int j = 0; j < 200; j++){
			//if(result_2[i][j] + result_1[i][j] < 20)
			//	fprintf(stdout, "0\t");
			//else
				fprintf(stdout, "%ld\t", result_1[i][j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "\n");
	for(int i = 0; i < 100; i++){
		fprintf(stdout, "Y_%d\t", i);
		for(int j = 0; j < 200; j++){
			//if(result_2[i][j] + result_1[i][j] < 20)
			//	fprintf(stdout, "0\t");
			//else
				fprintf(stdout, "%ld\t", result_2[i][j]);
		}
		fprintf(stdout, "\n");
	}
	return 0;
}

int signalPattenAnalysis(int argc, char *argv[]){
	//parameters
	char * ref_fn = argv[1];
	char * vcf_fn_in = argv[2];//separate by ','
	char * input_bam_fn = argv[3];
	char * ucsc_mappibility_binary_data_dump_path = argv[4];
	int ucsc_mappibility_key_value = atoi(argv[5]);
	int edge_length = atoi(argv[6]);

	char * input_centromeres_bed_fn = argv[7];

	//open vcf files
	BCF_FILE vcf_read;//vcf for read
	VCF_open_read(&vcf_read, vcf_fn_in);//open for read
	//open bam file
	Bam_file c_b;
	bam_file_open(input_bam_fn, ref_fn, NULL, &c_b);
	bam_hdr_t* bam_header = c_b._hdr;
	//VCF info buffs
	int SV_CHR_ID = 0;
	int SV_POS = 0;
	int SV_END = 0;
	int SV_length = 0;
	char *c_sv_type = (char *)malloc(1000);
	bcf_hdr_t * vcf_header = vcf_read.header;
	faidx_t *c_ref_idx = reference_index_load(ref_fn);

	//bam info buffs

	std::vector<REF_BASIC_DATA> UCSC_mappibility_ref_l;
	{
		std::vector<std::string> centromeres_bed_l;
		std::vector<std::string> item_value;
		//input_centromeres_bed_fn
		load_string_list_from_file(input_centromeres_bed_fn, centromeres_bed_l);
		//load reference index

		int N_seq = faidx_nseq(c_ref_idx);
		N_seq = MIN(24, N_seq);
		for(int i = 0; i < N_seq; i++){
			char fn[1024];
			sprintf(fn, "%d_dump_file.bin", i);
			//mappibility file
			UCSC_mappibility_ref_l.emplace_back();
			UCSC_mappibility_ref_l.back().ref_length =
					vector_load_bin(ucsc_mappibility_binary_data_dump_path, fn, (void **)&(UCSC_mappibility_ref_l.back().mappbility_value));
			//centromeres file
			split_string(item_value, centromeres_bed_l[i].c_str(), "\t");
			//BED begin and end
			UCSC_mappibility_ref_l.back().centromeres_region_bg = atoi(item_value[1].c_str());
			UCSC_mappibility_ref_l.back().centromeres_region_ed = atoi(item_value[2].c_str());
		}
	}

	//load all vcf data
	do
	{
		bcf1_t *c_r = &( vcf_read.r);
		//unpack the vcf data to get the Filter
		bcf_unpack(c_r, BCF_UN_INFO);
		//vcf filters
		if(*c_r->d.flt != 0)//filter: PASS
			continue;
		if(c_r->rid > 23)
			continue;
		//vcf get data
		SV_CHR_ID = c_r->rid;
		SV_POS = c_r->pos;
		vcf_get_sv_END(vcf_read.header, c_r, &SV_END);
		vcf_get_sv_LENGTH(vcf_read.header, c_r, &SV_length);
		vcf_get_sv_type(vcf_read.header, c_r, c_sv_type);

		std::vector<R_region> SV_region_l;
		if(ABS_U(SV_POS, SV_END) < 10){
			int SV_POS_middle = (SV_POS + SV_END) / 2;
			setR_Region(SV_region_l, SV_CHR_ID, SV_POS_middle - edge_length, SV_POS_middle + edge_length);
		}else{
			setR_Region(SV_region_l, SV_CHR_ID, SV_POS - edge_length, SV_POS + edge_length);
			setR_Region(SV_region_l, SV_CHR_ID, SV_END - edge_length, SV_END + edge_length);
		}
		int region_ID = 0;
		for(R_region & analysis_region:SV_region_l){
			region_ID++;
			//set analysis region
			resetRegion_ID(&c_b, &analysis_region);
			//mappability information

			uint8_t * mv = UCSC_mappibility_ref_l[analysis_region.chr_ID].mappbility_value;
			float total_mappibility = 0;
			for(int i = analysis_region.st_pos; i < analysis_region.ed_pos; i++){
				total_mappibility += mv[i];
			}
			float average_mappibility = (float)total_mappibility/ucsc_mappibility_key_value/(analysis_region.ed_pos - analysis_region.st_pos + 1);

			//reset calculators
			uint32_t readNum = 0;
			uint64_t total_readLen = 0;
			//clip signals
			uint32_t total_clip_len = 0;
			uint32_t readNum_clip = 0;
			uint32_t readNum_clip_0_10 = 0;
			uint32_t readNum_clip_11_50 = 0;
			uint32_t readNum_clip_over_50 = 0;
			//DRP signals
			uint32_t DRP_read_num = 0;//both forward:

			uint32_t insertSizeFF_read_num = 0;//both forward:
			uint32_t insertSizeRR_read_num = 0;//both reverse:
			uint32_t insertSizeRF_read_num = 0;//the read of smaller position is reverse, the other is forward

			uint32_t insertSizeOver1K_read_num = 0;//the read of insert size > 1023
			uint32_t insertSizeOver100K_read_num = 0;//the read of insert size > 10000
			uint32_t mate_unmapped_read_num = 0;
			uint32_t mate_different_chromsome_read_num = 0;
			//NM signals
			//SNP and INDEL
			uint32_t total_INDEL_len = 0;
			uint32_t total_NM_len = 0;
			uint32_t NM_over2_read_num = 0;
			uint32_t NM_over10_read_num = 0;

			//MAPQ
			uint32_t total_MAPQ = 0;
			uint32_t MAPQ_less_than30_read_num = 0;
			uint32_t MAPQ_less_than10_read_num = 0;
			uint32_t MAPQ_is_0_read_num = 0;

			// analysis the reads
			while (bam_next(&c_b)) {
				bam1_t *b = &(c_b._brec);
				//basic
				{
					readNum++;
					total_readLen += b->core.l_qseq;
				}
				//clip analysis:
				{
					int soft_clip_len_left; int soft_clip_len_right;
					bam_has_SH_cigar(b, &soft_clip_len_left, &soft_clip_len_right);
					int read_total_clip_len = soft_clip_len_left + soft_clip_len_right;
					total_clip_len += read_total_clip_len;
					if(read_total_clip_len != 0)
						readNum_clip++;
					if(read_total_clip_len == 0){/*Do nothing*/}
					else if(read_total_clip_len < 11)
						readNum_clip_0_10++;
					else if(read_total_clip_len < 51)
						readNum_clip_11_50++;
					else
						readNum_clip_over_50++;
				}

				//DRP signals
				{
					//read pair orientation
					bool read_direction = bam_is_fwd_strand(b);
					bool mate_direction = bam_is_mate_fwd_strand(b);
					if(b->core.pos > b->core.mpos) std::swap(read_direction, mate_direction);
					//mate read:
					bool mate_mapped = !bam_is_mate_unmapped(b);
					if(!mate_mapped) 	mate_unmapped_read_num ++;

					if(mate_mapped){
						if		(read_direction == FORWARD && mate_direction == FORWARD) insertSizeFF_read_num++;
						else if (read_direction == REVERSE && mate_direction == FORWARD) insertSizeRF_read_num++;
						else if (read_direction == REVERSE && mate_direction == REVERSE) insertSizeRR_read_num++;
					}
					//Trans
					bool mateInSame_chromsome = b->core.tid == b->core.mtid;
					if(!mateInSame_chromsome)	mate_different_chromsome_read_num ++;
					//insert size
					int32_t insertSizeOri = b->core.isize;
					uint32_t insertSize = ABS(insertSizeOri);

					if(insertSize > 1000) insertSizeOver1K_read_num++;
					if(insertSize > 100000) insertSizeOver100K_read_num++;

					//normal reads
					if(mate_mapped && read_direction == FORWARD && mate_direction == REVERSE && mateInSame_chromsome && insertSize <= 1000)
						{ /*DO NOTHING*/}
					else
						DRP_read_num++;
				}

				//SNP and INDEL
				{
					int INDEL_len = 0;
					int NM_len = 0;
					bam_get_INDEL_NM(b,&INDEL_len, &NM_len);

					total_INDEL_len += INDEL_len;
					total_NM_len += NM_len;

					if(NM_len > 2)
						NM_over2_read_num++;
					if(NM_len > 10)
						NM_over10_read_num++;
				}

				//MAPQ
				{
					uint8_t mapq = b->core.qual;
					total_MAPQ += mapq;
					if(mapq < 30)	MAPQ_less_than30_read_num++;
					if(mapq < 10)	MAPQ_less_than10_read_num++;
					if(mapq == 0)	MAPQ_is_0_read_num++;
				}
			}

			bool withIncentromeres_region = false;
			if(UCSC_mappibility_ref_l[SV_CHR_ID].centromeres_region_bg < SV_POS &&
					UCSC_mappibility_ref_l[SV_CHR_ID].centromeres_region_ed > SV_POS){
				withIncentromeres_region = true;
			}

			// out put information for that region
			//title:
			//SV INFO
			printf("!1\t%d_%d_%s_%d\t", SV_CHR_ID, SV_POS, c_sv_type, region_ID);
			printf("!2\t%d\t%d\t%d\t%d\t%s\t%d\t", SV_CHR_ID, SV_POS, SV_END, SV_length, c_sv_type, region_ID);
			//REGION
			printf("!3\t%d\t%d\t%f\t%d\t", analysis_region.st_pos, analysis_region.ed_pos, average_mappibility, withIncentromeres_region);
			//BAISC
			float depth = total_readLen/( analysis_region.ed_pos - analysis_region.st_pos + 148);
			printf("!4\t%d\t%ld\t%f\t", readNum, total_readLen, depth);
			if(readNum == 0) readNum = 1;
			//clip analysis
			printf("!5\t%d\t%d\t%d\t%d\t%d\t",total_clip_len ,readNum_clip ,readNum_clip_0_10, readNum_clip_11_50, readNum_clip_over_50);
			//DRP signals
			printf("!6\t%d\t%d\t%d\t%d\t", DRP_read_num, insertSizeFF_read_num, insertSizeRR_read_num, insertSizeRF_read_num);
			printf("!7\t%d\t%d\t%d\t%d\t", insertSizeOver1K_read_num, insertSizeOver100K_read_num,
					mate_unmapped_read_num, mate_different_chromsome_read_num);
			//SNP and INDEL
			printf("!8\t%d\t%d\t%d\t%d\t", total_INDEL_len, total_NM_len, NM_over2_read_num, NM_over10_read_num);
			//MAPQ
			printf("!9\t%d\t%d\t%d\t%d\t", total_MAPQ, MAPQ_less_than30_read_num, MAPQ_less_than10_read_num, MAPQ_is_0_read_num);
			printf("\n");
		}
	}while(VCF_next(&vcf_read));//read one

	//close files
	bam_file_close(&c_b);
	bcf_close(vcf_read.file);
	return 0;
}

#define MAX_LINE_LENGTH_1 100000
int UCSC_mappibility_binary_dump(int argc, char *argv[]){
	//parameters
	char * ref_fn = argv[1];
	char * ucsc_mappibility_fn_in = argv[2];
	int ucsc_mappibility_key_value = atoi(argv[3]);
	char * ucsc_mappibility_binary_data_dump_path = argv[4];

	//load reference index
	faidx_t *c_ref_idx = reference_index_load(ref_fn);
	int N_seq = faidx_nseq(c_ref_idx);
	N_seq = MIN(26, N_seq);
	//m-alloc for all data
	std::vector<REF_BASIC_DATA> UCSC_mappibility_ref_l;
	std::map<std::string , int > map_2_id;
	for(int i = 0; i < N_seq; i++){
		const char * chrName = faidx_iseq(c_ref_idx, i);
		int chrLength = faidx_seq_len(c_ref_idx, chrName);
		UCSC_mappibility_ref_l.emplace_back(i, chrLength, chrName);
		std::string chrName_str(chrName);
		map_2_id[chrName_str] = i;
	}

	//load data:

	{
		//get file names
		char *temp = new char[MAX_LINE_LENGTH_1];//100000
		std::ifstream ucsc_mappibility_file(ucsc_mappibility_fn_in);
		std::vector<std::string> item_value;
		std::vector<std::string> item_value_inner;
		int st_POS = 0;
		uint8_t * cur_mappbility_value = NULL;
		int store_index = 0;
		int block_count = 0;
		while(true){
			ucsc_mappibility_file.getline(temp, MAX_LINE_LENGTH_1);
			if(*temp == 0)	break;
			//analysis data:
			if(strlen(temp) > 10){
				if(block_count++ % 50 == 0)
					fprintf(stderr, "%s\n", temp);
				split_string(item_value, temp, " ");
				split_string(item_value_inner, item_value[1].c_str(), "=");
				std::string chrNAME = item_value_inner[1];

				std::map<std::string, int>::iterator it = map_2_id.find(chrNAME);
				cur_mappbility_value = UCSC_mappibility_ref_l[it->second].mappbility_value;
				split_string(item_value_inner, item_value[2].c_str(), "=");
				st_POS = atoi(item_value_inner[1].c_str());
				store_index = 0;
			}else{
				float map_value_float = atof(temp);
				int map_value_int = map_value_float*ucsc_mappibility_key_value;
				cur_mappbility_value[st_POS + store_index] = map_value_int;
				store_index ++;
			}
		}
		ucsc_mappibility_file.close();
	}

	//dump results
	{
		for(int i = 0; i < N_seq; i++){
			//file name
			char fn[1024];
			sprintf(fn, "%d_dump_file.bin", i);
			vector_dump_bin(ucsc_mappibility_binary_data_dump_path, fn, UCSC_mappibility_ref_l[i].mappbility_value, UCSC_mappibility_ref_l[i].ref_length);
		}
	}
	return 0;
}

int get_ref_region(int argc, char *argv[]){

  //load reference string;
  faidx_t *fai = fai_load(argv[optind]);//load reference
  int true_region_load_len = 0;
  char *ref_seq = fai_fetch(fai, argv[optind + 1], &true_region_load_len);
  fprintf(stdout, "%s\n", ref_seq);

  fai_destroy(fai);

  return 0;
}

int VENN_generate(int argc, char *argv[]){
	int number[31];
	int mode = atoi(argv[1]);
	int TYPE_N = 0;
	int TOOL_N = mode + 3;
	if(mode == 0)//3 tools
		TYPE_N = 7;
	if(mode == 1)//4 tools
		TYPE_N = 15;
	else if(mode == 2)
		TYPE_N = 31;
	
	for(int i= 0; i < TYPE_N; i++)
		number[i] = atoi(argv[i+2]);
	
	int data_idx = 0;
	std::vector<int> data[5];

	bool ADD[5] = {0};
	for (int J = 1; J <= TYPE_N; ++J)
	{	
		for (int K = 0; K < TOOL_N; ++K){
			if(((J >> K) & 0x1) == 1){
				for(int i = 0; i< number[J-1]; i++){
					data[TOOL_N-K-1].emplace_back(data_idx + i);
				}
			}
		}
		data_idx += number[J-1];
	}

	if(mode == 0){
		for(int i:data[0])	fprintf(stdout, "%d\tT2T\n", i);
		for(int i:data[1])	fprintf(stdout, "%d\tGCSV\n", i);
		for(int i:data[2])	fprintf(stdout, "%d\t1kgp_phase4\n", i);
	}
	else if(mode == 1){
		for(int i:data[0])	fprintf(stdout, "%d\tGCSV\n", i);
		for(int i:data[1])	fprintf(stdout, "%d\t1kgp_phase4\n", i);
		for(int i:data[2])	fprintf(stdout, "%d\tManta\n", i);
		for(int i:data[3])	fprintf(stdout, "%d\tINSurVeyor\n", i);
	}else if(mode == 2){
		for(int i:data[0])	fprintf(stdout, "%d\tT2T\n", i);
		for(int i:data[1])	fprintf(stdout, "%d\tGCSV\n", i);
		for(int i:data[2])	fprintf(stdout, "%d\t1kgp_phase4\n", i);
		for(int i:data[3])	fprintf(stdout, "%d\tManta\n", i);
		for(int i:data[4])	fprintf(stdout, "%d\tINSurVeyor\n", i);
	}

	return 0;
}


int N_count_ANA(const std::string &word){
	int N = 0;
	const char *s = word.c_str();
	int n = word.size();
	for(int i = 0; i < n; i++){
		if(s[i] == 'N')
			N++;
	}
	return N;
}

int fa_stat_full(int argc, char *argv[]){

	char * fa_fn_in = argv[1];//separate by ','
	int chrID = atoi(argv[2]);

	int REGION_LEN = atoi(argv[3]);
	int REGION_STEP_LEN = atoi(argv[4]);

	//FILE * bed_file = xopen(bed_fn_out, "w");
	//default parameters
	//int REGION_LEN = 500;
	//int REGION_STEP_LEN = 250;
	int MIN_K_size = 5;//4^5=
	int MAX_K_size = 200;
	int K_step_LEN = 1;

	//test file exist
	FILE* try_open = xopen(fa_fn_in, "r");
	fclose(try_open);
	faidx_t * c_ref_idx = reference_index_load(fa_fn_in);
	int N_seq = faidx_nseq(c_ref_idx);
	N_seq = MIN(26, N_seq);

	//header line
	fprintf(stdout, "%d\t%d\t%d\n", REGION_LEN, REGION_STEP_LEN, 0);
	std::set<std::string> kmerCountingSet;
	for(int i = chrID; i < chrID+1; i++){
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
					if (N_count_ANA(word) > 0 )
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

struct COMP_SV_infomation{
	std::string SV_info_ori;
	int pos_begin;
	std::vector<int8_t> tools_cover;
};

//------------------------------MAIN--------------------------------------//
int analysis_main(int argc, char *argv[])
{
	COMMAND_HANDLER ch;
	ch.set_main_command("analysis");
	ch.add_function("vntr_analysis", "", vntr_analysis);

	return ch.run(argc, argv);
	//not used functions
	ch.add_function("SV_in_simple_repeat", "", SV_in_simple_repeat);
	ch.add_function("bed_in_repeat_masker_region", "", bed_in_repeat_masker_region);

	ch.add_function("get_ref_region", "get_ref_region hg38.fa chr1:100000-2000000", get_ref_region);
	ch.add_function("ref_split", "Split reference file by chr_ID [input.fa]", ref_split);
	ch.add_function("getReverseStr", "Get Reverse Str of a DNA string [string (ACGT.....)]", getReverseStr);
	ch.add_function("genomic_contest_AnalysisNGS", "", genomic_contest_AnalysisNGS);
	ch.add_function("genomic_contest_AnalysisCCS", "", genomic_contest_AnalysisCCS);

	ch.add_function("fa_stat_full", "", fa_stat_full);
	ch.add_function("vcf_add_alt_string", "", vcf_add_alt_string);
	//ch.add_function("randomGenerateSV", "Generate 20000 random SV [sam_header_fn.sam, ref.fa, int_rand_seed]", randomGenerateSV);
	//ch.add_function("ROC_PR", "", analysis_ROC_PR);
	//ch.add_function("signalPattenAnalysis", "", signalPattenAnalysis);
	//ch.add_function("UCSC_mappibility_binary_dump", "", UCSC_mappibility_binary_dump);
	//ch.add_function("heatMapData", "", heatMapData);
	//ch.add_function("heatMapDataSum", "", heatMapDataSum);
	//ch.add_function("VENN_generate", "", VENN_generate);
	ch.add_function("combine_sort_vcf", "Combine vcf records from multiple files", combine_sort_vcf);
	ch.add_function("isize_count", "Count ISIZE for a bam file [input.bam]", isize_count);
	ch.add_function("bam2Fastq", "Convert bam into fastq [input.bam, output.fq]", bam2Fastq);
	ch.add_function("bamDump", "Dump first N record in a bam file [input.bam, ref.fa, output.bam, N]", bamDump);
	ch.add_function("gz_head", "Get the N characters from offset P of a XXX.gz file.[input.gz, N, P]", gz_head);
	ch.add_help_msg_back("zlib only support SEEK_SET, and read from begin P is at most 200M, not support from-where = 2: SEEK_END");
	ch.add_function("read_ACGT_analysis", "Analysis the acgt distribution of cram read", read_ACGT_analysis);
	ch.add_function("vcf_dump", "Filter and dump the SV items in VCF file using 'sample_ID'， ‘SV_TYPE‘ or ’chrID’", vcf_dump);
	ch.add_help_msg_back("[in_fn, out_fn, sample_ID, SV_TYPE, chrID], set 'ALL' for options the get all sample or SV type");

	//not used functions
	//	 if (argc <= 1) 							analysis_usage();
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
	//else if (strcmp(argv[1],"dump_ref_by_region") == 0) dump_ref_by_region(argv[2], argv[3]);
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
	//else 		 						 			analysis_usage();

	return 0;

}

