/*
 * Genome_context_analysis.hpp
 *
 *  Created on: 2025-4-30
 *      Author: fenghe
 */

#ifndef SVCALLING_CORE_ANALYSIS_GENOME_CONTEXT_ANALYSIS_HPP_
#define SVCALLING_CORE_ANALYSIS_GENOME_CONTEXT_ANALYSIS_HPP_


#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
#include "cpp_lib/cpp_utils.hpp"
extern "C"
{
#include "clib/utils.h"
#include "clib/bam_file.h"
#include "clib/vcf_lib.h"
}


struct Genomic_context_Analysis_handler{
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

	void get_block_ID_by_pos(int & id1, int &id2, int pos, int ana_block_length, int ana_block_step_size){
		id1 = (pos - ana_block_length + ana_block_step_size)/ana_block_step_size;
		id1 = MAX(0,id1);
		id2 = pos/ana_block_step_size;
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
		int MIN_GCA_sv_len = 30;
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
						if(p->length >= (unsigned)MIN_GCA_sv_len){
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
						if(p->length >= (unsigned)MIN_GCA_sv_len){
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

	//#define ANALYSIS_BLOCK_SIZE 500
	#define MIN_COMPACT_LEN 8
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

	void setR_Region(std::vector<R_region>& SV_region_l, int rid, int st, int end){
		SV_region_l.emplace_back();
		SV_region_l.back().chr_ID = rid;
		SV_region_l.back().st_pos = st;
		SV_region_l.back().ed_pos = end;
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

	//random generate 10000 deletion + 10000 insertion + 500 duplication
		bool pass_compact_filter(uint8_t *s, const int len) {
		int compact_str_len = 1; //AAAACCCT --> ACT
		for (int i = 1; i < len; i++)
			if (s[i - 1] != s[i])
				compact_str_len++;
		if (compact_str_len < MIN_COMPACT_LEN)
			return false;
		return true;
	}

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
};



#endif /* SVCALLING_CORE_ANALYSIS_GENOME_CONTEXT_ANALYSIS_HPP_ */
