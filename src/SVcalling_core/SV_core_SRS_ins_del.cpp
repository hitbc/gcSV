/*
 * SV_core_SRS_ins_del.cpp
 *
 *  Created on: 2020-4-29
 *      Author: fenghe
 */

#include <SVcalling_core/SV_core.hpp>
#include "../cpp_lib/cpp_utils.hpp"
#include "../SVcalling_core/Contig_polishing.hpp"
extern "C"
{
#include "../clib/utils.h"
#include "../clib/vcf_lib.h"
}
//******************************************filters************************************************/
#define MIN_MAPQ 20
static bool pass_mapQ_filter(bam1_t *br) {
	if (br->core.qual >= MIN_MAPQ)
		return true;
	if (bam_mate_map_qual(br) >= MIN_MAPQ && br->core.tid == br->core.mtid
			&& ABS(br->core.isize) < 20000)
		return true;
	return false;
}

//tools
inline bool is_mate_fwd(uint16_t flag) {
	return (!((flag & BAM_MATE_STRAND) != 0));
}

bool inline pos_within(int pos, int r_st, int r_ed) {
	if (pos >= r_st && pos <= r_ed)
		return true;
	return false;
}

void debug_print_flag(uint8_t flag) {
	std::cerr << "BAM_PAIRED " << (flag & 0x001) << "\t" << "BAM_PROPER_PAIR "
			<< (flag & 0x002) << "\t" << "BAM_UNMAPPED " << (flag & 0x004)
			<< "\t" << "BAM_MATE_UNMAPPED " << (flag & 0x008) << "\t"
			<< "BAM_STRAND " << (flag & 0x010) << "\t" << "BAM_MATE_STRAND "
			<< (flag & 0x020) << "\t" << "BAM_FIRST_READ " << (flag & 0x040)
			<< "\t" << "BAM_SECOND_READ " << (flag & 0x080) << "\t"
			<< "BAM_SECONDARY " << (flag & 0x100) << "\t" << "BAM_FILTER "
			<< (flag & 0x200) << "\t" << "BAM_DUPLICATE " << (flag & 0x400)
			<< "\t" << "BAM_SUPPLEMENTARY " << (flag & 0x800) << "\t" << "\t";
}

//*************************************Handle SA signals**********************************************/
void SV_CALLING_Handler::SRS_storeClipSignals(bool isClipAtRight, uint32_t ori_pos, uint8_t read_mapq) {
	//SVE c_sve;
	RefRegion ori_r(ref_handler->get_chr_ID(), ori_pos, ori_pos);
	if(isClipAtRight)
		ori_r.ed_pos += 10;
	else
		ori_r.st_pos -= 10;
	//when unmapped
	//if (rst.empty()) {	//clip not mapped, insertion
	RST::T isBegin = (isClipAtRight) ? RST::BEGIN :RST::END;
	SRS_sve[SIG::SH][SV::INS].emplace_back(isBegin, true, MIN(15, read_mapq), ori_r, ori_r);
}

#define MIN_MISMATCH_QUAL 20
#define MINI_KMER_LEN 10
extern uint64_t kmerMask[33];

void SV_CALLING_Handler::SRS_storeMismatchSignals(bam1_t *br, READ_record &c_r) {
	//s1: get high quality NM number:(filter)
	uint8_t * tseq = ref_handler->getRefStr(br->core.pos);
	uint8_t * qseq = SRS_read.getReadStr(c_r);
	uint8_t * qqual = SRS_read.getQualStr(c_r);

	mismatch_position.clear();

	int seq_i = 0;
	int ref_i = 0;

	uint32_t* bam_cigar = bam_get_cigar(br);
	uint32_t n_cigar = br->core.n_cigar;
	for (uint i = 0; i < n_cigar; ++i)
	{
		int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		switch (c_type){
		case CIGAR_MATCH: case CIGAR_SEQ_MATCH:
			for(int i = 0; i < c_size; i++, seq_i++, ref_i++)
				if((qseq[seq_i] != tseq[ref_i]) && (qqual[seq_i] >= MIN_MISMATCH_QUAL))
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
			return; //not handle reads with hard/soft clip
			break;
		default:	break;
		}
	}
	//check independence event number
	int independence_event_number = 1;
	int mismatch_position_size = mismatch_position.size();
	for(int i = 0; i < mismatch_position_size - 1 ; i++)
		if(mismatch_position[i] + 1 != mismatch_position[i + 1] && mismatch_position[i] != mismatch_position[i + 1])
			independence_event_number ++;

	if(independence_event_number < 3)
		return;

	if(false){
		int NM = 0;
		bam_get_num_tag(br, "NM", &NM);
		char * MD_string = bam_get_string_tag(br, "MD");
		for(int i = 0; i < br->core.l_qseq; i++){
			fprintf(stderr, "%c", qqual[i] + '#');
		}
		fprintf(stderr, " br->core.pos %d independence_event_number %d, MD_string %s mismatch_position is: ", br->core.pos, independence_event_number, MD_string);
		for(int i:mismatch_position){
			fprintf(stderr, "%d ", i);
		}
		for (uint i = 0; i < n_cigar; ++i)
		{
			int c_type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
			int c_size = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
			fprintf(stderr, "%d%c", c_size, "NMIDSSHXXXX"[c_type]);
		}

		fprintf(stderr, "\n");
	}

	int global_clip_point;
	//M1
	if(false){
		int min_kmer_idx = -1;
		//s2: get position of minimizer
		uint32_t kmer = bit2_nextKmer_init(qseq, MINI_KMER_LEN);
		uint32_t min_kmer  = MAX_uint32_t;
		int kmer_number = br->core.l_qseq - MINI_KMER_LEN + 1;
		uint64_t MASK = kmerMask[MINI_KMER_LEN];
		for(int i = 0; i < kmer_number; i++){
			kmer     = bit2_nextKmerMASK( qseq + i, kmer, MINI_KMER_LEN);
			if(min_kmer > kmer){
				min_kmer = kmer;
				min_kmer_idx = i;
			}
		}
		global_clip_point = min_kmer_idx + br->core.pos;
	}
	//M2
	{
		global_clip_point =  br->core.pos - (br->core.pos % 200);
	}

//	//S3: decide clip point and direction
//	int mis_before_minimizer = 0;
//	int total_mis_size = mismatch_position.size();
//	for(;mis_before_minimizer < total_mis_size && mismatch_position[mis_before_minimizer] < min_kmer_idx; mis_before_minimizer++);
//	bool clip_after_minimizer = (mis_before_minimizer + mis_before_minimizer < total_mis_size);

	RefRegion ori_r(ref_handler->get_chr_ID(), global_clip_point, global_clip_point + 200);
	//if(clip_after_minimizer)
	//	ori_r.ed_pos += 10;
	//else
	//	ori_r.st_pos -= 10;
	//RST::T isBegin = (clip_after_minimizer) ? RST::BEGIN :RST::END;
	RST::T isBegin =  RST::BEGIN;
	SRS_sve[SIG::SH][SV::INS].emplace_back(isBegin, true, MIN(3, br->core.qual), ori_r, ori_r);
}

#define MIN_CLIP_LEN 1
void SV_CALLING_Handler::SRS_handleSASignal(READ_record &c_r, bam1_t *br) {
	//uint8_t clip_string[MAX_SH_SIGNAL_STR_LEN + 1];//to store clip string
	//PART1： get left clip string
	if (c_r.soft_left >= MIN_CLIP_LEN) {
		uint32_t ori_pos = c_r.position;
		SRS_storeClipSignals(false, ori_pos, br->core.qual);
	}
	//PART1： get right clip string
	if (c_r.soft_right >= MIN_CLIP_LEN) {
		uint32_t ori_pos = c_r.position + c_r.read_l - c_r.soft_right - c_r.soft_left;
		SRS_storeClipSignals(true, ori_pos, br->core.qual);
	}
	int minNM = 3;
	if(c_r.NM_NUM >= minNM && c_r.soft_left == 0 && c_r.soft_right == 0){
		SRS_storeMismatchSignals(br, c_r);
	}
}

//combine signals in the try list, and get the max possible break point
int SV_CALLING_Handler::SRS_getTopPossibilityIdx(int r_min, int r_max, SVE_L &l, bool isR1, bool is_forward,
		float &max_poss, std::vector<float> &dis_bp_percent) {
	//clear possibility
	int dis_bp_percent_size = dis_bp_percent.size();
	int p_r_size = r_max - r_min + 2;
	p_r_size = MIN(5000, p_r_size);
	for (int i = 0; i < p_r_size; i++)
		possibility_r[i] = 0;

	if(is_forward){
		for (auto idx : cmb_try_list) {	//calculation
			int region_st = ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)) - r_min;
			//fprintf(stderr, " %d %d-%d \n", idx,  ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)),  ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)));
			for (int i = 0; i < dis_bp_percent_size; i++)
				possibility_r[region_st + i] += (dis_bp_percent[i]);
		}
	}else{
		for (auto idx : cmb_try_list) {	//calculation
			int region_ed = ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)) - r_min;
			//fprintf(stderr, " %d %d-%d \n", idx,  ((isR1)?(l[idx].r1.st_pos):(l[idx].r2.st_pos)),  ((isR1)?(l[idx].r1.ed_pos):(l[idx].r2.ed_pos)));
			for (int i = 0; i < dis_bp_percent_size; i++){
				//xassert(region_ed - i >= 0, "");
				if(region_ed - i < 0) continue;
				possibility_r[region_ed - i] += (dis_bp_percent[i]);
			}
		}
	}

	float max_possibility = 0;
	int max_index = 0;
	//debug code
	for (int i = 0; i < p_r_size; i++){
		//debug code:
		if (max_possibility < possibility_r[i]) {
			max_possibility = possibility_r[i];
			max_index = i;
		}
	}
	max_poss = possibility_r[max_index];
	return (r_min + max_index);
}

//combined signals to begin/end SV type
void SV_CALLING_Handler::SRS_single_type_sve_combine(bool print_log, SVE_L &l, int min_score_cutoff, SIG::T sigT, SV::T svt)//method 2:
		{
	float MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_SH;
	int MAX_ACCEPT_REGION;
	int BP_REGION;
	uint MIN_ACC_READ;
	std::vector<float> *breakpoint_distribution;
	if (sigT == SIG::DR) {
		MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_DR;
		MAX_ACCEPT_REGION = DR_bp_distribution.size();
		BP_REGION = 200;
		breakpoint_distribution = &(DR_bp_distribution);
		MIN_ACC_READ = 3;
	} else if (sigT == SIG::SH) {
		MIN_ACCEPT_POSS = MIN_ACCEPT_POSSIBILITY_SH;
		MAX_ACCEPT_REGION = 8;
		BP_REGION = 200;
		breakpoint_distribution = &(SH_bp_distribution);
		MIN_ACC_READ = 3;
	}else if (sigT == SIG::LRS) {
		//todo:::
		MIN_ACCEPT_POSS = 0.001;//todo::
		MAX_ACCEPT_REGION = 100;//todo::
		BP_REGION = 200;//todo::
		breakpoint_distribution = &(SH_bp_distribution);//todo::
		MIN_ACC_READ = 1;
	}

	std::sort(l.begin(), l.end(), SVE::cmp_by_position);

	auto sve_bg = l.begin();
	auto sve_ed = l.end();	//for each sve
	//for each sve
	cmb_store_tmp.clear();
	for (auto sve = sve_bg; sve < sve_ed; sve++) {
		if (sve->info.is_solid == RST::UNKNOWN)
			continue;	//already combined sve
		cmb_try_list.clear();
		RST::T is_solid = sve->info.is_solid;
		int max_score = 0;
		int r1_min = sve->r1.st_pos, r1_max = r1_min + 1;
		int r2_min = sve->r2.st_pos, r2_max = r2_min + 1;
		for (auto sve_try = sve;
				sve_try < sve_ed && sve_try->r1.region_overlap(r1_min, r1_max);
				sve_try++) {	//use sve as main, than try to combine
			if (is_solid != sve_try->info.is_solid)
				continue;	//condition 1
			if (sve_try->r2.chr_ID != sve->r2.chr_ID)
				continue;	//condition 3
			if (!sve_try->r2.region_overlap(r2_min, r2_max))
				continue;	//condition 3
			int try_idx = sve_try - sve_bg;
			cmb_try_list.emplace_back(try_idx);
			r1_min = MIN(r1_min, sve_try->r1.st_pos);
			r1_max = MAX(r1_max, sve_try->r1.ed_pos);
			r2_min = MIN(r2_min, sve_try->r2.st_pos);
			r2_max = MAX(r2_max, sve_try->r2.ed_pos);
			max_score += sve_try->info.getScore();
			sve_try->info.is_solid = RST::UNKNOWN;
		}
		//try to combine
		float max_possibility_r1, max_possibility_r2;
		//at least 1 LRS signals
		if(sigT == SIG::LRS){
			//store new sve into
			int64_t sum_r1 = 0, sum_r2 = 0;
			for(int i:cmb_try_list){
				sum_r1 += l[i].r1.st_pos + l[i].r1.ed_pos;
				sum_r2 += l[i].r2.st_pos + l[i].r2.ed_pos;
			}
			int64_t ave_r1 = sum_r1/(cmb_try_list.size())/2;
			int64_t ave_r2 = sum_r2/(cmb_try_list.size())/2;
			cmb_store_tmp.emplace_back(sve->r1.chr_ID, ave_r1, sve->r2.chr_ID,
					ave_r2, cmb_try_list.size(), 500, is_solid, BP_REGION);
		}
		//at least 3 normal clip signal or 8 multiple SNPs signals
		else if (cmb_try_list.size() >= MIN_ACC_READ
				&& max_score > sig_para.SVE_combine_min_score_step1
				&& (r1_max - r1_min < 5000) && (r2_max - r2_min < 5000)){
//
			bool search_forward = true;//handling R1 and is forward, for deletion, INV_1 and TRA
			if(svt == SV::INV_2) search_forward = false; //handling R1 and is reverse, for INV_2
			int r1_bp_position = SRS_getTopPossibilityIdx(r1_min, r1_max, l, true, search_forward,
					max_possibility_r1, *breakpoint_distribution);
			int max_accecpt_r1 = r1_bp_position;
			int min_accecpt_r1 = max_accecpt_r1 - MAX_ACCEPT_REGION;

			//handling R2 and reverse, for deletion
			search_forward = false;//handling R2 and reverse, for deletion, TRA and INV2
			if(svt == SV::INV_1) search_forward = true; //handling R2 and forward, for INV2
			int r2_bp_position = SRS_getTopPossibilityIdx(r2_min, r2_max, l, false, search_forward,
					max_possibility_r2, *breakpoint_distribution);
			int min_accecpt_r2 = r2_bp_position;
			int max_accecpt_r2 = min_accecpt_r2 + MAX_ACCEPT_REGION;

			if (max_possibility_r1 < MIN_ACCEPT_POSS
					&& max_possibility_r2 < MIN_ACCEPT_POSS)
				continue;
			int sve_n = 0;
			for (auto idx : cmb_try_list) {
				bool region_check_pass = false;
				if(svt != SV::INS){
					region_check_pass = (pos_within(l[idx].r1.st_pos, min_accecpt_r1, max_accecpt_r1)
							&& pos_within(l[idx].r2.ed_pos, min_accecpt_r2,	max_accecpt_r2));
				}
				else{
					region_check_pass = (is_solid == RST::BEGIN)? pos_within(l[idx].r1.st_pos, min_accecpt_r1, max_accecpt_r1) :
							pos_within(l[idx].r2.ed_pos, min_accecpt_r2, max_accecpt_r2);
				}
				if (region_check_pass)
					sve_n++;
				else
					l[idx].info.is_solid = is_solid;
			}
			//get new score:
			int score = (MAX(max_possibility_r1, max_possibility_r2))*2 / MIN_ACCEPT_POSS;
			if(score < min_score_cutoff) continue;
			//store new sve into
			cmb_store_tmp.emplace_back(sve->r1.chr_ID, r1_bp_position, sve->r2.chr_ID,
					r2_bp_position, sve_n, score, is_solid, BP_REGION);
		}
	}
	l.swap(cmb_store_tmp);
}

#define MIN_REALIGNMENT_SCORE 50
void SV_CALLING_Handler::SRS_handleUMSignal(bam1_t *br) {

	uint8_t *query = UMQueryBuff;
	const int read_len = br->core.l_qseq;
	if (get_bam_seq_bin(0, read_len, query, br) == false)//get string failed, return
		return;
	if (!BAM_handler::pass_compact_filter(query, read_len))// || !BAM_handler::passComplexFilter(query, read_len))//filter for string
		return;
	//store to unmapped read list
	SRS_read.storeReadUM(br, query);
}

//just equal : !bam_is_DR_signal(); but is more understandable
bool bamIsNormalDrPair(bam1_t *br, int i_size_min, int i_size_MAX) {
	bool direction = bam_is_fwd_strand(br);
	bool mate_dir = bam_is_mate_fwd_strand(br);
	int ABSisize = ABS(br->core.isize);

	if (br->core.tid == br->core.mtid && direction != mate_dir
			&& ABS(ABSisize) < i_size_MAX && ABS(ABSisize) > i_size_MAX) {
		if (direction == FORWARD && br->core.pos <= br->core.mpos)
			return true;	//normal condition 1
		if (direction == REVERSE && br->core.pos <= br->core.mpos)
			return true;	//normal condition 2
	}
	return false;
}

bool bam_aligned_analysis(bam1_t *b, int *clip_left, int *clip_right, int *gap_mismatch_inside){
	if(b->core.n_cigar == 0)//CIGAR un-available
		return false;
	int gap_number = 0;
	uint32_t* bam_cigar = bam_get_cigar(b);
	for (uint i = 0; i < b->core.n_cigar; ++i){
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		if(type == CIGAR_INSERT || type == CIGAR_DELETE){
			gap_number += (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		}
	}
	int NM = 0;
	bam_get_num_tag(b, "NM", &NM);
	*gap_mismatch_inside = MAX(gap_number, NM);
	//NM = MIS + INS + DEL
	int begin_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
	int end_type   = (int)(1 + (bam_cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK));
	*clip_left = 0; *clip_right = 0;
	if(	begin_type == CIGAR_SOFT_CLIP || begin_type == CIGAR_HARD_CLIP)
		*clip_left = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
	if(	end_type == CIGAR_SOFT_CLIP ||	end_type == CIGAR_HARD_CLIP)
		*clip_right = (bam_cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);

	return (*clip_left > 0 || *clip_right > 0 || *gap_mismatch_inside >= 2);
}

void SV_CALLING_Handler::SRS_GET_SIGNALS_FROM_READs(){
		Bam_file *c_b = &SRS_read.file;
		R_region region;
		ref_handler->get_cur_region()->toR_region(region);
		resetRegion_ID(c_b, &region);	//reset region
		read_counter = 0;

		int clip_left, clip_right, gap_mismatch_inside; //the results of analysis alignment
		rdc.set_dc_st(region.st_pos);
		rdc.clear_read_depth_list();

		SRS_read.tr_loader.clear();

		while (bam_next(c_b)) {
			bam1_t *br = &(c_b->_brec);
//			if(br->core.pos >= 108732977 && br->core.pos <= 108733329){
//				fprintf(stderr, " " );
//			}
			if(bam_is_duplicate(br))
				continue;
			rdc.add_read_depth_item(read_counter, br->core.pos);			//depth counter:
			//debug code:
			read_counter++;
			if(read_counter % 10000 == 0) fprintf(stderr, "%d\n", read_counter);
			//basic filter
			if(bam_is_secondary(br))
				continue;
			if(bam_is_supplementary(br))
				continue;
			//WARNING: the unmapped reads will always not pass the MAPQ filter, therefore it should be before the mapQ filter
			if(bam_is_unmapped(br)){
				SRS_handleUMSignal(br);
				continue;
			}//unmapped read
			if (!pass_mapQ_filter(br))
				continue;
			//add new signal type:
			//SR signal
			if (bam_aligned_analysis(br, &clip_left, &clip_right, &gap_mismatch_inside)) {
//				if(br->core.pos >= 108732977 && br->core.pos <= 108733329){
//					fprintf(stderr, "%d %d %d %d\n ",br->core.pos, clip_left, clip_right, gap_mismatch_inside);
//				}
				if(SRS_read.storeReadSR(br, clip_left, clip_right, gap_mismatch_inside))//store signal reads
					SRS_handleSASignal(SRS_read.read_list.back(), br);
			}
			//DR signal
			int middle_size = br->core.l_qseq - clip_left - clip_right;	//length of match bases
			if (bam_is_DR_signal(br, sig_para.insert_size_min, sig_para.insert_size_max)){
				SRS_handleDRSignal(&(br->core), middle_size);
			}

			//other reads pair not normal paired
			if(br->core.mtid != br->core.tid || ABS(br->core.pos - br->core.mpos) > 200000){
				SRS_read.tr_loader.trans_list.emplace_back(br);
			}
		}

        fprintf(stderr, "Signal process Done\n" );
        if(false){
            fprintf(stderr, "r:read.read_list.show\n");
            for(READ_record & r:SRS_read.read_list){
                r.show();
            }
        }
        if(false){
            fprintf(stderr, "show all single signals\n");
            std::cerr << SIG::STR[SIG::SH] << std::endl;
            for (int svt = 0; svt < SV::LEN; svt++) {
                SVE_L &svl = SRS_sve[SIG::SH][svt];
                if (svl.size() == 0)
                    continue;
                std::cerr << SV::STR[svt] << std::endl;
                for(unsigned int i = 0; i < svl.size(); i++){
                    svl[i].setSignalType((SV::T)svt, SIG::SH);
                    svl[i].show();
                }
            }
        }
		return;
	}

void SV_CALLING_Handler::SRS_load_realignment_read_from_file(){
	Bam_file *c_b = &SRS_read.re_alignment_file;
	if(c_b->_hdr == NULL)
		return;
	R_region region;
	ref_handler->get_cur_region()->toR_region(region);
	resetRegion_ID(c_b, &region);	//reset region
	while (bam_next(c_b)) {
		bam1_t *br = &(c_b->_brec);
		if(bam_is_duplicate(br))
			continue;
		uint8_t *query = UMQueryBuff;
		const int read_len = br->core.l_qseq;
		if (get_bam_seq_bin(0, read_len, query, br) == false)//get string failed, return
			continue;
		if (!BAM_handler::pass_compact_filter(query, read_len))// || !BAM_handler::passComplexFilter(query, read_len))//filter for string
			continue;
		//store to unmapped read list
		//get the position for store UM signal:
		int32_t read_store_pos = br->core.pos;
		//reading the SV tag of original reads
		char* SV_tag_char = bam_get_string_tag(br, "SV");
		if(SV_tag_char == NULL)
			xassert(0, "Without SV Tag!");
		//SV:Z:11131_9_26939007_517_INS_svim
		split_string(realn_sv_field_split, SV_tag_char, "_");
		if(realn_sv_field_split[4].compare("INS") == 0){
			int32_t SV_pos = atoi(realn_sv_field_split[2].c_str());
			read_store_pos = SV_pos;
		}

		SRS_read.storeReadCore(Read_type::UM, br->core.l_qseq, query, bam_get_qual(br), 0, 0, br->core.flag, read_store_pos, 0, 0, 0);
	}
	//sort UM reads signals
	std::sort(SRS_read.um_read_list.begin(), SRS_read.um_read_list.end(), READ_record::cmp_position);
}

void SV_CALLING_Handler::SRS_handleDRSignal(bam1_core_t *core, int middle_size) {
	bool dir = (((core->flag & BAM_STRAND) == 0)), mate_dir = ((core->flag & BAM_MATE_STRAND) == 0);
	int absIsize = (int) ABS(core->isize);
	SV::T t = SV::UNKNOWN;
	RST::T is_begin = RST::UNKNOWN;
	if (core->tid == core->mtid && dir != mate_dir && absIsize < sig_para.max_del_dup_length) {//insertion/ tandem duplication/ deletion
		is_begin = (dir == FORWARD) ? RST::BEGIN : RST::END;
		bool normalOri = ((dir == FORWARD && core->pos <= core->mpos) || (dir == REVERSE && core->pos >= core->mpos));
		if (normalOri){
			if(absIsize > sig_para.insert_size_max)			t = SV::DEL;
			else if(absIsize < sig_para.insert_size_min)	t = SV::DUP;
			else 											t = SV::UNKNOWN;
		} else 												t = SV::DUP;
	} else if (core->tid == core->mtid && dir == mate_dir) {
		is_begin = (core->isize > 0) ? RST::BEGIN : RST::END;
		t = (dir == FORWARD) ? SV::INV_1 : SV::INV_2;
	} else if (core->tid != core->mtid) {
		is_begin = (core->tid < core->mtid) ? RST::BEGIN : RST::END;
		t = (dir != mate_dir) ? SV::TRA : SV::TRA_INV;
	}
	if (t != SV::UNKNOWN)
		SRS_sve[SIG::DR][t].emplace_back(core, dir, mate_dir,
				sig_para.insert_region_len, middle_size, is_begin);
}

//**********************************MAIN*********************************************************/

//combine begin and end of each signal types, to be solid or unstable; 	//when it is solid SV, is begin set to be "SV_SOLID"
void SRS_sve_begin_end_combine(bool print_log, SVE_L & l, int MIN_SOLID_S, int MIN_READ_NUM){
	std::sort(l.begin(), l.end(), SVE::cmp_by_position);
	int store_index = 0;
	auto sve_ed = l.end();
	for(auto sve = l.begin();sve != sve_ed; sve++){
		if(sve->info.is_solid == RST::UNKNOWN) continue;
		//search pair:
		for(auto sve_try = sve + 1; sve_try < sve_ed && (sve->r1.region_overlap(sve_try->r1)); sve_try++){//use sve as main, than try to combine
			if(sve->info.is_solid == sve_try->info.is_solid) continue;//condition 1: only combine begin+end
			if(!sve->r2.region_overlap(sve_try->r2))	continue;//condition 3
			sve->r1.Combine(sve_try->r1, true);
			sve->r2.Combine(sve_try->r2, true);
			sve->info.is_solid = RST::SOLID;
			sve->info.combine_score(sve_try->info);
			sve_try->info.is_solid = RST::UNKNOWN;
			break;
		}
		bool pass_filter = true;
		if(pass_filter && sve->info.is_solid  < RST::SOLID && (sve->info.getScore() < MIN_SOLID_S*2 || sve->info.getReadNum() < MIN_READ_NUM*2)) pass_filter = false;
		if(pass_filter && sve->info.is_solid == RST::SOLID && (sve->info.getScore() < MIN_SOLID_S || sve->info.getReadNum() < MIN_READ_NUM)) pass_filter = false;//condion 1: solid but skip TODO:: 10% of depth

		if(!pass_filter){
			if(print_log) std::cerr << "Not pass filter:\t" << *sve;
			continue;
		}
		else{
			if(print_log) std::cerr << "Pass filter:\t" << *sve;
		}

		l[store_index++] = sve[0];
	}
	l.resize(store_index);
}

void SV_CALLING_Handler::SRS_combine_duplication(SVE_L & l) {
	auto sve_bg = l.begin();
	auto sve_ed = l.end();	//for each sve
	//for each sve
	cmb_store_tmp.clear();
	for (auto sve = sve_bg; sve < sve_ed; sve++) {
		if (sve->info.is_solid == RST::UNKNOWN)
			continue;	//already combined sve
		cmb_try_list.clear();
		for (auto sve_try = sve + 1; sve_try < sve_ed && sve->r1.region_overlap(sve_try->r1); sve_try++) {	//use sve as main, than try to combine
			if(sve->SVE_combine_UNION(*sve_try)){
				sve_try->info.is_solid = RST::UNKNOWN;
			}
		}
		cmb_store_tmp.emplace_back(*sve);
	}
	l.swap(cmb_store_tmp);
}


void SV_CALLING_Handler::SRS_CLUSTERING_AND_COMBINE_signals() {
	bool print_log = true;
	//combine the DR signals, signals combined in two steps, firstly, signals will be clustered by overlap, then, the break points will be calculated using those signals
	std::cerr << SIG::STR[SIG::DR] << std::endl;
	for (int svt = 0; svt < SV::LEN; svt++) {
		SVE_L &svl = SRS_sve[SIG::DR][svt];
		if (svl.size() == 0)
			continue;
		std::cerr << SV::STR[svt] << std::endl;
		SRS_single_type_sve_combine(print_log, svl, 2, SIG::DR, (SV::T)svt);//10% * read depth

		if (svt == SV::DEL) //DEL
			SRS_sve_begin_end_combine(print_log, svl, sig_para.SVE_MIN_SOLID_SCORE, sig_para.SVE_MIN_READ_NUM);
		else //DUP and BND
			SRS_sve_begin_end_combine(print_log, svl, sig_para.SVE_MIN_SOLID_SCORE*1.5, sig_para.SVE_MIN_READ_NUM*1.5);
		for(unsigned int i = 0; i < svl.size(); i++)
			svl[i].setSignalType((SV::T)svt, SIG::DR);
	}

	//combined the SH signals, signals combined in two steps
	// firstly, signals will be clustered by overlap, then, the break points will be calculated using those signals
	std::cerr << SIG::STR[SIG::SH] << std::endl;
	for (int svt = 0; svt < SV::LEN; svt++) {
		SVE_L &svl = SRS_sve[SIG::SH][svt];
		if (svl.size() == 0)
			continue;
		std::cerr << SV::STR[svt] << std::endl;
		SRS_single_type_sve_combine(print_log, svl, 2, SIG::SH, (SV::T)svt);
		SRS_sve_begin_end_combine(print_log, svl, sig_para.SVE_MIN_SOLID_SCORE, sig_para.SVE_MIN_READ_NUM);

		for(unsigned int i = 0; i < svl.size(); i++)
			svl[i].setSignalType((SV::T)svt, SIG::SH);
		//delete duplication:
		SRS_combine_duplication(svl);
	}
	//std::sort(sve[SIG::SH][SV::INS].begin(), sve[SIG::SH][SV::INS].end(), SVE::cmp_by_position);

	//add repeat region as signals
	//load signal region by reference information
	R_region region;
	ref_handler->get_cur_region()->toR_region(region);
	rlr->search(region.chr_ID, region.st_pos, region.ed_pos, rlr_result);
	for(RefLocalRepeatItemR &rr : rlr_result){
		int middle_pos;	int range;
		rlr->get_middle(rr, middle_pos, range);
		SRS_sveVNTR.emplace_back(rr.chrID, middle_pos, rr.chrID, middle_pos, 99, 99, RST::SOLID, range);
		SRS_sveVNTR.back().setSignalType(SV::INS, SIG::SH);
		if(print_log) std::cerr << "Repeat region adding:\t" << SRS_sveVNTR.back();
	}
	std::sort(SRS_sveVNTR.begin(), SRS_sveVNTR.end(), SVE::cmp_by_position);
	SRS_combine_duplication(SRS_sveVNTR);
}
void printContig(FILE * output, int ori_read_number, AssemblyContig & contig, int contig_ID, std::vector<ASS_reads_info> &ass_read_list){
	//basic part
	if(true){
		fprintf(output,
				"\ncontig_ID: [%d] "
				"word length: [%d] "
				"CONTIG size: [%ld]\n "	"CONTIG seq:  [%s]\n "
				"supportReads [%ld] "
				"ending_reason: [%d %d]"
				"new_support_read[%d] "
				,
				contig_ID,
				contig.kmerLength, contig.seq.size(), contig.seq.c_str(), contig.supportReads.size(), contig.ending_reason[0], contig.ending_reason[1], contig.new_support_read);
	}
	//supplementary information
	if(true){
		fprintf(output,
				"contig_begin_offset: [%d] "
				"seedCount: [%d] "
				"rejectReads[%ld] "
				,
				contig.ass_begin_offset_in_contig,
				contig.seedReadCount,  contig.rejectReads.size());
		fprintf(output, "supportReads: ");
		for (auto &r : contig.supportReads){
			if((int)r < ori_read_number){
				fprintf(output, "[ %d %d %d]\t", r,
						ass_read_list[r].read_in_ref_offset,
						ass_read_list[r].read_list_index);
			}
		}
	}

	fprintf(output, "\n");
}
//
void SV_CALLING_Handler::SRS_get_suggention_alignment_position_list(AssemblyContig & contig, int ori_read_number, std::vector<std::string> &read_list, std::vector<ASS_reads_info> &ass_read_list, FILE* log_f ){

	const char *contig_seq = contig.seq.c_str();
	int contig_seq_len = contig.seq.size();
	//step1: get all the read positions at assembled CONTIG
	//set read position for each actions
	contig.remove_read_set.clear();
	for (auto &ca : contig.actions)
		if(ca.read_ID < ori_read_number)
			ca.set_read_pos(read_list[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.kmerLength, contig.remove_read_set, ass_read_list[ca.read_ID].read_in_ref_offset);

	//step2: find the suggestion alignment position of CONTIG at reference from each read position
	suggest_st_pos_map.clear();
	if(false){
		for(auto &s : contig.remove_read_set){
			fprintf(stderr, "remove_read_set %d\n", s);
		}
	}
	//get contig coverage
	if((int)contig_depth.size() < contig_seq_len)	contig_depth.resize(contig_seq_len);
	memset(&(contig_depth[0]), 0, contig_seq_len*sizeof(uint16_t));

	for (auto &ca : contig.actions){
		ca.wrong_base = 999;
		if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
		if(ass_read_list[ca.read_ID].signal_type != Read_type::SR && ass_read_list[ca.read_ID].signal_type != Read_type::TL) continue;
		if(contig.remove_read_set.find(ca.read_ID) != contig.remove_read_set.end()) continue;
		//set coverage:
		const char* read_seq = read_list[ca.read_ID].c_str();
		int st_pos_ref = ca.position_in_contig - contig.ass_begin_offset_in_contig - ca.position_read;
		int ed_pos_ref = st_pos_ref + read_list[ca.read_ID].size();
		int st_pos_read = 0; if(st_pos_ref < 0){	st_pos_read -= st_pos_ref;st_pos_ref = 0;}
		ed_pos_ref = MIN(contig_seq_len, ed_pos_ref);
		ca.wrong_base = 0;

		for(int i = st_pos_ref; i < ed_pos_ref && ca.wrong_base <= MAX_WRONG_BASE; i++, st_pos_read++){
			if(read_seq[st_pos_read] != 'N' && contig_seq[i] != read_seq[st_pos_read])
				ca.wrong_base++;
		}

		if(ass_read_list[ca.read_ID].signal_type == Read_type::SR){
			if(ca.wrong_base <= MAX_WRONG_BASE){
				st_pos_read = 0;
				for(int i = st_pos_ref; i < ed_pos_ref; i++, st_pos_read++){
					if(contig_seq[i] == read_seq[st_pos_read]) contig_depth[i] ++;
				}
				std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
				if(it != suggest_st_pos_map.end())	it->second += 2;
				else suggest_st_pos_map[ca.suggest_contig_offset_in_ref] = 2;
			}else{
				std::map<int, int >::iterator it = suggest_st_pos_map.find(ca.suggest_contig_offset_in_ref);
				if(it != suggest_st_pos_map.end())	it->second--;
			}
		}
	}
	//if no result is there, simple selected the first read position
	if(suggest_st_pos_map.empty())	suggest_st_pos_map[contig.actions[0].suggest_contig_offset_in_ref] = 2;

	//step3: find the max suggest position
	//get max suggest alignment position
	//try simple merge:
	int max_suggention = 0; int max_count = 0;
	for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++)
		if(max_count < it->second) {max_suggention = it->first; max_count = it->second;}

	//step4: find the suggestion list
	//if a suggestion covers most of read, [unique_alignment_pos] is true and use it as alignment suggestion, otherwise store all suggestions in a list [suggent_pos_list]
	// when no result or result > 1: [unique_alignment_pos] = false
	//suggest_coverage only use when suggent_pos_list.size() > 1;
	contig.SUGGEST_pos_list.clear();
	if(false) fprintf(stderr, "suggest_st_pos_map BG\n");
	for(std::map<int, int >::value_type &sug : suggest_st_pos_map){
		if(false) fprintf(stderr, "suggest_st_pos_map:  sug.first %d , sug.second %d \n", sug.first, sug.second);
		if(sug.second >= 4){
			contig.SUGGEST_pos_list.emplace_back(sug.first);
			SUGGEST_POS_LIST_ITEM & sug = contig.SUGGEST_pos_list.back();
			int c_suggest_pos = sug.suggest_pos;
			for(auto &ca : contig.actions){
				if(ca.read_ID < ori_read_number && ca.suggest_contig_offset_in_ref == c_suggest_pos)
					sug.add_read_start_pos(ca.position_in_contig - ca.position_read, ca.wrong_base);
			}
			sug.count_ave_read_start();
			if(false) fprintf(stderr, "sug.suggest_pos %d sug.high_wrong_base_read_number: sug.low_wrong_base_read_number %d %d \n", sug.suggest_pos, sug.high_wrong_base_read_number, sug.low_wrong_base_read_number);
			if(sug.low_wrong_base_read_number < 1 || sug.high_wrong_base_read_number > sug.low_wrong_base_read_number * 4)//only 0 support reads or less than 20% right supported reads
				contig.SUGGEST_pos_list.erase(contig.SUGGEST_pos_list.end() - 1);
		}
	}

	std::sort(contig.SUGGEST_pos_list.begin(), contig.SUGGEST_pos_list.end(), SUGGEST_POS_LIST_ITEM::cmp_by_read_position);
	//std::sort(SUGGEST_pos_list.begin(), SUGGEST_pos_list.end(), SUGGEST_POS_LIST_ITEM::cmp_by_read_count);

	//print suggestion alignment range information
	if(false){
		//all possible ranges
		fprintf(log_f, "All suggestions:\t");
		for(std::map<int, int >::iterator it = suggest_st_pos_map.begin(); it != suggest_st_pos_map.end(); it++) fprintf(log_f, "[SUG: %d CN: %d]\t", it->first, it->second);
		fprintf(log_f, "\n Used suggestions: \t");
		for(SUGGEST_POS_LIST_ITEM & sug : contig.SUGGEST_pos_list) sug.printf(log_f);
		//possible true ranges
		if(contig.SUGGEST_pos_list.empty()) 			{fprintf(log_f, "\nNO suggestion\n"); }
		else if(contig.SUGGEST_pos_list.size() == 1) 	{fprintf(log_f, "\nUNIQUE: [MAX_SUG: %d]\n", max_suggention);}
		else									{fprintf(log_f, "\nMULTY\n"); }
		if(true){
			for (auto &ca : contig.actions)
				if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
				{
					ass_read_list[ca.read_ID].print(log_f, 0);
					ca.print(log_f, true);
				}
			fprintf(log_f, "\n");
		}
	}
}

uint8_t charToDna5n[] =
{
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,//4 the last but one
    /*   		 A     C           G                    N */
    /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,//'Z'
    /*          	      T */
    /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
    /*   		 a     c           g                    n */
    /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*          	      t */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

void store_bin_from_char(const char * contig_seq, int contig_seq_len, std::vector<uint8_t> &bin_contig){
	xassert(nullptr != contig_seq, "");
	//store bin contig
	bin_contig.resize(contig_seq_len);
	for (int i = 0; i < contig_seq_len; ++i)
		bin_contig[i] = charToDna5n[(uint8_t)contig_seq[i]];
}

void store_bin_contig(std::string &contig_string, std::vector<uint8_t> &bin_contig){
	store_bin_from_char(contig_string.c_str(), contig_string.size(), bin_contig);
}

//using "bin_contig" to store contig sequence
void SV_CALLING_Handler::SRS_alignment_and_get_var(bool print_log, int ab_idx, int contig_ID, int suggest_ref_st_pos, int contig_seq_len, int ref_region_length){
	//alignment contig into reference
	int suggest_ref_ed_pos = suggest_ref_st_pos + contig_seq_len + region_addition_load;
	int suggest_contig_st_pos = 0;
	if(suggest_ref_st_pos < 0){	suggest_contig_st_pos = (- suggest_ref_st_pos);	suggest_ref_st_pos = 0;}
	suggest_ref_st_pos = MAX(0, suggest_ref_st_pos);//get 15 addition alignment region at begin
	suggest_ref_ed_pos += 60; suggest_ref_ed_pos = MIN(suggest_ref_ed_pos, ref_region_length);//get 60 additional alignment region at end
	bool contig_out_range = (suggest_ref_ed_pos < suggest_ref_st_pos + 20 || suggest_contig_st_pos > contig_seq_len); //when true region length less then 20 bp, it is out of range
	//if(output_assembly_contig) fprintf(detail_output, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");
	if(!contig_out_range){
		//print alignment range
		ca.align_non_splice_default_ref(&(bin_contig[0]) + suggest_contig_st_pos, contig_seq_len - suggest_contig_st_pos,
				suggest_ref_st_pos, suggest_ref_ed_pos, 400);
		//if(output_assembly_contig) print_contig_aln_ori_cigar(ca.ez, detail_output, suggest_contig_st_pos);
		int ref_adj_size = ca.adjustCIGAR();
		suggest_ref_st_pos += ref_adj_size;
		fprintf(stderr, "[aln_ref_region %d ~ %d ] [aln contig st %d ] %s", suggest_ref_st_pos , suggest_ref_ed_pos, suggest_contig_st_pos, contig_out_range?("[ contig_out_range]"):"");

		if(print_log){
			ca.printf_alignment_detail(stderr, suggest_ref_st_pos, &contig_depth[0] + suggest_contig_st_pos, 0);
			fprintf(stderr, "ref_pos_enough_match_base: %d\t", 0);
		}
		//contig alignment filter:
		int match_size;
		 int contig_NM = ca.get_contig_NM(suggest_ref_st_pos, match_size);
		 fprintf(stderr, "Contig filter: contig_NM %d match_size %d \n", contig_NM, match_size);
		 if(contig_NM * 20 > match_size){
			 fprintf(stderr, "Contig filter: Too mant NM: contig_NM %d error rate: %f%%(over 5%)\n", contig_NM, (float)contig_NM/match_size*100);
			 return;
		 }

		ca.get_canditate_SVs(print_log, SV_result_TMP_buff, MIN_sv_len, suggest_ref_st_pos, suggest_SV_length, region_ref_global_position);
	}
}

bool SV_CALLING_Handler::SRS_SV_region_depth_filter(SVE & c_sve){
	float r1_max_depth, r2_max_depth;
	float r1_depth = rdc.get_ave_depth(c_sve.r1.st_pos, c_sve.r1.ed_pos, &r1_max_depth);
	float r2_depth = rdc.get_ave_depth(c_sve.r2.st_pos, c_sve.r2.ed_pos, &r2_max_depth);
	fprintf(stderr, "r1:AVE(max) depth: [%f, (%f)], r2:AVE(max) depth: [%f, (%f)] ", r1_depth, r1_max_depth, r2_depth, r2_max_depth);
	std::cerr << "r1: " <<  c_sve.r1 << ", r2:" << c_sve.r2 << std::endl;
	if(r1_max_depth > 15 || r2_max_depth > 15){
		fprintf(stderr, "Assembly will be skipped because over depth in this region\n");
		return false;
	}
	if(r1_depth > 5 || r2_depth > 5){
		fprintf(stderr, "Assembly will be skipped because over depth in this region\n");
		return false;
	}
	return true;
}

void RunLengthCompect(std::string &S, std::vector<uint8_t> &count){
	std::string t;
	std::swap(t, S);
	if(t.empty()){
		return;
	}
	char c = t[0];
	S.push_back(c);
	count.emplace_back(1);
	for(uint i = 1; i < t.size(); i++){
		if(t[i] != c){
			S.push_back(t[i]);
			c = t[i];
			count.emplace_back(1);
		}
		else
			count.back() += 1;
	}
}

void show_read_string_RLC(const char * title,  std::string & reads_str, ASS_reads_info & ai, int rid){
	if(reads_str.size() != 0){
		fprintf(stderr, "%s [rid: %d, full begin %d, full end %d ai.soft_left %6d , ai.soft_right %6d , len %ld ]:%s\n", title, rid,
				ai.LRS_full_begin,ai.LRS_full_end,
				ai.soft_left, ai.soft_right,
				reads_str.size(), reads_str.c_str());
	}
}

void RunLengthCompect_All(bool run_length_all_read, bool print_log, Ass_Block &ass_block, std::vector< std::string> &reads_str, std::vector<std::vector<uint8_t>> &RunLength_count_read){
	if(run_length_all_read){
		//run-length compact all reads:
		uint r_size = reads_str.size();
		RunLength_count_read.resize(r_size);
		for(uint rid = 0; rid < r_size; rid++){
			if(print_log) show_read_string_RLC("ORI", reads_str[rid], ass_block.ass_read_list[rid], rid);
			RunLengthCompect(reads_str[rid], RunLength_count_read[rid]);
			//show:
			if(print_log){
				show_read_string_RLC("CMP", reads_str[rid], ass_block.ass_read_list[rid], rid);
				fprintf(stderr, "CNT:");
				for(uint8_t i: RunLength_count_read[rid])
					fprintf(stderr, "%c", i + '0');
				fprintf(stderr, "\n");
			}
		}
	}else{
		//run-length compact all reads:
		if(print_log){
			for(uint rid = 0; rid <  reads_str.size(); rid++)
				show_read_string_RLC("ORI", reads_str[rid], ass_block.ass_read_list[rid], rid);
		}
	}
}

void RunLengthCompect_BIN(uint8_t *from, int len, std::vector<uint8_t> &to){
	to.clear();
	if(len == 0){
		return;
	}
	uint8_t c = from[0];
	to.push_back(c);
	for(uint i = 1; i < len; i++){
		if(from[i] != c){
			to.push_back(from[i]);
			c = from[i];
		}
	}
}

//load NOT only signal reads, but ALL reads in this region(not only signal reads in this region)
void SV_CALLING_Handler::LRS_assembly_load_read(bool print_log, RefRegion &main, Ass_Block & a_b, int min_LRS_read_length, BAM_handler *cur_LRS_read){//buff to store contig
	RefRegion &r = main;
	//load all reads iin this regions:
	r.show();
	//combine the LRS signals with NGS signals
	//load part of LRS reads
	int read_idx = -1;
	BAM_handler &cur = *cur_LRS_read;
	bool is_full_begin = true;
	bool is_full_end = true;
	for(READ_record &rr :cur.read_list){
		std::string s;
		read_idx++;
		cur.load_read_LRS(s, rr, r.st_pos, r.ed_pos, is_full_begin, is_full_end);
		if(!s.empty() && s.size() >= min_LRS_read_length)
			a_b.add_read_LRS(s, read_idx, r.chr_ID, r.st_pos, is_full_begin, is_full_end, rr.soft_left, rr.soft_right);
	}
}

#define MAX_CONTIG_LEN 3000000
struct ASM_CMB_ITEM{
	int region_AT_read;
	int ref_span_pos;
	int direction;
};

struct ASM_CMB_LIST{
	std::string r_name;
	int primary_chrID;
	int primary_pos;
	std::string s_left;
	std::string s_right;
	std::vector<ASM_CMB_ITEM> l_left;
	std::vector<ASM_CMB_ITEM> l_right;
	void store_new(std::string &r_name, int primary_chrID, int primary_pos){
		this->r_name = r_name;
		this->primary_chrID = primary_chrID;
		this->primary_pos = primary_pos;
		l_left.clear();
		l_right.clear();
	}

	void store_one_data(int region_AT_read, int ref_span, bool is_left_side, std::string & s, int direction){
		if(region_AT_read == -1)
			return;
		std::vector<ASM_CMB_ITEM> * l_cur = &(l_left);
		if(false == is_left_side)   l_cur = &(l_right);
		(*l_cur).emplace_back();
		(*l_cur).back().region_AT_read = region_AT_read;
		(*l_cur).back().ref_span_pos = ref_span;
		(*l_cur).back().direction = direction;
		if(true == is_left_side)	this->s_left = s;
		else 						this->s_right = s;
	}

	void get_max_span(int & max_idx_left, int & max_idx_right){
		max_idx_left = -1;
		max_idx_right = -1;
		int MAX_span = 0;
		int MAX_read_len = 0;
		for(uint i = 0; i < l_left.size(); i++){
			for(uint j = 0; j < l_right.size(); j++){
				int ref_span_len = l_right[j].ref_span_pos - l_left[i].ref_span_pos;
				int read_len = l_right[j].region_AT_read - l_left[i].region_AT_read;
				if(read_len > MAX_CONTIG_LEN) continue;
				if(read_len < 0) continue;
				if(l_right[j].direction != l_left[i].direction) continue;

				if(MAX_span < ref_span_len){
					MAX_span = ref_span_len;
					max_idx_left = i;
					max_idx_right = j;
					MAX_read_len = read_len;
				}
			}
		}

		fprintf(stderr, "max_idx_left %d, max_idx_right %d MAX_span %d , MAX_read_len  %d \n", max_idx_left, max_idx_right, MAX_span, MAX_read_len);
	}

	bool pairing_and_get_string(Bam_file *c_b, std::string &s_store, int &region_AT_read_left, int &region_AT_read_right, int &ref_span_left, int &ref_span_right){
		//try to get the MAX span:
		s_store.clear();
		int max_idx_left; int max_idx_right;
		get_max_span(max_idx_left, max_idx_right);
		if(max_idx_left == -1 || max_idx_right == -1){
			fprintf(stderr, "FAIL, NO contig generated \n");
			return false;
		}
		region_AT_read_left = l_left[max_idx_left].region_AT_read;
		region_AT_read_right = l_right[max_idx_right].region_AT_read;

		ref_span_left = l_left[max_idx_left].ref_span_pos;
		ref_span_right = l_right[max_idx_right].ref_span_pos;

		int region_direction = l_left[max_idx_left].direction;
		int load_length = region_AT_read_right - region_AT_read_left;

		if(load_length <= 0 || load_length > MAX_CONTIG_LEN)//AT most 3M
			return false;
		//load the string
		R_region region;
		region.chr_ID = primary_chrID;
		region.st_pos = primary_pos + 20;
		region.ed_pos = primary_pos + 21;
		//load the original read string from the primary
		resetRegion_ID(c_b, &region);	//reset region
		while (bam_next(c_b)) {
			bam1_t *br = &(c_b->_brec);
			if(r_name.compare(bam_get_qname(br)) != 0)	continue; //name check:
			if(bam_is_supplementary(br) || bam_is_secondary(br))	continue; //primary check
			uint8_t *bam_seq =  bam_get_seq(br);
			s_store.resize(region_AT_read_right - region_AT_read_left);
			int direction_primary = bam_is_fwd_strand(br);
			int read_len = br->core.l_qseq;
			if(direction_primary == region_direction)
			{
				for(int i = region_AT_read_left; i < region_AT_read_right; i++){
				  uint8_t b_c = bam_seqi(bam_seq, i);
				  char c_c = 'N';
				  if(b_c < 10) c_c = "NACNGNNNTNN"[b_c];
				  s_store[i - region_AT_read_left] = (c_c);//store binary to acgt
				}
				//fprintf(stderr, "%s \n", s_store.c_str());
			}
			else
			{
				for(int i = region_AT_read_left; i < region_AT_read_right; i++){
				  uint8_t b_c = bam_seqi(bam_seq, read_len - i -1);
				  char c_c = 'N';
				  if(b_c < 10) c_c = "NTGNCNNNANN"[b_c];
				  s_store[i - region_AT_read_left] = (c_c);//store binary to acgt
				}
				//fprintf(stderr, "%s \n", s_store.c_str());
			}

			break;
		}
		return true;
	}

	bool is_single_read(){
		return (l_left.size() == 1 && l_right.size() == 1 && s_left == s_right);
	}

};
//load NOT only signal reads, but ALL reads in this region(not only signal reads in this region)
void SV_CALLING_Handler::LRS_assembly_load_read_ASM(bool print_log, RefRegion &main, Ass_Block & a_b, int min_LRS_read_length, BAM_handler *cur_LRS_read){//buff to store contig
	RefRegion &r = main;
	//load all reads iin this regions:
	r.show();
	//combine the LRS signals with NGS signals
	//load part of LRS reads
	int read_idx = -1;
	BAM_handler &cur = *cur_LRS_read;
	bool is_full_begin = true;
	bool is_full_end = true;
	int region_AT_read_st;
	int region_AT_read_ed;

	struct read_mapping_in_ref{
		int pos;
		int end;
	};

	std::vector<read_mapping_in_ref> rm_l;
	std::map<std::string, ASM_CMB_LIST> asm_combine_l;

	for(READ_record &rr :cur.read_list){
		std::string s;
		read_idx++;
		int ref_span_st;
		int ref_span_ed;
		cur.load_read_LRS_ASM(s, rr, r.st_pos, r.ed_pos,
				region_AT_read_st, region_AT_read_ed, ref_span_st, ref_span_ed,
				is_full_begin, is_full_end);
		if(!s.empty()){
			fprintf(stderr, "rname %s PRIMARY_MAPPING_POS:%d:%d, "
					"Current_region_AT_ref [%d-%d ,len %d]"
					"Current_region_AT_read [%d-%d ,len %d]"
					"Current_alignment_total_ref_span [%d-%d, len %d] "
					"Region_2_ref_span_edge [%d-%d] "
					"Current_region is_full_begin-end [%d-%d] "
					"Mapping Direction %d\n",
					rr.r_name.c_str(), rr.primary_chrID, rr.primary_pos,
					r.st_pos,r.ed_pos, r.ed_pos - r.st_pos,
					region_AT_read_st, region_AT_read_ed, region_AT_read_ed -region_AT_read_st,
					ref_span_st, ref_span_ed, ref_span_ed - ref_span_st,
					r.st_pos - ref_span_st, ref_span_ed - r.ed_pos,
					is_full_begin, is_full_end,
					rr.direction
					);

			//store to list: and try to combine
			//load the main strings
			std::map<std::string, ASM_CMB_LIST>::iterator it = asm_combine_l.find(rr.r_name);
			if(it != asm_combine_l.end()){//find pairs
				it->second.store_one_data(region_AT_read_st, ref_span_st, true, s, rr.direction);  //store the left:
				it->second.store_one_data(region_AT_read_ed, ref_span_ed, false, s, rr.direction); //store the right:
			}else{//store a new one:
				ASM_CMB_LIST as;
				as.store_new(rr.r_name, rr.primary_chrID, rr.primary_pos);
				as.store_one_data(region_AT_read_st, ref_span_st, true, s, rr.direction);  //store the left:
				as.store_one_data(region_AT_read_ed, ref_span_ed, false, s, rr.direction); //store the right:
				asm_combine_l[rr.r_name] = as;
			}
		}
	}
	//try to get the full cover contigs
	std::map<std::string, ASM_CMB_LIST>::iterator it = asm_combine_l.begin();
	for(;it!= asm_combine_l.end() ;it++){
		if(it->second.l_left.size() == 0 || it->second.l_right.size() == 0){
			/* When one side is zero, DO NOTHING*/
		}
		else if(it->second.is_single_read()){
		    a_b.add_read_LRS(it->second.s_left, -1, r.chr_ID, r.st_pos, 1, 1, -1, -1);
			rm_l.emplace_back();
			rm_l.back().pos = it->second.l_left[0].ref_span_pos;
			rm_l.back().end = it->second.l_right[0].ref_span_pos;
		}else{
			//try to get the MAX span:
			int region_AT_read_left; int region_AT_read_right; int ref_span_left; int ref_span_right;
			std::string s_store;
			//fprintf(stderr, "Cur Hand %s\n", it->first.c_str());
			it->second.pairing_and_get_string(&LRS_read.file, s_store, region_AT_read_left, region_AT_read_right, ref_span_left, ref_span_right);
			if(s_store.empty() )
				continue;
			a_b.add_read_LRS(s_store, read_idx, r.chr_ID, r.st_pos, 1, 1, -1, -1);
			rm_l.emplace_back();
			rm_l.back().pos = ref_span_left;
			rm_l.back().end = ref_span_right;

			fprintf(stderr, "[AFTER COMBINE] rname %s PRIMARY_MAPPING_POS:%d:%d, "
					"Current_region_AT_ref [%d-%d ,len %d]"
					"Current_region_AT_read [%d-%d ,len %d]"
					"Current_alignment_total_ref_span [%d-%d, len %d]\n",
					it->second.r_name.c_str(), it->second.primary_chrID, it->second.primary_pos,
					r.st_pos,r.ed_pos, r.ed_pos - r.st_pos,
					region_AT_read_left, region_AT_read_right, region_AT_read_right - region_AT_read_left,
					rm_l.back().pos, rm_l.back().end, rm_l.back().end - rm_l.back().pos);
		}
	}
	//remove results when more then one:
	if(a_b.ass_read_list.size() > 1){
		uint load_read_n = a_b.ass_read_list.size();
		int middle_pos = (r.st_pos + r.ed_pos) /2;
		int max_idx = -1;
		int max_dis = 0;
		for(uint i = 0 ; i < load_read_n; i++ ){
			int a = ABS_U(rm_l[i].pos, middle_pos);
			int b = ABS_U(rm_l[i].end, middle_pos);
			//int dis_clip = MIN(a, b);
			int dis_clip = a+b;
			if(max_dis < dis_clip){
				max_dis = dis_clip;
				max_idx = i;
			}else if(max_dis == dis_clip){
				fprintf(stderr, "WARNING::LRS_assembly_load_read_ASM::SAME\n");
			}
		}
		std::swap(a_b.ass_read_list[max_idx], a_b.ass_read_list[0]);
		a_b.ass_read_list.resize(1);
		std::swap(a_b.reads[max_idx], a_b.reads[0]);
		a_b.reads.resize(1);
	}
}
//load NOT only signal reads, but ALL reads in this region(not only signal reads in this region)
//1: load All SRS reads in the region
//2: load all unmapped SRS reads in this region
//3: load all TL SRS reads in this region
//4: load all LRS reads in this region

bool SV_CALLING_Handler::assembly_load_read(bool print_log,
		RefRegion r1, RefRegion r2,
		RefRegion main,
		RefRegion supp, int minQUAL, uint MAX_load_reads, bool usingLRS, bool usingNGS, int min_LRS_read_length, BAM_handler *cur_LRS_read){//buff to store contig
//am->clear();//clear
	RefRegion r = *(ref_handler->get_cur_region());
	if(!(r1.region_overlap(r) && r2.region_overlap(r)))
		return false;
	int region_st_pos = r.st_pos;
	fprintf(stderr, "in assembly_load_read_ALL1\n");
	//step1: get SV regions
	//am->reads.reserve(MAX_load_read_number*4);//malloc for reads list
	//for one SV region, search all reads from position: [st_pos - edge_len] to [ed_pos + edge_len]
	bool with_supp = true;
	if(main.region_overlap(supp)){ main.Combine(supp, true); with_supp = false; }
	//step2: load reference for contig aligner
	//step3: load reads
	//load reads:
	ass_block.clear();
	//load read signals S1: load from SR read list
	//if(with_supp) read.search_reads(Read_type::SR, ass_block_r1, supp.st_pos, supp.ed_pos, MAX_load_read_number, region_st_pos);

	if(usingLRS){
		//load reads LRS: load reads from LRS
		if(preset_platform == PRESET_ASM){
			//loading read for ASM dataset
			LRS_assembly_load_read_ASM(print_log, main, ass_block, min_LRS_read_length, cur_LRS_read);
		}else{
			LRS_assembly_load_read(print_log, main, ass_block, min_LRS_read_length, cur_LRS_read);
		}
	}

	if(usingNGS){
		for(int mode = 0; mode < 2; mode++){
			Bam_file *c_b = &SRS_read.file;
			R_region region;
			if(mode == 0)			main.toR_region(region);
			else if(with_supp)		supp.toR_region(region);
			else					break;

			resetRegion_ID(c_b, &region);	//reset region
			while (bam_next(c_b) && ass_block.reads.size() < MAX_load_reads) {
				bam1_t *br = &(c_b->_brec);
				if(bam_is_duplicate(br))				continue;
				if(bam_is_secondary(br))				continue;
				if(bam_is_supplementary(br))			continue;
				//bam_aligned_analysis(br, &clip_left, &clip_right, &gap_mismatch_inside);
				ass_block.add_read_bin_compact(Read_type::SR, bam_get_seq(br), bam_get_qual(br), br->core.l_qseq, 0, 0, 0 ,0, 0, minQUAL);
				//fprintf(stderr, "%d %d %s\n", br->core.tid, br->core.pos, ass_block_r1.reads.back().c_str());
			}
		}
		if(ass_block.reads.size() >= MAX_load_reads){
			ass_block.clear();
			return false;
		}

		int UM_st_pos = main.st_pos - 1000; UM_st_pos = MAX(UM_st_pos, 0);
		int UM_ed_pos = main.ed_pos + 1000;
		//load read signals S2: load from UM read list
		int MAX_load_read_number = 1000;
		SRS_read.search_reads(Read_type::UM, ass_block, UM_st_pos, UM_ed_pos, MAX_load_read_number, region_st_pos);

		if(print_log) fprintf(stderr, "S4\n");

		char log_title[1024]; sprintf(log_title, "[Load TR signals]");
		double __cpu_time = cputime(); double __real_time = realtime();

		//load read signals S3: load from TR read list
		TL_read_number = SRS_read.tr_loader.load_read(print_log, main.chr_ID, UM_st_pos, UM_ed_pos, ass_block, ref_handler, MAX_load_read_number);

		if(print_log) fprintf(stderr, "S5\n");

		fprintf(stderr, "\n%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);

	}

	return true;
}

//load ONLY signal reads(CLIP/DRP/MIS-MATCH)
bool SV_CALLING_Handler::SRS_assembly_load_read(bool print_log, SVE &sv, RefRegion main, RefRegion supp){//buff to store contig
	fprintf(stderr, "in assembly_load_read\n");
	//am->clear();//clear
	RefRegion r = *(ref_handler->get_cur_region());
	if(sv.isInRegion(r) == false)return false;//check
	int region_st_pos = r.st_pos;

	//step1: get SV regions
	int MAX_load_read_number = 1000;
	//am->reads.reserve(MAX_load_read_number*4);//malloc for reads list
	//for one SV region, search all reads from position: [st_pos - edge_len] to [ed_pos + edge_len]
	bool with_supp = true;
	if(main.region_overlap(supp)){ main.Combine(supp, true); with_supp = false; }
	//step2: load reference for contig aligner
	//step3: load reads
	//load reads:
	ass_block.clear();
	//load read signals S1: load from SR read list

	if(print_log) fprintf(stderr, "S1\n");

	SRS_read.search_reads(Read_type::SR, ass_block, main.st_pos, main.ed_pos, MAX_load_read_number, region_st_pos);

	if(print_log) fprintf(stderr, "S2\n");

	if(with_supp) SRS_read.search_reads(Read_type::SR, ass_block, supp.st_pos, supp.ed_pos, MAX_load_read_number, region_st_pos);

	if(print_log) fprintf(stderr, "S3\n");

	int UM_st_pos = main.st_pos - 1000; UM_st_pos = MAX(UM_st_pos, 0);
	int UM_ed_pos = main.ed_pos + 1000;
	//load read signals S2: load from UM read list
	SRS_read.search_reads(Read_type::UM, ass_block, UM_st_pos, UM_ed_pos, MAX_load_read_number, region_st_pos);

	if(print_log) fprintf(stderr, "S4\n");

	char log_title[1024]; sprintf(log_title, "[Load TR signals]");
	double __cpu_time = cputime(); double __real_time = realtime();

	//load read signals S3: load from TR read list
	TL_read_number = SRS_read.tr_loader.load_read(print_log, main.chr_ID, UM_st_pos, UM_ed_pos,ass_block, ref_handler, MAX_load_read_number);

	if(print_log) fprintf(stderr, "S5\n");

	fprintf(stderr, "\n%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);

	return true;
}

//return a list:
void SV_CALLING_Handler::SRS_getSuggestSVlength(AssemblyContig &contig){
	int sum_deletions = 0; int deletion_num = 0;
	int sum_insertions = 0; int insertion_num = 0;
	suggest_SV_length.clear();
	int suggest_list_size = contig.SUGGEST_pos_list.size();
	for(int i = 1; i < suggest_list_size; i++){
		int suggest_SV_len = contig.SUGGEST_pos_list[i - 1].suggest_pos - contig.SUGGEST_pos_list[i].suggest_pos;
		if(suggest_SV_len < 0){
			deletion_num ++; sum_deletions += suggest_SV_len;
		}else{
			insertion_num ++; sum_insertions += suggest_SV_len;
		}
		suggest_SV_length.emplace_back(suggest_SV_len);
	}
	if(deletion_num > 1)suggest_SV_length.emplace_back(sum_deletions);
	if(insertion_num > 1)suggest_SV_length.emplace_back(sum_insertions);
	int sum_total = sum_deletions + sum_insertions;
	int sum_number_total = deletion_num + insertion_num;
	if(sum_number_total > 1)suggest_SV_length.emplace_back(sum_total);
}

struct contig_repeat_info{
	contig_repeat_info(	int begin_){
		begin = begin_; repeat_times = 1;
	}
	contig_repeat_info(){
		begin = 0; repeat_times = 1;
	}
	int begin;
	int repeat_times;
};

//void ref_search_repeat_core(int wordLength, uint8_t *tseq, int tlen, std::map<std::string, contig_repeat_info> &ref_kmer_set){
//	std::string ref;
//	ref.resize(tlen);
//	for(int i = 0; i < tlen; i++)
//		ref[i] = ("ACGTN"[tseq[i]]);
//
//	ref_kmer_set.clear();
//	unsigned kmer_number = ref.size() - wordLength + 1;
//	for (unsigned j = 0; j < kmer_number; ++j) {
//		const std::string word(ref.substr(j, wordLength));
//		auto it = ref_kmer_set.find(word);
//		if(it == ref_kmer_set.end()){
//			ref_kmer_set.insert(std::make_pair(word, contig_repeat_info(j)));
//		}else{
//			it->second.repeat_times++;
//			it->second.begin = j;
//		}
//	}
//}

//#define REPEAT_WORD_LEN 31
void repeat_search_contig_core(	bool print_log, int REPEAT_WORD_LEN, std::string &contig, std::vector<char> &repeat_pos, std::map<std::string, contig_repeat_info> &kmer_set_buff){
	repeat_pos.resize(contig.size());
	for(uint32_t i = 0; i < contig.size(); i++)
		repeat_pos[i] = 0;

	int last_dis = -1;
	unsigned kmer_number = contig.size() - REPEAT_WORD_LEN + 1;
	for (unsigned j = 0; j < kmer_number; ++j) {
		const std::string word(contig.substr(j, REPEAT_WORD_LEN));//get kmer
		auto it = kmer_set_buff.find(word);
		if(it == kmer_set_buff.end()){
			contig_repeat_info cr(j);
			kmer_set_buff[word] = cr;
			last_dis = -1;
		}else{
			//fprintf(stderr, "%d %d %d %s\n", j, it->second.begin, j - it->second.begin, word.c_str());
			int set_char;
			if(last_dis == (int)(j - it->second.begin))
				set_char = 1;
			else
				set_char = 0;
			last_dis = j - it->second.begin;

			//when a kmer repeat 2 times in contig
			if(it->second.repeat_times == 1){
				for(int i = 0; i < REPEAT_WORD_LEN; i++)
					repeat_pos[i+it->second.begin] = set_char;
			}
			it->second.repeat_times++;
			it->second.begin = j;
			repeat_pos[j] = set_char;
			for(int i = 0; i < REPEAT_WORD_LEN; i++)
				repeat_pos[i+j] = set_char;
		}
	}
	if(print_log){
		for(uint32_t i = 0; i < contig.size(); i++){
			if(repeat_pos[i] == 0)
				fprintf(stderr,"-");
			else
				fprintf(stderr,"1");
		}
		fprintf(stderr,"\n");
	}
}

void repeat_search_contig(
	bool print_log,
		std::string &contig, std::vector<char> &repeat_pos, std::map<std::string, contig_repeat_info> &kmer_set_buff){
	repeat_pos.resize(contig.size());
	for(uint32_t i = 0; i < contig.size(); i++)
		repeat_pos[i] = 0;
	std::vector<char> repeat_pos_TMP;
	repeat_search_contig_core(print_log, 15, contig, repeat_pos_TMP, kmer_set_buff);
	for(uint32_t i = 0; i < contig.size(); i++)
		repeat_pos[i] |= repeat_pos_TMP[i];
	repeat_search_contig_core(print_log, 31, contig, repeat_pos_TMP, kmer_set_buff);
	for(uint32_t i = 0; i < contig.size(); i++)
		repeat_pos[i] |= repeat_pos_TMP[i];
	repeat_search_contig_core(print_log, 55, contig, repeat_pos_TMP, kmer_set_buff);
	for(uint32_t i = 0; i < contig.size(); i++)
		repeat_pos[i] |= repeat_pos_TMP[i];
}

bool repeat_search_contig_core_M1(int wordLength, std::string &c, std::vector<int> &repeat_pos, std::map<std::string, contig_repeat_info> &ref_kmer_set){
	std::map<std::string, contig_repeat_info> kmer_set;

	repeat_pos.resize(c.size());
	for(uint32_t i = 0; i < c.size(); i++){
		repeat_pos[i] = 0;
	}

	unsigned kmer_number = c.size() - wordLength + 1;
	bool find_repeat = false;
	for (unsigned j = 0; j < kmer_number; ++j) {
		const std::string word(c.substr(j, wordLength));
		auto it = kmer_set.find(word);
		if(it == kmer_set.end()){
			kmer_set.insert(std::make_pair(word, contig_repeat_info(j)));
		}else{
			//when a kmer repeat 2 times in contig
			if(it->second.repeat_times == 1){
				for(int i = 0; i < wordLength; i++){repeat_pos[i+it->second.begin] = 1;	}
			}
			it->second.repeat_times++;
			it->second.begin = j;
			for(int i = 0; i < wordLength; i++){repeat_pos[i+j] = 1;}
		}
	}
	//debug
	//if(kmer_set.size() < kmer_number){ for(auto it = kmer_set.begin(); it != kmer_set.end();it++){ if(it->second.repeat_times > 1)	fprintf(stderr,"Find repeat of wordLength %d at [%s] [BEGIN: %d REP:%d END:%d]\n", wordLength, it->first.c_str(), it->second.begin, it->second.repeat_times, it->second.end);}}	else fprintf(stderr,"No repeat found at wordLength %d\n", wordLength);

	for(uint32_t i = 0; i < c.size(); i++){
		if(repeat_pos[i] == 0){
			fprintf(stderr,"0");
		}
		else{
			fprintf(stderr,"1");
		}
	}
	fprintf(stderr," @ wordLength %d\n", wordLength);
	return find_repeat;
}


//the check whether ENOUGH reads cover the both EDGE of repeat
bool SV_CALLING_Handler::SRS_read_cover_repeat_filter(
		bool print_log,
		AssemblyContig &contig,
		uint32_t SV_size_ori,
		std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list, int ori_read_number){
	//when new SVs is found
	//search repeat for this contig
	bool withRepeat = false;
	std::vector<char> repeat_pos_INS_type;//the repeat position: used when check INSERTIONs
	std::vector<char> repeat_pos_DEL_type;//the repeat position: used when check DELETIONs
	std::map<std::string, contig_repeat_info> repeat_kmer_set_buff;
//			if(ref_kmer_set.empty()){
	//get ref string for this contig:
	//CORE:

	bool with_DEL = false;
	bool with_INS = false;
	for(uint SV_ID = SV_size_ori; SV_ID < SV_result_TMP_buff.size(); SV_ID++){
		int SV_type_int = SV_result_TMP_buff[SV_ID].get_SV_TYPE_int();
		if(1 == SV_type_int)//0 is INS and 1 is DEL
			with_DEL = true;
		else
			with_INS = true;
	}

	if(with_DEL){
		repeat_kmer_set_buff.clear();
		if(print_log) fprintf(stderr, "Search repeat in the new contig(DEL):\n");
		//get repeat MASK for DEL
		uint8_t * ref_bin = ca.get_ref() + ca.last_ref_st_pos;
		uint ref_len = ca.last_ref_end_pos - ca.last_ref_st_pos + 1;
		std::string ref_string;	ref_string.resize(ref_len);
		for(uint i = 0; i < ref_len; i++)
			ref_string[i] = "ACGT"[ref_bin[i]];
		repeat_search_contig(print_log, ref_string, repeat_pos_DEL_type, repeat_kmer_set_buff);
		//get CIGAR:
		uint32_t* bam_cigar = NULL;
		int n_cigar = ca.get_cigar(&bam_cigar);
		std::vector<char> repeat_pos_tmp;

		int seq_i = 0;
		for(int cigar_ID = 0;cigar_ID < n_cigar; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			switch(type){
			case 0:	for(int i = 0; i < cigar_len; i++, seq_i++) {repeat_pos_tmp.emplace_back(repeat_pos_DEL_type[seq_i]);}	break;//M
			case 1:	for(int i = 0; i < cigar_len; i++) {repeat_pos_tmp.emplace_back(0);}	break;//I, print 0
			case 2:	seq_i += cigar_len;	break;//D, print nothing
			case 3:	for(int i = 0; i < cigar_len; i++) {repeat_pos_tmp.emplace_back(0);}	break;//N, print 0
			case 4:	for(int i = 0; i < cigar_len; i++) {repeat_pos_tmp.emplace_back(0);}	break;//S, print 0
			default: break;//do nothing
			}
		}
		std::swap(repeat_pos_tmp, repeat_pos_DEL_type);
		if(print_log){
			for(uint32_t i = 0; i < repeat_pos_DEL_type.size(); i++){
				if(repeat_pos_DEL_type[i] == 0)
					fprintf(stderr,"-");
				else if(repeat_pos_DEL_type[i] == 1)
					fprintf(stderr,"1");
				else
					fprintf(stderr,"%c", repeat_pos_DEL_type[i]);
			}
			fprintf(stderr,"\n");
		}
	}

	if(with_INS){
		repeat_kmer_set_buff.clear();
		fprintf(stderr, "Search repeat in the new contig(INS):\n");
		//get repeat MASK for INS
		repeat_search_contig(print_log, contig.seq, repeat_pos_INS_type, repeat_kmer_set_buff);
		if(print_log){
			for(uint32_t i = 0; i < repeat_pos_INS_type.size(); i++){
				if(repeat_pos_INS_type[i] == 0)
					fprintf(stderr,"-");
				else if(repeat_pos_INS_type[i] == 1)
					fprintf(stderr,"1");
				else
					fprintf(stderr,"%c", repeat_pos_INS_type[i]);
			}
			fprintf(stderr,"\n");
		}
	}
	std::vector<int> support_depth_INS;
	std::vector<int> support_depth_DEL;

	//get support_depth
	for(int SVTYPE = 0; SVTYPE < 2; SVTYPE++){
		if(SVTYPE == 0 && with_INS == false)
			continue;
		if(SVTYPE == 1 && with_DEL == false)
			continue;

		std::vector<int> *support_depth = &support_depth_INS;
		std::vector<char> *repeat_pos = &repeat_pos_INS_type;

		if(SVTYPE == 0){//0 is INS and 1 is DEL
			fprintf(stderr, "Get DEPTH(INS):\n");

		}else{
			support_depth = &support_depth_DEL;
			repeat_pos = &repeat_pos_DEL_type;
			fprintf(stderr, "Get DEPTH(DEL):\n");
		}
		support_depth->resize(contig.seq.size());
		for(uint32_t i = 0; i < contig.seq.size(); i++){
			support_depth[0][i] = 0;
		}
		//read check
		for (auto &ca : contig.actions){
			if(ca.read_ID >= ori_read_number || ca.isAdd == false) continue;
			if(ass_read_list[ca.read_ID].signal_type == Read_type::SR && ca.wrong_base >= 1) continue;

			int contig_pos = -contig.ass_begin_offset_in_contig + ca.position_in_contig - ca.position_read;

			if(print_log){
				for(int i = 0; i < contig_pos ;i++)fprintf(stderr, " ");
				fprintf(stderr, "%s\t", read_list[ca.read_ID].c_str());
				ass_read_list[ca.read_ID].print(stderr, 0);
				fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
				ca.print(stderr, true);
			}

			int read_st_contig = contig_pos; read_st_contig = MAX(read_st_contig, 0);
			int read_ed_contig = contig_pos + read_list[ca.read_ID].size();
			read_st_contig += 3;
			read_ed_contig -= 3;
			read_ed_contig = std::min(read_ed_contig, (int)contig.seq.size() - 1);

			//check begin:
			int check_begin = 0;
			int check_end = 0;
			//search read support region
			for(check_begin = read_st_contig; check_begin < read_ed_contig; check_begin++){	if(repeat_pos[0][check_begin] == 0) break;	}
			for(check_end = read_ed_contig - 1; check_end >= read_st_contig; check_end--) {	if(repeat_pos[0][check_end] == 0)   break; }
			//add supported region depth
			for(int i = check_begin; i < check_end; i++)
				support_depth[0][i] ++;
		}
		if(print_log){
			for(uint32_t i = 0; i < support_depth->size(); i++){
				fprintf(stderr, "%c", '0' + support_depth[0][i]);
			}
			fprintf(stderr,"\n");
		}
	}

	uint32_t SV_final_size = SV_size_ori;
	while(SV_final_size < SV_result_TMP_buff.size()){
		bool SV_erase = false;
		auto &c_sv = SV_result_TMP_buff[SV_final_size];

		int SV_contig_bg = 0; int SV_contig_ed = 0; int SV_type_int = 0;
		c_sv.get_sv_ST_EN_in_contig(SV_contig_bg, SV_contig_ed, SV_type_int);
		int unsupport_base_number = 0;
		if(SV_type_int == 0){ //INS
			SV_contig_bg = MAX(SV_contig_bg - 10, 0);
			SV_contig_ed = MIN(SV_contig_ed + 10, (int)support_depth_INS.size());
			for(int i = SV_contig_bg; i < SV_contig_ed; i++){
				if(support_depth_INS[i] < 2)//at least 2 read support
					unsupport_base_number ++;
			}
			if(unsupport_base_number > 20)
				SV_erase = true;

		}else if(SV_type_int == 1){ //DEL
			SV_contig_bg = MAX(SV_contig_bg - 10, 0);
			SV_contig_ed = MIN(SV_contig_ed + 10, (int)support_depth_DEL.size());
			for(int i = SV_contig_bg; i < SV_contig_ed; i++){
				if(support_depth_DEL[i] < 2)
					unsupport_base_number ++;
			}
			if(unsupport_base_number > 0)
				SV_erase = true;
		}
		fprintf(stderr, "Unsupport_base_number %d [%d %d]\t", unsupport_base_number,SV_contig_bg -10,SV_contig_ed+10);

		if(SV_erase){
			withRepeat = true;
			fprintf(stderr, "SV was erased by repeat filter\n");
			if(print_log){
				SV_result_TMP_buff[SV_final_size].print(stderr, SRS_read.file._hdr, ref_handler);
			}
			SV_result_TMP_buff.erase(SV_result_TMP_buff.begin() + SV_final_size);
		}else{
			fprintf(stderr, "SV was KEPT by repeat filter\n");
			if(print_log){
				SV_result_TMP_buff[SV_final_size].print(stderr, SRS_read.file._hdr, ref_handler);
			}
			SV_final_size ++;
		}
	}
	return withRepeat;
}

void SV_CALLING_Handler::SRS_force_calling_MEI_ins(bool print_log, SVE &sv, std::vector<AssemblyContig> &contigs, int main_st_pos){//buff to store contig
	fprintf(stderr, "Trigger the force calling ALU process\n");
	int contig_ID = -1;
	//int ref_is_alu = -1;//-1: unknown; 0: false; 1: true
	for (AssemblyContig &contig : contigs) {
		contig_ID++;
		if(contig.contig_is_filtered) continue;
		if(contig.ending_reason[0] != 0 ||  contig.ending_reason[1] != 0) continue;
		int alu_kmer_idx;
		int search_pos = alu_handler->kmer_check(contig.seq, alu_kmer_idx);
		if(alu_kmer_idx != -1){
			fprintf(stderr, "contig_ID %d , contig.seq.c_str() %s alu_kmer_idx %d, search_pos %d, ", contig_ID, contig.seq.c_str(), alu_kmer_idx, search_pos);
			int ALU_DIR;//0:AA; 1: TT
			int suggest_BP_in_contig = alu_handler->get_suggest_BP_in_contig(alu_kmer_idx, search_pos, ALU_DIR);
			fprintf(stderr, "ALU_DIR %d, suggest_BP_in_contig %d\n", ALU_DIR, suggest_BP_in_contig);

			if(!contig.SUGGEST_pos_list.empty()){
				int contig_suggest_pos_in_ref = 0;
				if(ALU_DIR == 1) contig_suggest_pos_in_ref = contig.SUGGEST_pos_list.back().suggest_pos;   //TT
				else	 		 contig_suggest_pos_in_ref = contig.SUGGEST_pos_list.begin()->suggest_pos; //AA
				int suggest_ALU_position_in_ref = contig_suggest_pos_in_ref + suggest_BP_in_contig;

				//generate demo SVs and cigar:
				std::vector<char> ref_buff; ref_buff.emplace_back('N'); ref_buff.emplace_back(0);
				std::vector<char> alt_buff;
				std::vector<uint8_t> contig_modified;
				int ALU_len = 320 + ((suggest_ALU_position_in_ref%15)-7);//set as an random number
				uint32_t CIGAR_v[3];
				int contig_pos_in_local_region;

				//generate ALT string:
				{
					int ALU_beg_in_contig = (ALU_DIR == 1)?(suggest_BP_in_contig - ALU_len):suggest_BP_in_contig;
					int ALU_end_in_contig = ALU_beg_in_contig + ALU_len;
					//S1
					{//generate the CONTIG for ALU
						//generate the CONTIG for ALU: head part
						if(ALU_beg_in_contig < 50){
							int N_str_pre = 50 - ALU_beg_in_contig;
							for(int i = 0; i < N_str_pre;i++)
								contig_modified.emplace_back(4);
							ALU_beg_in_contig += N_str_pre;
							ALU_end_in_contig += N_str_pre;
							contig_suggest_pos_in_ref -= N_str_pre;
						}
						//generate the CONTIG for ALU: middle part
						const char * char_contig = contig.seq.c_str();
						for(uint32_t i = 0; i < contig.seq.size(); i++){
							uint8_t bin_c;
							switch (char_contig[i]) {
								case 'A': bin_c = 0; break;
								case 'C': bin_c = 1; break;
								case 'G': bin_c = 2; break;
								case 'T': bin_c = 3; break;
								default:  bin_c = 4; break;
							}
							contig_modified.emplace_back(bin_c);
						}
						////generate the CONTIG for ALU: tail part
						int N_str_suf = ALU_end_in_contig + 50 - (int)contig_modified.size();
						if(N_str_suf > 0){
							for(int i = 0; i < N_str_suf;i++)
								contig_modified.emplace_back(4);
						}
					}
					//S2.1
					{//contig_in_local_ref: the contig position in the local refe
						if(ALU_DIR == 1)	contig_pos_in_local_region = contig_suggest_pos_in_ref- main_st_pos + ALU_len;
						else				contig_pos_in_local_region = contig_suggest_pos_in_ref- main_st_pos;
					}
					//S2
					{//ADD cigar:
						CIGAR_v[0] = (ALU_beg_in_contig << BAM_CIGAR_SHIFT) + 0;//Match
						CIGAR_v[1] = (ALU_len << BAM_CIGAR_SHIFT) + 1;//Insertion
						CIGAR_v[2] = ((contig_modified.size() - ALU_end_in_contig) << BAM_CIGAR_SHIFT) + 0;//Match
					}
					//S3
					{//generate ALT string:
						for(int i = ALU_beg_in_contig; i < ALU_end_in_contig; i++){
							alt_buff.emplace_back("ACGTN"[contig_modified[i]]);
						}
						alt_buff.emplace_back(0);
					}
					if(print_log){
						//modified contig_string:
						fprintf(stderr, "ALU_beg_in_contig %d %d %s \n", ALU_beg_in_contig, ALU_end_in_contig, contig.seq.c_str());
						fprintf(stderr, "suggest_ALU_position_in_ref %d %s\n", suggest_ALU_position_in_ref, &(alt_buff[0]));
					}
					//S4
					{//filter: check whether the ALU string is SAME as reference, is true, the variant will be deleted
						//align the ALU: position
						int alu_align_suggest_ref_st_pos = suggest_ALU_position_in_ref - main_st_pos - ALU_len - 80;
						fprintf(stderr, "ALIGN_BIGIN, alu_align_suggest_ref_st_pos %d, ca.get_tlen() %d\n", alu_align_suggest_ref_st_pos, ca.get_tlen());
						//alignment
						alu_align_suggest_ref_st_pos = MAX(alu_align_suggest_ref_st_pos, 0);
						uint8_t *INS_string_bin = &(contig_modified[0]) + ALU_beg_in_contig;
						int INS_string_len = ALU_end_in_contig - ALU_beg_in_contig;
						ca.align_non_splice_default_ref(INS_string_bin, INS_string_len,	alu_align_suggest_ref_st_pos, ca.get_tlen(), 400);
						int ref_adj_size = ca.adjustCIGAR();
						alu_align_suggest_ref_st_pos += ref_adj_size;
						fprintf(stderr, "ALIGN_END, alu_align_suggest_ref_st_pos %d\n", alu_align_suggest_ref_st_pos);
						//count the 'match' base number
						int ref_same_base = ca.check_match_base(alu_align_suggest_ref_st_pos);
						//count the 'N' base number
						int N_base = 0;
						for(int i = 0;i < INS_string_len; i++ )
							if(INS_string_bin[i] > 3)
								N_base++;
						if(print_log){
							ca.print_X_E(stderr, alu_align_suggest_ref_st_pos);
							fprintf(stderr, "alu_align_suggest_ref_st_pos %d ref_same_base %d, N_base %d, INS_string_len %d\n", alu_align_suggest_ref_st_pos, ref_same_base, N_base, INS_string_len);
						}
						//get the bases uniquely belonged to ALU
						int ALU_BASE = INS_string_len - ref_same_base - N_base;
						if(ALU_BASE < 80){
							if(print_log)fprintf(stderr, "ALU_BASE %d, SKIP\n", ALU_BASE);
							continue;
						}
					}
					//S5 alu string check
					{
						uint8_t *alu_bin_st = &(contig_modified[0]) + ALU_beg_in_contig;
						int ALU_align_len = ALU_end_in_contig - ALU_beg_in_contig;
						bool pass_ins_filter = alu_handler->ALU_INS_check(alu_bin_st, ALU_align_len, ALU_DIR);
						if(!pass_ins_filter)
							continue;
					}
				}
				//position check
				{
					bool pass_alu_position_check = true;
					//INS position must be in the middle part of the SV region, not be the edge
					if(suggest_ALU_position_in_ref < sv.r1.st_pos + 100 || suggest_ALU_position_in_ref > sv.r2.ed_pos - 100)
						pass_alu_position_check = false;
					fprintf(stderr, "ALU_region position check %s [%d %d] %d \n", (pass_alu_position_check)?"PASS":"FAIL",sv.r1.st_pos, sv.r2.ed_pos, suggest_ALU_position_in_ref);
					if(!pass_alu_position_check)
						continue;
				}

				//store the final results
				NOVA_SV_FINAL_RST_item::add_to_vector(SV_result_TMP_buff, sv.r1.chr_ID, suggest_ALU_position_in_ref,
						"INS", &(ref_buff[0]), &(alt_buff[0]),	alt_buff.size() - 1,
						&(contig_modified[0]), contig_modified.size(),
						CIGAR_v, 3, 1, 1, contig_pos_in_local_region, region_ref_global_position);
				SV_result_TMP_buff.back().set_SVLEN_is_Imprecise(); //SV is Imprecise
				SV_result_TMP_buff.back().set_force_ALU();// SV is ALU(force-calling)
				//store ONLY one final result
				break;
			}
		}
	}
}

//assembly based on simple greedy-assembly
bool SV_CALLING_Handler::SRS_assembly_variations_MEI_AND_LONG_SV(SVE &sv, bool call_ALU_TL, bool &withRepeat){
	bool print_log = false;
	int edge_len = sig_para.MaxReadLen;

	RefRegion main(sv.r1.chr_ID, sv.r1.st_pos - edge_len, sv.r1.ed_pos);
	RefRegion supp(sv.r2.chr_ID, sv.r2.st_pos - edge_len, sv.r2.ed_pos);

	main.st_pos = MAX(main.st_pos, 0);
	supp.st_pos = MAX(supp.st_pos, 0);

	//set ori_reference
	int ref_region_length;
	if(supp.ed_pos - main.st_pos < 50000)
		ref_region_length = supp.ed_pos - main.st_pos;
	else
		ref_region_length = main.getLen();
	region_ref_global_position = main.st_pos;
	ca.setRef(ref_handler->getRefStr(region_ref_global_position), ref_region_length, ref_handler->get_cur_region()->chr_ID, main.st_pos);

	//show reference string
	if(print_log){ fprintf(stderr, "C_REF_STRING: [%d-%d, len: %d]\n", main.st_pos, main.st_pos + ref_region_length, ref_region_length);
	for(int i = 0; i < ref_region_length; i++)
		fprintf(stderr, "%c", "ACGT"[ ca.getTseq(i)]); }

	if(print_log) fprintf(stderr, "\nBegin assembly_load_read\n");

	if( !SRS_assembly_load_read(print_log, sv, main, supp)) return false;

	if(false){
		Bam_file *c_b = &SRS_read.file;
		R_region region;
		ref_handler->get_cur_region()->toR_region(region);
		resetRegion_ID(c_b, &region);	//reset region
	}

	if(print_log) fprintf(stderr, "\nEND assembly_load_read\n");

	if(print_log){ fprintf(stderr, "after assembly_load_read: C_REF_STRING: [%d-%d, len: %d]\n", main.st_pos, main.st_pos + ref_region_length, ref_region_length);}

	if(print_log){//show reads list
		fprintf(stderr, "Current read list:\n");
		int read_num = ass_block.ass_read_list.size();
		int first_read_ID = (read_num > 0)?ass_block.ass_read_list[0].read_list_index:0;
		for(int i = 0; i < read_num; i++){
			ass_block.ass_read_list[i].print(stderr, first_read_ID);
			std::cerr << ass_block.reads[i] << std::endl;
		}
	}

	ass_block.run_assembly(am);
	if(am->reachMaxKmerLength) withRepeat = true;
	if(ass_block.contigs.empty()){  fprintf(stderr, "No assembly results\n"); return false; }

	std::vector<std::string> &read_list = ass_block.reads;
	std::vector<ASS_reads_info> &ass_read_list = ass_block.ass_read_list;
	int ori_read_number = ass_read_list.size();
	int contig_ID = -1;
	std::vector<AssemblyContig> &contigs = ass_block.contigs;
	FILE * detail_output = stderr;
	//handle each contig
	SV_result_TMP_buff.clear();
	//if(print_log){ fprintf(stderr, "T1 assembly_load_read: C_REF_STRING: [%d-%d, len: %d]\n", main.st_pos, main.st_pos + ref_region_length, ref_region_length);}
	//call variant using single contig
	withRepeat = false;
	for (AssemblyContig &contig : contigs) {
		//fprintf(stderr, "contig.kmerLength %d \n", contig.kmerLength);
		contig_ID++;
		if(contig.ending_reason[0] != 0 ||  contig.ending_reason[1] != 0) {contig.contig_is_filtered = true; continue; withRepeat = true;}
		if(print_log) printContig(detail_output, ori_read_number, contig, contig_ID, ass_read_list);
		if(contig.seq.size() < 100)  {contig.contig_is_filtered = true; continue;}
		bool support_read_filter_fail = (contig_ID != 0 && (contig.new_support_read <= 2 && contig.kmerLength < 100));
		if(print_log && support_read_filter_fail)  fprintf(detail_output, "This CONTIG is discarded, reason: 'MIN_NEW_SUPPORT_READ'\n");
		if(support_read_filter_fail)  {contig.contig_is_filtered = true; continue;}

		//get suggestion alignment position
		SRS_get_suggention_alignment_position_list(contig, ori_read_number, read_list, ass_read_list, stderr);

		//output actions
		std::sort(contig.actions.begin(), contig.actions.end(), AssemblyReadAction::cmp_by_position);
		if(print_log){
			for (auto &ca : contig.actions){
				if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
				{
					int contig_pos = -contig.ass_begin_offset_in_contig + ca.position_in_contig;
					if(ca.isAdd){
						for(int i = 0; i < contig_pos - ca.position_read ;i++)fprintf(stderr, " ");
						fprintf(stderr, "%s\t", read_list[ca.read_ID].c_str());
					}
					ass_read_list[ca.read_ID].print(stderr, 0);
					fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
					ca.print(stderr, true);
				}
			}
		}

		if(contig.SUGGEST_pos_list.size() <= 1){
			fprintf(stderr, "No SV found in this contig, skip realignment.\n");
			continue;
		}

		if(true){
			//for each suggestion position:
			store_bin_contig(contig.seq, bin_contig);
			SRS_getSuggestSVlength(contig);
			if(print_log){
				fprintf(stderr, "Suggest SV length list for this contig:\t");
				for(int sug_len: suggest_SV_length)
					fprintf(stderr, "%d\t", sug_len);
				fprintf(stderr, "\n");
			}

			uint32_t SV_size_ori = SV_result_TMP_buff.size();
			if(!contig.SUGGEST_pos_list.empty()){
				if(print_log){ fprintf(stderr, "SUGGEST_pos_list[0].suggest_pos [%d] main.st_pos %d\n", contig.SUGGEST_pos_list[0].suggest_pos, main.st_pos);}
				SRS_alignment_and_get_var(print_log, 0, contig_ID, contig.SUGGEST_pos_list[0].suggest_pos - main.st_pos, contig.seq.size(), ref_region_length);
			}

			if(SV_result_TMP_buff.size() != SV_size_ori){
				withRepeat |= SRS_read_cover_repeat_filter(print_log, contig, SV_size_ori, ass_read_list, read_list, ori_read_number);
			}
		}
	}
	//tran-ins process:
	if(call_ALU_TL){
		fprintf(stderr, "Trigger the contig-combine process\n");
		//combined for repeat tail
		SRS_combine_repeat_tail_of_contigs(print_log, ori_read_number, ass_read_list, read_list, contigs);
		if(SV_result_TMP_buff.empty()){
			fprintf(stderr, "Trigger the tran_based_ins calling process\n");
			SRS_handle_tran_based_ins(print_log, ori_read_number, ass_read_list, read_list, contigs);
		}
		//force calling ALU process
		if(SV_result_TMP_buff.empty() && force_calling_ALU)
			SRS_force_calling_MEI_ins(print_log, sv, contigs, main.st_pos);
	}

	//output candidate SVs
	uint32_t candidateSVSize = SV_result_TMP_buff.size();
	if(candidateSVSize > 0){//select best result and remove duplication results
		NOVA_SV_FINAL_RST_item::resultFilter(SV_result_TMP_buff, MIN_sv_len);
		for(uint32_t i = 0; i < candidateSVSize; i++){
			SV_result_TMP_buff[i].writeVCF_final(&vcfBuffStr, SRS_read.file._hdr, NULL);
			fprintf(stderr, "%s: %s", (SV_result_TMP_buff[i].will_be_output_to_vcf)?"With_VCF":"Without VCF", vcfBuffStr.s);

			if(SV_result_TMP_buff[i].will_be_output_to_vcf){
				result_SRS_INDEL.emplace_back();
				std::swap(SV_result_TMP_buff[i], result_SRS_INDEL.back());
			}
		}
	}
	return true;
}
//
//void debug_code_(std::vector<NOVA_SV_FINAL_RST_item> & result_sv_l, const char * vcf_fn, RefHandler *ref){
//	result_sv_l.clear();
//	BCF_FILE vcf_r;//vcf for read
//	VCF_open_read(&vcf_r, vcf_fn);//open for read
//
//	char *c_sv_type = (char *)malloc(1000);
//	bcf_hdr_t *header = vcf_r.header;
//	int *support_read_number = (int*)xcalloc(3,4);
//	int32_t SV_LEN = 0;
//	RefRegion * r = ref->get_cur_region();
//	char *GT[1]; char GT_1[3]; GT[0] = GT_1;
//	do{//read one
//		bcf1_t *c_r = &( vcf_r.r);
//		if(r->chr_ID != c_r->rid || c_r->pos < r->st_pos || c_r->pos > r->ed_pos)
//			continue;
//		//unpack the vcf data to get the alt string
//		bcf_unpack(c_r, BCF_UN_STR);
//		vcf_get_sv_GT(vcf_r.header, c_r, GT);
//		vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
//		vcf_get_sv_LENGTH(vcf_r.header, c_r, &SV_LEN);
//		vcf_get_sv_SR(vcf_r.header, c_r, support_read_number);
//		bool is_vntr = vcf_get_sv_flag(header, c_r, "VNTR");
//		fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t\n", GT_1, is_vntr, c_r->rid , c_r->pos, SV_LEN, c_sv_type, c_r->d.allele[0], c_r->d.allele[1]);// chrID+st+len
//		NOVA_SV_FINAL_RST_item::add_to_vector(result_sv_l, c_r->rid, c_r->pos + 1, c_sv_type, c_r->d.allele[0], c_r->d.allele[1],
//				0,  0, 0, NULL, 0, 0, 0, 0, r->st_pos);
//		result_sv_l.back().isVNTRcallerRst = is_vntr;
//		if(*c_r->d.flt == 0)//filter: PASS
//			result_sv_l.back().QUAL_presetVNTR = 30;
//		else
//			result_sv_l.back().QUAL_presetVNTR = 0;
//		result_sv_l.back().set_supp_read_NGS(support_read_number[0], support_read_number[1], support_read_number[2]);
//		int GT_final = 0;
//		if(GT[0][0] == 4) GT_final++;
//		if(GT[0][1] == 4) GT_final++;
//		result_sv_l.back().setGenotype_directly(GT_final);
//		//store data into:
//	}while(VCF_next(&vcf_r));
//	//close
//	bcf_close(vcf_r.file);
//}

//assembly based on complex DBG-based assembly
bool SV_CALLING_Handler::SRS_assembly_variations_INS_DEL_VNTR(SVE &sv, bool using_LRS, bool using_SRS){//buff to store contig
	SV_result_TMP_buff.clear();
	//parameters
	bool print_log = false;
	int edge_len = sig_para.MaxReadLen;
	int refAdditionLoad = 250;
	int readAdditionLoad = 250;
	//main handlers
	dbgHandler->clear();
	//sv.r1.st_pos = sv.r2.st_pos = 54446887;
	//sv.r1.ed_pos = sv.r2.ed_pos = 54447887;
	//set core region
	RefRegion mainN(sv.r1.chr_ID, sv.r1.st_pos - edge_len, sv.r1.ed_pos);
	RefRegion suppN(sv.r2.chr_ID, sv.r2.st_pos - edge_len, sv.r2.ed_pos);

	mainN.st_pos = MAX(mainN.st_pos, 0);
	suppN.st_pos = MAX(suppN.st_pos, 0);

	//set reference
	RefRegion loadRef;
	{
		//set ori_reference
		bool regionNearby = (suppN.ed_pos - mainN.st_pos < 10000 && suppN.chr_ID == mainN.chr_ID);
		int ref_region_length = (regionNearby)?(suppN.ed_pos - mainN.st_pos):mainN.getLen();
		loadRef.chr_ID = sv.r1.chr_ID;
		loadRef.st_pos = mainN.st_pos - refAdditionLoad;
		loadRef.st_pos = MAX(loadRef.st_pos, 0);
		loadRef.ed_pos = mainN.st_pos + ref_region_length + refAdditionLoad;
		ca.setRef(ref_handler->getRefStr(loadRef.st_pos), loadRef.getLen(), loadRef.chr_ID, loadRef.st_pos);
		dbgHandler->refRegion = loadRef;
		dbgHandler->addRef(print_log, ref_handler->getRefStr(loadRef.st_pos), loadRef.getLen(), &ca);
		//RunLengthCompect(dbgHandler->reference[0]);
	}
	//load reads
	{
		char log_title[1024]; sprintf(log_title, "[VNTR load reads]");
		double __cpu_time = cputime(); double __real_time = realtime();
		RefRegion mainNLoadRead(sv.r1.chr_ID, mainN.st_pos - readAdditionLoad, mainN.ed_pos + readAdditionLoad);
		RefRegion suppNLoadRead(sv.r1.chr_ID, suppN.st_pos - readAdditionLoad, suppN.ed_pos + readAdditionLoad);
		fprintf(stderr, "Loading reads with loop!\n");
		mainNLoadRead.print(stderr);
		suppNLoadRead.print(stderr);
		assembly_load_read(print_log, sv.r1,sv.r2, mainNLoadRead, suppNLoadRead, 3, 60000,using_LRS, using_SRS,0, &LRS_read);
		fprintf(stderr, "TL_read_number %d\n", TL_read_number);
		std::swap(ass_block.reads, dbgHandler->reads_str);
		//store reads type:
		uint read_size = ass_block.ass_read_list.size();
		dbgHandler->read_info.resize(read_size);
		for(uint i = 0; i < read_size; i++){
			if(ass_block.ass_read_list[i].signal_type == Read_type::LRS)
				dbgHandler->read_info[i].is_LRS_read = true;
			else
				dbgHandler->read_info[i].is_LRS_read = false;
		}
		fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
	}

	//error correction for reads
	{
		char log_title[1024]; sprintf(log_title, "[VNTR err corr]");
		double __cpu_time = cputime(); double __real_time = realtime();
		if(false) dbgHandler->showAllRead();
		for(int KMER_len = 65; KMER_len >= 65; KMER_len-=30){
			dbgHandler->kmerCountingTable.clear();
			dbgHandler->kmerCounting(KMER_len, false, true);
			dbgHandler->errorCorrection(KMER_len, 5, 1);
		}
		fprintf(stderr, "SHOW READs debug point\n");
		if(print_log) dbgHandler->showAllRead();
		fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
	}

	//debug:130
	int minDepthNGS = 3;
	int kmerLenMax = 120;
	int kmerLenMin = 60;
	int kmerLenStep = 20;
	//60: 130- --- 60*2/15=8; 110: 60*4/15=16; 90: 60*6/15=24; 70: 60*8/15=32; 50: 60*10/15=40;
	//30: 130- --- 30*2/15=4; 110: 30*4/15=08; 90: 30*6/15=12; 70: 30*8/15=16; 50: 30*10/15=20;
	for(int KMER_len =kmerLenMin  ; KMER_len <=kmerLenMax; KMER_len+=kmerLenStep){
		//minDepth += 1;
		fprintf(stderr, "++KMER_len %d\n", KMER_len);
		dbgHandler->kmerCountingTable.clear();
		//get the counting of all KMERs
		dbgHandler->kmerCounting(KMER_len, true, true);
		//dbgH->kmerCountingTest(KMER_len);
		//using the KMERs with enough depth(>=minDepth) to build DBG
		dbgHandler->buildDBG(minDepthNGS, 0);
		if(print_log)
			dbgHandler->showDBG_graph();
		//show path for all reads
		if(print_log)
			dbgHandler->prFilter.showStringPath("REF",dbgHandler->reference[0], KMER_len, dbgHandler->kmerIdx, dbgHandler->unitig_l, print_log);
		//set path for all reads
		dbgHandler->prFilter.generateLink(dbgHandler->reads_str, KMER_len, dbgHandler->kmerIdx, dbgHandler->unitig_l, print_log);
		if(print_log)dbgHandler->prFilter.show();

		bool withCPXNode = dbgHandler->getDBGPathM1(print_log);
		//generate VCFs
		for(KeyNodePath &k :dbgHandler->pc_l_final){
			k.analysisAndGenerateSV_M2(SV_result_TMP_buff, dbgHandler->id2info, dbgHandler->reads_str, dbgHandler->abH, dbgHandler->curKmerLen,
					dbgHandler->reference[0],	dbgHandler->ca, dbgHandler->refRegion, print_log, dbgHandler->unitig_l, &(dbgHandler->gra), &SRS_read.file);
		}

		//final filters
		for(NOVA_SV_FINAL_RST_item & sv: SV_result_TMP_buff)
			sv.printSimple(stderr);
		fprintf(stderr, "Cur KmerLen is %d, show all SVs END(suggestFinalSVLength)\n", dbgHandler->curKmerLen);
		if(!withCPXNode){
			fprintf(stderr, "All node is not complex( %d), other kmer is skipped\n", dbgHandler->curKmerLen);
			break;
		}
		else{
			fprintf(stderr, "With complex path( %d), higher kmer length is used\n", dbgHandler->curKmerLen);
		}
	}

	//output candidate SVs
	bool withDEL = false;
	bool withINS = false;
	uint32_t candidateSVSize = SV_result_TMP_buff.size();
	if(candidateSVSize > 0){//select best result and remove duplication results
		//NOVA_SV_FINAL_RST_item::resultFilter(SV_result_TMP_buff, MIN_sv_len);
		for(uint32_t i = 0; i < candidateSVSize; i++){
			if(SV_result_TMP_buff[i].get_SV_TYPE_int() == 0) withINS = true;
			if(SV_result_TMP_buff[i].get_SV_TYPE_int() == 1) withDEL = true;
		}
	}

	if(candidateSVSize > 0 &&(withDEL == false || withINS == false)){//select best result and remove duplication results
		//NOVA_SV_FINAL_RST_item::resultFilter(SV_result_TMP_buff, MIN_sv_len);
		for(uint32_t i = 0; i < candidateSVSize; i++){
			//filters:
			result_SRS_INDEL.emplace_back();
			std::swap(SV_result_TMP_buff[i], result_SRS_INDEL.back());
		}
	}
	return true;
}

void SV_CALLING_Handler::LRS_contig_candidate_generate(bool print_log, int minDepthTGS, int KMER_len){//buff to store contig
	target_path_idx_l.clear();
	updated_targets_by_current_read.clear();
	special_branch_handler.clear();
	int MIN_score_final = 100;
	if(PRESET_ERR == preset_platform){
		MIN_score_final = -300;
	}

	//contig candidate generating
	for(int mode = 0; mode < 2; mode ++){
		if(mode == 1){
			special_branch_handler.special_branch_check(print_log, minDepthTGS);
			fprintf(stderr, "SP check: special_branch_handler.need_reclassify is %d minDepthTGS %d\n", special_branch_handler.need_reclassify, minDepthTGS);
			if(special_branch_handler.need_reclassify)
				target_path_idx_l.clear();
			else
				break;
		}

		for(READ_MATCH_NUM &rmn : read_path_handler.read_match_num_l){
			int read_id = rmn.read_id;
			//if(36 == read_id)		continue;
			int read_length = dbgHandler->reads_str[read_id].size();
			bool LRS_full_begin = ass_block.ass_read_list[read_id].LRS_full_begin;
			bool LRS_full_end = ass_block.ass_read_list[read_id].LRS_full_end;
			bool is_clip_read = ass_block.ass_read_list[read_id].is_clip_read();
			bool is_clip_AT_RIGHT_read = ass_block.ass_read_list[read_id].is_clip_right();
			bool is_clip_AT_LEFT_read = ass_block.ass_read_list[read_id].is_clip_left();
			if(is_clip_AT_RIGHT_read && is_clip_AT_LEFT_read)
				continue;//SKIP when both end CLIP, method cannot handler this type reads right now
			int clip_target_mode = ass_block.ass_read_list[read_id].get_clip_mode();

			fprintf(stderr, "Current handle read id: %d %d full_bg %d full_ed %d clip %d [mode %d right %d left %d] \n", read_id, read_length, LRS_full_begin, LRS_full_end, is_clip_read, clip_target_mode, is_clip_AT_RIGHT_read, is_clip_AT_LEFT_read);
			if(print_log) {
				read_path_handler.read_path_list[read_id].show_combine_read_path(read_id, read_length);
			}
			updated_targets_by_current_read.clear();
			for(uint cur_target_ID = 0; cur_target_ID < target_path_idx_l.size(); cur_target_ID++){
				TARGET_path_and_index & tpi = target_path_idx_l[cur_target_ID];
				bool is_other_snp_reject_read_target_pair = false;
				if(snp_sh.is_with_data()){
					is_other_snp_reject_read_target_pair = snp_sh.is_reject_read_target_pair(read_id, tpi.support_read_l);
					int length_diff = ABS_U(tpi.target_total_len + KMER_len, read_length);
					if(print_log) fprintf(stderr, "length_diff %d is_other_snp_reject_read_target_pair %d\n", length_diff, is_other_snp_reject_read_target_pair);
					if(length_diff < 15)
						is_other_snp_reject_read_target_pair = false;
				}
				if(print_log){
					fprintf(stderr, "Current used target id: %d:\n", cur_target_ID);
					tpi.show_the_target_path();
				}
				if(is_other_snp_reject_read_target_pair == true){
					if(print_log)
						fprintf(stderr, "The read is not used because reject by SNP checker\n");
				}
				else if(tpi.clip_target_mode != clip_target_mode){
					//tpi。
					if(print_log)
						fprintf(stderr, "The read is not used because different CLIP condition %d %d \n", tpi.clip_target_mode, clip_target_mode);
				}
				else{
					tpi.SDP(print_log, read_path_handler.read_path_list[read_id].combine_path, LRS_full_begin, LRS_full_end, read_length);//already check:
					if(print_log) fprintf(stderr, "[SDP log:] [r_id:%d;t_id:%d]The final max score is %d @\n", read_id, cur_target_ID, tpi.max_score_final);
					if(print_log) tpi.showSDP_result(&(read_path_handler.read_path_list[read_id].combine_path[0]), read_id, cur_target_ID);
					bool is_absorb_by_target =
							tpi.SDP_update_and_branch_condition_check(print_log, MIN_score_final, read_path_handler.read_path_list[read_id].combine_path,
									read_id, LRS_full_begin, LRS_full_end, is_clip_read,
									dbgHandler->unitig_l, special_branch_handler.long_branch_check_list, special_branch_handler.special_branch_check_list);
					if(is_absorb_by_target){
						updated_targets_by_current_read[cur_target_ID] = tpi.max_score_final;
					}
				}
			}
			if(updated_targets_by_current_read.empty()){
				fprintf(stderr, "The read[%d] is treat as a new branch [target ID %ld]\n", read_id, target_path_idx_l.size());
				target_path_idx_l.emplace_back();
				target_path_idx_l.back().build(read_path_handler.read_path_list[read_id].combine_path, read_id, dbgHandler->unitig_l,
						target_path_idx_l.size() - 1, ass_block.ass_read_list[read_id].get_clip_mode());
			}else{
				if(print_log || true){
					fprintf(stderr, "The read[%d] is used by following branch: ", read_id);
					std::map<int, int>::iterator it = updated_targets_by_current_read.begin();
					for(; it != updated_targets_by_current_read.end(); it++)
						fprintf(stderr, "target: %d:score %d\t", it->first, it->second);
					fprintf(stderr, "\n");
				}
				if(updated_targets_by_current_read.size() > 1){
					//select the max score node:
					int max_score = -1;
					int max_score_target = -1;
					std::map<int, int>::iterator it = updated_targets_by_current_read.begin();
					for(; it != updated_targets_by_current_read.end(); it++){
						if(max_score < it->second){
							max_score = it->second;
							max_score_target = it->first;
						}
					}
					fprintf(stderr, "The read[%d] has max_score_target %d\n ", read_id, max_score_target);
					it = updated_targets_by_current_read.begin();
					for(; it != updated_targets_by_current_read.end(); it++){
						int target_id = it->first;
						target_path_idx_l[target_id].support_read_l.back().share_target_num = updated_targets_by_current_read.size();
						if(target_id != max_score_target)
							target_path_idx_l[target_id].support_read_l.back().is_best_target = false;
					}
				}
			}
		}
		fprintf(stderr, "END of mode %d\n", mode);
	}
}

bool SV_CALLING_Handler::LRS_assembly_variations(bool print_log, BAM_handler *cur_LRS_read,
		LRS_SIG &sig, RefRegion &loadRef, std::string &full_ref_str, std::vector<uint8_t> &bin_ref){//buff to store contig
	print_log = false;
	RefRegion r;
	r.chr_ID = sig.chrID;
	r.st_pos = sig.POS;
	r.ed_pos = sig.END;
	SV_result_TMP_buff.clear();
	//parameters
	int edge_len = sig_para.MaxReadLen;
	int refAdditionLoad = 500;
	int readAdditionLoad = 500;
	//main handlers
	dbgHandler->clear();
	//sv.r1.st_pos = sv.r2.st_pos = 54446887;
	//sv.r1.ed_pos = sv.r2.ed_pos = 54447887;
	//set core region
	RefRegion mainN(r.chr_ID, r.st_pos - edge_len, r.ed_pos);
	RefRegion suppN(r.chr_ID, r.st_pos - edge_len, r.ed_pos);

	mainN.st_pos = MAX(mainN.st_pos, 0);
	suppN.st_pos = MAX(suppN.st_pos, 0);

	fprintf(stderr, "\n\n[BEGIN]SV region handler begin:\t");
	//10775964-10776605
//	if(mainN.st_pos == 10775964){
//		fprintf(stderr, "\n");
//	}
	mainN.print(stderr);
	fprintf(stderr, "\n");

	bool regionNearby = (suppN.ed_pos - mainN.st_pos < 10000 && suppN.chr_ID == mainN.chr_ID);
	int ref_region_length = (regionNearby)?(suppN.ed_pos - mainN.st_pos):mainN.getLen();

	bool run_length_all_read = false;

	//for bin reference
	std::vector<uint8_t> RunLength_count_ref;
	std::vector<uint8_t> cmp_bin_ref;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S1:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load the reference~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	{
		//set ori_reference
		loadRef.chr_ID = r.chr_ID;
		loadRef.st_pos = mainN.st_pos - refAdditionLoad;
		loadRef.st_pos = MAX(loadRef.st_pos, 0);
		loadRef.ed_pos = mainN.st_pos + ref_region_length + refAdditionLoad;

		//load full reference:
		dbgHandler->refRegion = loadRef;
		dbgHandler->addRef(print_log, ref_handler->getRefStr(loadRef.st_pos), loadRef.getLen(), NULL);
		full_ref_str = dbgHandler->reference[0];//copy:
		store_bin_contig(full_ref_str, bin_ref);
		if(run_length_all_read){
			RunLengthCompect_BIN(ref_handler->getRefStr(loadRef.st_pos), loadRef.getLen(), cmp_bin_ref);
			ca.setRef(&(cmp_bin_ref[0]), cmp_bin_ref.size(), loadRef.chr_ID, loadRef.st_pos);
		}
		else{
			ca.setRef(ref_handler->getRefStr(loadRef.st_pos), loadRef.getLen(), loadRef.chr_ID, loadRef.st_pos);
		}

		ca.setZdrop(5000, 5000);

		if(run_length_all_read){
			if(print_log){
				fprintf(stderr, "REF ORI:%s\n", dbgHandler->reference[0].c_str());
			}
			RunLengthCompect(dbgHandler->reference[0], RunLength_count_ref);
			//show:
			if(print_log){
				if(dbgHandler->reference[0].size() != 0){
					fprintf(stderr, "REF CMP:%s\nREF CNT:", dbgHandler->reference[0].c_str());
					for(uint8_t i: RunLength_count_ref)
						fprintf(stderr, "%c", i + '0');
					fprintf(stderr, "\n");
				}
			}
		}else{
			if(print_log){
				fprintf(stderr, "REF ORI:%s\n", dbgHandler->reference[0].c_str());
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S2:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<std::vector<uint8_t>> RunLength_count_read;
	{
		char log_title[1024]; sprintf(log_title, "[VNTR load reads]");
		double __cpu_time = cputime(); double __real_time = realtime();
		RefRegion mainNLoadRead(r.chr_ID, mainN.st_pos - readAdditionLoad, mainN.ed_pos + readAdditionLoad);
		RefRegion suppNLoadRead(r.chr_ID, suppN.st_pos - readAdditionLoad, suppN.ed_pos + readAdditionLoad);
		fprintf(stderr, "Loading reads with loop!\t"); mainNLoadRead.print(stderr);	suppNLoadRead.print(stderr);
		int min_LRS_read_length = readAdditionLoad*2 - 200;//1800
		//int min_LRS_read_length = 1500;
		assembly_load_read(print_log, r,r, mainNLoadRead, suppNLoadRead, 3, 60000, true, false,  min_LRS_read_length, cur_LRS_read);
 		fprintf(stderr, "TL_read_number %d\n", TL_read_number);
		std::swap(ass_block.reads, dbgHandler->reads_str);
		//store reads type:
		if(dbgHandler->reads_str.size() > 100){
			//the depth is too high, skip
			fprintf(stderr, "Skip this SV because the read depth is too high![read: %ld]\n", ass_block.ass_read_list.size());
			return false;
		}
		uint read_size = ass_block.ass_read_list.size();
		dbgHandler->read_info.resize(read_size);
		for(uint i = 0; i < read_size; i++){
			if(ass_block.ass_read_list[i].signal_type == Read_type::LRS)
				dbgHandler->read_info[i].is_LRS_read = true;
			else
				dbgHandler->read_info[i].is_LRS_read = false;
		}
		fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
		RunLengthCompect_All(run_length_all_read, print_log, ass_block, dbgHandler->reads_str, RunLength_count_read);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S3:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~KMER counting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//debug:130
	int minDepth_SRS = 3;
	int minDepth_LRS = 0;
	if(preset_platform == PRESET_ONT_Q20){
		minDepth_LRS = dbgHandler->reads_str.size() / 6;
		if( dbgHandler->reads_str.size() < 10)
			minDepth_LRS = (dbgHandler->reads_str.size()+2) / 3;
	}else if(preset_platform == PRESET_ERR){
		minDepth_LRS = (dbgHandler->reads_str.size())/4;
	}else{
		minDepth_LRS = dbgHandler->reads_str.size() / 6;
	}
	int KMER_len = 20;
	if(preset_platform == PRESET_ERR)
		KMER_len = 10;
	//reads kmer counting
	{
		//minDepth += 1;
		fprintf(stderr, "++KMER_len %d\n", KMER_len);
		dbgHandler->kmerCountingTable.clear();
		//get the counting of all KMERs
		dbgHandler->kmerCounting(KMER_len, true, false);
		//dbgH->kmerCountingTest(KMER_len);
		//using the KMERs with enough depth(>=minDepth) to build DBG
		dbgHandler->buildDBG(minDepth_SRS, minDepth_LRS);

		if(print_log)
			dbgHandler->showDBG_graph();
		if(dbgHandler->unitig_l.empty()){
			fprintf(stderr, "Skip this SV because the unitig_l is empty\n");
			return false;
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S4:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~get the path of reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	{
		read_path_handler.clear();
		read_path_handler.get_the_read_path(print_log, dbgHandler->prFilter, KMER_len, dbgHandler->kmerIdx, dbgHandler->unitig_l,
				dbgHandler->reads_str);
		//int number_full_cover = 0;
		for(uint read_id = 0; read_id < read_path_handler.read_match_num_l.size() ;read_id++){
			ASS_reads_info &ai = ass_block.ass_read_list[read_id];
			bool is_full_cover = (ai.LRS_full_begin && ai.LRS_full_end) ;
			read_path_handler.read_match_num_l[read_id].is_full_cover = is_full_cover;
			read_path_handler.read_match_num_l[read_id].is_clip_read = ai.is_clip_read();
			read_path_handler.read_match_num_l[read_id].is_clip_left = ai.is_clip_left();
			read_path_handler.read_match_num_l[read_id].is_clip_right = ai.is_clip_right();
		}
		//search, SDP all reads to target
		if(read_path_handler.read_path_list.empty()){
			fprintf(stderr, "Skip this SV because the read_path_list is empty\n");
			return false;
		}
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S5:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate the contig candidate using LRS reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//set, init and clear
	if(snp_sh.is_with_data()){
		std::vector<int> used_read_list;
		for(uint i = 0; i < ass_block.ass_read_list.size(); i++){
			used_read_list.emplace_back(ass_block.ass_read_list[i].read_list_index);
		}
		if(false) snp_sh.SNP_read_h.show_all_data();
		snp_sh.analysis(used_read_list);
	}
	region_ID_global ++;
	read_path_handler.sort_read_path_by_match_number();
	int min_depth_new_branch = MAX(minDepth_LRS, 2);
	LRS_contig_candidate_generate(print_log, min_depth_new_branch, KMER_len);
	read_path_handler.restore_sort_read_path_by_id();
	for(uint contig_ID = 0; contig_ID < target_path_idx_l.size(); contig_ID++){
		target_path_idx_l[contig_ID].get_the_contig_string(print_log, dbgHandler->unitig_l, dbgHandler->reads_str);			//get the contig strings
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S6:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~contig selecting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	contig_selecter.contig_basic_count(print_log, target_path_idx_l, read_path_handler.read_match_num_l);
	if(print_log || true){
		fprintf(stderr, "Before contig select: \n\n");
		for(auto & s :contig_selecter.select_list){
			s.show(true);
			//target_path_idx_l[s.target_ID].show_read_list();
			//fprintf(stderr, "%s\n\n", target_path_idx_l[s.target_ID].contig.c_str());
		}
	}

	{
		//checking if the contigs need combination??
		int MIN_CLIP_READ_N = 3;
		clip_contig_combine_Handler.clear();
		fprintf(stderr, "Before contig select: \n\n");
		for(auto & s :contig_selecter.select_list){
			if(s.clip_right_read_number >= MIN_CLIP_READ_N) clip_contig_combine_Handler.add_right_contig_n();
			if(s.clip_left_read_number >= MIN_CLIP_READ_N) clip_contig_combine_Handler.add_left_contig_n();
			if(s.support_read_number >= MIN_CLIP_READ_N) clip_contig_combine_Handler.add_total_contig_n();
		}

		//contig combination:
		if(clip_contig_combine_Handler.is_contig_need_combine()){
			fprintf(stderr, "Clip contigs combination begin:\n");
			//contig combination:
			//collect the lest
			for(auto & s :contig_selecter.select_list){
				if(s.clip_right_read_number >= MIN_CLIP_READ_N)  	clip_contig_combine_Handler.add_CLIP_AT_right_contig(target_path_idx_l[s.target_ID].contig, s.target_ID);
				else if(s.clip_left_read_number >= MIN_CLIP_READ_N) clip_contig_combine_Handler.add_CLIP_AT_left_contig(target_path_idx_l[s.target_ID].contig, s.target_ID);
				else if(s.full_cover_read_n >= MIN_CLIP_READ_N) 	clip_contig_combine_Handler.add_full_contig(target_path_idx_l[s.target_ID].contig, s.target_ID);
			}
			clip_contig_combine_Handler.combine();
			//replace the original target sequence:
			for(auto & cc: clip_contig_combine_Handler.combine_contigs){
				//set the original contig to useless:
				//,	cc.target_ID_part2;
				fprintf(stderr, "[%d %d] %s\n",	cc.target_ID_part1,	cc.target_ID_part2,	cc.combine_s.c_str());
				std::vector<SUPPORT_READ_LIST_ITEM> &p1_l = target_path_idx_l[cc.target_ID_part1].support_read_l;
				std::vector<SUPPORT_READ_LIST_ITEM> &p2_l = target_path_idx_l[cc.target_ID_part2].support_read_l;
				std::vector<SUPPORT_READ_LIST_ITEM> sr_c_l;
				sr_c_l = p1_l;
				sr_c_l.insert(sr_c_l.end(), p2_l.begin(), p2_l.end());
				//p2_l.clear();
				target_path_idx_l[cc.target_ID_part1].contig.clear();
				target_path_idx_l[cc.target_ID_part2].contig.clear();
				//store the combined contig as new target
				target_path_idx_l.emplace_back();
				std::swap(target_path_idx_l.back().contig, cc.combine_s);
				std::swap(target_path_idx_l.back().support_read_l, sr_c_l);
				target_path_idx_l.back().target_id = target_path_idx_l.size() - 1;
			}
			for(uint contig_ID = 0; contig_ID < target_path_idx_l.size(); contig_ID++){
				if(target_path_idx_l[contig_ID].contig.empty())
					target_path_idx_l[contig_ID].support_read_l.clear();
			}
			//recount all the target:
			contig_selecter.contig_basic_count(print_log, target_path_idx_l, read_path_handler.read_match_num_l);
			if(print_log || true){
				fprintf(stderr, "COMBINE: Before contig select: \n\n");
				for(auto & s :contig_selecter.select_list){
					s.show(false);
					//target_path_idx_l[s.target_ID].show_read_list();
					fprintf(stderr, "%s\n", target_path_idx_l[s.target_ID].contig.c_str());
				}
			}
		}
	}

	//select the contigs
	if(contig_selecter.select_list.size() == 0 || false == contig_selecter.select_contig(print_log, preset_platform)){
		fprintf(stderr, "Skip this SV because the no contig generated\n");
		return false;
	}
	//show the contigs
	fprintf(stderr, "\nAfter contig select: \n\n");
	for(uint haplotype_ID = 0; haplotype_ID < contig_selecter.select_list.size(); haplotype_ID++){
		int cur_target_ID = contig_selecter.select_list[haplotype_ID].target_ID;
		//target_path_idx_l[cur_target_ID].get_the_contig_string(print_log, dbgHandler->unitig_l, dbgHandler->reads_str);			//get the contig strings
		contig_selecter.select_list[haplotype_ID].show(true);
		if(print_log){
			target_path_idx_l[cur_target_ID].show_the_final_result(print_log, cur_target_ID, dbgHandler->unitig_l);			//show all the results:
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S6:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~candidate polishing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if(SRS_BAM_F != NULL && false){
		for(int haplotype_ID = 0; haplotype_ID < contig_selecter.use_contig_num; haplotype_ID++){
			int cur_target_ID = contig_selecter.select_list[haplotype_ID].target_ID;
			fprintf(stderr, "[candidate polishing]: current haplotype_ID %d contig ID %d\n", haplotype_ID, cur_target_ID);
			//contig polishing using NGS reads
			R_region NGS_load_region;
			int ref_region_length = (regionNearby)?(suppN.ed_pos - mainN.st_pos):mainN.getLen();
			int NGS_additional_load = 300;
			NGS_load_region.chr_ID = mainN.chr_ID;
			NGS_load_region.st_pos = MAX(mainN.st_pos - NGS_additional_load, 0);
			NGS_load_region.ed_pos = mainN.st_pos + ref_region_length + NGS_additional_load;
			//todo::
			store_bin_contig(target_path_idx_l[cur_target_ID].contig, bin_contig);
			contig_polsihing_handler.run(bin_contig, NGS_load_region);
		}
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~S7:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SV generating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fprintf(stderr, "Successfully call contig.\n");
	return true;
}

void getGT_string(int final_genotype, std::string & GT_STR){
	GT_STR.clear();
	switch (final_genotype) {
		case 0:	GT_STR.append("0/0"); break;
		case 1:	GT_STR.append("0/1"); break;
		case 2:	GT_STR.append("1/1"); break;
		case 3:	GT_STR.append("0|0"); break;
		case 4:	GT_STR.append("0|1"); break;
		case 5:	GT_STR.append("1|0"); break;
		case 6:	GT_STR.append("1|1"); break;
		default: GT_STR.append("0/0"); break;
	}
}

void SV_CALLING_Handler::LRS_SV_generating_Germline(bool print_log, LRS_SIG &sig, RefRegion &loadRef, std::string &full_ref_str, std::vector<uint8_t> &bin_ref,
		 std::vector<NOVA_SV_FINAL_RST_item> & SV_result_final,
		 std::vector<TARGET_path_and_index> &target_path_idx_l, CONTIG_SELECTER &contig_selecter){
	for(int haplotype_ID = 0; haplotype_ID < contig_selecter.use_contig_num; haplotype_ID++){
		int cur_target_ID = contig_selecter.select_list[haplotype_ID].target_ID;
		int supp_read_n, not_supp_read_n;
		contig_selecter.get_haplotype_support_read_n(haplotype_ID, supp_read_n, not_supp_read_n);
		if(supp_read_n == 0)
			continue;
		//store the bin contig
		std::string & strContig = target_path_idx_l[cur_target_ID].contig;
		store_bin_contig(strContig, bin_contig);

		//loading the reference for contig-realignment
		RefRegion r_p1;
		RefRegion r_p2;

		std::string full_ref_string;
		std::string re_alignment_ref_string;
		std::vector<uint8_t> re_alignment_bin_ref_string;

		fprintf(stderr, "Handle hap %d \n\n", haplotype_ID);
		//try to set the reference string and reference region:
		bool is_full_cover_contig = contig_selecter.select_list[haplotype_ID].is_full_cover_contig(3);
		bool is_long_range_contig = contig_selecter.select_list[haplotype_ID].is_long_range_contig() && (!is_full_cover_contig);
		bool is_left_clip_contig =  contig_selecter.select_list[haplotype_ID].is_left_clip_contig(2);
		bool is_right_clip_contig =  contig_selecter.select_list[haplotype_ID].is_right_clip_contig(2);
		//???
		if((is_left_clip_contig || is_right_clip_contig || is_long_range_contig) &&
				contig_selecter.select_list[haplotype_ID].is_single_support_contig() && bin_contig.size() > 10000){
			continue;
		}
		contig_align_region_Handler.load_alignment_reference(
				is_long_range_contig, is_left_clip_contig, is_right_clip_contig,
				loadRef, bin_ref, full_ref_str,
				bin_contig,
				r_p1, r_p2, re_alignment_ref_string, re_alignment_bin_ref_string, full_ref_string);

		//special handle of combine multy-region
		//print_log = true;
		bool is_replacement = (sig.type == LRS_SIG_X);
		LRS_contig_realignment(
				print_log, &ca,//background
				strContig, bin_contig,//input: contig
				r_p1, r_p2, re_alignment_ref_string, re_alignment_bin_ref_string, full_ref_string, //input: reference
				supp_read_n, not_supp_read_n, 0, haplotype_ID, contig_selecter.use_contig_num,//info
				SV_result_TMP_buff, is_replacement);//output
		//store result to the final SV list
		for(uint i = 0; i < SV_result_TMP_buff.size(); i++){
			SV_result_final.emplace_back();
			std::swap( SV_result_final.back(), SV_result_TMP_buff[i]);
		}
	}
}

void SV_CALLING_Handler::LRS_SV_generating_Somatic(bool print_log, LRS_SIG &sig,
		RefRegion &loadRef, std::string &full_ref_str,
		std::vector<uint8_t> &bin_ref,
		std::vector<NOVA_SV_FINAL_RST_item> &SV_result_final,
		std::vector<TARGET_path_and_index> &target_path_idx_l_Tumer,
		CONTIG_SELECTER &contig_selecter_Tumer,
		std::vector<TARGET_path_and_index> &target_path_idx_l_Normal,
		CONTIG_SELECTER &contig_selecter_Normal) {
	int MIN_READ_N_TUMOR = 3;
	int MIN_READ_N_NORMAL = 2;
	std::vector<uint8_t> bin_tc; std::vector<uint8_t> bin_nc;
	//std::vector<uint8_t> bin_tc_cmp; std::vector<uint8_t> bin_nc_cmp;
	//SELECT: THE SAME CONTIG:
	//for ALL contigs in Tumor BAM, check if it is one of contig in Normal
	//(1) whether the length is similar ? if not the SAME( += 15%+ ),they are not the same;
	//(2) whether the sequence is similar ? using KWS2 to detect whether is the SAME

	if(true){
		fprintf(stderr, "Begin Somatic SV calling\n");
		fprintf(stderr, "ALL tumor contig used:\n\n");
		//show:
		for (uint tumor_haplotype_ID = 0; tumor_haplotype_ID < contig_selecter_Tumer.select_list.size(); tumor_haplotype_ID++) {
			//skip;
			int supp_read_n, not_supp_read_n;
			contig_selecter_Tumer.get_haplotype_support_read_n(tumor_haplotype_ID, supp_read_n, not_supp_read_n);
			if (supp_read_n == 0) continue;
			CONTIG_SELECT &tc = contig_selecter_Tumer.select_list[tumor_haplotype_ID];
			if (tc.support_read_number < MIN_READ_N_TUMOR)	continue;
			tc.show(true);

		}
		fprintf(stderr, "ALL normal contig used:\n\n");
		for (uint normal_haplotype_ID = 0; normal_haplotype_ID < contig_selecter_Normal.select_list.size(); normal_haplotype_ID++) {
			CONTIG_SELECT &nc = contig_selecter_Normal.select_list[normal_haplotype_ID];
			if (nc.support_read_number < MIN_READ_N_NORMAL)	continue;
			nc.show(true);
		}
		fprintf(stderr, "END SHOW\n");
	}

	for (uint tumor_haplotype_ID = 0; tumor_haplotype_ID < contig_selecter_Tumer.select_list.size(); tumor_haplotype_ID++) {
		//skip;
		int supp_read_n, not_supp_read_n;
		contig_selecter_Tumer.get_haplotype_support_read_n(tumor_haplotype_ID, supp_read_n, not_supp_read_n);
		//
		if (supp_read_n == 0) continue;
		CONTIG_SELECT &tc = contig_selecter_Tumer.select_list[tumor_haplotype_ID];
		if (tc.support_read_number < MIN_READ_N_TUMOR)	continue;
		fprintf(stderr, "Find same contig in normal for tumor contig:\n");
		tc.show(true);
		bool tumer_is_full_cover = tc.is_full_cover_contig(MIN_READ_N_TUMOR);
		bool tumer_is_long_range_contig = tc.is_long_range_contig();
		if(tumer_is_full_cover == false && tumer_is_long_range_contig == false)
			continue;

		bool is_same_contig_with_normal = false;
		std::string & full_contig_tc = target_path_idx_l_Tumer[tc.target_ID].contig;
		store_bin_contig(full_contig_tc, bin_tc);
		//RunLengthCompect_BIN(&(bin_tc[0]), bin_tc.size(), bin_tc_cmp);

		//find the similar contig of tumor in the normal contigs
		if(tumer_is_full_cover){
			//find the similar contig of tumor in the normal contigs
			for (uint normal_haplotype_ID = 0; normal_haplotype_ID < contig_selecter_Normal.select_list.size(); normal_haplotype_ID++) {
				CONTIG_SELECT &nc = contig_selecter_Normal.select_list[normal_haplotype_ID];
				if (nc.support_read_number < MIN_READ_N_NORMAL)continue;
				if(false == nc.is_full_cover_contig(MIN_READ_N_NORMAL))	continue;//skip when the normal contig is not full covered
				//size similarL
				int tumor_len = tc.contig_len;
				int tumor_cmp_len = tc.contig_cmp_len;
				//size similarL
				int normal_len = nc.contig_len;
				int normal_cmp_len = nc.contig_cmp_len;
				//length similarity check....
				int max_diff_contig_base_len = MIN(tumor_len, normal_len) *0.15;
				max_diff_contig_base_len = MAX(max_diff_contig_base_len, 50);
				if(ABS_U(tumor_len, normal_len) >= max_diff_contig_base_len || ABS_U(tumor_cmp_len, normal_cmp_len) >= max_diff_contig_base_len){
					fprintf(stderr, "SKIP: LENGTH DIFF %d %d \n", tumor_len, normal_len);
					continue;
				}
				store_bin_contig(target_path_idx_l_Normal[nc.target_ID].contig, bin_nc);
				//RunLengthCompect_BIN(&(bin_nc[0]), bin_nc.size(), bin_nc_cmp);
				//xassert(tumor_cmp_len == (int)bin_tc_cmp.size(), "");
				//xassert(normal_cmp_len == (int)bin_nc_cmp.size(), "");
				int tumor_aln_len = tumor_len;
				int normal_aln_len = normal_len;

				int zdrop = MAX(tumor_aln_len, normal_aln_len);
				int bandwidth = MAX(tumor_aln_len, normal_aln_len) *0.3;
				somatic_same_ca.setZdrop(zdrop, bandwidth);
				somatic_same_ca.align_non_splice_normal_size(bin_nc, bin_tc);
				int tc_mapped_len, nc_mapped_len, match_base, gap_INS, gap_DEL;
				somatic_same_ca.analysis_cigar_M2(nc_mapped_len, tc_mapped_len, match_base, 50, gap_INS, gap_DEL);
				fprintf(stderr, "T_N_Contig_alignment_summary: MATCH %d gap_INS %d gap_DEL %d \n",  match_base, gap_INS, gap_DEL);
				//check: wrong alignment:
				if(tc_mapped_len != tumor_aln_len || nc_mapped_len != normal_aln_len ){
					fprintf(stderr, "SKIP: Wrong ALN\n");
					continue;
				}
				if((float)match_base*1.1 < MIN(tumor_aln_len, normal_aln_len)){
					fprintf(stderr, "SKIP: match_base too low, CIGAR is: \n");
					somatic_same_ca.printCIGAR(stderr);
					continue;
				}
				fprintf(stderr, "IS_same_contig_with_normal\n");
				is_same_contig_with_normal = true;
				break;
				//seq sim check......
				//store the bin contig
			}
		}else if(tumer_is_long_range_contig){//the tumor contig is not full covered
			for (uint normal_haplotype_ID = 0; normal_haplotype_ID < contig_selecter_Normal.select_list.size(); normal_haplotype_ID++) {
				CONTIG_SELECT &nc = contig_selecter_Normal.select_list[normal_haplotype_ID];
				if (nc.support_read_number < MIN_READ_N_NORMAL)continue;
				//if(false == nc.is_long_range_contig())	continue;//skip when the normal contig is full covered
				if((nc.is_left_clip_contig(MIN_READ_N_NORMAL) && tc.is_left_clip_contig(MIN_READ_N_TUMOR)) ||
						(nc.is_right_clip_contig(MIN_READ_N_NORMAL) && tc.is_right_clip_contig(MIN_READ_N_TUMOR))){
					is_same_contig_with_normal = true;
					break;
				}
			}
		}else{
			xassert(0 , "DATA NEED SKIP\n");
		}
		int tumer_flag = (is_same_contig_with_normal == true)?SOMATIC_TUMOR_NORMAL_SHARE:SOMATIC_TUMOR_UNIQ;
		if(false == is_same_contig_with_normal)	fprintf(stderr, "LOG: find tumor uniq contig\n");
		else									fprintf(stderr, "LOG: find tumor share contig\n");

		//set to unknown when normal contig is no data:
		if(contig_selecter_Normal.select_list.empty()) tumer_flag = SOMATIC_UNKNOWN;

		if(SOMATIC_TUMOR_UNIQ != tumer_flag){
			continue;
		}

		//loading the reference for contig-realignment
		RefRegion r_p1; RefRegion r_p2;
		std::string full_ref_string;
		std::string re_alignment_ref_string;
		std::vector<uint8_t> re_alignment_bin_ref_string;
		fprintf(stderr, "Handle hap %d \n\n", tumor_haplotype_ID);
		//try to set the reference string and reference region:
		bool is_full_cover_contig = tc.is_full_cover_contig(MIN_READ_N_TUMOR);
		bool is_long_range_contig = tc.is_long_range_contig() && (!is_full_cover_contig);
		bool is_left_clip_contig = tc.is_left_clip_contig(2);
		bool is_right_clip_contig = tc.is_right_clip_contig(2);
		if ((is_left_clip_contig || is_right_clip_contig || is_long_range_contig) &&
				tc.is_single_support_contig() && bin_contig.size() > 10000)
			continue;
		contig_align_region_Handler.load_alignment_reference( is_long_range_contig, is_left_clip_contig, is_right_clip_contig, loadRef, bin_ref, full_ref_str, bin_tc, r_p1, r_p2, re_alignment_ref_string,	re_alignment_bin_ref_string, full_ref_string);
		bool is_replacement = (sig.type == LRS_SIG_X);
		LRS_contig_realignment(print_log,
				&ca, //background, using default aligner
				full_contig_tc, bin_tc, //input: contig
				r_p1, r_p2, re_alignment_ref_string,
				re_alignment_bin_ref_string,
				full_ref_string, //input: reference
				supp_read_n, not_supp_read_n, 0, tumor_haplotype_ID,
				contig_selecter_Tumer.use_contig_num, //info
				SV_result_TMP_buff, is_replacement); //output
		//store result to the final SV list
		for (uint i = 0; i < SV_result_TMP_buff.size(); i++) {
			SV_result_final.emplace_back();
			std::swap(SV_result_final.back(), SV_result_TMP_buff[i]);
			//todo:: store the tumor/normal flags
			SV_result_final.back().tumor_flag.set(tumer_flag);
		}
	}
}
bool SV_CALLING_Handler::LRS_contig_realignment(bool print_log, Contig_String_aligner *ca_p,//background
		std::string & full_contig, std::vector<uint8_t> &bin_contig,//input, contig
		RefRegion & ref_r1, RefRegion & ref_r2, //
		std::string & re_alignment_ref_str, std::vector<uint8_t> & re_alignment_bin_ref_string, //input, reference
		std::string & full_reference_string,
		int supp_read_n, int not_supp_read_n, int unknown_read_n, int haplotype_ID, int total_haplotype_number,  //info
		std::vector<NOVA_SV_FINAL_RST_item> & SV_output_buff, bool is_replacement//output
		){//buff to store contig
	SV_output_buff.clear();//clear buff:
	uint32_t rl = re_alignment_ref_str.size();
	if(rl == 0)
		return false;
	uint32_t cl = bin_contig.size();
	ca_p->chr_ID = ref_r1.chr_ID;
	ca_p->global_ref_pos = ref_r1.st_pos;

	bool is_aln_super_length = ((ABS_U(rl, cl) > 1.5*(MIN(rl, cl))));
	int addition_load = ABS_U(rl, cl) * 1.2;
	if(rl > 10000 && cl > 10000){
		is_aln_super_length = true;
		int addition_load2 = MAX(rl, cl) * 0.2;
		addition_load = MAX(addition_load, addition_load2);
	}
	bool right_aln = true;
	if(is_aln_super_length) right_aln = ca_p->align_non_splice_super_long(bin_contig,  re_alignment_bin_ref_string, addition_load);  /*LONG*/
	else								ca_p->align_non_splice_normal_size(bin_contig, re_alignment_bin_ref_string); 				/*middle-size*/

	if(!right_aln){
		fprintf(stderr, "KSW2 align failed, return false!\n");
		return false;
	}

	if(print_log){
		fprintf(stderr, "CIGAR before adjust:\n");
		ca_p->printCIGAR(stderr);
		ca_p->show_score();
	}
	//check alignment:
	{
		int contig_mapped_len, ref_mapped_len, match_base;
		ca_p->analysis_cigar(contig_mapped_len, ref_mapped_len, match_base);
		fprintf(stderr, "ALN: contig_len %d, ref_len %d\n", contig_mapped_len, ref_mapped_len);
		if(is_aln_super_length == false && (contig_mapped_len + 500 < (int)bin_contig.size() || ref_mapped_len + 500 < (int)re_alignment_bin_ref_string.size())){
			fprintf(stderr, "Basic_ALN failed, reaigned using long-align mode\n");
			right_aln = ca_p->align_non_splice_super_long(bin_contig,  re_alignment_bin_ref_string, addition_load);  /*LONG*/
			if(!right_aln){
				fprintf(stderr, "KSW2 align failed, return false!\n");
				return false;
			}
		}
	}

	int contig_Pos_Bgn_LOCAL_region = 0;
	contig_Pos_Bgn_LOCAL_region += ca_p->adjustCIGAR();

	if(print_log) ca_p->printf_alignment_detail(stderr, contig_Pos_Bgn_LOCAL_region, NULL, 0);
	ca_p->printCIGAR(stderr);

	if(!ref_r1.same_as(ref_r2)){
		fprintf(stderr, "Long region cigar adjusting begin: \n");
		//cigar_adjusting
		uint32_t* bam_cigar;
		int cigar_len = ca_p->get_cigar(&bam_cigar);
		int ref_idx = contig_Pos_Bgn_LOCAL_region;

		int middle_position = ref_r1.getLen();
		int additional_deletion_size = ref_r2.ed_pos - ref_r1.st_pos - re_alignment_bin_ref_string.size();
		bool successfully_modification = false;

		for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
			int cigar_len =	bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
			//try test:
			if(ref_idx < middle_position && ref_idx + cigar_len > middle_position){
				if(2 == type){//deletion
					int cigar_len_new = cigar_len + additional_deletion_size;
					bam_cigar[cigar_ID] = (cigar_len_new << BAM_CIGAR_SHIFT) + type;
					successfully_modification = true;
					break;
				}
			}
			switch(type){
			case 0: case 2: case 3: case 4: case 7: case 8: /*M,D,N,S,X,=*/ ref_idx += cigar_len; break;
			case 1:	/*DO nothing*/;	break;//I, print nothing
			default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len);
			}
		}
		if(false == successfully_modification){
			fprintf(stderr, "Long region cigar adjusting [FAIL]: \n");
			return false;
		}
		else{
			fprintf(stderr, "Long region cigar adjusting [SUCCESSFUL]: \n");
			ca_p->printCIGAR(stderr);
		}
	}

	//store the results
	//store SVs as result vector, then convert to SVE_SVs
	std::vector<SV_REGION_LRS> sv_r;
	//truSV_num is the number of SVs that not shorter than MIN_sv_len
	int truSV_num = ca_p->get_canditate_SVs_TGS(false, MIN_sv_len, contig_Pos_Bgn_LOCAL_region, sv_r, output_small_var,
			full_reference_string, full_contig);
	if(truSV_num != 0 || is_replacement){
		//generate SVs
		int SV_in_HAP_ID = -1;
		//S1: ref position
		uint true_ref_position = 0;
		int refPosBgn = contig_Pos_Bgn_LOCAL_region;
//		if(run_length_all_read){ for(int i = 0; i < refPosBgn; i++) true_ref_position += RunLength_count_ref[i];}
//		else{ true_ref_position += refPosBgn; }

		true_ref_position += refPosBgn;

		for(SV_REGION_LRS & r: sv_r){
			if(print_log) r.show();
			SV_in_HAP_ID++;
			//recover from the running-length compact:
			//for deletion:
			std::string ref_str;
			std::string alt_str;
			if(r.ref_begin < 1 || r.ref_begin - 1 >= (int)full_reference_string.size())
				continue;

			alt_str += full_reference_string[r.ref_begin - 1];
			ref_str += full_reference_string[r.ref_begin - 1];

			uint32_t* cigar;
			int cigar_num = ca_p->get_cigar(&cigar);

			if(r.ref_begin != r.ref_end){//deletion
				//run_length_recover_ref(r.ref_begin, r.ref_end, RunLength_count_ref, reference_str, ref_str);
				ref_str += full_reference_string.substr(r.ref_begin, r.ref_end - r.ref_begin);
				//xassert(r.ref_end <= full_ref_str.size(),"");

				if(!ref_str.empty() && !alt_str.empty()){
					if(ref_str.size() - alt_str.size() < 10){
						NOVA_SV_FINAL_RST_item::add_to_vector(SV_output_buff, ca_p->chr_ID, r.ref_begin + ca_p->global_ref_pos,
								"INDEL", &(ref_str[0]), &(alt_str[0]),
								alt_str.size() - ref_str.size(), NULL, 0, NULL, 0, 0, 0, 0, 0);
					}
					else{
						NOVA_SV_FINAL_RST_item::add_to_vector(SV_output_buff, ca_p->chr_ID, r.ref_begin + ca_p->global_ref_pos,
								"DEL", &(ref_str[0]), &(alt_str[0]),
								alt_str.size() - ref_str.size(), &(bin_contig[0]), bin_contig.size(),
								cigar, cigar_num, r.cigar_idx, r.cigar_idx, contig_Pos_Bgn_LOCAL_region, ref_r1.st_pos + true_ref_position);
					}
					SV_output_buff.back().set_LRS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, total_haplotype_number, SV_in_HAP_ID, sv_r.size(), r.contig_begin);
				}
			}
			else if(r.contig_begin != r.contig_end){//insertion
				alt_str += full_contig.substr(r.contig_begin, r.contig_end - r.contig_begin);
				if(!ref_str.empty() && !alt_str.empty()){
					//run_length_recover_contig(print_log, r.contig_begin, r.contig_end, ch->support_read_l, RunLength_count_read, reads, alt_str);
					if(alt_str.size() - ref_str.size() < 10){
						NOVA_SV_FINAL_RST_item::add_to_vector(SV_output_buff, ca_p->chr_ID, r.ref_begin + ca_p->global_ref_pos, "INDEL", &(ref_str[0]), &(alt_str[0]),
								alt_str.size() - ref_str.size(), NULL, 0, NULL, 0, 0, 0, 0, 0);
					}
					else{
						NOVA_SV_FINAL_RST_item::add_to_vector(SV_output_buff, ca_p->chr_ID, r.ref_begin + ca_p->global_ref_pos, "INS", &(ref_str[0]), &(alt_str[0]),
								alt_str.size() - ref_str.size(), &(bin_contig[0]), bin_contig.size(),
								cigar, cigar_num, r.cigar_idx, r.cigar_idx, contig_Pos_Bgn_LOCAL_region, ref_r1.st_pos + true_ref_position);
					}
					SV_output_buff.back().set_LRS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID,
							total_haplotype_number, SV_in_HAP_ID, sv_r.size(), r.contig_begin);
				}
			}
			//S3: insertion string
			else{//SNPs
				bool skip_var = false;
				if(r.contig_begin < 0 || r.contig_begin >= (int)full_contig.size())
					skip_var = true;
				if(r.ref_begin < 0 || r.ref_begin >= (int)full_reference_string.size())
					skip_var = true;
				if(skip_var == false){
					alt_str.clear(); alt_str += full_contig[r.contig_begin];
					ref_str.clear(); ref_str += full_reference_string[r.ref_begin];
					if(!ref_str.empty() && !alt_str.empty()){
						//run_length_recover_contig(print_log, r.contig_begin, r.contig_end, ch->support_read_l, RunLength_count_read, reads, alt_str);
						NOVA_SV_FINAL_RST_item::add_to_vector(SV_output_buff, ca_p->chr_ID, r.ref_begin + ca_p->global_ref_pos + 1, "SNP", &(ref_str[0]), &(alt_str[0]),
								alt_str.size() - ref_str.size(), NULL, 0, NULL, 0, 0, 0, 0, 0);
						SV_output_buff.back().set_LRS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, total_haplotype_number, SV_in_HAP_ID, sv_r.size(),r.contig_begin);
					}
				}
			}
		}
	}
	//remove the replacement-SV support only by one reads
	if(supp_read_n <= 1){
		//CSV check
		int SV_number = 0;
		for(uint i = 0; i < SV_output_buff.size(); i++){
			int SV_len = SV_output_buff[i].SV_length;
			if(SV_len >= MIN_sv_len || -SV_len >= MIN_sv_len ){
				SV_number ++;
			}
		}
		if(SV_number > 5)
			SV_output_buff.clear();
		else if(SV_number > 1){
			std::vector<int> SV_output_IDX;
			for(uint i = 0; i < SV_output_buff.size(); i++){
				int SV_len = SV_output_buff[i].SV_length;
				if(SV_len >= MIN_sv_len || -SV_len >= MIN_sv_len ){
					SV_output_IDX.emplace_back(i);
				}
			}

			for(uint i = 1; i < SV_output_IDX.size(); i++){
				NOVA_SV_FINAL_RST_item & pre_s = SV_output_buff[SV_output_IDX[i - 1]];
				NOVA_SV_FINAL_RST_item & cur_s = SV_output_buff[SV_output_IDX[i]];
				int l1 = pre_s.SV_length; int st1 = pre_s.st_pos; int ed1 = st1 + pre_s.ref.size();
				int l2 = cur_s.SV_length; int st2 = cur_s.st_pos; int ed2 = st2 + cur_s.ref.size();
				if(l1 * l2 < 0 && (ABS_U(st1, st2) < 10 || ABS_U(ed1, ed2) < 10) && ABS(l1 + l2) != 0){
					SV_output_buff.clear();
					break;
				}
			}
		}
	}
	if(supp_read_n == 1 && SV_output_buff.size() >= 1){
		for(NOVA_SV_FINAL_RST_item & n :SV_output_buff){
			if(n.alt.empty())
				continue;
			if(n.SV_length < MIN_sv_len && n.SV_length > -MIN_sv_len )
				continue;

			std::string & s = n.alt;
	        extern uint64_t kmerMask[33];
	        int CONTIG_IDX_KMER_LEN = 10;//4^13=2^26=64M
	        uint64_t CONTIG_IDX_KMER_MASK = kmerMask[CONTIG_IDX_KMER_LEN];
	        int kmer_number = s.size() - CONTIG_IDX_KMER_LEN + 1;

			uint64_t MASK = CONTIG_IDX_KMER_MASK;
			//long AAA strings add, additional add kmer < 0.8
			std::vector<uint8_t> bin_s;
			store_bin_contig(s, bin_s);

			uint8_t * buff_bin = &(bin_s[0]);
			std::map<int32_t, int>::iterator it;
			std::map<int32_t, int> contig_kmer_index;
			std::map<int32_t, int> contig_kmer_index_tmp;

			uint64_t kmer = bit2_nextKmer_init(buff_bin, CONTIG_IDX_KMER_LEN);
			contig_kmer_index.clear();

			for(int i = 0; i < kmer_number; i++){
				kmer = bit2_nextKmerMASK( buff_bin + i, kmer, CONTIG_IDX_KMER_LEN);
				it = contig_kmer_index.find(kmer);
				if(it!=contig_kmer_index.end())
					it->second ++;
				else
					contig_kmer_index[kmer] = 1;
			}
			fprintf(stderr, "  ");

			int total_simple = 0;
			for(it = contig_kmer_index.begin(); it != contig_kmer_index.end(); it++){
				if(it->second >= 6){
					contig_kmer_index_tmp[it->first] = it->second;
					total_simple += it->second;
				}
			}
			std::swap(contig_kmer_index_tmp, contig_kmer_index); contig_kmer_index_tmp.clear();
			float rate = 0;
			if(kmer_number-CONTIG_IDX_KMER_LEN > 0){
				rate = (float)total_simple/(kmer_number-CONTIG_IDX_KMER_LEN);
				fprintf(stderr, "total_simple %d, kmer_number %d, rate %f\n", total_simple, kmer_number, rate);
			}

			if(rate > 0.4){
				//
				int contig_bg = n.LRS_INFO.SV_hap_position;
				int contig_ed = n.LRS_INFO.SV_hap_position + s.size();

				std::string s_before = full_contig.substr(MAX(0,contig_bg - 100), 100);
				std::string s_after = full_contig.substr(contig_ed, 100);
				std::string around = s_before + s_after;
				fprintf(stderr, "contig_bg %d contig_ed  %d \n full_contig %s\n, s_after %s\n, s_before %s \n s%s\n",
						contig_bg, contig_ed,
						full_contig.c_str(), s_after.c_str(), s_before.c_str(), s.c_str());
				std::vector<uint8_t> bin_around;
				store_bin_contig(around, bin_around);
				uint8_t * buff_bin = &(bin_around[0]);
				int kmer_number_around = bin_around.size() - CONTIG_IDX_KMER_LEN + 1;
				uint64_t kmer = bit2_nextKmer_init(buff_bin, CONTIG_IDX_KMER_LEN);
				bool kmer_found_around = false;
				for(int i = 0; i < kmer_number_around; i++){
					kmer = bit2_nextKmerMASK( buff_bin + i, kmer, CONTIG_IDX_KMER_LEN);
					it = contig_kmer_index.find(kmer);
					if(it!=contig_kmer_index.end()){
						kmer_found_around = true;
						break;
					}
				}
				fprintf(stderr, "kmer_found_around %d\n", kmer_found_around);
				if(kmer_found_around == false){
					SV_output_buff.clear();
					break;
				}
			}
		}
	}

	return true;
}

void ref_load_char_2_bin(bool cur_ref_is_forward, int load_len, std::vector<uint8_t> & to_bin_string, char* cur_load_ref){
	to_bin_string.resize(load_len);
	uint8_t *bin_ref = &(to_bin_string[0]);

	if(cur_ref_is_forward){
		for(int i = 0; i < load_len; i++){
			switch(cur_load_ref[i]){
			case 'A': case 'a': bin_ref[i] = 0; break;
			case 'C': case 'c': bin_ref[i] = 1; break;
			case 'G': case 'g': bin_ref[i] = 2; break;
			case 'T': case 't': bin_ref[i] = 3; break;
			case 'n': case 'N': bin_ref[i] = 4; break;
			}
		}
	}else{
		for(int i = 0; i < load_len; i++){
			switch( cur_load_ref[load_len - i - 1]){
			case 'A': case 'a': bin_ref[i] = 3; break;
			case 'C': case 'c': bin_ref[i] = 2; break;
			case 'G': case 'g': bin_ref[i] = 1; break;
			case 'T': case 't': bin_ref[i] = 0; break;
			case 'n': case 'N': bin_ref[i] = 4; break;
			}
		}
	}
}

