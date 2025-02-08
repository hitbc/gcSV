/*
 * ReadHandler.cpp
 *
 *  Created on: 2020年5月4日
 *      Author: fenghe
 */


#include "../SVcalling_core/ReadHandler.hpp"

#define MIN_COMPACT_LEN 8
bool BAM_handler::pass_compact_filter(uint8_t *s, const int len) {
	int compact_str_len = 1; //AAAACCCT --> ACT
	for (int i = 1; i < len; i++)
		if (s[i - 1] != s[i])
			compact_str_len++;
	if (compact_str_len < MIN_COMPACT_LEN)
		return false;
	return true;
}

//when a clip string clip to an "AAAAAAAA..." tail, return false
bool BAM_handler::clip_AAA_TailFilter(uint8_t *read_bin, const int read_len, int soft_left, int soft_right) {
	int middle_len = read_len - soft_left - soft_right;
	if(middle_len < 20) return false;
	//left check:
	if(soft_left > 0){
		int base_ACGT_left[4] = { 0 };
		uint8_t * st_base = read_bin + soft_left;
		for (int i = 0; i < 20; i++)
			base_ACGT_left[st_base[i]]++;
		if(base_ACGT_left[0] > 18 || base_ACGT_left[3] > 18)
			return false;
	}
	if(soft_right > 0){
		//right check:
		int base_ACGT_right[4] = { 0 };
		uint8_t * st_base = read_bin + read_len - soft_right - 20;
		for (int i = 0; i < 20; i++)
			base_ACGT_right[st_base[i]]++;
		if(base_ACGT_right[0] > 18 || base_ACGT_right[3] > 18)
			return false;
	}
	return true;
}

float get_ave_quality_value(int st_pos, int end_pos, bam1_t* _bp)
{
	uint32_t total_qual = 0;
	uint8_t* bam_quality = bam_get_qual(_bp);
	end_pos = MIN(end_pos, _bp->core.l_qseq);
	if(st_pos > end_pos){ return -100; }
	for(int i = st_pos; i <= end_pos; i++)
		total_qual += bam_quality[i];
	return (float)total_qual / (end_pos - st_pos + 1);
}

////when a clip string has low quality, return false
bool BAM_handler::clip_low_quality_Filter(bam1_t *br, int read_len, int soft_left, int soft_right) {
	float mapped_part_ave_quality_value = get_ave_quality_value(soft_left, read_len - soft_right, br);
	float clip_part_ave_quality_value_left;
	float clip_part_ave_quality_value_right;
	if(soft_left > 0){
		clip_part_ave_quality_value_left = get_ave_quality_value(0, soft_left, br);
		if(clip_part_ave_quality_value_left < 0.6 * mapped_part_ave_quality_value)
			soft_left = 0;
	}
	if(soft_right > 0){
		clip_part_ave_quality_value_right = get_ave_quality_value(read_len - soft_right - 1, read_len - 1, br);
		if(clip_part_ave_quality_value_right < 0.6 * mapped_part_ave_quality_value)
			soft_right = 0;
	}
	if(soft_left > 0 || soft_right > 0)
		return true;
	return false;
}

bool BAM_handler::storeReadSR(bam1_t *br, int soft_left, int soft_right, int gap_mismatch_inside){
	const int read_len = br->core.l_qseq;
	if(get_bam_seq_bin(0, read_len, storeReadBuff, br) == false) return false;//get string failed, return
	if (!pass_compact_filter(storeReadBuff, read_len))//
		return false;
	if (!clip_AAA_TailFilter(storeReadBuff, read_len, soft_left, soft_right))
		return false;
	if(soft_left > 0 || soft_right > 0){
		if (!clip_low_quality_Filter(br, read_len, soft_left, soft_right))
			return false;
	}
	storeReadCore(Read_type::SR, br->core.l_qseq, storeReadBuff, bam_get_qual(br), br->core.n_cigar,
			bam_get_cigar(br), br->core.flag, br->core.pos, soft_left, soft_right, gap_mismatch_inside);
	return true;
}

void BAM_handler::storeReadUM(bam1_t *br, uint8_t *query){
	storeReadCore(Read_type::UM, br->core.l_qseq, query, bam_get_qual(br), 0, 0, br->core.flag, br->core.pos, 0, 0, 0);
}

static void cigar2Path(const uint32_t* bam_cigar, const unsigned n_cigar, path_segment* path){
	for (unsigned int i = 0; i < n_cigar; ++i){
		path[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		path[i].type   = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
	}
}

//core function to store NGS reads
void BAM_handler::storeCCS_ReadCore(int readLen, uint8_t *seq, uint8_t * qual, int cigarLen, uint32_t* bam_cigar, int32_t pos){//check
	kv_resize(uint8_t, base, base.n + readLen);
	kv_resize(uint8_t, quality, quality.n + readLen);
	kv_resize(path_segment, cigar, cigar.n + cigarLen);

	//store string
	memcpy(base.a + base.n, seq, readLen);
	if(qual != nullptr)
		memcpy(quality.a + quality.n, qual, readLen);
	else
		memset(quality.a + quality.n, 0, readLen);
	//store data::
	READ_LIST *rl =&read_list;
	cigar2Path(bam_cigar, cigarLen, cigar.a + cigar.n);
	// the direction of the read, and will be direction of mate when reads are unmapped
	rl->emplace_back(0, 0, pos, readLen, cigarLen, base.n, quality.n, cigar.n, 0, 0, 0);//new a read record
	//show data:
	if(false){
		//seq
		fprintf(stderr, "CCS_read@POS:%d readLen %d, cigarLen %d\n", pos, readLen, cigarLen);
		for(int i = 0; i < readLen; i++)
			fprintf(stderr, "%c", "ACGTNNNN"[base.a[base.n + i]]);
		fprintf(stderr, "\n");
		for(int i = 0; i < readLen; i++)
			fprintf(stderr, "%c", '!' + quality.a[quality.n + i]);
		fprintf(stderr, "\n");
		for(int i = 0; i < cigarLen; i++){
			fprintf(stderr, "%d%c",cigar.a[cigar.n + i].length, pathSTR[cigar.a[cigar.n + i].type]);
		}
		fprintf(stderr, "\n");
	}
	//reset Ns; <must be after the [store data] >
	base.n += readLen;
	quality.n += readLen;
	cigar.n += cigarLen;
}

//core to store NGS reads
void BAM_handler::storeReadCore(//check
		Read_type::T t, //SR or UM
		int readLen, uint8_t *seq, uint8_t * qual,
		int cigarLen, uint32_t* bam_cigar,
		uint16_t flag, int32_t pos,
		int soft_left, int soft_right, int gap_mismatch_inside){
	if(t == Read_type::DR) return;
	if(t == Read_type::UM) cigarLen = 0; //not store cigar for UM reads
	//resize
	kv_resize(uint8_t, base, base.n + readLen);
	kv_resize(uint8_t, quality, quality.n + readLen);
	kv_resize(path_segment, cigar, cigar.n + cigarLen);
	//store read seq for SR signals
	bool stored_data_reverse = false;
	uint8_t mate_direction = ((flag & BAM_MATE_STRAND) == 0);
	if(t == Read_type::UM && (((flag & BAM_MATE_STRAND) == 0) == FORWARD))
		stored_data_reverse = true;
	//store string
	memcpy(base.a + base.n, seq, readLen);
	if(qual != nullptr)
		memcpy(quality.a + quality.n, qual, readLen);
	else
		memset(quality.a + quality.n, 0, readLen);

	if(stored_data_reverse){
		 getReverseStr(base.a + base.n, readLen);
		 if(qual != nullptr)
			 getReverseStr_qual((uint8_t*)quality.a + quality.n, readLen);
	}
	//store data::
	if(t == Read_type::SR){
		READ_LIST *rl =&read_list;
		cigar2Path(bam_cigar, cigarLen, cigar.a + cigar.n);
		// the direction of the read, and will be direction of mate when reads are unmapped
		uint8_t direction = ((flag & BAM_STRAND) == 0);
		rl->emplace_back(flag, direction, pos, readLen, cigarLen, base.n, quality.n, cigar.n, soft_left, soft_right, gap_mismatch_inside);//new a read record
	}
	else{
		READ_LIST *rl =&um_read_list;
		rl->emplace_back(flag, 1 - mate_direction, pos, readLen, -1, base.n, quality.n, cigar.n, -1, -1, -1);//new a read record
	}

	//reset Ns; <must be after the [store data] >
	base.n += readLen;
	quality.n += readLen;
	cigar.n += cigarLen;
}

bool BAM_handler::isSameHeader(const BAM_handler &B) const {
	bam_hdr_t* h1 = file._hdr;
	bam_hdr_t* h2 = B.file._hdr;
	if (h1->n_targets != h2->n_targets)								return false;
	for (int32_t i = (0); i < h1->n_targets; ++i) {
		if (h1->target_len[i] != h2->target_len[i])					return false;
		if (0 != strcmp(h1->target_name[i], h2->target_name[i]))	return false;
	}
	return true;
}
