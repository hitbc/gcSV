#include "bam_file.h"
#include "htslib/sam.h"
#include "kvec.h"
#include <stdlib.h>
#include <string.h>

void parse_bam_region(char* region, char* chrom, int32_t* begin_pos, int32_t* end_pos)
{
	char* chrom_;
	int32_t begin_pos_ = 0;
	int32_t end_pos_ = 0;

	// make first split:
	char region_cp[512];
	strcpy(region_cp, region);
	char* afterChrom = strchr(region, ':');
    if (NULL != afterChrom)
    {
    	chrom_ = region;
    	afterChrom[0] = '\0';
    	xassert(afterChrom[1] != '\0', "");
    	afterChrom++;
    }

    bool isWholeChrom = (NULL == afterChrom);
    if (!isWholeChrom)
    {
    	// make second split
    	char *tokens = NULL;
    	tokens = strtok(afterChrom, "-") - 1;
    	begin_pos_ = strtoul(tokens,NULL, 10);
    	tokens = strtok(NULL, "\0");
    	end_pos_ = strtoul(tokens,NULL, 10);
    	//"Can't parse begin and end positions from bam_region"
    }
    if(begin_pos_ < 0 || begin_pos_ > end_pos_)
    	isWholeChrom = true;
    xassert(chrom_ != NULL, "Can't parse contig name from bam_region ");
	chrom     = chrom_;
    if (isWholeChrom)
    {
    	*begin_pos = 0;
    	*end_pos   = MAX_uint32_t;
    }
    else
    {
    	*begin_pos = begin_pos_;
    	*end_pos   = end_pos_;
    }
}

void R_region2region_string(
	bam_hdr_t* header,
	char* region_str,
	R_region *region)
{
	sprintf(region_str, "%s:%d-%d",header->target_name[region->chr_ID], region->st_pos, region->ed_pos);
}

char* strTailModify(char* to, char *pre, char tail_char){
	strcpy(to, pre);
	to[strlen(to) - 1] = tail_char;
	return to;
}

void region_string2R_region(
    bam_hdr_t* header,
	char* region_str,
	R_region *region)
{
  xassert(NULL != header, "");
  xassert(NULL != region_str, "");
  char * chrom = NULL;
  parse_bam_region(region_str, chrom, &(region->st_pos), &(region->ed_pos));

  region->chr_ID = bam_name2id(header, chrom);
  xassert_3(region->chr_ID >= 0, "Contig [%s] from bam_region [%s] not found in BAM/CRAM header",chrom, region_str, NULL);
  region->ed_pos = MIN(region->ed_pos, header->target_len[region->chr_ID]);
}

void bam_load_index(Bam_file * bf)
{
  if (NULL != bf->_hidx)//return directly when header file exist
	  return;

  char *bam_cram_file_name = bf->_hfp->fn;

  // check whether index files exist
  char index_name[1024];
  if(
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".bai"))) &&
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".csi"))) &&
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".crai")))&&
		  (!fexist( strTailModify(index_name, bam_cram_file_name, 'i')))
  ){
	  //build index for cram/bam files
	  fprintf(stderr, "BAM/CRAM index is not available for file %s, now building new index file for it.\n", bam_cram_file_name);
	  bam_build_index(bam_cram_file_name);
	  fprintf(stderr, "End for building BAM/CRAM index for file %s\n", bam_cram_file_name);
  }

  bf->_hidx = sam_index_load(bf->_hfp, bam_cram_file_name);
  xassert(NULL != bf->_hidx, "BAM/CRAM index is not available for file");
}

int bam_build_index(const char *fn)
{
	return sam_index_build3(fn, NULL, 0, 4);
}

void resetRegion_ID(Bam_file * bf, R_region *region)
{
  if (NULL != bf->_hitr)
	  hts_itr_destroy(bf->_hitr);
  bam_load_index(bf);

  xassert(region->chr_ID >= 0, "Invalid region specified for BAM/CRAM file");

  bf->_hitr = sam_itr_queryi(bf->_hidx, region->chr_ID, region->st_pos, region->ed_pos);
  xassert(bf->_hitr != NULL, "Failed to fetch region specified for BAM/CRAM file");
  bf->_is_region = true;
  //bf->_region.clear();//todo::

  bf->_is_record_set = false;
  bf->_record_no     = 0;
}

void resetRegion_char(Bam_file * bf, char* region_str)
{
	R_region region;
	region_string2R_region(bf->_hdr, region_str, &region);
    resetRegion_ID(bf,  &region);
    bf->_region = region_str;
}

void bam_file_open(const char* filename, const char* referenceFilename, char* region, Bam_file * bf)
{
	xassert(filename != NULL, "Can't initialize bam_streamer with empty filename\n");
	xassert(*filename != '\0', "Can't initialize bam_streamer with empty filename\n");
	//try open:
	FILE* try_open = xopen(filename, "r");
	fclose(try_open);
	memset(bf, 0, sizeof(Bam_file));
	bf->_hfp = hts_open(filename, "rb");
	xassert(bf->_hfp != NULL, "Failed to open SAM/BAM/CRAM file for reading\n");

	//set reference file
	if (NULL != referenceFilename)
	{
		char referenceFilenameIndex[1024];
		strcpy(referenceFilenameIndex, referenceFilename);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(bf->_hfp, referenceFilenameIndex);
		xassert(ret == 0, "Failed to use reference for BAM/CRAM file");
	}
	bf->_hdr = sam_hdr_read(bf->_hfp);
	xassert(bf->_hdr != NULL, "Failed to parse header from SAM/BAM/CRAM file");

	//set region
	if (NULL == region)
	{
		// setup to read the whole BAM file by default if resetRegion() is not called:
		if (bf->_hdr->n_targets)
			// parse any contig name so that header->hash is created; ignore returned tid value, so doesn't matter if fake name exists
		bam_name2id(bf->_hdr, "fake_name");
	}
	else
		// read a specific region of the bam file:
		resetRegion_char(bf, region);
}

//read a record and stored in "_brec"
bool bam_next(Bam_file * bf)
{
	if (NULL == bf->_hfp)
		return false;
	int ret;
	if (NULL == bf->_hitr)
	{
		ret = sam_read1(bf->_hfp, bf->_hdr, &(bf->_brec));
		// Semi-documented sam_read1 API: -1 is expected read failure at end of stream, any other negative value
		xassert(ret >= -1, "Unexpected return value from htslib sam_read1 function while attempting to read BAM/CRAM file:\n");
	}
	else
	{
		ret = sam_itr_next(bf->_hfp, bf->_hitr, &(bf->_brec));
		xassert(ret >= -1, "Unexpected return value from htslib sam_read1 function while attempting to read BAM/CRAM file:\n");
	}
	bf->_is_record_set = (ret >= 0);
	if (bf->_is_record_set) bf->_record_no++;

	return bf->_is_record_set;
}

//read a record and stored in "r"
bool bam_next_dump(Bam_file * bf, bam1_t *r){
	if (NULL == bf->_hfp)
		return false;
	int ret = (NULL == bf->_hitr)?sam_read1(bf->_hfp, bf->_hdr, r):sam_itr_next(bf->_hfp, bf->_hitr, r);
	xassert(ret >= -1, "Unexpected return value from htslib sam_read1 function while attempting to read BAM/CRAM file:\n");
	bf->_is_record_set = (ret >= 0);
	if (bf->_is_record_set) bf->_record_no++;
	return bf->_is_record_set;
}

void bam_file_close(Bam_file * bf){
	hts_close(bf->_hfp);
	bam_hdr_destroy(bf->_hdr);
	hts_idx_destroy(bf->_hidx);
	hts_itr_destroy(bf->_hitr);
	if(bf->_stream_name) free(bf->_stream_name);
	if(bf->_region) free(bf->_region);
	memset(bf, 0, sizeof(Bam_file));
}

void bam2fastqWrite(FILE* output_file, bam1_t *br, char * seq_buff, uint8_t * qual_buff){
	const int read_len = br->core.l_qseq;
	get_bam_seq(0, read_len, seq_buff, br);//store in binary format
	get_bam_quality_str(0, read_len, qual_buff, br);
	if(!bam_is_unmapped(br) && !bam_is_fwd_strand(br)){//reverse the string and qual when bam is reversed
		getReverseStr_char(seq_buff, read_len);
		getReverseStr_qual(qual_buff, read_len);
	}
	fprintf(output_file, ""
			"@%s\n"
			"%s\n"
			"+%s\n"
			"%s\n",
			bam_qname(br), seq_buff, bam_qname(br), qual_buff);
}

void bam2fastqWrite_additional_str(FILE* output_file, bam1_t *br, char * seq_buff, uint8_t * qual_buff, char* additionalStr1, char* additionalStr2){
	const int read_len = br->core.l_qseq;
	get_bam_seq(0, read_len, seq_buff, br);//store in binary format
	get_bam_quality_str(0, read_len, qual_buff, br);
	if(!bam_is_unmapped(br) && !bam_is_fwd_strand(br)){//reverse the string and qual when bam is reversed
		getReverseStr_char(seq_buff, read_len);
		getReverseStr_qual(qual_buff, read_len);
	}
	fprintf(output_file, ""
			"@%s %s%s\n"
			"%s\n"
			"+\n"
			"%s\n",
			bam_qname(br), additionalStr1, additionalStr2,
			seq_buff,
			qual_buff);
}

//all defined in faidx.h, //use follow functions to handle fna file and fai files
int reference_index_build(const char *ref_file_name)
{
	return fai_build(ref_file_name);
}

//Build fai index for fna file
void reference_index_free(faidx_t *fai)
{
	fai_destroy(fai);
}

faidx_t *reference_index_load(const char *fn)
{
	  // check whether index files exist
	  const char index_name[1024];
	  if(!fexist( strcmb(index_name, fn, ".fai")))
	  {
		  //build index for fna files
		  fprintf(stderr, "Index is not available for file %s, now building new index file for it.\n", fn);
		  reference_index_build(index_name);
		  fprintf(stderr, "End for building index for file %s\n", index_name);
	  }

	  faidx_t *index = fai_load(fn);
	  xassert(NULL != index, "FAI index is not available");
	  return index;
}

//get a reference sequence from REF, need to free result after using this function
// 	@param  fai:  Pointer to the faidx_t struct
// 	@param  region:  Region in the format "chr2:20,000-30,000"
// 	@param  len(is a return value):  Length of the region; -2 if seq not present, -1 general error
// 	@return      Pointer to the sequence; `NULL` on failure
char *get_reference_region(const faidx_t *fai, const char *region, int *len)
{
	 return fai_fetch(fai, region, len);
}

char *get_reference_region_by_chrID(Bam_file * bf, const faidx_t *fai, R_region *region, int *len)
{
	char region_str[1024];
	R_region2region_string(bf->_hdr, region_str, region);
	return fai_fetch(fai, region_str, len);
}

uint8_t getReverseStrEXC[16] = {3, 2, 1, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
void getReverseStr(uint8_t * s, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half; i++){
		uint8_t tmp = s[i];
		s[i] = getReverseStrEXC[s[len - 1 - i]];
		s[len - 1 - i] = getReverseStrEXC[tmp];
	}
	if(len & 0x1){
		s[len_half] = getReverseStrEXC[s[len_half]];
	}
}

void getReverseStr_4bit(uint8_t * s, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half + 1; i++){
		int ri = len - 1 - i;
		uint8_t tmp1 = getReverseStrEXC[bam_seqi(s, i)];
		uint8_t tmp2 = getReverseStrEXC[bam_seqi(s, ri)];
		s[i >>1] = (s[ i>>1] & (0xf << ((( i)&1)<<2))) + (tmp2 << ((~( i)&1)<<2));
		s[ri>>1] = (s[ri>>1] & (0xf << (((ri)&1)<<2))) + (tmp1 << ((~(ri)&1)<<2));
	}
}

char getReverseChar(char c){
	switch(c){
	case 'A': case 'a': return 'T';
	case 'C': case 'c': return 'G';
	case 'G': case 'g': return 'C';
	case 'T': case 't': return 'A';
	}
	return 'N';
}

void getReverseStr_char(char * s, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half; i++){
		char tmp = s[i];
		s[i] = getReverseChar(s[len - 1 - i]);
		s[len - 1 - i] = getReverseChar(tmp);
	}
	if(len & 0x1){
		s[len_half] = getReverseChar(s[len_half]);
	}
}

void getReverseStr_qual(uint8_t * q, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half + 1; i++){
		int ri = len - 1 - i;
		uint8_t tmp = q[i ];
		q[i ] = q[ri];
		q[ri] = tmp;
	}
}

void getReverseStr_qual_char(char * q, int len){
	int len_half = (len >> 1);
	for(int i = 0; i < len_half + 1; i++){
		int ri = len - 1 - i;
		char tmp = q[i];
		q[i ] = q[ri];
		q[ri] = tmp;
	}
}

/***********************************************
 *											   *
 * 					BAM handler                *
 *											   *
 **********************************************/

//when  _bp == NULL, show header, otherwise show record
void print_bam_record(FILE * log, bam1_t* br)
{
	if(br == NULL)
	{
		fprintf(log,
			//basics
				"q name "				"chr_ID "				"position "				"direction "				"mapQ "
			//mate
				"is_mate_unmapped "				"mate_chrom_ID "				"mate_position "				"mate_direction "				"insert_size "
			//flags:
				"secondary "				"supplementary "				"passing_QC "				"properly_alignment "				"PCR_duplication "
			//strings:
				"seq "				"CIGAR "
			//SA tag:
				"chr_ID "				"pos "				"CIGAR "				"MAPQ "
				);
		return;
	}

	//basics: 	q name/chr_ID/position/direction/mapQ
	fprintf(log,
			"qname:%s\t"			"tid:%d\t"			"pos:%d\t"			"D:%d\t"					"BF:%d\t"			"mapQ:%d\t",
			bam_qname(br),			br->core.tid,		br->core.pos,		bam_is_fwd_strand(br),		bam_is_first(br),	br->core.qual
	);

	//mate: 	is_mate_unmapped/mate chrom_ID/mate position/mate direction/insert size/
	fprintf(log,
			"MUM:%d "					"mtid:%d\t"			"mpos:%d\t"			"MD:%d "						"MD1:%d "													"isize:%d\t",
			bam_is_mate_unmapped(br),	br->core.mtid,		br->core.mpos,		bam_is_mate_fwd_strand(br),		bam_is_mate_fwd_strand(br) + bam_is_fwd_strand(br),			br->core.isize
	);

	//flags:	secondary/supplementary/passing QC/properly alignment/PCR duplication
	fprintf(log,
			"SEC:%d "					"SUP:%d "					"QC:%d "				"PAIR:%d "					"DUP:%d ",
			bam_is_secondary(br),		bam_is_supplementary(br),	bam_is_filter(br),		bam_is_proper_pair(br),		bam_is_duplicate(br)
	);

	//strings: 	seq/CIGAR
	char q_seq[1024];
	get_bam_seq(0, br->core.l_qseq, q_seq, br);
	fprintf(log, "%s ", q_seq);
	print_cigar(br, log);

	//SA tag:	chr_ID/pos/CIGAR/MAPQ
	char* SA_tag_char = bam_get_string_tag(br, "SA");
	if(SA_tag_char != NULL)
		fprintf(log, "%s ", SA_tag_char);

	fprintf(log, "\n");
}


/// Get string AUX field, return NULL if field is not found, or field if not a string
///
/// \param[in] tag AUX field tag. This is a char array of length two, null term is not required
///
/// example tag: static const char smtag[] = {'S','M'};
///
char* bam_get_string_tag(bam1_t* _bp, const char* tag)
{
  // retrieve the BAM tag
  uint8_t* pTag = bam_aux_get(_bp, tag);
  if (!pTag) return NULL;

  // skip tags that are not encoded as a null-terminated string
  if (pTag[0] != 'Z') return NULL;
  ++pTag;

  return (char*)pTag;
}

bool bam_is_int_code(char c)
{
	switch (c) {
	case 'c':
	case 's':
	case 'i':
	case 'C':
	case 'S':
	case 'I':
		return true;
	default:
		return false;
	}
}

//return true when successfully get the number tag
bool bam_get_num_tag(bam1_t* _bp, char* tag, int32_t* num)
{
  // retrieve the BAM tag
  uint8_t* pTag = bam_aux_get(_bp, tag);
  if (!pTag) return false;

  // skip tags that are not encoded as integers
  if (!bam_is_int_code(pTag[0]))
	  return false;
  *num = bam_aux2i(pTag);

  return true;
}

//return true when with SA tag
bool bam_get_sa_tag(bam1_t* b, SA_tag *sa, bam_hdr_t *h)
{
	char* SA_tag_char = bam_get_string_tag(b, "SA");
	if(SA_tag_char == NULL)
		return false;
	char SA_tag_char_cp[200], *char_p = SA_tag_char_cp;
	strcpy(SA_tag_char_cp, SA_tag_char);
	while(*char_p++)
		if(*char_p == ',')
			*char_p = ' ';

	char chr_ID[200];
	char direction;
	sscanf(SA_tag_char_cp, "%s %d %c %s %*d %d", chr_ID, &(sa->pos), &direction, sa->cigar, &(sa->mapq));
	sa->chrom_ID = bam_name2id(h, chr_ID);
	sa->direction = (direction == '+')?FORWARD:REVERSE;
	sa->next = 0;
	SA_cigar_break_point(sa);
	if(sa->break_point_left == 0 && sa->break_point_right == 0)
		return false;
	return true;
}

inline uint8_t* bam_qname(bam1_t* _bp)					{	return _bp->data;}
//inline void bam_set_qname(bam1_t* _bp, char* name)	{	edit_bam_qname(name, *_bp);}//todo::
inline bool bam_is_paired(bam1_t* _bp)				{	return ((_bp->core.flag & BAM_PAIRED) != 0);}
inline bool bam_is_proper_pair(bam1_t* _bp) 		{	return ((_bp->core.flag & BAM_PROPER_PAIR) != 0);}
inline bool bam_is_unmapped(bam1_t* _bp) 			{	return ((_bp->core.flag & BAM_UNMAPPED) != 0);}
inline bool bam_is_mate_unmapped(bam1_t* _bp) 		{	return ((_bp->core.flag & BAM_MATE_UNMAPPED) != 0);}
inline bool bam_is_fwd_strand(bam1_t* _bp) 			{	return (((_bp->core.flag & BAM_STRAND) == 0));}
inline bool bam_is_mate_fwd_strand(bam1_t* _bp) 	{	return (((_bp->core.flag & BAM_MATE_STRAND) == 0));}
inline bool bam_is_first(bam1_t* _bp) 				{	return ((_bp->core.flag & BAM_FIRST_READ) != 0);}
inline bool bam_is_second(bam1_t* _bp) 				{	return ((_bp->core.flag & BAM_SECOND_READ) != 0);}
inline bool bam_is_secondary(bam1_t* _bp) 			{	return ((_bp->core.flag & BAM_SECONDARY) != 0);}
inline bool bam_is_filter(bam1_t* _bp) 				{	return ((_bp->core.flag & BAM_FILTER) != 0);}
inline bool bam_is_duplicate(bam1_t* _bp) 				{	return ((_bp->core.flag & BAM_DUPLICATE) != 0);}
inline bool bam_is_supplementary(bam1_t* _bp) 		{	return ((_bp->core.flag & BAM_SUPPLEMENTARY) != 0);}

bool bam_is_DR_signal(bam1_t* br, int i_size_min, int i_size_MAX)
{
	if(bam_is_mate_unmapped(br) || bam_is_unmapped(br))//condition: unmapped
		return false;
	//if(!bam_is_proper_pair(br))
	//	return true;
	if(br->core.tid != br->core.mtid)
		return true;
	bool direction = bam_is_fwd_strand(br);
	if(direction == bam_is_mate_fwd_strand(br))
		return true;

	int insert_size = ABS(br->core.isize);
	if(insert_size <= i_size_min || insert_size >= i_size_MAX)
		return true;
	if(direction == FORWARD && br->core.pos > br->core.mpos)
		return true;
	if(direction == REVERSE && br->core.pos < br->core.mpos)
		return true;

	return false;
}

//just equal : !bam_is_DR_signal(); but is more understandable
bool bamIsNormalDrPair(bam1_t* br, int i_size_min, int i_size_MAX)
{
	bool direction = bam_is_fwd_strand(br);
	bool mate_dir = bam_is_mate_fwd_strand(br);
	int ABSisize = ABS(br->core.isize);

	if(br->core.tid == br->core.mtid && direction != mate_dir && ABSisize < i_size_MAX && ABSisize > i_size_MAX){
		if(direction == FORWARD && br->core.pos <= br->core.mpos) return true;//normal condition 1
		if(direction == REVERSE && br->core.pos <= br->core.mpos) return true;//normal condition 2
	}
	return false;
}

//not directly set flags to be 1, but set 0 to be 1, and set 1 tp be 0
inline void bam_set_paired(bam1_t* _bp) 		{	_bp->core.flag ^= BAM_PAIRED;}
inline void bam_set_filtered(bam1_t* _bp) 		{	_bp->core.flag ^= BAM_FILTER;}
inline void bam_set_unmapped(bam1_t* _bp) 		{	_bp->core.flag ^= BAM_UNMAPPED;}
inline void bam_set_mate_unmapped(bam1_t* _bp)	{	_bp->core.flag ^= BAM_MATE_UNMAPPED;}
inline void bam_set_rvs_strand(bam1_t* _bp)		{	_bp->core.flag ^= BAM_STRAND;}
inline void bam_set_mate_rvs_strand(bam1_t* _bp){	_bp->core.flag ^= BAM_MATE_STRAND;}
inline void bam_set_duplicate(bam1_t* _bp)		{	_bp->core.flag ^= BAM_DUPLICATE;}
inline void bam_set_first(bam1_t* _bp) 			{	_bp->core.flag ^= BAM_FIRST_READ;}
inline void bam_set_second(bam1_t* _bp) 		{	_bp->core.flag ^= BAM_SECOND_READ;}
inline void bam_set_secondary(bam1_t* _bp)		{	_bp->core.flag ^= BAM_SECONDARY;}
inline void bam_set_supplementary(bam1_t* _bp) 	{	_bp->core.flag ^= BAM_SUPPLEMENTARY;}

int bam_read_no(bam1_t* _bp) 		{	return ((bam_is_second(_bp) && (!bam_is_first(_bp))) ? 2 : 1);}
int bam_target_id(bam1_t* _bp)  	{	return _bp->core.tid;}
int bam_mate_target_id(bam1_t* _bp) {	return _bp->core.mtid;}
int bam_pos(bam1_t* _bp)  			{	return (_bp->core.pos + 1);}
int bam_mate_pos(bam1_t* _bp)  		{	return (_bp->core.mpos + 1);}
uint8_t bam_map_qual(bam1_t* _bp) 	{	return _bp->core.qual;}
uint8_t bam_mate_map_qual(bam1_t* _bp){	int32_t num = 0; bam_get_num_tag(_bp, "MQ", &num); return (uint8_t)num; }


bool bam_is_chimeric(bam1_t* _bp)
{
	return ((bam_target_id(_bp) != bam_mate_target_id(_bp)) &&
			(bam_target_id(_bp) >= 0) && (bam_mate_target_id(_bp) >= 0));
}
/// \brief Test if this read contains an 'SA' tag, used to annotate split read alignments
///
/// \return True if the 'SA' tag is found
bool bam_isSASplit(bam1_t* _bp)
{
	static char satag[] = { 'S', 'A' };
	return (NULL != bam_get_string_tag(_bp, satag));
}

/// \brief Test if this read contains an 'MC' tag, containing mate cigar alignment information
///
/// \return True if the 'MC' tag is found
bool bam_hasMateCigar(bam1_t* _bp)
{
	static char satag[] = { 'M', 'C' };
	return (NULL != bam_get_string_tag(_bp, satag));
}

//1 for A, 2 for C, 4 for G,
// 8 for T and 15 for N.
void get_bam_seq(int st_pos, int end_pos, char * seq, bam1_t* _bp)
{
	uint8_t *bam_seq = bam_get_seq(_bp);
	end_pos = MIN(end_pos, _bp->core.l_qseq);

	while(st_pos < end_pos)
	{
		uint8_t c_c = bam_seqi(bam_seq, st_pos);
		switch(c_c)
		{
		case 1: *seq++ ='A'; break;
		case 2: *seq++ ='C'; break;
		case 4: *seq++ ='G'; break;
		case 8: *seq++ ='T'; break;
		case 15: *seq++ ='N'; break;
		default: fprintf(stderr, "Wrong base!");
		}
		st_pos++;
	}
	*seq = 0;
}

//0~3: ACGT
//when contain N or other chars, return false
bool get_bam_seq_bin(int st_pos, int end_pos, uint8_t * seq, bam1_t* _bp)
{
	uint8_t *bam_seq = bam_get_seq(_bp);
	end_pos = MIN(end_pos, _bp->core.l_qseq);
	int n_count = 0;
	while(st_pos < end_pos)
	{
		uint8_t c_c = bam_seqi(bam_seq, st_pos);
		switch(c_c)
		{
		case 1: *seq++ =0; break;
		case 2: *seq++ =1; break;
		case 4: *seq++ =2; break;
		case 8: *seq++ =3; break;
		default: *seq++ =0; n_count ++; break;
		}
		st_pos++;
	}
	if(n_count > 10)
		return false;
	return true;
}

void get_bam_quality_str(int st_pos, int end_pos, uint8_t * quality, bam1_t* _bp)
{
	uint8_t* bam_quality = bam_get_qual(_bp) + st_pos;
	end_pos = MIN(end_pos, _bp->core.l_qseq);
	while(st_pos < end_pos)
	{
		*quality++ = (*(bam_quality++)) + 33;
		st_pos++;
	}
	*quality = 0;
}

int get_bam_low_quality_num(int st_pos, int end_pos, uint8_t min_quality, bam1_t* _bp)
{
	int low_quality_number = 0;
	uint8_t* bam_quality = bam_get_qual(_bp);
	end_pos = MIN(end_pos, _bp->core.l_qseq);
	while(st_pos < end_pos){
		if((*(bam_quality++)) < min_quality)
			low_quality_number++;
		st_pos++;
	}
	return low_quality_number;
}

//-------------------------------MAPQ-----------------------
unsigned bam_alt_map_qual(bam1_t* _bp, const char* tag)
{
	uint8_t* alt_ptr = bam_aux_get(_bp, tag);
	if((NULL != alt_ptr) && bam_is_int_code(alt_ptr[0]))
	{
		const int alt_map = bam_aux2i(alt_ptr);
		xassert_3(alt_map >= 0, "Unexpected negative value in optional BAM/CRAM tag: [%s]\n", tag, NULL, NULL);
		return alt_map;
	}
	else
	{
		uint8_t mapQ = bam_map_qual(_bp);
		return mapQ;
	}
}

/// return single read mapping score if it exists,
/// else return MAPQ:
unsigned bam_se_map_qual(bam1_t* _bp)
{
	static const char smtag[] = { 'S', 'M' };
	return bam_alt_map_qual(_bp, smtag);
}

int32_t bam_template_size(bam1_t* _bp)
{
	return _bp->core.isize;
}

const uint32_t* bam_raw_cigar(bam1_t* _bp){
	return bam_get_cigar(_bp);
}
unsigned bam_n_cigar(bam1_t* _bp){
	return _bp->core.n_cigar;
}

unsigned bam_read_size(bam1_t* _bp){
	return _bp->core.l_qseq;
}

//bam_seq bam_get_bam_read(bam1_t* _bp)todo:
//{
//	return bam_seq(bam_get_seq(_bp), read_size());
//}

const uint8_t* bam_qual(bam1_t* _bp){ return bam_get_qual(_bp);}
void bam_set_target_id(bam1_t* _bp, int32_t tid) {
	if (tid < -1)
		tid = -1;
	_bp->core.tid = tid;
}

// read should be null terminated, qual should already have offset removed:todo
//
//void bam_set_readqual(const char* read, const uint8_t* init_qual) {
//	edit_bam_read_and_quality(read, init_qual, *_bp);
//}

bool bam_empty(bam1_t* _bp){
	xassert(NULL != _bp, "");
	return (_bp->l_data == 0);
}

void bam_freeBam(bam1_t* _bp)
{
	if (NULL != _bp) {
		if (NULL != _bp->data)
			free(_bp->data);
		free(_bp);
	}
}

bool isReadFilteredCore(bam1_t* b)
{
  if (bam_is_filter(b))
    return true;
  else if (bam_is_duplicate(b))
    return true;
  // supplementary reads without SA tag
  else if (bam_is_supplementary(b) && (!bam_isSASplit(b)))
    return true;
  else
  {
    // hack to work with bwamem '-M' formatting,
    // keep secondary reads when they contain an SA tag
    if (bam_is_secondary(b))
    {
      if (!bam_isSASplit(b)) return true;
    }
  }
  return false;
}

bool isReadUnmappedOrFilteredCore(bam1_t* b)
{
  if (isReadFilteredCore(b))
	  return true;
  return bam_is_unmapped(b);
}


/// \brief Test if the read is supplemental, using a more liberal community criteria to define
/// 'supplemental' compared to that from the BAM spec.
///
/// Reads are considered supplemental if either:
/// 1. The 'supplemental' bit is set in the bam record.
/// 2. The 'secondary' bit is set in the bam record and the record contains an 'SA' tag.
///
/// The second condition supports the common workaround typified by bwamem's '-M' option,
/// which allows split reads to be added to the alignment without creating BAM's which could break
/// on older tools.
///
/// \return True if this read is treated as supplemental
bool isNonStrictSupplement(bam1_t* b)
{
  if (bam_is_supplementary(b)) return true;
  if (!bam_is_secondary(b)) return false;
  return bam_isSASplit(b);
}

bool bam_is_mapped_pair(bam1_t* b)
{
  if (!bam_is_paired(b)) return false;
  if (bam_is_unmapped(b) || bam_is_mate_unmapped(b)) return false;
  return true;
}

bool is_mapped_chrom_pair(bam1_t* b)
{
  if (!bam_is_mapped_pair(b)) return false;
  if (b->core.tid != b->core.mtid) return false;
  return true;
}

/*********************************************************
 *
 *							CIGAR
 *
 ********************************************************/
inline char segment_type_to_cigar_code(const int id)
{
	switch (id)
	{
	case CIGAR_MATCH:
		return 'M';
	case CIGAR_INSERT:
		return 'I';
	case CIGAR_DELETE:
		return 'D';
	case CIGAR_SKIP:
		return 'N';
	case CIGAR_SOFT_CLIP:
		return 'S';
	case CIGAR_HARD_CLIP:
		return 'H';
	case CIGAR_PAD:
		return 'P';
	case CIGAR_SEQ_MATCH:
		return '=';
	case CIGAR_SEQ_MISMATCH:
		return 'X';
	default:
		return 'X';
	}
}

inline int cigar_code_to_segment_type(const char c)
{
  switch (c) {
  case 'M':
    return CIGAR_MATCH;
  case 'I':
    return CIGAR_INSERT;
  case 'D':
    return CIGAR_DELETE;
  case 'N':
    return CIGAR_SKIP;
  case 'S':
    return CIGAR_SOFT_CLIP;
  case 'H':
    return CIGAR_HARD_CLIP;
  case 'P':
    return CIGAR_PAD;
  case '=':
    return CIGAR_SEQ_MATCH;
  case 'X':
    return CIGAR_SEQ_MISMATCH;
  default:
    return CIGAR_NONE;
  }
}

inline bool is_segment_type_read_length(const int id)
{
  switch (id)
  {
  case CIGAR_MATCH:
  case CIGAR_INSERT:
  case CIGAR_SOFT_CLIP:
  case CIGAR_SEQ_MATCH:
  case CIGAR_SEQ_MISMATCH:
    return true;
  default:
    return false;
  }
}

inline bool is_segment_type_ref_length(const int id)
{
  switch (id) {
  case CIGAR_MATCH:
  case CIGAR_DELETE:
  case CIGAR_SKIP:
  case CIGAR_SEQ_MATCH:
  case CIGAR_SEQ_MISMATCH:
    return true;
  default:
    return false;
  }
}

inline bool is_segment_align_match(const int id)
{
  switch (id) {
  case CIGAR_MATCH:
  case CIGAR_SEQ_MATCH:
  case CIGAR_SEQ_MISMATCH:
    return true;
  default:
    return false;
  }
}

inline bool is_segment_type_indel(const int id)
{
  switch (id) {
  case CIGAR_INSERT:
  case CIGAR_DELETE:
    return true;
  default:
    return false;
  }
}

void bam_cigar_to_apath(const uint32_t* bam_cigar, const unsigned n_cigar, path_t* apath)
{
  // this assertion isn't really required...
  //    assert(n_cigar>0);
	kv_resize(path_segment, *apath, n_cigar);
	for (int i = 0; i < n_cigar; ++i)
	{
		apath->a[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		apath->a[i].type   = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
  }
}

void get_cigar(bam1_t* b, path_t* apath)
{
	uint32_t* bam_cigar = bam_get_cigar(b);
	bam_cigar_to_apath(bam_cigar, b->core.n_cigar, apath);
}

void print_cigar(bam1_t* b, FILE *log_f)
{
	uint32_t* bam_cigar = bam_get_cigar(b);

	for (int i = 0; i < b->core.n_cigar; ++i)
	{
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		char c_type = segment_type_to_cigar_code(type);
		int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		fprintf(log_f,"%d%c",length,  c_type);
	}
	fprintf(log_f," ");
}

/*
 * FROM help in BWA-MEM
-B INT        penalty for a mismatch [4]
-O INT[,INT]  gap open penalties for deletions and insertions [6,6]
-E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
-L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
-U INT        penalty for an unpaired read pair [17]
*/

#define BWA_MEM_LIKE_match_score 1
#define BWA_MEM_LIKE_mismatch_penalty 2
#define BWA_MEM_LIKE_gap_open 4
#define BWA_MEM_LIKE_gap_extension 1
#define BWA_MEM_LIKE_clipping_penalty 1

uint32_t getScoreByCigar_BWA_MEM_LIKE(bam1_t* b){
	uint32_t score = 0;
	uint32_t* bam_cigar = bam_get_cigar(b);

	for (int i = 0; i < b->core.n_cigar; ++i)
	{
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		switch (type)
		{
		case CIGAR_MATCH:
			score += length; break;
		case CIGAR_INSERT:
		case CIGAR_DELETE:
			score -= (BWA_MEM_LIKE_gap_open + BWA_MEM_LIKE_gap_extension*(length - 1)); break;
		case CIGAR_SKIP:
			break;
		case CIGAR_SOFT_CLIP:
		case CIGAR_HARD_CLIP:
			score -= BWA_MEM_LIKE_clipping_penalty; break;
		case CIGAR_PAD:
			break;
		case CIGAR_SEQ_MATCH:
			score += length; break;
		case CIGAR_SEQ_MISMATCH:
			break;
		default:
			break;
		}
	}
	int NM = 0;
	bam_get_num_tag(b, "NM", &NM);
	if(score <= BWA_MEM_LIKE_mismatch_penalty*NM)
		return 0;
	else
		return score - BWA_MEM_LIKE_mismatch_penalty*NM;
}

void print_cigar_list(uint32_t* bam_cigar, int len)
{

	for (int i = 0; i < len; ++i)
	{
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		char c_type = segment_type_to_cigar_code(type);
		int length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
		printf("%d%c",length,  c_type);
	}
	printf(" ");
}

//CIGAR with soft/hard clip
//soft_left and right left are 0 when no clip
bool bam_has_SH_cigar(bam1_t* b, int *soft_left, int *soft_right)
{
	if(b->core.n_cigar == 0)//CIGAR un-available
		return false;
	uint32_t* bam_cigar = bam_get_cigar(b);

	int begin_type = (int)(1 + (bam_cigar[0] & BAM_CIGAR_MASK));
	int end_type   = (int)(1 + (bam_cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK));

	int has_SH = false;
	soft_left[0] = 0; soft_right[0] = 0;
	if(	begin_type == CIGAR_SOFT_CLIP || begin_type == CIGAR_HARD_CLIP)
	{
		soft_left[0] = (bam_cigar[0] >> BAM_CIGAR_SHIFT);
		has_SH = true;
	}
	if(	end_type == CIGAR_SOFT_CLIP ||	end_type == CIGAR_HARD_CLIP)
	{
		soft_right[0] = (bam_cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
		has_SH = true;
	}
	return has_SH;
}

//Get total INDEL and NMs in a read
void bam_get_INDEL_NM(bam1_t* b, int *INDEL_SIZE, int *NM_SIZE)
{
	int INDEL = 0;
	uint32_t* bam_cigar = bam_get_cigar(b);
	for(int i = 0; i < b->core.n_cigar; i++)
	{
		int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
		if(type == CIGAR_INSERT || type == CIGAR_DELETE)
			INDEL += (bam_cigar[i] >> BAM_CIGAR_SHIFT);
	}
	*INDEL_SIZE = INDEL;
	int NM = 0;
	bam_get_num_tag(b, "NM", &NM);
	*NM_SIZE = NM;
}

void apath_to_bam_cigar(path_t* apath, uint32_t* bam_cigar)
{
  const unsigned as = apath->n;
  for (int i = 0; i < as; ++i)
  {
    path_segment* ps = apath->a + i;
    xassert(ps->type != CIGAR_NONE, "");
    bam_cigar[i] = (ps->length << BAM_CIGAR_SHIFT | ((ps->type) - 1) );
  }
}

/// get insert size from bam record removing refskip (e.g. spliced) segments
int getFragSizeMinusSkip(bam1_t* b)
{
  int fragSize = ABS(bam_template_size(b));
  if (fragSize == 0)
	  return 0;

	uint32_t* bam_cigar = bam_raw_cigar(b); unsigned n_cigar = bam_n_cigar(b);
	for (int i = 0; i < n_cigar; ++i)
		 if ((int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK)) == CIGAR_SKIP)
			 fragSize -= (bam_cigar[i] >> BAM_CIGAR_SHIFT);

  if (fragSize <= 0)
	  fprintf(stderr, "Unexpected fragment size (%d), deduced from bam record: [%s]\n", fragSize, bam_get_qname(b));

  return fragSize;
}

bool hasRefSkip(bam1_t* b)
{
	uint32_t* bam_cigar = bam_raw_cigar(b); unsigned n_cigar = bam_n_cigar(b);
	for (int i = 0; i < n_cigar; ++i)
		 if ((int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK)) == CIGAR_SKIP)
			 return true;
	return false;
}

void char_cigar_to_path(char *string_cigar, path_t * path)
{
	int length = 0;
	path->n = 0;
	while((*string_cigar) != 0)
	{
		char c = *string_cigar;
		string_cigar++;
		if(c <= '9' && c >= '0')
			length = (length * 10) + c - '0';
		else
		{
			int type = cigar_code_to_segment_type(c);
			xassert(type != CIGAR_NONE, "Wrong Cigar!\n");
			path_segment seg = {type, length};
			xassert(path->n < path->m, "Path not enough space!\n");
			kv_push_2(path_segment, path, seg);
			length = 0;
		}
	}
}


///*********************************************FUNC read loader***************************************************//
//CIGAR adjust:
//adjust timing: before store sort item
//adjust SAM: when M number less 2 and followed by long D(ignore I and other cigar type), like: 25I2M55D121M; 60D148M; etc.
//adjust method: change small M into I, change length of D into 0, pos += (len_small_M + len_D)
//original: 25I2M55D121M (pos: 1913086) ----> 25I2I0D121M (pos:1913143)
int cigar_adjust(uint32_t *cigar_l_, uint32_t *cigar, bool add_blank, int delete_small_tail){
	FILE * output = stdout;
	if(*cigar_l_ == 0) return 0;
	int cigar_l = *cigar_l_;
	if(false){
		fprintf(output, "\noriginal cigar: \t");
		for(int i = 0; i < cigar_l; i++){
			fprintf(output, "%d%c", (cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(cigar[i] & BAM_CIGAR_MASK)]);
		}
		fprintf(output, "\n");
	}

	int M_len = 0; int stable_cigar_ID = 0;
	//search to find a stable long match:
	for(int cigar_ID = 0;cigar_ID < (int)cigar_l; cigar_ID++){
		int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
		if(type == 0){//M
			if(cigar_len > delete_small_tail){ stable_cigar_ID = cigar_ID; break; }
			M_len += cigar_len;
		}
	}
	int position_adjust = 0;
	int insertion_len = 0;
	if(stable_cigar_ID != 0){//need to be adjust at begin of alignment
		position_adjust += M_len;
		insertion_len += M_len;
		for(int cigar_ID = 0;cigar_ID < stable_cigar_ID; cigar_ID++){
			int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
			if(type == 1)	insertion_len += cigar_len;//INS
			else if(type == 2) position_adjust += cigar_len; //DEL
		}
		//store new cigar list:
		int new_cigar_len = 0;
		if(insertion_len != 0)
			cigar[new_cigar_len++] = (insertion_len << BAM_CIGAR_SHIFT) + 1;
		for(int cigar_ID = stable_cigar_ID; cigar_ID < (int)cigar_l; cigar_ID++)
			cigar[new_cigar_len++] = cigar[cigar_ID];
		cigar_l = new_cigar_len;
	}

	//adjust from tail of a cigar
	M_len = 0; stable_cigar_ID = 0;
	//search to find a stable long match:
	for(int cigar_ID = cigar_l - 1; cigar_ID >= 0; cigar_ID--){
		int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
		int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
		if(type == 0){//M
			if(M_len + cigar_len > delete_small_tail){ stable_cigar_ID = cigar_ID; break; }
			M_len += cigar_len;
		}
	}
	insertion_len = 0;
	if(stable_cigar_ID != cigar_l - 1){//need to be adjust at end of alignment
		insertion_len += M_len;
		for(int cigar_ID = cigar_l - 1;cigar_ID > stable_cigar_ID; cigar_ID--){
			int cigar_len =	cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
			int type = cigar[cigar_ID] & BAM_CIGAR_MASK;
			if(type == 1)	insertion_len += cigar_len;//INS
		}
		//store new cigar list:
		if(insertion_len != 0){
			cigar[stable_cigar_ID + 1] = (insertion_len << BAM_CIGAR_SHIFT) + 1;
			cigar_l = stable_cigar_ID + 2;
		}
		else
			cigar_l = stable_cigar_ID + 1;
	}

	if(false){
		fprintf(output, "AfterAdj cigar: \t");
		for(int i = 0; i < cigar_l; i++){
			fprintf(output, "%d%c", (cigar[i] >> BAM_CIGAR_SHIFT), "MIDNSHP=XB"[(cigar[i] & BAM_CIGAR_MASK)]);
		}
		fprintf(output, "\t new pos %d\n", position_adjust);
	}
	xassert((int)(*cigar_l_) >= cigar_l, "");
	if(add_blank){
		for(int i = cigar_l; i < (int)(*cigar_l_); i++){
			cigar[i] = 0;
		}
	}
	else
		*cigar_l_ = cigar_l;
	return position_adjust;
}


bool readFilteredAlignment(bam1_t* b)
{
	bool isMatched = false;
	bool isSkip = false;
	bool isClipped = false;
	bool reverse = !bam_is_fwd_strand(b);
	uint32_t* bam_cigar = bam_raw_cigar(b); unsigned n_cigar = bam_n_cigar(b);
	for (int i = 0; i < n_cigar; ++i)
	{
		uint32_t cigar = reverse?bam_cigar[n_cigar - i - 1]:bam_cigar[i];
		int len = (cigar >> BAM_CIGAR_SHIFT);
		int type = (int)(1 + (cigar & BAM_CIGAR_MASK));

		if (is_segment_align_match(type)) {
			if (isClipped) return true;
			isMatched = true;
		}
		else if (type == CIGAR_SKIP) {
			if (isSkip) return true;
			isSkip = true;
		}
		else if (type == CIGAR_SOFT_CLIP)
			isClipped = true;
		else
			return true;
	}
	return (!isMatched);
}

//input: SA position/read position; CIGAR in char string;
//output: possible break point position
//left break point position is +0, right break point is +M/+D
void SA_cigar_break_point(SA_tag *sa)
{
	//string to path
	path_segment path_a[100];
	path_t path = {0};
	path.m = 100; path.a = path_a;
	char_cigar_to_path(sa->cigar, &path);

	sa->soft_sa_left = sa->soft_sa_right = sa->break_point_left = sa->break_point_right = 0;
	//init
	if(path.n == 0)//CIGAR un-available
		return;

	int begin_type = path_a[0].type;
	int end_type   = path_a[path.n - 1].type;

	//handle left
	if(	begin_type == CIGAR_SOFT_CLIP || begin_type == CIGAR_HARD_CLIP)
	{
		sa->soft_sa_left = path_a[0].length;
		sa->break_point_left = sa->pos;
	}

	//handle right
	if(	end_type == CIGAR_SOFT_CLIP ||	end_type == CIGAR_HARD_CLIP)
	{
		int match_delete_len = 0;
		for(int i = 0; i < path.n; i++)
			if(path_a[i].type == CIGAR_MATCH || path_a[i].type == CIGAR_DELETE)
				match_delete_len += path_a[i].length;
		sa->break_point_right = sa->pos + match_delete_len;
		sa->soft_sa_right = path_a[path.n - 1].length;
	}
}
