#ifndef BAM_FILE_H
#define BAM_FILE_H

#include "utils.h"
#include "kvec.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"



/**********************************************	*
 *											   	*
 * 		BAM FILE handler (I/O and index)   		*
 *											   	*
 ************************************************/
typedef struct Bam_file{
	bool _is_record_set;//set to 1 when _brec has data
	htsFile* _hfp;//the file of BAM/CRAM
	bam_hdr_t* _hdr;//header for BAM/CRAM file
	hts_idx_t* _hidx;//index for bam/cram file
	hts_itr_t* _hitr;//Iterator for bam/cram file
	bam1_t _brec;//current BAM record

	// track for debug only:
	unsigned _record_no;
	char * _stream_name;
	bool _is_region;
	char* _region;
}Bam_file;

kvec_T(char*, strList)
kvec_T(Bam_file, bam_List )

/// @param filename CRAM/BAM input file
/// @param referenceFilename: Corresponding reference file. NULL can be given here, but many CRAM files cannot be read in this case.
/// @param region :Restrict the stream to iterate through a specific region. If 'region' is NULL, the stream is configured to iterate through the entire alignment file.
void bam_file_open(const char* filename, const char* referenceFilename, char* region, Bam_file * bf);
void open_BAM_files(strList * BAM_file_name, char * ref_file_name, bam_List* bam_list);
void bam_file_close(Bam_file * bf);
//when read one data successfully, return 1, otherwise return 0
bool bam_next(Bam_file * bf);
bool bam_next_dump(Bam_file * bf, bam1_t *r);
void bam2fastqWrite(FILE* output_file, bam1_t *br, char * seq_buff, uint8_t * qual_buff);
void bam2fastqWrite_additional_str(FILE* output_file, bam1_t *br, char * seq_buff, uint8_t * qual_buff, char* additionalStr1, char* additionalStr2);
void bamFileReset(Bam_file *bf);

void bam_load_index(Bam_file * bf);
int  bam_build_index(const char * fn);

void getReverseStr(uint8_t * s, int len);
void getReverseStr_4bit(uint8_t * s, int len);
void getReverseStr_char(char * s, int len);
void getReverseStr_qual(uint8_t * q, int len);
void getReverseStr_qual_char(char * q, int len);

//bam read/write
//STEP1: htsFile *hts_open(const char *fn, const char *mode);//open bam file, bf->_hfp = hts_open(filename, "rb");
//STEP2: bam_hdr_t *sam_hdr_read(htsFile *fp);//read header
//STEP4: int sam_hdr_write(htsFile *fp, const bam_hdr_t *h) HTS_RESULT_USED;//write header
//STEP3: int sam_read1(htsFile *fp, bam_hdr_t *h, bam1_t *b) HTS_RESULT_USED;//read record
//STEP5: int sam_write1(htsFile *fp, const bam_hdr_t *h, const bam1_t *b) HTS_RESULT_USED;//write record
//STEP6: int hts_close(htsFile *fp);//close file

/// Set new region to iterate over.
/// @param region: string in format: "chromName:beginPos-endPos", cannot be NULL
/// @param tid: reference contig id
/// @param beginPos endPos: start and end position (zero-indexed, closed) of reference
//region to ID/ ID to region

typedef struct R_region{//break point candidate
	uint16_t	chr_ID;
	int 		st_pos;
	int 		ed_pos;
}R_region;
void R_region2region_string(bam_hdr_t* header,	char* region_str, R_region *region);
void region_string2R_region(bam_hdr_t* header,	char* region_str, R_region *region);
void resetRegion_char(Bam_file * bf, char* region);
void resetRegion_ID(Bam_file * bf, R_region *region);

/********************************************************
 *														*
 * 		   reference file handler: SAMTOOLS index		*
 *														*
 ********************************************************/
//all defined in faidx.h, //use follow functions to handle fna file and fai files
int reference_index_build(const char *ref_file_name);//Build fai index for fna file
void reference_index_free(faidx_t *fai);//Free fai index for fna file, and close reference file
faidx_t *reference_index_load(const char *fn);//load index (build new index when without index file), open reference file
//get a reference sequence from REF, need to free result after using this function
// 	@param  fai:  Pointer to the faidx_t struct
// 	@param  region:  Region in the format "chr2:20,000-30,000"
// 	@param  len(is a return value):  Length of the region; -2 if seq not present, -1 general error
// 	@return      Pointer to the sequence; `NULL` on failure
char *get_reference_region(const faidx_t *fai, const char *region, int *len);
char *get_reference_region_by_chrID(Bam_file * bf, const faidx_t *fai, R_region *region, int *len);

/***********************************************
 *											   *
 * 			BAM record handler         		   *
 *											   *
 **********************************************/

#define  BAM_PAIRED         0x001
#define  BAM_PROPER_PAIR    0x002
#define  BAM_UNMAPPED       0x004
#define  BAM_MATE_UNMAPPED  0x008
#define  BAM_STRAND  		0x010
#define  BAM_MATE_STRAND    0x020
#define  BAM_FIRST_READ     0x040
#define  BAM_SECOND_READ    0x080
#define  BAM_SECONDARY      0x100
#define  BAM_FILTER         0x200
#define  BAM_DUPLICATE      0x400
#define  BAM_SUPPLEMENTARY  0x800


//-----------------------------------------get tags-----------------------------------------
//when  _bp == NULL, show header, otherwise show record
void print_bam_record(FILE * log, bam1_t* br);
char* bam_get_string_tag(bam1_t* _bp, const char* tag);
bool bam_get_num_tag(bam1_t* _bp, char* tag, int32_t* num);
//-----------------------------------------from sam.h-----------------------------------------
//uint8_t* bam_get_data(b)  ((b)->data)
//uint8_t* bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
//uint8_t* bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
//uint8_t* bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
//uint8_t* bam_get_l_aux(b) ((b)->l_data - ((b)->core.n_cigar<<2) - (b)->core.l_qname - (b)->core.l_qseq - (((b)->core.l_qseq + 1)>>1))
//uint8_t  bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)
#define bam_get_data(b) ((b)->data)
void get_bam_seq(int st_pos, int end_pos, char * seq, bam1_t* _bp);
bool get_bam_seq_bin(int st_pos, int end_pos, uint8_t * seq, bam1_t* _bp);//store in binary format
void get_bam_quality_str(int st_pos, int end_pos, uint8_t * quality, bam1_t* _bp);
int get_bam_low_quality_num(int st_pos, int end_pos, uint8_t min_quality, bam1_t* _bp);
//-----------------------------------------bool values-----------------------------------------
// void bam_set_qname(bam1_t* _bp, char* name)	{	edit_bam_qname(name, *_bp);}//todo::
bool bam_is_paired(bam1_t* _bp);
bool bam_is_proper_pair(bam1_t* _bp);
bool bam_is_unmapped(bam1_t* _bp);
bool bam_is_mate_unmapped(bam1_t* _bp);
bool bam_is_fwd_strand(bam1_t* _bp);
bool bam_is_mate_fwd_strand(bam1_t* _bp);
bool bam_is_duplicate(bam1_t* _bp);
bool bam_is_filter(bam1_t* _bp);
bool bam_is_first(bam1_t* _bp);
bool bam_is_second(bam1_t* _bp);
bool bam_is_secondary(bam1_t* _bp);
bool bam_is_supplementary(bam1_t* _bp);
bool bam_is_DR_signal(bam1_t* br, int i_size_min, int i_size_MAX);//discordant read pairs
bool bam_is_UM_signal(bam1_t* _bp);//mate read unmapped pairs
//-----------------------------------------set bool values-----------------------------------------
void bam_set_paired(bam1_t* _bp);
void bam_set_filtered(bam1_t* _bp);
void bam_set_unmapped(bam1_t* _bp);
void bam_set_mate_unmapped(bam1_t* _bp);
void bam_set_rvs_strand(bam1_t* _bp);
void bam_set_mate_rvs_strand(bam1_t* _bp);
void bam_set_duplicate(bam1_t* _bp);
void bam_set_first(bam1_t* _bp);
void bam_set_second(bam1_t* _bp);
void bam_set_secondary(bam1_t* _bp);
void bam_set_supplementary(bam1_t* _bp);
bool bam_is_chimeric(bam1_t* _bp);
bool bam_hasMateCigar(bam1_t* _bp);/// \brief Test if this read contains an 'MC' tag, containing mate cigar alignment information, return True if the 'MC' tag is found
bool bam_isNonStrictSupplement(bam1_t* _bp);
//-----------------------------------------int values----------------------------
int bam_read_no(bam1_t* _bp);
int bam_target_id(bam1_t* _bp);
int bam_mate_target_id(bam1_t* _bp);
int bam_pos(bam1_t* _bp);
int bam_mate_pos(bam1_t* _bp);
uint8_t bam_map_qual(bam1_t* _bp);
uint8_t bam_mate_map_qual(bam1_t* _bp);
int bam_seq_len(bam1_t* _bp);
//-------------------------------MAPQ-----------------------
unsigned bam_alt_map_qual(bam1_t* _bp, const char* tag);
/// return single read mapping score if it exists,
/// else return MAPQ:
unsigned bam_se_map_qual(bam1_t* _bp);
int32_t bam_template_size(bam1_t* _bp);
unsigned bam_read_size(bam1_t* _bp);
const uint8_t* bam_qual(bam1_t* _bp);
void bam_set_target_id(bam1_t* _bp, int32_t tid);
bool bam_empty(bam1_t* _bp);
void bam_freeBam(bam1_t* _bp);
bool isReadFilteredCore(bam1_t* b);
bool isReadUnmappedOrFilteredCore(bam1_t* b);
bool isNonStrictSupplement(bam1_t* b);
bool bam_is_mapped_pair(bam1_t* b);//return true when read and mate all mapped and they are paired
bool is_mapped_chrom_pair(bam1_t* b);

/************************************************
 *											   	*
 * 			SA tag handler         				*
 *											   	*
 ************************************************/

typedef struct SA_tag{
	int chrom_ID;
	int pos;
	int direction;
	int mapq;
	char cigar[100];
	int next;//== 0 when no sa signal,
	int break_point_left;
	int break_point_right;
	int soft_sa_left;
	int soft_sa_right;
}SA_tag;
bool bam_get_sa_tag(bam1_t* b, SA_tag *sa, bam_hdr_t *h);//return null when no SA tag
uint8_t* bam_qname(bam1_t* _bp);//string values
bool bam_isSASplit(bam1_t* _bp);/// \return True if the 'SA' tag is found, brief Test if this read contains an 'SA' tag, used to annotate split read alignments

/********************************************************
 *														*
 *				   CIGAR string handler					*
 *														*
 ********************************************************/


enum align_t { 	CIGAR_NONE, CIGAR_MATCH, CIGAR_INSERT, CIGAR_DELETE, CIGAR_SKIP, CIGAR_SOFT_CLIP,
				CIGAR_HARD_CLIP, CIGAR_PAD, CIGAR_SEQ_MATCH, CIGAR_SEQ_MISMATCH };

//#define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
typedef struct path_segment{
  int  type;
  unsigned length;
}path_segment;

#define pathSTR "N=IDXSHP=X"

kvec_T(path_segment, path_t);
const uint32_t* bam_raw_cigar(bam1_t* _bp);
unsigned 		bam_n_cigar(bam1_t* _bp);
void 			get_cigar(bam1_t* b, path_t* apath);
void 			print_cigar(bam1_t* b, FILE *log_f);

uint32_t getScoreByCigar_BWA_MEM_LIKE(bam1_t* b);

void 			print_cigar_list(uint32_t* bam_cigar, int len);
void bam_cigar_to_apath(const uint32_t* bam_cigar, const unsigned n_cigar, path_t* apath);
char segment_type_to_cigar_code(const int id);
bool bam_has_SH_cigar(bam1_t* b, int *soft_left, int *soft_right);
void bam_get_INDEL_NM(bam1_t* b, int *INDEL_SIZE, int *NM_SIZE);
int cigar_code_to_segment_type(const char c);
bool is_segment_type_read_length(const int id);
bool is_segment_type_ref_length(const int id);
bool is_segment_align_match(const int id);
bool is_segment_type_indel(const int id);
int getFragSizeMinusSkip(bam1_t* b);
bool hasRefSkip(bam1_t* b);
//convert a string cigar to path format
void char_cigar_to_path(char *string_cigar, path_t * path);
int cigar_adjust(uint32_t *cigar_l_, uint32_t *cigar, bool add_blank, int delete_small_tail);
bool readFilteredAlignment(bam1_t* b);

//input: SA position/read position; CIGAR in char string;
//output: possible break point position
//left break point position is +0, right break point is +M/+D
void SA_cigar_break_point(SA_tag *sa);

//modify bam data


#endif
