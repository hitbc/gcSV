/*
 * vcf_file.h
 *
 *  Created on: 2019年12月28日
 *      Author: fenghe
 */

#ifndef LIB_VCF_FILE_H_
#define LIB_VCF_FILE_H_

#include "htslib/vcf.h"
#include "utils.h"
#include "htslib/sam.h"

typedef struct
{
    htsFile *file; //file of VCF
    bcf_hdr_t *header;//header file of VCF
    bcf1_t 		r; //current record
    char fname[1024]; //name of vcf file
    bool _is_record_set ;
}
BCF_FILE;


int VCF_next(BCF_FILE *vcf);
int VCF_next_dump(BCF_FILE *vcf, bcf1_t *r);
int VCF_open(BCF_FILE *vcf, const char *fn);
int VCF_open_read(BCF_FILE *vcf, const char *fn);
int VCF_open_write(BCF_FILE *vcf, const char *fn, bool is_compression);
void VCF_close(BCF_FILE *vcf);

//tags:info
int vcf_get_sample(bcf_hdr_t *hdr, bcf1_t *line, char *dst);//tag for sample
int vcf_get_sv_type(bcf_hdr_t *hdr, bcf1_t *line, char *svtype);
int vcf_get_sv_END(bcf_hdr_t *hdr, bcf1_t *line, int *end);
int vcf_get_sv_LENGTH(bcf_hdr_t *hdr, bcf1_t *line, int *length);
int vcf_get_sv_flag(bcf_hdr_t *hdr, bcf1_t *line, const char * tag_name);
int vcf_get_sv_GT(bcf_hdr_t *hdr, bcf1_t *line, char **GT);
int vcf_get_sv_SR(bcf_hdr_t *hdr, bcf1_t *line, int *SR);

//tags: format
int vcf_get_genome_type(bcf_hdr_t *hdr, bcf1_t *line);
int vcf_get_genotype(bcf_hdr_t *hdr, bcf1_t *line);

//generate a new vcf/bcf header

void out_put_vcf_header(const char * ref_fn, bam_hdr_t* bam_header, BCF_FILE *vcf_w);
#endif /* LIB_VCF_FILE_H_ */
