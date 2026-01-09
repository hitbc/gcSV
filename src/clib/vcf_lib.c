/*
 * vcf_file.c
 *
 *  Created on: 2019-12-28
 *      Author: fenghe
 */

#include <time.h>

#include "desc.h"
#include "vcf_lib.h"

int VCF_next(BCF_FILE *vcf)
{
	if(vcf->header == NULL)
		return false;
	if(bcf_read(vcf->file, vcf->header, &(vcf->r)) == -1)
		return false;
	return true;
}

int VCF_next_dump(BCF_FILE *vcf, bcf1_t *r)
{
	if(vcf->header == NULL)
		return false;
	if(bcf_read(vcf->file, vcf->header, r) == -1)
		return false;
	return true;
}

int VCF_open(BCF_FILE *vcf, const char *fn)
{
	memset(vcf,0,sizeof(BCF_FILE));
	vcf->file = 		hts_open(fn, "r");
	vcf->header = 	bcf_hdr_read(vcf->file);
	strcpy(vcf->fname, fn);
	if(vcf->header == NULL)
		return false;
	return true;
}

int VCF_open_read(BCF_FILE *vcf, const char *fn)
{
	memset(vcf,0,sizeof(BCF_FILE));
	vcf->file = 		hts_open(fn, "r");
	vcf->header = 	bcf_hdr_read(vcf->file);
	strcpy(vcf->fname, fn);
	if(vcf->header == NULL)
		return false;
	if(bcf_read(vcf->file, vcf->header, &(vcf->r)) == -1)
		return false;
	return true;
}

int VCF_open_write(BCF_FILE *vcf, const char *fn, bool is_compression)
{
	memset(vcf,0,sizeof(BCF_FILE));
	if(is_compression)
		vcf->file = hts_open(fn, "wb");//open for write/BCF
	else
		vcf->file = hts_open(fn, "w");//open for write/VCF
	return true;
}

void VCF_close(BCF_FILE *vcf){
	hts_close(vcf->file);
}

//---------------------------CVF INFO-------------------------------------------//
int vcf_get_sample(bcf_hdr_t *hdr, bcf1_t *line, char *sample_name)
{
	int ndst = 1000;
	return bcf_get_info_values(hdr, line, "SAMPLE", (void**)&sample_name, &ndst, BCF_HT_STR);
}

int vcf_get_sv_type(bcf_hdr_t *hdr, bcf1_t *line, char *svtype)
{
	int ndst = 1000;
	return bcf_get_info_values(hdr, line, "SVTYPE", (void**)&svtype, &ndst, BCF_HT_STR);
}

int vcf_get_sv_END(bcf_hdr_t *hdr, bcf1_t *line, int *end)
{
	int ndst = 2;
	void *END_p = end;
	int r = bcf_get_info_values(hdr, line, "END", (void**)(&END_p), &ndst, BCF_HT_INT);
	return r;
}

int vcf_get_sv_LENGTH(bcf_hdr_t *hdr, bcf1_t *line, int *length)
{
	int ndst = 2;
	void *LEN_p = length;
	int r = bcf_get_info_values(hdr, line, "SVLEN", (void**)(&LEN_p), &ndst, BCF_HT_INT);
	return r;
}

int vcf_get_sv_flag(bcf_hdr_t *hdr, bcf1_t *line, const char * tag_name)
{
	int ndst = 100;
	int r = bcf_get_info_flag(hdr,line,tag_name,NULL,&ndst);
	return r;
}

int vcf_get_sv_GT(bcf_hdr_t *hdr, bcf1_t *line, char **GT)
{
	int ndst = 9;
	char ***GT_p = &GT;
	return bcf_get_format_string(hdr, line, "GT", GT_p, &ndst);
}

//support read number
int vcf_get_sv_SR(bcf_hdr_t *hdr, bcf1_t *line, int *SR)
{
	int ndst = 100;
	int **SR_p = &SR;
	return bcf_get_format_int32(hdr, line, "SR", SR_p, &ndst);
}

//------------------------------VCF FORMAT-------------------------------------------//
int vcf_get_genotype_quality(bcf_hdr_t *hdr, bcf1_t *line)
{
	int ndst = 0; int *dst = NULL;
	if ( bcf_get_format_int32(hdr, line, "GQ", &dst, &ndst) > 0 )
		for (int i=0; i<bcf_hdr_nsamples(hdr); i++)
			fprintf(stderr, "%d-", dst[i]);
//	free(dst[0]);
	free(dst);
	fprintf(stderr, "\n");
	return 0;
}

int vcf_get_genotype(bcf_hdr_t *hdr, bcf1_t *line)
{
	int ndst = 0; char**dst = NULL;
	if ( bcf_get_format_string(hdr, line, "GT", &dst, &ndst) > 0 )
		for (int i=0; i<bcf_hdr_nsamples(hdr); i++)
			fprintf(stderr, "%s-", dst[i]);
//	free(dst[0]);
	free(dst);
	fprintf(stderr, "\n");
	return 0;
}

//---------------------------------BEMO-----------------------------------------------//
void DEMO_VCF_READ_WRITE(int nfiles, char ** fn_l)
{
	bool is_compression = false;//BCF or VCF
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, fn_l[0]);//open for read
	BCF_FILE vcf_w;//vcf for write
	VCF_open_write(&vcf_w, fn_l[1], is_compression);
	bcf_hdr_write(vcf_w.file, vcf_r.header);
	while(VCF_next(&vcf_r))//read one
	{
		bcf1_t *c_r = &( vcf_r.r);
		bcf_write(vcf_w.file, vcf_r.header, c_r);
	}
	//close
	bcf_close(vcf_r.file);
	bcf_close(vcf_w.file);
}

bcf_hdr_t *get_vcf_header(const char * ref_fn, bam_hdr_t* header)
{
	char header_line[4096];
	bcf_hdr_t *hearder = bcf_hdr_init("w");
	//BASIC PART
	time_t c_time; struct tm*p;
	time(&c_time);
	p = gmtime(&c_time);
	//date
	sprintf(header_line, "##fileDate=%d%d%d\n",1990+p->tm_year, 1 + p->tm_mon, p->tm_mday); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##source=%s%s\n",PACKAGE_NAME, PACKAGE_VERSION); bcf_hdr_append(hearder, header_line);//software
	sprintf(header_line, "##reference=file://%s\n",ref_fn); bcf_hdr_append(hearder, header_line);
	///INFO PART
	sprintf(header_line, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"sample_id from dbVar submission; every call must have SAMPLE\">\n"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical homology at event breakpoints\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical homology at event breakpoints\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description=\"Length of insertion\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of insertion\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description=\"Known left side of insertion for an insertion of unknown length\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description=\"Known right side of insertion for an insertion of unknown length\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##INFO=<ID=JUNCTION_QUAL,Number=1,Type=Integer,Description=\"If the SV junction is part of an EVENT (ie. a multi-adjacency variant), this field provides the QUAL value for the adjacency in question only\">"); bcf_hdr_append(hearder, header_line);
	///FORMAT PART
	sprintf(header_line, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter, 'PASS' indicates that all filters have passed for this sample\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">"); bcf_hdr_append(hearder, header_line);
	///FILTER PART

	sprintf(header_line, "##FILTER=<ID=Ploidy,Description=\"For DEL & DUP variants, the genotypes of overlapping variants (with similar size) are inconsistent with diploid expectation\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=MaxDepth,Description=\"Depth is greater than 3x the median chromosome depth near one or both variant breakends\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=MaxMQ0Frac,Description=\"For a small variant (<1000 bases), the fraction of reads in all samples with MAPQ0 around either breakend exceeds 0.4\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=NoPairSupport,Description=\"For variants significantly larger than the paired read fragment size, no paired reads support the alternate allele in any sample.\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=MinQUAL,Description=\"QUAL score is less than 20\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=SampleFT,Description=\"No sample passes all the sample-level filters (at the field FORMAT/FT)\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=MinGQ,Description=\"GQ score is less than 15 (filter applied at sample level)\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##FILTER=<ID=HomRef,Description=\"homozygous reference call (filter applied at sample level)\">"); bcf_hdr_append(hearder, header_line);
	//ALT PART
	sprintf(header_line, "##ALT=<ID=DEL,Description=\"Deletion\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##ALT=<ID=INS,Description=\"Insertion\">"); bcf_hdr_append(hearder, header_line);
	sprintf(header_line, "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">"); bcf_hdr_append(hearder, header_line);
	//CONTIG PART
	for(int i = 0; i < header->n_targets; i++)
	{
		sprintf(header_line, "##contig=<ID=%s,length=%d>", header->target_name[i], header->target_len[i]);
		bcf_hdr_append(hearder, header_line);
	}
	sprintf(header_line, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	UNION"); bcf_hdr_append(hearder, header_line);//TODO::
	return hearder;
}

void out_put_vcf_header(const char * ref_fn, bam_hdr_t* bam_header, BCF_FILE *vcf_w)
{
	vcf_w->header = get_vcf_header(ref_fn, bam_header);
	bcf_hdr_write(vcf_w->file, vcf_w->header);
}

