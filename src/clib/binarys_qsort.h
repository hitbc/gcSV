#ifndef _BINARYS_QSORT_H
#define _BINARYS_QSORT_H
// #include "read_seeding.h"

#define UNPIPATH_OFF_K20

#ifdef UNPIPATH_OFF_K20
int multi_binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset, int64_t seed_binary[], int8_t k_r);
#else
int multi_binsearch_offset(uint32_t x, uint32_t v[], int64_t n, uint32_t offset, int64_t seed_binary[], int8_t k_r);
#endif


int binsearch_range(uint64_t key, uint32_t *v, int64_t n, int64_t *range, int8_t k_off);

//seed_id_r = binsearch_interval_unipath64(kmer_pos_uni, buffer_seqf, result_seqf);
int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n);

int64_t binsearch_interval_unipath(uint32_t x, uint32_t v[], uint32_t n);

#endif
