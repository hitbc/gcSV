
#ifndef L_UTILS_H
#define L_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/time.h>

#ifdef __GNUC__
// Tell GCC to validate printf format string and args
#define ATTRIBUTE(list) __attribute__ (list)
#else
#define ATTRIBUTE(list)
#endif

//used to open file bigger than 4G, it`s used to instead fopen64
//it should be defined before stdio.h, the value are used in that file
#undef   _FILE_OFFSET_BITS
#define   _FILE_OFFSET_BITS   64

//BOOL
#ifndef __cplusplus
 #ifndef bool
  typedef int bool;
 #endif

 #ifndef true
  #define true 1
 #endif

 #ifndef false
  #define false 0
 #endif
#endif

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define ABS(a) (((a) > 0)?(a): (- (a)))
#define ABS_U(a,b) (((a) > (b))?((a) - (b)): ((b) - (a)))
#define SWAP(a, b) 							\
		{									\
			__typeof__(a) c = a; 			\
			a = b;							\
			b = c;							\
		}

#define FORWARD 1
#define REVERSE 0

#define MAX_uint8_t 	0xff
#define MAX_uint16_t	0xffff
#define MAX_uint32_t 	0xffffffff
#define MAX_uint64_t 	0xffffffffffffffff
#define MAX_int8_t		0x7f
#define MAX_int16_t		0x7fff
#define MAX_int32t		0x7fffffff
#define MAX_int64_t		0x7fffffffffffffff

#define MASK_41 0x1ffffffffff
#define MASK_40 0xffffffffff
#define MASK_39 0x7fffffffff
#define MASK_38 0x3fffffffff
#define MASK_37 0x1fffffffff
#define MASK_36 0xfffffffff
#define MASK_35 0x7ffffffff
#define MASK_34 0x3ffffffff
#define MASK_33 0x1ffffffff
#define MASK_32 0xffffffff
#define MASK_31 0x7fffffff
#define MASK_30 0x3fffffff


#define one_eighth_GigaByte 0x8000000
#define quarter_GigaByte	0x10000000
#define half_GigaByte 		0x20000000
#define one_GigaByte 		0x40000000
#define two_GigaByte 		0x80000000
#define four_GigaByte 		0x100000000
#define eight_GigaByte 		0x200000000
#define sixteen_GigaByte 	0x400000000

#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

#define xopen(fn, mode) 		err_xopen_core	(__func__, __LINE__, fn, mode)
#define xzopen(fn, mode) 		err_xzopen_core	(__func__, __LINE__, fn, mode)
#define xreopen(fn, mode, fp) 	err_xreopen_core(__func__, fn, mode, fp)

#define xcalloc(count,eltsize) 	err_xcalloc_core(__func__,__LINE__,count,eltsize)
#define xcalloc_t(type, count) 	(type*)err_xcalloc_core(__func__,__LINE__,count,sizeof(type))
#define xmalloc(size)			err_xmalloc_core(__func__,__LINE__,size)
#define xmalloc_t(type, count)	(type*) err_xmalloc_core(__func__,__LINE__,count*sizeof(type))
#define xrealloc(ptr,newSize)	err_xrealloc_core(__func__,__LINE__,ptr,newSize)
#define xfree(p) 				{free(p);p = NULL;}
#define xassert(cond, msg) 		if ((cond) == 0) _err_fatal_simple_core(__func__,__LINE__, msg)
#define xassert_3(cond, msg, m1, m2, m3) if ((cond) == 0) _err_fatal_core(__func__,__LINE__, msg, m1, m2, m3)
#define xread(data, size, count, stream) 	xassert(count == fread(data, size, count, stream),"[xREAD] Wrong in read file, data not enough.")
#define xget_file(dirPath, postfix, opentype) get_file_(dirPath, PACKAGE_NAME, postfix, opentype)

//reallocate memory for a buff, at least 20 new block are allocated
//equal: void BUFF_REALLOC(TYPE * data, uint64_t old_length, uint64_t new_length){};
#define BUFF_REALLOC(p, L_max, L_new) 	\
	if(L_new > L_max)					\
	{									\
		L_max =  L_new + 20; 			\
		p = realloc(p, L_max*(sizeof(__typeof__(*p))));\
	}

#define FUNC_GET_TIME(code_block, sum_time, isInUse) 	\
	struct timeval tv1, tv2;							\
	if(isInUse)											\
		gettimeofday(&tv1,NULL);						\
	code_block;											\
	if(isInUse)											\
	{													\
		gettimeofday(&tv2,NULL);						\
		sum_time += (float)(tv2.tv_sec - tv1.tv_sec) + 	\
		(float)(tv2.tv_usec - tv1.tv_usec)/1000000;		\
	}													\

#define FUNC_GET_TIME_P(code_block, string, isInUse) 		\
	do														\
	{														\
		struct timeval tv1, tv2;							\
		double sum_time;									\
		if(isInUse)											\
			gettimeofday(&tv1,NULL);						\
		code_block;											\
		if(isInUse)											\
		{													\
			gettimeofday(&tv2,NULL);						\
			sum_time = (float)(tv2.tv_sec - tv1.tv_sec) + 	\
			(float)(tv2.tv_usec - tv1.tv_usec)/1000000;		\
			fprintf(stderr, "%s:[%f]\n", string, sum_time);	\
		}													\
	}														\
  	while(0) 												\


/*2 bit bin char to KMER functions
*	each base pair use 1 byte, but only the last 2 bit are used
*	bit2 mean that all base pair are store in "A - 0; C - 1; G - 2; T - 3".
* 	no matter RC or FC, pre/next kmer is the kmer after offset --/++
* 	s of FC kmer is defined as the start pointer of a string, kmer is [s, s + L_kmer - 1]
* 	s of RC kmer is defined as the end   pointer of a string, kmer is [s - L_kmer + 1, s]
* Notice:
* 	use: extern uint64_t kmerMask[33]; before using 'bit2_nextKmer'
*/
#define bit2_nextKmer_init(		s, __L_kmer)   			(binchar2Kmer(  s, __L_kmer) >>2)
#define bit2_preKmer_init(		s, __L_kmer)   			(binchar2Kmer(  s, __L_kmer) <<2)
#define bit2_nextKmerRC_init(	s, __L_kmer) 			(binchar2KmerRC(s, __L_kmer) <<2)

//faster function: predefined mask and move,
//to use those function, you should predefine MASK and MOVE
//#define uint64_t MASK = kmerMask[_kmer];
//#define uint8_t  MOVE = ((_kmer << 1) - 2);
#define bit2_nextKmer(  	s, old_kmer, __L_kmer) 	(((old_kmer) << 2 	| (s)[__L_kmer - 1]) & kmerMask[__L_kmer])
#define bit2_nextKmerMASK(  s, old_kmer, __L_kmer)	((((old_kmer) << 2) | (s)[__L_kmer - 1]) & MASK)
#define bit2_nextKmerRC(	s, old_kmer, __L_kmer) 	((((old_kmer) >> 2) | ((uint64_t)(3 - 	(s)[0]) << ((__L_kmer << 1) - 2))))
#define bit2_nextKmerRCMOVE(s, old_kmer)           	((((old_kmer) >> 2) | ((uint64_t)(3 - 	(s)[0]) << MOVE)))
#define bit2_preKmerMOVE(	s, old_kmer) 			((((old_kmer) >> 2) | ((uint64_t)(		(s)[0]) << MOVE)))

typedef struct {
	uint64_t x, y;
} pair64_t;

//for string
#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

//for kseq
typedef struct __kstream_t {
	unsigned char *buf;
	int begin, end, is_eof;
	gzFile f;
} kstream_t;

typedef struct {
	kstring_t name, comment, seq, qual;
	int last_char;
	kstream_t *f;
} kseq_t;

#ifdef __cplusplus
extern "C" {
#endif
	//file
	bool fexist(const char *filename);
	void fexist_check(const char *filename);
	//err
	void err_fatal(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void err_fatal_core(const char *header, const char *fmt, ...) ATTRIBUTE((noreturn));
	void _err_fatal_simple(const char *func, const char *msg) ATTRIBUTE((noreturn));
	void _err_fatal_simple_core(const char *func, const int line, const char *msg)ATTRIBUTE((noreturn));
	void _err_fatal_core(const char *func, const int line, const char *msg, const char *m1, const char *m2, const char *m3)ATTRIBUTE((noreturn));
	//safety basic I/O
	FILE *err_xopen_core(const char *func,  const int line, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func,  const int line, const char *fn, const char *mode);

    size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
	size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream);

	int err_gzread(gzFile file, void *ptr, unsigned int len);
	int err_gzwrite(gzFile file, void *ptr, unsigned int len);
	int err_fseek(FILE *stream, long offset, int whence);
#define err_rewind(FP) err_fseek((FP), 0, SEEK_SET)
	long err_ftell(FILE *stream);
	long err_fsize(FILE *stream);
	int err_fprintf(FILE *stream, const char *format, ...)
        ATTRIBUTE((format(printf, 2, 3)));
	int err_printf(const char *format, ...)
        ATTRIBUTE((format(printf, 1, 2)));
	int err_fputc(int c, FILE *stream);
#define err_putchar(C) err_fputc((C), stdout)
	int err_fputs(const char *s, FILE *stream);
	int err_puts(const char *s);
	int err_fflush(FILE *stream);
	int err_fclose(FILE *stream);
	int err_gzclose(gzFile file);
	void xmkdir(const char *dirPath);
	void xrm(const char *fileName);
	FILE* get_file_(const char *dirPath, const char *package_name, const char *postfix, const char *opentype);

	//safety basic m-allocation
	void* err_xcalloc_core(const char *func, const int line, size_t count,size_t eltsize);
	void* err_xmalloc_core(const char *func, const int line, size_t size);
	void *err_xrealloc_core(const char *func, const int line, void *ptr, size_t newsize);
	//void  xfree(void **p);

	//time
	double cputime();
	double realtime();
	double realduration(struct timeval start);
	long peak_memory(void);
	char* strcmb(char* to, char *pre, char * suf);
	bool hasEnding(char * fullString, char * ending);
	#define kstring_init(v) ((v).l = (v).m = 0, (v).s = 0)
	#define kstring_initp(v) ((v)->l = (v)->m = 0, (v)->s = 0)

	void ks_introsort_64 (size_t n, uint64_t *a);
	void ks_introsort_128(size_t n, pair64_t *a);

	typedef int (*compar_fn) (const void *, const void *);
	void ksort_stable(void *base, size_t nmemb, size_t size, compar_fn cmp);
	void ksort_stable_mt(void *base, size_t nmemb, size_t size, compar_fn cmp, int n_thread);

	//kseq
	kseq_t *kseq_init(gzFile fd);
	int kseq_read(kseq_t *seq);
	void kseq_destroy(kseq_t *ks);
	void kseq_rewind(kseq_t *ks);
	//ks
	kstream_t *ks_init(gzFile f);
	void ks_destroy(kstream_t *ks);
	//kmer
	uint64_t binchar2Kmer(uint8_t *s, uint8_t L_kmer);
	uint64_t binchar2KmerRC(uint8_t *s, uint8_t L_kmer);
	uint64_t char2Kmer(const char *str_kmer, uint8_t L_kmer, uint8_t *Bit);
	uint64_t char2KmerRC(const char *str_kmer, uint8_t L_kmer, uint8_t *Bit);
	uint64_t getRcKmer(uint64_t kmer, uint8_t l_kmer);
	//64 bit hash methods
	uint64_t hash64_1(uint64_t key);
	uint64_t hash64_2(uint64_t key);
	//alignment NW
	int simpleNW(uint8_t *q, int q_len, uint8_t *t, int t_len);
	//simple NW algorithm with extension mode, return score; and return max_q and max_t as max extend length of ref and read
	int simpleNW_ext(uint8_t *q, int q_len, uint8_t *t, int t_len, int *max_q, int *max_t, bool is_forward);

	uint64_t get_random_seed_by_time();
	void simple_rand(uint64_t *rand_next);

#ifdef __cplusplus
}
#endif

#endif

//USAGE
#ifdef KSEQ_MAIN
#include <stdio.h>
int main()
{
	gzFile fp = gzopen("abcdef", "r");
	kseq_t *seqs = kseq_init(fp);
	while (kseq_read(seqs)>=0)
	{
		//do some thing for each sequence
	}
	kseq_destroy(seqs);
	return 0;
}
#endif

#ifdef KSTRING_MAIN
#include <stdio.h>
int main()
{
	kstring_t *s;
	s = (kstring_t*)calloc(1, sizeof(kstring_t));
	ksprintf(s, "abcdefg: %d", 100);
	printf("%s\n", s->s);
	free(s);
	return 0;
}
#endif
