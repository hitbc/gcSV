#define _GNU_SOURCE
#define FSYNC_ON_FLUSH

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>
#ifdef FSYNC_ON_FLUSH
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <sys/resource.h>
#include "utils.h"
#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))

/********************
 * System utilities *
 ********************/

bool fexist(const char *filename)
{
	if(access(filename, 0) == 0)
		return true;
	else
		return false;
}

void fexist_check(const char *filename){
	FILE* f = xopen(filename,"r");
	fclose(f);
}

FILE *err_xopen_core(const char *func, const int line, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		err_fatal(func, "line:[%d] fail to open file '%s' : %s", line, fn, strerror(errno));
	}
	return fp;
}

gzFile err_xzopen_core(const char *func,  const int line, const char *fn, const char *mode)
{
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		fprintf(stderr, "using stdin as input\n");
		/* According to zlib.h, this is the only reason gzdopen can fail */
		if (!fp) err_fatal(func, "Out of memory");
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		err_fatal(func, "line:[%d] fail to open file '%s' : %s", line, fn, errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp)
{
	if (freopen(fn, mode, fp) == 0) {
		err_fatal(func, "fail to open file '%s' : %s", fn, strerror(errno));
	}
	return fp;
}

void xmkdir(const char *dirPath){
	char mkdir[200] = "mkdir ";
	strcat(mkdir, dirPath);
	if(0 == system(mkdir)) fprintf(stderr,"\nmkdir error!\n");
}

void xrm(const char *fileName){
	char remove[200] = "rm ";
	strcat(remove, fileName);
	int rst = system(remove);
	if(rst < 0)
		fprintf(stderr, "remove file error[%s]\n", fileName);
}

//open file: dirPath/package_name.postfix, retuen the file handle
FILE* get_file_(const char *dirPath, const char *package_name, const char *postfix, const char *opentype){
	kstring_t filePath = {0};
	kputs(dirPath, &filePath);

	if (filePath.s[filePath.l - 1] != '/')
		kputc('/', &filePath);
	kputs(package_name, &filePath);
	kputs(postfix, &filePath);

	FILE * file = xopen(filePath.s, opentype); free(filePath.s);
	return file;
}

inline void* err_xcalloc_core(const char *func, const int line, size_t count,size_t eltsize)
{
	void* p = calloc(count,eltsize);
	if(NULL == p)
		err_fatal(func, "Out of memory in line '%d'", func, line);
	return p;
}

inline void* err_xmalloc_core(const char *func, const int line, size_t size)
{
	void *p = malloc(size);
	if(NULL == p)
		err_fatal(func, "Out of memory in line '%d'", func, line);
	return p;
}

inline void *err_xrealloc_core(const char *func, const int line, void *ptr, size_t newsize)
{
	void *p = realloc(ptr, newsize);
	if(NULL == p)
		err_fatal(func, "Out of memory in line '%d'", func, line);
	return p;
}

/*Abandoned
void xfree(void **p)
{
	free(*p);
	*p = NULL;
}
*/

void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	va_end(args);
	exit(EXIT_FAILURE);
}

void err_fatal_core(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}

void _err_fatal_simple(const char *func, const char *msg)
{
	fprintf(stderr, "[%s] %s\n", func, msg);
	exit(EXIT_FAILURE);
}

//set msgN tobe 0 when no messages
//at most 3 additional messages
//fprintf(stderr, msg, m1, m2, m3)
void _err_fatal_core(const char *func, const int line, const char *msg, const char *m1, const char *m2, const char *m3)
{
	fprintf(stderr, "Function [%s] Abort, line [%d], message:[", func, line);
	if(m1 == NULL)
		fprintf(stderr, msg, "");
	else if(m2 == NULL)		fprintf(stderr, msg, m1);
	else if(m3 == NULL)		fprintf(stderr, msg, m1, m2);
	else					fprintf(stderr, msg, m1, m2, m3);
	fprintf(stderr,"]\n");
	abort();
}

void _err_fatal_simple_core(const char *func, const int line, const char *msg)
{
	fprintf(stderr, "[%s] %s Abort, line [%d]!\n", func, msg, line);
	abort();
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb) 
		_err_fatal_simple("fwrite", strerror(errno));
	return ret;
}

size_t err_fread_noeof(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
	{
		_err_fatal_simple("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
	}
	return ret;
}

int err_gzread(gzFile file, void *ptr, unsigned int len)
{
	int ret = gzread(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

int err_gzwrite(gzFile file, void *ptr, unsigned int len)
{

	int ret = gzwrite(file, ptr, len);

	if (ret < 0)
	{
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		_err_fatal_simple("gzread", Z_ERRNO == errnum ? strerror(errno) : msg);
	}

	return ret;
}

//SEEK_SET, SEEK_CUR, or SEEK_END
int err_fseek(FILE *stream, long offset, int whence)
{
	int ret = fseek(stream, offset, whence);
	if (0 != ret)
	{
		_err_fatal_simple("fseek", strerror(errno));
	}
	return ret;
}

long err_ftell(FILE *stream)
{
	long ret = ftell(stream);
	if (-1 == ret)
	{
		_err_fatal_simple("ftell", strerror(errno));
	}
	return ret;
}

long err_fsize(FILE *stream)
{
	long old = err_ftell(stream);
	err_fseek(stream, -1, SEEK_END);
	long rst = err_ftell(stream) + 1;
	//reset
	err_fseek(stream, old, SEEK_SET);
	return rst;
}

int err_printf(const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stdout, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf(stdout)", strerror(saveErrno));
	return done;
}

int err_fprintf(FILE *stream, const char *format, ...) 
{
	va_list arg;
	int done;
	va_start(arg, format);
	done = vfprintf(stream, format, arg);
	int saveErrno = errno;
	va_end(arg);
	if (done < 0) _err_fatal_simple("vfprintf", strerror(saveErrno));
	return done;
}

int err_fputc(int c, FILE *stream)
{
	int ret = putc(c, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputc", strerror(errno));
	}

	return ret;
}

int err_fputs(const char *s, FILE *stream)
{
	int ret = fputs(s, stream);
	if (EOF == ret)
	{
		_err_fatal_simple("fputs", strerror(errno));
	}

	return ret;
}

int err_puts(const char *s)
{
	int ret = puts(s);
	if (EOF == ret)
	{
		_err_fatal_simple("puts", strerror(errno));
	}

	return ret;
}

int err_fflush(FILE *stream) 
{
    int ret = fflush(stream);
    if (ret != 0) _err_fatal_simple("fflush", strerror(errno));

#ifdef FSYNC_ON_FLUSH
	/* Calling fflush() ensures that all the data has made it to the
	   kernel buffers, but this may not be sufficient for remote filesystems
	   (e.g. NFS, lustre) as an error may still occur while the kernel
	   is copying the buffered data to the file server.  To be sure of
	   catching these errors, we need to call fsync() on the file
	   descriptor, but only if it is a regular file.  */
	{
		struct stat sbuf;
		if (0 != fstat(fileno(stream), &sbuf))
			_err_fatal_simple("fstat", strerror(errno));
		
		if (S_ISREG(sbuf.st_mode))
		{
			if (0 != fsync(fileno(stream)))
				_err_fatal_simple("fsync", strerror(errno));
		}
	}
#endif
    return ret;
}

int err_fclose(FILE *stream) 
{
	int ret = fclose(stream);
	if (ret != 0) _err_fatal_simple("fclose", strerror(errno));
	return ret;
}

int err_gzclose(gzFile file)
{
	int ret = gzclose(file);
	if (Z_OK != ret)
	{
		_err_fatal_simple("gzclose", Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}

	return ret;
}

//return 1 if big endian, else return 0
int big_or_small_endian(){
	int nNum = 0x12345678;
	char chData = *(char*)(&nNum);
	if (chData == 0x12)	return 1;
	else				return 0;
}

/*********
 * Timer *
 *********/

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double realduration(struct timeval start)
{
	struct timeval end;
	gettimeofday(&end, NULL);
	end.tv_usec -= start.tv_usec;
	end.tv_sec 	-= start.tv_sec;
	if (end.tv_usec < 0) {
		end.tv_sec--;
		end.tv_usec += 1000000;
	}
	double seconds = end.tv_usec;
	seconds /= 1e6;
	seconds += end.tv_sec;
	return seconds;
}

//return peak memory, in Kbp
long peak_memory(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss;
}

/**********
 * Sort *
 **********/

void ksort_stable_step(void *base, size_t nmemb, size_t size, __compar_fn_t cmp, size_t step)
{
	size_t total_size = nmemb*size;
	char small_data_buff[1000];
	void *base_ = (total_size <= 1000)?small_data_buff:xmalloc(total_size);
	void *from = base, *to = base_;
	while (step < total_size)
	{
		size_t two_times_step = step << 1;
		size_t bgn1;
		for(bgn1 = 0; bgn1 + step < total_size; bgn1 += two_times_step)
		{
			size_t end1 = bgn1 + step;
			size_t bgn2 = end1;
			size_t end2 = MIN(bgn2 + step, total_size);
			{//merge two
				size_t p1 = bgn1, p2 = bgn2, p = bgn1;
				while (p1 < end1 && p2 < end2)
				{
					int cmp_rst = cmp(from + p1, from + p2);
					if(cmp_rst >= 0)
					{//mem_cpy:m1 + p1 -> base + p, size:size;
						char * cp_to = to + p, * cp_from = from + p1;
						size_t width = size;
						p  += size;
						p1 += size;
						while(width--)
							*cp_to++ = *cp_from++;
					}
					else
					{//mem_cpy:m2 + p2 -> base + p, size:size;
						char * cp_to = to + p, *cp_from = from + p2;
						size_t width = size;
						p += size;
						p2 += size;
						while(width--)
							*cp_to++ = *cp_from++;
					}
				}
				if(p1 < end1)
				{//mem_cpy:m1 + p1 -> base + p;
					char * cp_to = to + p, *cp_from = from + p1;
					size_t width = (end1 - p1);
					while(width--)
						*cp_to++ = *cp_from++;
				}
				else
				{//p2 < len2 //mem_cpy:m1 + p1 -> base + p;
					char * cp_to = to + p, *cp_from = from + p2;
					size_t width = (end2 - p2);
					while(width--)
						*cp_to++ = *cp_from++;
				}
			}
		}
		//copy the last part
		if(bgn1 < total_size)
		{
			char * cp_to = to + bgn1, *cp_from = from + bgn1;
			size_t width = (total_size - bgn1);
			while(width--)
				*cp_to++ = *cp_from++;
		}
		step <<= 1;
		void * swap = from;
		from = to, to = swap;
	}
	if(from != base)
		memcpy(base, base_, total_size);
	if(total_size > 1000)
		free(base_);
}

void ksort_stable(void *base, size_t nmemb, size_t size, __compar_fn_t cmp)
{
	ksort_stable_step(base, nmemb, size, cmp, size);
}

typedef struct
{
	void *base;
	size_t nmemb;
	size_t size;
	__compar_fn_t cmp;
}STABLE_SORT_DATA;

static void *ksort_stable_worker(void *data)
{
	STABLE_SORT_DATA *d = data;
	ksort_stable_step(d->base, d->nmemb, d->size, d->cmp, d->size);
	return NULL;
}

void ksort_stable_mt(void *base, size_t nmemb, size_t size, __compar_fn_t cmp, int n_thread)
{
	//distributing works
	pthread_t *tid = xcalloc_t(pthread_t, n_thread);
	STABLE_SORT_DATA *data = xcalloc_t(STABLE_SORT_DATA, n_thread);
	size_t nmemb_t = nmemb / n_thread;
	for (int i = 0; i < n_thread; ++i)
	{
		data[i].base = base + i * nmemb_t * size;
		data[i].nmemb = nmemb_t;
		data[i].size = size;
		data[i].cmp = cmp;
	}
	data[n_thread - 1].nmemb += nmemb - nmemb_t * n_thread;
	for (int i = 0; i < n_thread; ++i)
		pthread_create(&tid[i], 0, ksort_stable_worker, &data[i]);
	for (int i = 0; i < n_thread; ++i)
		pthread_join(tid[i], 0);
	//combine result
	ksort_stable_step(base, nmemb, size, cmp, nmemb_t * size);
}

/**
 * Stable bucket sort, using for sorting extremely large data-set,
 * use no additional memory, please use main function :radix_sort.
 * when use, you should define:
 * rstype_t and rskey as data type and function to get key and sizeof key
 *
 * Example1: when sort int list:
 * #define rstype_t int
 * #define rskey(a) (a)
 * #define sizeof_key sizeof(int)
 *
 * Example2:when sort struct DATA{uint64_t key, int data}:
 * #define rstype_t DATA
 * #define rskey(a) (a->key)
 * #define sizeof_key sizeof(uint64_t)
 *  */
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8
#define rstype_t int
#define rskey(a) (a)
#define bit_sizeof_key 32//for int

typedef struct
{
	rstype_t *b, *e;
}rsbucket;

//rskey(*i) < rskey(*(i - 1))::cmp function
//stable insert sort
void __rs_insertsort(rstype_t *beg, rstype_t *end)
{
	rstype_t *i;
	for (i = beg + 1; i < end; ++i)//from begin to end
	{
		if (rskey(*i) < rskey(*(i - 1)))//compare function
		{
			rstype_t *j, tmp = *i;
			for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j)//insert position
				*j = *(j - 1);
			*j = tmp;
		}
	}
}

//bucket hash function: (key>>s)&0xff
//s:s should be the length of key - 8(as max as it can)
void __rs_sort(rstype_t *beg, rstype_t *end, int n_bits, int s)
{
	rstype_t *i;
	int size = 1 << n_bits;//256, bucket number
	int m = size - 1;//0xff, mask
	rsbucket *k;
	rsbucket b[1<<RS_MAX_BITS];//256 bucket
	rsbucket *be = b + size;//bucket end
	for (k = b; k != be; ++k)//from begin to end
		k->b = k->e = beg;
	for (i = beg; i != end; ++i)
		++b[rskey(*i)>>s&m].e;//count size of each bucket
	for (k = b + 1; k != be; ++k)
	{
		k->e += (k-1)->e - beg;//reset begin and end of each bucket
		k->b = (k-1)->e;
	}
	///bucket distribution
	for (k = b; k != be;)//for each bucket
	{
		if (k->b != k->e)
		{
			int hash = rskey(*k->b)>>s&m;//get one data in this bucket
			rsbucket *l = b + hash;
			if (l != k)//data not belong to this bucket
			{
				rstype_t tmp = *k->b, swap;
				do
				{
					swap = tmp;//exchange pointer
					tmp = *l->b;
					*l->b++ = swap;//store data in right position
					l = b + (rskey(tmp)>>s&m);//get new l
				} while (l != k);
				*k->b++ = tmp;//when data belong to this bucket, store
			}
			else//data belong to this bucket
				++k->b;
		}
		else//next bucket
			++k;
	}
	//reset bucket pointer
	for (b->b = beg, k = b + 1; k != be; ++k)
		k->b = (k-1)->e;
	//when number of data in one bucket is over 1, sort each bucket
	if (s)
	{
		s = s > n_bits? s - n_bits : 0;//reset s
		for (k = b; k != be; ++k)
			if (k->e - k->b > RS_MIN_SIZE && s != 0)
				__rs_sort(k->b, k->e, n_bits, s);
			else
				if (k->e - k->b > 1)
					__rs_insertsort(k->b, k->e);
	}
}

//bgn and end is the bgn and end pointer of
void radix_sort(rstype_t *beg, rstype_t *end)
{
	if (end - beg <= RS_MIN_SIZE)
		__rs_insertsort(beg, end);
	else
		__rs_sort(beg, end, RS_MAX_BITS, bit_sizeof_key - RS_MAX_BITS);
}

/**********
 * String *
 **********/

//to = pre + suf; return to
char* strcmb(char* to, char *pre, char * suf)
{
	strcpy(to, pre);
	strcat(to, suf);
	return to;
}

bool hasEnding(char * fullString, char * ending)
{
	int full_len = strlen(fullString);
	int end_len = strlen(ending);
	if (full_len < end_len)
		return false;
	if(strcmp(fullString + full_len - end_len, ending) == 0)
		return true;
	return false;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

//resize a string
inline void ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		s->m = size;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

//enlarge a string
void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m)
	{
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

//store a string, == memcat
int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

//store a string
inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

//store a char
inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

//store an int
inline int kputw(int c, kstring_t *s)
{
	char buf[16];
	int l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

//store an unsigned-int
inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

//store a long int
inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	long l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

void kstrcpy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

//equal sprintf, without float
void sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p)
	{
		if (*p == '%')
		{
			if (p > q)
				kstrcpy(s, q, p);
			++p;
			if (*p == 'd')
			{
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do
				{
					buf[l++] = x%10 + '0'; x /= 10;
				}
				while (x > 0);
				if (c < 0)
					buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i)
					s->s[s->l++] = buf[i];
			}
			else if (*p == 'u')
			{
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do
				{
					buf[l++] = x%10 + '0'; x /= 10;
				}
				while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i)
					s->s[s->l++] = buf[i];
			}
			else if (*p == 's')
			{
				char *r = va_arg(ap, char*);
				kstrcpy(s, r, r + strlen(r));
			}
			else if (*p == 'c')
			{
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			}
			else
				abort();
			q = p + 1;
		}
	}
	if (p > q)
		kstrcpy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

//for kseq

#define __kseq_bufsize 100000

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_MAX   1

static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append)
{
	if (dret) *dret = 0;
	str->l = append? str->l : 0;
	if (ks->begin >= ks->end && ks->is_eof) return -1;
	for (;;) {
		int i;
		if (ks->begin >= ks->end) {
			if (!ks->is_eof) {
				ks->begin = 0;
				ks->end = gzread(ks->f, ks->buf, __kseq_bufsize);
				if (ks->end < __kseq_bufsize) ks->is_eof = 1;
				if (ks->end == 0) break;
			} else break;
		}
		if (delimiter > KS_SEP_MAX) {
			for (i = ks->begin; i < ks->end; ++i)
				if (ks->buf[i] == delimiter) break;
		} else if (delimiter == KS_SEP_SPACE) {
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i])) break;
		} else if (delimiter == KS_SEP_TAB) {
			for (i = ks->begin; i < ks->end; ++i)
				if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break;
		} else i = 0; /* never come to here! */
		if (str->m - str->l < (size_t)i - ks->begin + 1) {
			str->m = str->l + (i - ks->begin) + 1;
			kroundup32(str->m);
			str->s = (char*)realloc(str->s, str->m);
		}
		memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
		str->l = str->l + (i - ks->begin);
		ks->begin = i + 1;
		if (i < ks->end) {
			if (dret) *dret = ks->buf[i];
			break;
		}
	}
	if (str->s == 0) {
		str->m = 1;
		str->s = (char*)calloc(1, 1);
	}
	str->s[str->l] = '\0';
	return str->l;
}

static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret)
{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

static inline int ks_getc(kstream_t *ks)
{
	if (ks->is_eof && ks->begin >= ks->end) return -1;
	if (ks->begin >= ks->end) {
		ks->begin = 0;
		ks->end = gzread(ks->f, ks->buf, __kseq_bufsize);
		if (ks->end < __kseq_bufsize) ks->is_eof = 1;
		if (ks->end == 0) return -1;
	}
	return (int)ks->buf[ks->begin++];
}

kstream_t *ks_init(gzFile f)
{
	kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
	ks->f = f;
	ks->buf = (unsigned char*)malloc(__kseq_bufsize);
	return ks;
}

void ks_destroy(kstream_t *ks)
{
	if (ks) {
		free(ks->buf);
		free(ks);
	}
}

kseq_t *kseq_init(gzFile fd)
{
	kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));
	s->f = ks_init(fd);
	return s;
}

inline void kseq_rewind(kseq_t *ks)
{
	ks->last_char = 0;
	ks->f->is_eof = ks->f->begin = ks->f->end = 0;
}

void kseq_destroy(kseq_t *ks)
{
	if (!ks) return;
	free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s);
	ks_destroy(ks->f);
	free(ks);
}

int kseq_read(kseq_t *seq)
{
	int c;
	kstream_t *ks = seq->f;
	if (seq->last_char == 0) { /* then jump to the next header line */
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
		if (c == -1) return -1; /* end of file */
		seq->last_char = c;
	} /* else: the first header char has been read in the previous call */
	seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */
	if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* normal exit: EOF */
	if (c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0); /* read FASTA/Q comment */
	if (seq->seq.s == NULL) { /* we can do this in the loop below, but that is slower */
		seq->seq.m = 256;
		seq->seq.s = (char*)malloc(seq->seq.m);
	}
	while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
		seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */
		ks_getuntil2(ks, '\n', &seq->seq, 0, 1); /* read the rest of the line */
	}
	if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
	if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */
		seq->seq.m = seq->seq.l + 2;
		kroundup32(seq->seq.m); /* rounded to the next closest 2^k */
		seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);
	}
	seq->seq.s[seq->seq.l] = 0;	/* null terminated string */
	if (c != '+') return seq->seq.l; /* FASTA */
	if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */
		seq->qual.m = seq->seq.m;
		seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);
	}
	while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
	if (c == -1) return -2; /* error: no quality string */
	while (ks_getuntil2(ks, '\n', &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l);
	seq->last_char = 0;	/* we have not come to the next header line */
	if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */
	return seq->seq.l;
}

uint64_t kmerMask_bit[65] = {
		0x0000000000000000,		0x0000000000000001,		0x0000000000000003,		0x0000000000000007,
		0x000000000000000f,		0x000000000000001f,		0x000000000000003f,		0x000000000000007f,
		0x00000000000000ff,		0x00000000000001ff,		0x00000000000003ff,		0x00000000000007ff,
		0x0000000000000fff,		0x0000000000001fff,		0x0000000000003fff,		0x0000000000007fff,
		0x000000000000ffff,		0x000000000001ffff,		0x000000000003ffff,		0x000000000007ffff,
		0x00000000000fffff,		0x00000000001fffff,		0x00000000003fffff,		0x00000000007fffff,
		0x0000000000ffffff,		0x0000000001ffffff,		0x0000000003ffffff,		0x0000000007ffffff,
		0x000000000fffffff,		0x000000001fffffff,		0x000000003fffffff,		0x000000007fffffff,
		0x00000000ffffffff,		0x00000001ffffffff,		0x00000003ffffffff,		0x00000007ffffffff,
		0x0000000fffffffff,		0x0000001fffffffff,		0x0000003fffffffff,		0x0000007fffffffff,
		0x000000ffffffffff,		0x000001ffffffffff,		0x000003ffffffffff,		0x000007ffffffffff,
		0x00000fffffffffff,		0x00001fffffffffff,		0x00003fffffffffff,		0x00007fffffffffff,
		0x0000ffffffffffff,		0x0001ffffffffffff,		0x0003ffffffffffff,		0x0007ffffffffffff,
		0x000fffffffffffff,		0x001fffffffffffff,		0x003fffffffffffff,		0x007fffffffffffff,
		0x00ffffffffffffff,		0x01ffffffffffffff,		0x03ffffffffffffff,		0x07ffffffffffffff,
		0x0fffffffffffffff,		0x1fffffffffffffff,		0x3fffffffffffffff,		0x7fffffffffffffff,
		0xffffffffffffffff,
};

//kmer values and functions
uint64_t kmerMask[33] = {
	0x0000000000000000,	0x0000000000000003,
	0x000000000000000f,	0x000000000000003f,
	0x00000000000000ff,	0x00000000000003ff,
	0x0000000000000fff,	0x0000000000003fff,
	0x000000000000ffff,	0x000000000003ffff,
	0x00000000000fffff,	0x00000000003fffff,
	0x0000000000ffffff,	0x0000000003ffffff,
	0x000000000fffffff,	0x000000003fffffff,
	0x00000000ffffffff,	0x00000003ffffffff,
	0x0000000fffffffff,	0x0000003fffffffff,
	0x000000ffffffffff,	0x000003ffffffffff,
	0x00000fffffffffff,	0x00003fffffffffff,
	0x0000ffffffffffff,	0x0003ffffffffffff,
	0x000fffffffffffff,	0x003fffffffffffff,
	0x00ffffffffffffff,	0x03ffffffffffffff,
	0x0fffffffffffffff,	0x3fffffffffffffff,
	0xffffffffffffffff,
};

uint64_t inline char2Kmer(char *str_kmer, uint8_t L_kmer, uint8_t *Bit){
	xassert(L_kmer <= 32 ,"Length of kmer should be no more than 32 base pair.");
	uint64_t value = 0;
	for(uint8_t i=0; i<L_kmer; ++i)
		value = (value<<2)|Bit[(int8_t)str_kmer[i]];
	return value;
}

uint64_t inline char2KmerRC(char *str_kmer, uint8_t L_kmer, uint8_t *Bit){
	xassert(L_kmer <= 32 ,"Length of kmer should be no more than 32 base pair.");
	uint64_t value = 0;
	for (uint8_t i=0; i<L_kmer;++i)
		value = (value<<2)|(Bit[(int8_t)str_kmer[L_kmer - i - 1]]^0x3);
	return value;
}

uint64_t inline getRcKmer(uint64_t kmer, uint8_t l_kmer)
{
	uint64_t rckv = -1;
	uint64_t mask = 0x3;
	for (uint8_t i = 0; i < l_kmer; ++i)
	{
		rckv = (rckv << 2) | (kmer & mask);
		kmer >>= 2;
	}
	return (~rckv);
}

uint64_t inline binchar2Kmer(uint8_t *s, uint8_t L_kmer){
	xassert(L_kmer <= 32 ,"Length of kmer should be no more than 32 base pair.");
	xassert(*s <= 4, "data should be bin char, value[0~3].");
	uint64_t kmer = 0;
	for(int i = 0; i < L_kmer; i++)
		kmer  = (kmer << 2) | s[i];
	return kmer;
}

uint64_t inline binchar2KmerRC(uint8_t *s, uint8_t L_kmer){
	xassert(L_kmer <= 32 ,"Length of kmer should be no more than 32 base pair.");
	xassert(*s < 4, "data should be bin char, value[0~3].");
	uint64_t kmer = 0;
	for(int i = 0; i >= 1 - L_kmer; --i)
		kmer  = (kmer << 2) | (3 - s[i]);
	return kmer;
}

//hash method 1
inline uint64_t hash64_1(uint64_t key)
{
	key = (~key + (key << 21));
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8));
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4));
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}

//hash method 2
inline uint64_t hash64_2(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

#define SNW_MATCH 1
#define SNW_MISMATCH -3
#define SNW_GAP -3
//implement a traditional NW algorithm, with score per base 1
int simpleNW(uint8_t *q, int q_len, uint8_t *t, int t_len)
{
	int M[q_len + 1][t_len + 1];
    M[0][0] = 0;	//init matrix
    for (int i = 1; i < q_len +  1; i++)
        M[i][0] = M[i-1][0] + (SNW_GAP);
    for (int j = 1; j < t_len + 1; j++)
        M[0][j] = M[0][j-1] + (SNW_GAP);
	//alignment
    for (int i = 1; i < q_len + 1; i++)
        for (int j = 1; j < t_len + 1; j++)
        {
            int scoreDiag = M[i-1][j-1] + ((q[i-1] == t[j-1])?SNW_MATCH:SNW_MISMATCH);
            int scoreLeft_up = MAX(M[i][j-1], M[i-1][j]) + SNW_GAP;
            M[i][j] = MAX(scoreDiag, scoreLeft_up);
        }
    return M[q_len][t_len];
}

//simple NW algorithm with extension mode
//forward == true:  from pos 0 to pos [len - 1]
//forward == false: from pos 0 to pos [- len + 1]
int simpleNW_ext(uint8_t *q, int q_len, uint8_t *t, int t_len, int *max_q, int *max_t, bool is_forward)
{
	int max_q_pos = 0; int max_t_pos = 0; int max_s = 0;
	int M[q_len + 1][t_len + 1];
    M[0][0] = 0;	//init matrix
    for (int i = 1; i < q_len +  1; i++)
        M[i][0] = M[i-1][0] + (SNW_GAP);
    for (int j = 1; j < t_len + 1; j++)
        M[0][j] = M[0][j-1] + (SNW_GAP);
	//alignment
    for (int i = 1; i < q_len + 1; i++)
        for (int j = 1; j < t_len + 1; j++)
        {
        	bool is_match = is_forward?(q[i-1] == t[j-1]):(q[-i+1] == t[-j+1]);
            int scoreDiag = M[i-1][j-1] + ((is_match)?SNW_MATCH:SNW_MISMATCH);
            int scoreLeft_up = MAX(M[i][j-1], M[i-1][j]) + SNW_GAP;
            M[i][j] = MAX(scoreDiag, scoreLeft_up);
            if(max_s < M[i][j])
            {
            	max_s = M[i][j]; max_q_pos = i; max_t_pos = j;
            }
        }
    *max_q = max_q_pos; *max_t = max_t_pos;
    return max_s;
}


uint64_t get_random_seed_by_time()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec * 1e6 + tp.tv_usec;
}

inline void simple_rand(uint64_t *rand_next) {
	*rand_next = ((*rand_next)>>2) * 1103515245 + 48561;
}
