/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:
#include "kvec.h"
int main() {
	kvec_t(int) array;
	kv_init(array);
	kv_push(int, array, 10); // append
	kv_a(int, array, 20) = 5; // dynamic
	kv_A(array, 20) = 4; // static
	kv_destroy(array);
	return 0;
}
*/

/*
  2008-09-22 (0.1.0):
	* The initial version.
*/

#ifndef AC_KVEC_H
#define AC_KVEC_H

#include <stdlib.h>

#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#define kvec_t(type) struct { size_t n, m; type *a; }
#define kvec_T(type, type_v_NAME) typedef struct { size_t n, m; type *a; } type_v_NAME;
#define kv_init(v) ((v).n = (v).m = 0, (v).a = 0)
#define kv_destroy(v) free((v).a)
#define kv_A(v, i) ((v).a[(i)])
#define kv_pop(v) ((v).a[--(v).n])
#define kv_size(v) ((v).n)
#define kv_max(v) ((v).m)

//type is the type of node; v is a vector (not a pointer); and s is an int
#define kv_resize(type, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			kv_roundup32((v).m); \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_copy(type, v1, v0) do { \
		if ((v1).m < (v0).n) kv_resize(type, (v1), (v0).n); \
		(v1).n = (v0).n; \
		memcpy((v1).a, (v0).a, sizeof(type) * (v0).n); \
	} while (0) \

#define kv_push(type, v, x) do { \
		if ((v).n == (v).m) { \
			(v).m = (v).m? (v).m<<1 : 10; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		(v).a[(v).n++] = (x); \
	} while (0 )

#define kv_push_2(type, v, x) do { \
		if ((v)->n == (v)->m) { \
			(v)->m = (v)->m? (v)->m<<1 : 10; \
			(v)->a = (type*)realloc((v)->a, sizeof(type) * (v)->m); \
		} \
		(v)->a[(v)->n++] = (x); \
	} while (0)

#define kv_pushp(type, v, p) do { \
		if ((v).n == (v).m) { \
			(v).m = (v).m? (v).m<<1 : 10; \
			(v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
		} \
		*(p) = &(v).a[(v).n++]; \
	} while (0)

//v is a pointer, while p is the original pointer
#define kv_pushp_2(type, v, p) do { \
		if ((v)->n == (v)->m) { \
			(v)->m = (v)->m? (v)->m<<1 : 10; \
			(v)->a = (type*)realloc((v)->a, sizeof(type) * (v)->m); \
		} \
		(p) = (v)->a + (v)->n++; \
	} while (0)

#define kv_reverse(type, v, start) do { \
		if ((v).m > 0 && (v).n > (start)) { \
			size_t __i, __end = (v).n - (start); \
			type *__a = (v).a + (start); \
			for (__i = 0; __i < __end>>1; ++__i) { \
				type __t = __a[__end - 1 - __i]; \
				__a[__end - 1 - __i] = __a[__i]; __a[__i] = __t; \
			} \
		} \
	} while (0)

//.---------------------------------define for hash table--------------------------------------//
//function: used to define a chain-format hash table; the table has 2^len
//			hash cells, the bottom of _len_ bases of _key_ will be used as
//			hash value.
// 1)type: type of data, and the type must implement at least two fields:
//		(1)	key, and it can be one of [uint64_t] or [uint32_t]
//		(2)	next, and it can be one of [uint64_t] or [uint32_t]
// 2)type_v_NAME : type of new hash table type
// 3)n: total number of data stored in hash table; m: total memory allocated for data
// 4)len: 	length of mask, the hash table will hash 2^_len blocks, 1M cells will be malloced
//			when _len == 20

// the lower [len] bit of the key will be used as hash key, if the lower bits of hash keys
//assembled together, use function to scatter them
#define re_hash_64(key) do{\
		uint64_4 __key__ = key;\
		__key__ = (~__key__ + (__key__ << 21));\
		__key__ = __key__ ^ __key__ >> 24;\
		__key__ = ((__key__ + (__key__ << 3)) + (__key__ << 8));\
		__key__ = __key__ ^ __key__ >> 14;\
		__key__ = ((__key__ + (__key__ << 2)) + (__key__ << 4));\
		__key__ = __key__ ^ __key__ >> 28;\
		__key__ = (__key__ + (__key__ << 31));\
		key = __key__;\
	} while (0)

#define hash_T(type, type_v_NAME) typedef struct { \
		size_t n, m; \
		type *a; \
		size_t block_size; \
		uint64_t MASK;\
}type_v_NAME;

//after define a hash table; you should init it: clean it and doing basic malloc
//v: a pointer to a hash table
#define hash_t_init(type, v, len) do{\
		(v)->MASK = ((0x1 << ((len))) - 1);\
		(v)->block_size = (0x1 << (len));\
		(v)->m = (v)->block_size + (v)->block_size;\
		(v)->n = (v)->block_size;\
		(v)->a = NULL; \
		(v)->a = (type*)realloc((v)->a, sizeof(type) * (v)->m); \
		for(uint64_t i = 0; i < (v)->block_size; i++) (v)->a[i].next = 0;\
	} while (0)

//clean a hash table, when you want to reuse a hash table to store new data.
#define hash_t_clean(v) do{\
		(v)->n = (v)->block_size;\
		for(uint64_t i = 0; i < (v)->block_size; i++) (v)->a[i].next = 0;\
	} while (0)

//free a hash table
#define hash_t_destroy(v) free((v)->a)

//store a new data to a hash table, when difference data have same key, they will be stored in difference place
//type: as defined in "hash_T":"type"
//v: a pointer to a hash table
//p: a pointer to a to be stored data
#define hash_t_push(type, v, p) do { \
		type *__p__ = (p);\
		if ((v)->n == (v)->m) { \
			(v)->m = (v)->m? (v)->m<<1 : 10; \
			(v)->a = (type*)realloc((v)->a, sizeof(type) * (v)->m); \
			 /*clean hash table*/\
			 for(uint64_t i = (v)->n; i < (v)->m; i++) (v)->a[i].next = 0; \
		} \
		uint64_t __index__ = (__p__->key & (v)->MASK), __new_index__;\
		while((__new_index__ = (v)->a[__index__].next) != 0){__index__ = __new_index__;}\
		(v)->a[(v)->n] = *__p__;\
		(v)->a[__index__].next = (v)->n;\
		(v)->a[(v)->n++].next = 0;\
	} while (0)

//search a data by key, when difference data have same key, only the first one will return
//return a pointer when results found, otherwise NULL
//key__: the key for search, it can be [uint64_t] or [uint32_t]
//rst: the pointer to a result
#define hash_t_search(v, key__, rst) do { \
		uint64_t __key__ = (key__);\
		uint64_t index = ((__key__) & (v)->MASK);\
		(rst) = NULL;\
		do{\
			if(__key__ ==  (v)->a[index].key){\
				(rst) = (v)->a + index;\
				break;\
			}\
		} while((index = (v)->a[index].next) != 0); \
	} while (0)

//search a data by key, but return the next data after the _pre data;
//return a pointer when results found, otherwise NULL
//_pre: the previous result search by [hash_t_search] or [hash_t_search_next]
#define hash_t_search_next(v, _pre, rst) do { \
		uint64_t __key__ = (_pre)->key;\
		uint64_t index = (_pre)->next;\
		if(index == 0) return NULL;\
		(rst) = NULL;\
		do{\
			if(__key__ ==  (v)->a[index].key){\
				(rst) = (v)->a + index;\
				break;\
			}\
		} while((index = (v)->a[index].next) != 0); \
	} while (0)

//#define EXAMPLE_OF_HASH_T
#ifdef EXAMPLE_OF_HASH_T

	//basic type
	#define HASH_T_LEN 10 // for 1K size hash table
	typedef struct MY_DATA_T{
			uint32_t key;//needed: for _len <= 32
			uint32_t next;//needed: for len <= 32
			uint32_t data;//custom: any type
	}MY_DATA_T;

	//define a hash table type
	hash_T(MY_DATA_T, MY_HASH_TABLE_T);
	//new and init a hash table
	MY_HASH_TABLE_T t;
	hash_t_init(MY_DATA_T, &t, HASH_T_LEN);
	//new some data and store： next part can be any thing
	MY_DATA_T d1 = {12, 0 , 1};
	MY_DATA_T d2 = {12 + 1024, 0 , 2};//with same hash value
	MY_DATA_T d3 = {14, 0 , 3};
	MY_DATA_T d4 = {12213323, 0 , 4};//with big key
	MY_DATA_T d5 = {12, 0 , 5};//with same key
	MY_DATA_T d6 = {12, 0 , 6};//with same key
	hash_t_push(MY_DATA_T,&t,&d1);
	hash_t_push(MY_DATA_T,&t,&d2);
	hash_t_push(MY_DATA_T,&t,&d3);
	hash_t_push(MY_DATA_T,&t,&d4);
	hash_t_push(MY_DATA_T,&t,&d5);
	hash_t_push(MY_DATA_T,&t,&d6);

	MY_DATA_T * rst = NULL;
	hash_t_search(&t, 12, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data); else printf("NULL\n");
	hash_t_search_next(&t, rst, rst);if(rst != NULL)  printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");
	hash_t_search_next(&t, rst, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");
	hash_t_search_next(&t, rst, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");
	hash_t_search(&t, 14, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");
	hash_t_search(&t, 12213323, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");
	hash_t_search(&t, 12 + 1024, rst); if(rst != NULL) printf("%d %d\n", rst->key, rst->data);else printf("NULL\n");

#endif




#endif
