/*
 * forceCalling.cpp
 *
 *  Created on: 2025年1月4日
 *      Author: fenghe
 */

extern "C"
{
#include "../clib/utils.h"
}

void showGT(int GT){
	switch (GT) {
		case 0:	fprintf(stderr, "./.");	break;
		case 1:	fprintf(stderr, "0/0");	break;
		case 2:	fprintf(stderr, "0/1");	break;
		case 3:	fprintf(stderr, "1/1");	break;
		default:
			break;
	}
}
