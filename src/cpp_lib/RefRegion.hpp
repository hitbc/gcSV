/*
 * RefRegion.hpp
 *
 *  Created on: 2020-4-29
 *      Author: fenghe
 */

#ifndef SRC_SIGNAL_REFREGION_HPP_
#define SRC_SIGNAL_REFREGION_HPP_
extern "C"
{
	#include "../clib/utils.h"
	#include "../clib/bam_file.h" //used to define R_region
}
#include<iostream>
#include<string>

struct RefRegion{//break point candidate
	RefRegion(){
		chr_ID = -1;
		st_pos = -1;
		ed_pos = -1;
	}

	RefRegion(uint16_t	chr_ID_, int st_pos_, int	ed_pos_)
		:chr_ID(chr_ID_), st_pos(st_pos_), ed_pos(ed_pos_){}

	void set(uint16_t	chr_ID_, int st_pos_, int	ed_pos_){
		chr_ID = (chr_ID_); st_pos = (st_pos_); ed_pos = (ed_pos_);
	}

	void copy(RefRegion &from){
		chr_ID = (from.chr_ID); st_pos = (from.st_pos); ed_pos = (from.ed_pos);
	}

	void set_empty(){ chr_ID = -1; st_pos = 0; ed_pos = 0;}
	bool is_empty(){ return (chr_ID == -1); }

	void inline Combine(RefRegion &R, bool isUnion){if(isUnion) Union(R); else Intersection(R);}
	bool region_overlap(RefRegion &R){
		if(chr_ID != R.chr_ID)
			return false;
		if(ed_pos < R.st_pos || R.ed_pos < st_pos)
			return false;
		return true;
	}
	bool pos_within_region(int pos){
		if(pos >= st_pos && pos <= ed_pos)
			return true;
		return false;
	}
	bool pos_within_same_chr_region(int chr, int pos){
		if(chr == chr_ID && pos >= st_pos && pos <= ed_pos)
			return true;
		return false;
	}
	bool region_overlap(int st_pos_, int ed_pos_){
		if(ed_pos < st_pos_ || ed_pos_ < st_pos)
			return false;
		return true;
	}

	int region_overlap_length(RefRegion &R){
		if(chr_ID != R.chr_ID)
			return 0;
		if(ed_pos < R.st_pos || R.ed_pos < st_pos)
			return 0;
		int overlap1 = ed_pos - R.st_pos + 1, overlap2 = R.ed_pos - st_pos + 1, overlap = MIN(overlap1, overlap2);
		int len1 = getLen(), len2 = R.getLen(), len = MAX(len1, len2);
		return overlap;
	}

	float region_overlap_percent(RefRegion &R){
		int overlap = region_overlap_length(R);
		int len1 = getLen(), len2 = R.getLen(), len = MAX(len1, len2);
		return (float)overlap/len;
	}
	//return true when this region is within the region R
	bool Within(RefRegion &R){
		if(chr_ID == R.chr_ID && ed_pos <= R.ed_pos && st_pos >= R.st_pos)
			return true;
		return false;
	}
	int inline getMiddle(){return (st_pos+ed_pos)/2;}
	bool isForwardRegion(RefRegion &R){
		if(chr_ID == R.chr_ID)
			return (st_pos < R.st_pos)?true:false;
		else
			return (chr_ID < R.chr_ID)?true:false;
	}
	friend std::ostream& operator<<(std::ostream& os, const RefRegion& r){
	  os << "[" << r.chr_ID << ":" << r.st_pos << "-" << r.ed_pos <<  ", len: " << r.ed_pos - r.st_pos << "]";
	  return os;
	}

	void show(){
		fprintf(stderr, "[%d:%d-%d,len %d]\n",chr_ID , st_pos , ed_pos , ed_pos - st_pos );
	}

	//std::string toString(bam_hdr_t* header){ return "" + header->target_name[chr_ID] + ":" + st_pos + "-" + ed_pos;	}
	void toR_region(R_region &region) const{
		 region.chr_ID = chr_ID; region.st_pos = st_pos; region.ed_pos = ed_pos;
	}

	void store_from_R_region(const R_region &region){
		 chr_ID = region.chr_ID; st_pos = region.st_pos; ed_pos = region.ed_pos;
	}

	inline int getLen() {return ed_pos - st_pos + 1;}
	static inline int cmp_by_pos(const RefRegion &a, const RefRegion &b){
		if(a.chr_ID != b.chr_ID)
			return a.chr_ID < b.chr_ID;
		else
			return a.st_pos < b.st_pos;
	}
	void print(FILE* o){
		fprintf(o, "[ %d:%d-%d ]:\t", chr_ID, st_pos, ed_pos);
	}

	bool same_as(RefRegion &R){
		return (chr_ID == R.chr_ID && ed_pos == R.ed_pos && st_pos == R.st_pos);
	}

	int	chr_ID;
	int st_pos;
	int ed_pos;

private:
	void inline Union(RefRegion &R){
		if(R.chr_ID != chr_ID) std::cerr << this << "and" << R << "Chr_ID not same." << std::endl;
		st_pos = MIN(st_pos, R.st_pos);ed_pos = MAX(ed_pos, R.ed_pos);
	}
	void inline Intersection(RefRegion &R)	{
		if(R.chr_ID != chr_ID) std::cerr << this << "and" << R << "Chr_ID not same." << std::endl;
		st_pos = MAX(st_pos, R.st_pos);ed_pos = MIN(ed_pos, R.ed_pos);
	}
};

#endif /* SRC_SIGNAL_REFREGION_HPP_ */
