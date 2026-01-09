/*
 * sve.hpp
 *
 *  Created on: 2020-4-28
 *      Author: fenghe
 */
#pragma once

extern "C"
{
	#include <clib/vcf_lib.h>
}

#include <vector>
#include<iostream>
#include "../cpp_lib/RefRegion.hpp"

namespace SIG{
	enum T					  { DR, SH, LEN, UNKNOWN, LRS};
	const std::string STR[] = {"DR", "SH", "LEN", "UNKNOWN", "TGS"};
}
namespace SV{
	enum T				      {  DEL,   DUP,   INS,   INV_1,   INV_2,   TRA,   TRA_INV, LEN, UNKNOWN};
	const std::string STR[] = { "DEL", "DUP", "INS", "INV_1", "INV_2", "TRA", "TRA_INV", "LEN", "UNKNOWN"};
}

namespace RST{
	enum T					  {  BEGIN,   END,   SOLID, UNKNOWN};
	const std::string STR[] = { "BEGIN", "END", "SOLID", "UNKNOWN"};
}

#define SAMPLE_UNKNOWN MAX_uint32_t

class SVE_core{
public:
	RefRegion 	r1;
	RefRegion 	r2;
	static inline int cmp_by_position(const SVE_core &a, const SVE_core &b){
		if(a.r1.chr_ID == b.r1.chr_ID)	return a.r1.st_pos < b.r1.st_pos;
		else							return a.r1.chr_ID < b.r1.chr_ID;
	}
	int getSvLen()const {	return ((r2.st_pos + r2.ed_pos) - (r1.st_pos + r1.ed_pos))/2; }
	bool isInRegion(RefRegion &block_region){
			if(r1.region_overlap(block_region) && r2.region_overlap(block_region))
				return true;
			return false;
	}
	friend std::ostream& operator<<(std::ostream& os, const SVE_core& r){
		os << 	"r1: " << r.r1 << "\t"
				"r2: " << r.r2 << "\t"
				"sv len: " << r.getSvLen() << "\t";
	  return os;
	}
	float overlapRate(SVE_core &s){
		RefRegion A(0,   r1.getMiddle(),   r2.getMiddle());
		RefRegion B(0, s.r1.getMiddle(), s.r2.getMiddle());
		return A.region_overlap_percent(B);
	}
};

struct SAMPLE_INFO{
	SAMPLE_INFO(RST::T is_begin_, uint8_t ori_at_r1_, uint16_t score_){
		signal_type = SIG::UNKNOWN;
		sample_ID = SAMPLE_UNKNOWN;
		ori_at_r1 = ori_at_r1_;
		is_solid = is_begin_;
		if(is_begin_ == RST::BEGIN){
			scoreB = score_; readNumB = 1;
			scoreE = 0;		 readNumE = 0;
		}
		else if(is_begin_ == RST::END){
			scoreB = 0; 	 readNumB = 0;
			scoreE = score_; readNumE = 1;
		}
		else if((is_begin_ == RST::UNKNOWN)){
			scoreB = 0; readNumB = 0;
			scoreE = 0; readNumE = 0;
		}
	}
	SAMPLE_INFO(RST::T is_begin_, uint8_t ori_at_r1_, uint16_t score_, int readNum){
		signal_type = SIG::UNKNOWN;
		sample_ID = SAMPLE_UNKNOWN;
		ori_at_r1 = ori_at_r1_;
		is_solid = is_begin_;
		if(is_begin_ == RST::BEGIN){
			scoreB = score_; readNumB = readNum;
			scoreE = 0;		 readNumE = 0;
		}
		else if(is_begin_ == RST::END){
			scoreB = 0; 	 readNumB = 0;
			scoreE = score_; readNumE = readNum;
		}
		else if((is_begin_ == RST::UNKNOWN)){
			scoreB = 0; readNumB = 0;
			scoreE = 0; readNumE = 0;
		}
		else if((is_begin_ == RST::SOLID)){
			scoreB = score_; readNumB = readNum;
			scoreE = score_; readNumE = readNum;
		}
	}
	SAMPLE_INFO(){}

	void combine_score(SAMPLE_INFO &i){
		if(i.is_solid == RST::BEGIN || i.is_solid == RST::SOLID){
			scoreB += i.scoreB; readNumB += i.readNumB;
			if(is_solid == RST::END) is_solid = RST::SOLID;//add solid
		}
		if(i.is_solid == RST::END || i.is_solid == RST::SOLID){
			scoreE += i.scoreE; readNumE += i.readNumE;
			if(is_solid == RST::BEGIN) is_solid = RST::SOLID;//add solid
		}
	}
	uint16_t getScore(){return scoreB + scoreE;}
	uint16_t getReadNum(){return readNumB + readNumE;}
	friend std::ostream& operator<<(std::ostream& os, const SAMPLE_INFO& r){
		std::string ori_at((r.ori_at_r1 == true)?"r1":"r2");
		std::string sampleID;	if(r.sample_ID == SAMPLE_UNKNOWN)sampleID = "UNKNOWN";else sampleID = r.sample_ID;
		os << 	"type_ID: " << SV::STR[r.type_ID] << "\t"
				"ori_at: " << ori_at  << "\t"
				"sample_ID: " << sampleID << "\t"
				"signal_type: " << SIG::STR[r.signal_type] << "\t"
				"is_solid: " << RST::STR[r.is_solid] << "\t"
				"scoreB: " << r.scoreB << "\t"
				"scoreE: " << r.scoreE << "\t"
				"readNumB: " << r.readNumB << "\t"
				"readNumE: " << r.readNumE << "\t";
	  return os;
	}

	void show(){
		std::string ori_at((ori_at_r1 == true)?"r1":"r2");
		std::string sampleID;	if(sample_ID == SAMPLE_UNKNOWN)sampleID = "UNKNOWN";else sampleID = sample_ID;
		fprintf(stderr, "type_ID: %s\t", SV::STR[type_ID].c_str());
		fprintf(stderr, "ori_at: %s\t", ori_at.c_str());
		fprintf(stderr, "sample_ID: %s\t", sampleID.c_str());
		fprintf(stderr, "signal_type: %s\t",SIG::STR[signal_type].c_str());
		fprintf(stderr, "is_solid: %s\t",RST::STR[is_solid].c_str());
		fprintf(stderr, "scoreB: %d\t",scoreB);
		fprintf(stderr, "scoreE: %d\t",scoreE);
		fprintf(stderr, "readNumB: %d\t",readNumB);
		fprintf(stderr, "readNumE: %d\t",readNumE);
	}

public:
	SV::T 		type_ID = SV::UNKNOWN;
	SIG::T 		signal_type;//Sh/Dr...
	RST::T  	is_solid;
	bool		ori_at_r1;//bool

	uint32_t 	sample_ID;
	uint16_t	scoreB;
	uint16_t	scoreE;
	uint16_t 	readNumB;
	uint16_t 	readNumE;
};

#define MIN_SCORE_SOLID 50
#define MIN_SCORE_SOLID_BE 15
#define MIN_SCORE_UNSTABLE 100

class SVE:public SVE_core{
public:
	SVE():info(RST::UNKNOWN, 0, 0){ r1.chr_ID = -1; r2.chr_ID = -1;}
	SVE(RST::T is_begin_, uint8_t ori_at_r1_, uint16_t score_, RefRegion &r1_ ,RefRegion &r2_)
	:info(is_begin_, ori_at_r1_, score_){ r1 = (r1_); r2 = (r2_);}//new for SH signal
	SVE(bam1_core_t *core, bool dir, bool m_dir, int insert_region_len, int middle_size, RST::T is_begin)
	:info(is_begin, true, MIN(15, core->qual)){//new for DR signal
		int this_region_st = core->pos + ((dir == FORWARD)?middle_size:(- insert_region_len));
		int mate_region_st = core->mpos + ((m_dir == FORWARD)?core->l_qseq:(- insert_region_len));//only use l_seq because we can`t get mate middle size directly
		RefRegion this_region(core->tid,  this_region_st, this_region_st + insert_region_len);
		RefRegion mate_region(core->mtid, mate_region_st, mate_region_st + insert_region_len);
		if(dir != m_dir){//for del and dup
			r1 = (dir == FORWARD)?this_region:mate_region;
			r2 = (dir == FORWARD)?mate_region:this_region;
		}else{//for reverse
			bool this_smaller = this_region.isForwardRegion(mate_region);
			r1 = (this_smaller)?this_region:mate_region;
			r2 = (this_smaller)?mate_region:this_region;
		}
	}
	SVE(int chr_ID_1, int top_r1, int chr_ID_2, int top_r2, int read_n, int score, RST::T is_begin, int accept_region)
		:info(is_begin, true, score, read_n){
		 r1 = RefRegion(chr_ID_1, top_r1 - accept_region, top_r1 + accept_region);
		 r2 = RefRegion(chr_ID_2, top_r2 - accept_region, top_r2 + accept_region);//store sve for combine step1
	}

	SVE(const SVE &b){
		memcpy(this, &b, sizeof(SVE));
	}

	void setTGS(int chr_ID_1, int top_r1, int chr_ID_2, int top_r2, int accept_region ){
		info.is_solid = RST::SOLID;
		info.scoreB = 999;//big enough
		info.scoreE = 999;
		 r1 = RefRegion(chr_ID_1, top_r1 - accept_region, top_r1 + accept_region);
		 r2 = RefRegion(chr_ID_2, top_r2 - accept_region, top_r2 + accept_region);//store sve for combine step1
	}

	bool SVE_combine_UNION(SVE & b){
		if(r1.region_overlap(b.r1) && r2.region_overlap(b.r2)){
			r1.Combine(b.r1, true);
			r2.Combine(b.r2, true);
			info.combine_score(b.info);
			return true;
		}
		return false;
	}

	//S4: filter for SVs
	bool vcfFilter(){
		if(info.is_solid == RST::SOLID){ // for stable results
			if(info.scoreB + info.scoreE < MIN_SCORE_SOLID) return false;
			if(info.scoreB < MIN_SCORE_SOLID_BE) return false;
			if(info.scoreE < MIN_SCORE_SOLID_BE) return false;
		}
		else if(info.scoreB + info.scoreE < MIN_SCORE_UNSTABLE) return false;// for unstable results
		return true;
	}

	void setSignalType(SV::T sv_type, SIG::T signal_type_){
		info.type_ID = sv_type;
		info.signal_type = signal_type_;
	}

	friend std::ostream& operator<<(std::ostream& os, const SVE& r){
		os << (SVE_core)r << r.info	<< std::endl;
	  return os;
	}

	void show(){
		info.show();
		fprintf(stderr, "r1: "); r1.print(stderr);
		fprintf(stderr, "r2: "); r2.print(stderr);
		fprintf(stderr, "\n");
	}

	SAMPLE_INFO info;

};

typedef std::vector<SVE> SVE_L;

void sve_combine_STEP1_SH(SVE_L & l, int min_score, bool is_union);
void sve_repeat_edge_det(SVE_L &l, uint8_t * ref);
