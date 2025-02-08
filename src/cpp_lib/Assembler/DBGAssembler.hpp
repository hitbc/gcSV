/*
 * loopDBG.hpp
 *
 *  Created on: 2024年1月6日
 *      Author: fenghe
 */

#ifndef CPP_LIB_ASSEMBLER_DBGASSEMBLER_HPP_
#define CPP_LIB_ASSEMBLER_DBGASSEMBLER_HPP_

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <string>
#include <vector>
#include "../../clib/utils.h"
#include "../../SVcalling_core/Contig_aligner.hpp"
#include "../RefRegion.hpp"
#include "math.h"
struct ReadPos{
	int readID;
	int readPos;
	ReadPos(int readID, int readPos){
		this->readID = readID;
		this->readPos = readPos;
	}
	void show(FILE * o){
		fprintf(o, "[%d:%d]\t", readID, readPos);
	}
};

struct KmerCounter{
	int read_Ngs_n;
	int read_Tgs_n;
	int refN;
	int ref_pos;
	std::vector<ReadPos> NGS_rp;
	std::vector<ReadPos> TGS_rp;
	KmerCounter(){
		NGS_rp.clear(); TGS_rp.clear();
		read_Ngs_n = 0; refN = 0; read_Tgs_n = 0;
		ref_pos = -1; // -1 : UNKNOWN -2 : repeat in REF
	}
	void add_new_kmer(int isRead, int ref_read_pos, int strIndex, bool storeReadInfo);
	void show(std::string &kmerStr){ fprintf(stderr, "%s: [ %d | %d ]\n", kmerStr.c_str(), read_Ngs_n, refN); }
	bool is_low_depth(int minDepthNGS, int minDepthTGS){return (read_Ngs_n < minDepthNGS && refN == 0 && read_Tgs_n < minDepthTGS);}
	int getReadDepth(int LOW_Depth, int minDepthTGS){
		int read_count = 0;
		if(read_Ngs_n >= LOW_Depth)
			read_count = read_Ngs_n;
		if(read_Tgs_n >= minDepthTGS)
			read_count += MAX(LOW_Depth, read_Tgs_n);
		return read_count;
	}
};

typedef std::unordered_map<std::string, KmerCounter> KmerCounterTable;

struct DBG_INFO_item{
	std::string kmer;
	KmerCounter * kcp;//the pointer to the original kmer counter
	int depth_read;
	int depth_ref;
	int kmerId;
	int unitigId = 0;
	int unitigPos = 0;
	int unitigRootId = 0;
	uint32_t in[4];
	uint32_t out[4];
	uint8_t IN_NUM = 0;
	uint8_t OUT_NUM = 0;

	uint8_t is_UNITIG_begin = false;
	uint8_t is_UNITIG_end = false;

	uint8_t isRepeat = false;
	int kmerPOS_ref;
	DBG_INFO_item(){
		kmerId = 0;
		depth_read = 0;
		depth_ref = 0;
		kmerPOS_ref = -1;
		kcp = NULL;
	}

	void setNew(int id, int LOW_Depth, int minDepthTGS, const std::string& kmer, KmerCounter * kcp){
		this->kcp = kcp;
		this->kmerPOS_ref = kcp->ref_pos;
		this->kmer = kmer;
		kmerId = id;
		this->depth_read = kcp->getReadDepth(LOW_Depth, minDepthTGS);
		this->depth_ref =  kcp->refN;
	}

	void setID(int kmer_Id){ kmerId = kmer_Id ; }

	void setInOut(bool isIn, int ID){
		if(isIn) in[IN_NUM++]  = ID;
		else     out[OUT_NUM++] = ID;
	}

	void showAction(const std::string &kmer, int AC_I){
		fprintf(stderr, "ACTION_%d: id: %d\t%s\n",AC_I, kmerId, kmer.c_str());
	}

	void showData();
	void showDataSimple();
	void showDataDetail();
};

struct UNITIG_INFO_item{
	int id;
	int totalDepthRead = 0;
	int totalDepthRef = 0;
	float averageKmerDepth_read = 0;//==((float)totalDepthRead)*100/kmer_list.size();
	float averageKmerDepth_read_head = 0;//==((float)totalDepthRead)*100/kmer_list.size();
	float averageKmerDepth_read_tail = 0;//==((float)totalDepthRead)*100/kmer_list.size();
	float averageKmerDepth_ref = 0;

	bool with_read = false;
	bool with_ref = false;

	std::vector<int> kmer_list;
	std::string unitig;
	std::vector<int> inBranchKmer;
	std::vector<int> outBranchKmer;
	std::vector<int> inBranchUni;
	std::vector<int> outBranchUni;
	int inBranch_n = 0;
	int outBranch_n = 0;
	int loopCurr = 0;
	int loopRoot = 0;

	int finalRootId = 0;
	bool isRepeat = false;

	int refPos = 0;

	UNITIG_INFO_item(int uid){ id = uid; }
	void setAverageKmerDepth(std::vector<DBG_INFO_item *> &id2info);
	void setWithRead(int headNodeId, int endNodeId, std::vector<DBG_INFO_item *>& id2info, int MIN_DEPTH);
	int get_ref_pos(std::vector<DBG_INFO_item *>& id2info);
	void setString(std::vector<DBG_INFO_item *>& id2info);
	void setInOutKmer(std::vector<DBG_INFO_item *>& id2info);
	void setBranchUni(std::vector<DBG_INFO_item *>& id2info);
	void show(std::vector<DBG_INFO_item *>& id2info);
	std::vector<ReadPos> & getAllReadTGS(std::vector<DBG_INFO_item *>& id2info){
		 return id2info[kmer_list[0]]->kcp->TGS_rp;
	}
	std::vector<ReadPos> & getAllReadTGSTail(std::vector<DBG_INFO_item *>& id2info){
		 return id2info[kmer_list.back()]->kcp->TGS_rp;
	}

};

#define REPEAT_DEBUG false

struct KEY_UNITIG{
	UNITIG_INFO_item* p;
	int ref_pos;
	static inline int cmp_by_ref_pos(const KEY_UNITIG &a, const KEY_UNITIG &b){
		return a.ref_pos < b.ref_pos;
	}
};

struct PathNodeCount{
	UNITIG_INFO_item *up;
	int count;
};

int getBinBase(char c);

struct Path{
	std::vector<UNITIG_INFO_item *> unitigPL;//the UNITIG list for the SVs
	std::vector<PathNodeCount> pnc;//the count for the list of SVs
	int bgnUID;//the begin UID
	int endUID;//the end UID
	int SV_len;//suggest SV length
	int is_final_path;//suggest SV length
	int repeatRstNum;
	int type;//==1: repeat; ==2: linear
	uint64_t pathKey;
	int ID_SUM;
	double UNIQ_ID_product;
	int UNIQ_ID_SUM;
	//
	std::vector<ReadPos> support_read_l;

	bool isSAMEwith(Path &c);
	void setNodeCount(std::map<UNITIG_INFO_item *, int > &pathNodeCountBuff);
	void setSVLen(std::vector<DBG_INFO_item *> &id2info);

	static inline int cmp_by_uid(const Path &a, const Path &b){
		if(a.bgnUID != b.bgnUID)					return a.bgnUID < b.bgnUID;
		if(a.endUID != b.endUID)					return a.endUID < b.endUID;
		if(a.unitigPL.size() != b.unitigPL.size())	return a.unitigPL.size() < b.unitigPL.size();
		if(a.SV_len != b.SV_len)					return a.SV_len < b.SV_len;
		if(a.ID_SUM != b.ID_SUM)					return a.ID_SUM < b.ID_SUM;
		if(a.pathKey != b.pathKey)					return a.pathKey < b.pathKey;
		return a.type < b.type;
	}

	static inline int cmp_by_idSUM(const Path &a, const Path &b){
		if(a.bgnUID != b.bgnUID)					return a.bgnUID < b.bgnUID;
		if(a.endUID != b.endUID)					return a.endUID < b.endUID;
		if(a.UNIQ_ID_product != b.UNIQ_ID_product)	return a.UNIQ_ID_product < b.UNIQ_ID_product;
		if(a.UNIQ_ID_SUM != b.UNIQ_ID_SUM)			return a.UNIQ_ID_SUM < b.UNIQ_ID_SUM;
		if(a.SV_len != b.SV_len)					return a.SV_len < b.SV_len;
		return a.type < b.type;
	}

	static inline int cmp_by_SV_LEN(const Path &a, const Path &b){
		return a.SV_len < b.SV_len;
	}

	static inline int cmp_by_support_read_number(const Path &a, const Path &b){
		return a.support_read_l.size() > b.support_read_l.size();
	}

	bool isSAMENodeSet(Path &c);
	void showSim(std::vector<DBG_INFO_item *> &id2info);
	void show(std::vector<DBG_INFO_item *> &id2info);
	void getContig(std::vector<DBG_INFO_item *> &id2info, std::string & contig){
		bool isHead = true;
		for(UNITIG_INFO_item* up: unitigPL){
			for(uint i = 0; i < up->kmer_list.size(); i++){
				if(isHead){
					contig = id2info[up->kmer_list[i]]->kmer;
					isHead = false;
				}else
					contig += id2info[up->kmer_list[i]]->kmer.back();
			}
		}
	}

	void getReadPosition_to_contig(std::vector<DBG_INFO_item *> &id2info, std::vector<ReadPos> & rst){
		//the first UNITIG
		UNITIG_INFO_item* up = unitigPL[0];
		//the last KMER in the first UNITIG
		uint kmer_size = up->kmer_list.size();
		DBG_INFO_item * k = id2info[up->kmer_list[kmer_size - 1]];
		//get the read positions
		rst = k->kcp->TGS_rp;
		//
		for(ReadPos &p: rst){
			p.readPos -= (kmer_size - 1);
		}
	}

	void showContig(std::vector<DBG_INFO_item *> &id2info, std::string & contig){
		bool isHead = true;
		fprintf(stderr, "Show contig\t %s\n", contig.c_str());
		fprintf(stderr, "Show contig\t");
		int uid = 0;
		isHead = true;
		for(UNITIG_INFO_item* up: unitigPL){
			char c = 'A' + uid; uid++;
			for(uint i = 0; i < up->kmer_list.size(); i++){
				if(isHead){
					for(uint j = 0; j < id2info[up->kmer_list[i]]->kmer.size(); j++){
						fprintf(stderr, "S");
					}
					isHead = false;
				}else
					fprintf(stderr, "%c", c);
			}
		}
		fprintf(stderr, "\n");
	}

	void showReference(int refPosBgn,  int refPosEnd, std::string & reference_str){
		int refLoadLen = refPosEnd - refPosBgn;
		fprintf(stderr, "Show Reference [%d~%d, len %d]\t%s\n", refPosBgn, refPosEnd, refLoadLen, reference_str.substr(refPosBgn, refLoadLen).c_str());
	}

	void getReferenceRegion(int curKmerLen, int &refPosBgn,  int &refPosEnd){
		refPosBgn = unitigPL[0]->refPos;
		if(refPosBgn < 0)
			refPosBgn = 0;
		refPosEnd = unitigPL.back()->refPos + unitigPL.back()->kmer_list.size() + curKmerLen - 1;
	}
};

struct NodeAnalysisINFO{
	int currentCountDiff = 0;
	int basic_Count = -1;
	int step_Count = -1;
	float depth = 0;
	UNITIG_INFO_item * up;
	void show(){
		fprintf(stderr, " NODE_INFO ID %5d,"
				" basic_Count %5d , step_Count %5d ,"
				" depth %.1f\n",
				up->id, basic_Count, step_Count, depth);
	}
};

struct A_B_para_suggest{
	int A; int B;
	A_B_para_suggest(int A, int B){ this->A = A; this->B = B; }
};

struct AB_SuggsetItem{
	std::vector<int> diff3_NodeID;
	std::vector<float> diff3;
	float averageDepthBasic;
	int A; int B;//
	float diffSUM;
	bool suggsetSuccess;
	static inline int cmp_by_DiffSum(const AB_SuggsetItem &a, const AB_SuggsetItem &b){
		return a.diffSUM < b.diffSUM;
	}
	void show(){
		fprintf(stderr, "suggestFinalSVLength: Final suggested diffSUM_MIN %.2f, @ A %d B %d suggsetSuccess %d averageDepthBasic %.2f\n",
				diffSUM, A, B, suggsetSuccess, averageDepthBasic);
		for(int i: diff3_NodeID)
			fprintf(stderr, "%6d ", i);
		fprintf(stderr, "\n");
		for(float i: diff3)
			fprintf(stderr, "%6.2f ", i);
		fprintf(stderr, "\n\n");
	}
};

struct AB_SuggsetHandler{

	std::vector<A_B_para_suggest> ab_sug_Diploid;
	std::vector<AB_SuggsetItem> dr_l;
	std::vector<A_B_para_suggest> ab_sug_Haploid;
	bool isDiploid;

	bool isFit(){
		return dr_l[0].suggsetSuccess;
	}
	//only using repeat node
	AB_SuggsetItem * getMaxSuggsetM2(std::vector<NodeAnalysisINFO> &nfl, bool isDiploid, bool usingHeadNode,
			int bgnUID, int endUID, int adjustRate, float maxError){
		float totalDepth = 0;
		int totalBaseCount = 0;
		int totalStepCount = 0;
		for(NodeAnalysisINFO &na: nfl){
			if(na.step_Count != 0){
				totalBaseCount += na.basic_Count;
				totalStepCount += na.step_Count;
				totalDepth += na.depth*(1+adjustRate*na.step_Count);
			}
		}
		if(nfl.size() >= 6){
			usingHeadNode = false;
		}
		if(totalBaseCount == 0 ||totalStepCount == 0){
			fprintf(stderr, "The depth-based caller is not used because insufficient signals");
			return NULL;
		}

		this->isDiploid = isDiploid;
		std::vector<A_B_para_suggest> * suggsetMode = NULL;
		suggsetMode = &(ab_sug_Diploid);//(isDiploid)?&(ab_sug_Diploid):&(ab_sug_Haploid);
		dr_l.clear();
		for(A_B_para_suggest &abs: *suggsetMode){
			int A = abs.A; int B = abs.B;
			float averageDepthBasic = totalDepth/(A*totalBaseCount + B*totalStepCount);
			float diffSUM = 0;
			dr_l.emplace_back();
			dr_l.back().A = A; dr_l.back().B = B;
			for(NodeAnalysisINFO &na: nfl){
				if(na.step_Count == 0) continue;
				if(!usingHeadNode && (na.up->id == bgnUID || na.up->id == endUID))		continue;
				float finalDiffRate = 0;
				if(A*na.basic_Count + B*na.step_Count == 0){
					finalDiffRate = (na.depth < 1.5)?0:3;
				}else{
					float diffRate1 = (A*na.basic_Count + B*na.step_Count)*averageDepthBasic/(na.depth*(1+adjustRate*na.step_Count));
					float diffRate2 = 1/diffRate1;
					float diffRate1_ABS = ABS_U(diffRate1, 1);
					float diffRate2_ABS = ABS_U(diffRate2, 1);
					finalDiffRate = MAX(diffRate1_ABS, diffRate2_ABS);
				}
				dr_l.back().diff3_NodeID.emplace_back(na.up->id);
				dr_l.back().diff3.emplace_back(finalDiffRate);
				diffSUM += finalDiffRate;
				dr_l.back().averageDepthBasic = averageDepthBasic;
			}
			dr_l.back().diffSUM = diffSUM;
			dr_l.back().suggsetSuccess = (dr_l.back().diffSUM/nfl.size()>maxError)?false:true;
		}
		std::sort(dr_l.begin(), dr_l.end(), AB_SuggsetItem::cmp_by_DiffSum);
		if(dr_l[0].suggsetSuccess == false){
			fprintf(stderr, "suggestFinalSVLength FAIL::diffSUM HIGH\n");
		}
		//suggest MAX
		dr_l[0].show();
		dr_l[1].show();

		bool isSimilar = false;
		if(dr_l[1].diffSUM < 1.5*dr_l[0].diffSUM){
			isSimilar = true;
			fprintf(stderr, "suggestFinalSVLength::similar diffSUM(M1), %.2f\n", dr_l[1].diffSUM/(dr_l[0].diffSUM+0.001));
		}else if(dr_l[1].diffSUM - dr_l[0].diffSUM < 0.1 && dr_l[0].diffSUM > 0.1){
			isSimilar = true;
			fprintf(stderr, "suggestFinalSVLength::similar diffSUM(M2), %.2f\n", dr_l[1].diffSUM-(dr_l[0].diffSUM));
		}

		if(isSimilar){
			if(dr_l[0].A == 1 && dr_l[0].B == 1 && dr_l[1].A == 2 && dr_l[1].B == 1 ){

			}else if(dr_l[0].A == 2 && dr_l[0].B == 1 && dr_l[1].A == 1 && dr_l[1].B == 1 ){
				std::swap(dr_l[1], dr_l[0]);
			}
			else{
				dr_l[0].suggsetSuccess = false;
			}
		}


		//dr_l[2].show();
		return &(dr_l[0]);
	}

	AB_SuggsetItem * getMaxSuggset(std::vector<NodeAnalysisINFO> &nfl, bool isDiploid, bool usingHeadNode,
			int bgnUID, int endUID, int adjustRate, float maxError){
		float totalDepth = 0;
		int totalBaseCount = 0;
		int totalStepCount = 0;
		for(NodeAnalysisINFO &na: nfl){
			totalBaseCount += na.basic_Count;
			totalStepCount += na.step_Count;
			totalDepth += na.depth*(1+adjustRate*na.step_Count);
		}
		if(nfl.size() >= 6){
			usingHeadNode = false;
		}

		this->isDiploid = isDiploid;
		std::vector<A_B_para_suggest> * suggsetMode = NULL;
		suggsetMode = (isDiploid)?&(ab_sug_Diploid):&(ab_sug_Haploid);
		dr_l.clear();
		for(A_B_para_suggest &abs: *suggsetMode){
			int A = abs.A; int B = abs.B;
			float averageDepthBasic = totalDepth/(A*totalBaseCount + B*totalStepCount);
			float diffSUM = 0;
			dr_l.emplace_back();
			dr_l.back().A = A; dr_l.back().B = B;
			for(NodeAnalysisINFO &na: nfl){
				if(!usingHeadNode && (na.up->id == bgnUID || na.up->id == endUID))		continue;
				float finalDiffRate = 0;
				if(A*na.basic_Count + B*na.step_Count == 0){
					finalDiffRate = (na.depth <= 2)?0:3;
				}else{
					float diffRate1 = (A*na.basic_Count + B*na.step_Count)*averageDepthBasic/(na.depth*(1+adjustRate*na.step_Count));
					float diffRate2 = 1/diffRate1;
					float diffRate1_ABS = ABS_U(diffRate1, 1);
					float diffRate2_ABS = ABS_U(diffRate2, 1);
					finalDiffRate = MAX(diffRate1_ABS, diffRate2_ABS);
				}
				dr_l.back().diff3_NodeID.emplace_back(na.up->id);
				dr_l.back().diff3.emplace_back(finalDiffRate);
				diffSUM += finalDiffRate;
			}
			dr_l.back().diffSUM = diffSUM;
			dr_l.back().suggsetSuccess = (dr_l.back().diffSUM/nfl.size()>maxError)?false:true;
		}
		std::sort(dr_l.begin(), dr_l.end(), AB_SuggsetItem::cmp_by_DiffSum);
		if(dr_l[0].suggsetSuccess == false){
			fprintf(stderr, "suggestFinalSVLength FAIL diffSUM HIGH\n");
		}
		//suggest MAX
		dr_l[0].show();
		//dr_l[1].show();
		//dr_l[2].show();
		return &(dr_l[0]);
	}

	bool getSVLen(AB_SuggsetItem * bsetS, int basicSVLength, int stepSVLength, int &suggestFinalSVLength, bool &isHOM, bool using_Linear){
		if(bsetS == NULL)
			return false;
		fprintf(stderr, "getSVLen begin: basicSVLength %d, stepSVLength %d \n", basicSVLength, stepSVLength);
		bool isFitRst = false;
		if(bsetS->suggsetSuccess){
			isFitRst = true;
		}
		if(bsetS->B == 0 || bsetS->A >3){//1:0
			if(using_Linear){
				suggestFinalSVLength = basicSVLength;
			}else{
				suggestFinalSVLength = 0;
				isFitRst = false;
			}
		}
		if(bsetS->A == 2 && bsetS->B == 1){//2:1
			//there is too much wrong when handle 2:1 deletion, there for 2:1 is used as
			//same as the 1:1
//				suggestFinalSVLength = basicSVLength + stepSVLength;
			if(basicSVLength == 0)//2:1 insertion
				suggestFinalSVLength = basicSVLength + stepSVLength;
			else if(basicSVLength + stepSVLength == 0){
				suggestFinalSVLength = basicSVLength;
//					bsetS->suggsetSuccess = (bsetS->diffSUM/bsetS->diff3.size()>0.2)?false:true;
//					if(!bsetS->suggsetSuccess)
//						isFitRst = false;
				//2:1 deletion
			}else{
				suggestFinalSVLength = basicSVLength + stepSVLength;
			}
			isHOM = false;
		}else if(bsetS->A == 1 && bsetS->B == 1){//1:1
			suggestFinalSVLength = basicSVLength + stepSVLength;
		}else if(bsetS->A == 1 && bsetS->B <= 4){//1:?(2,3,4)
			suggestFinalSVLength = basicSVLength + bsetS->B*stepSVLength;
		}else{
			isFitRst = false;
		}
		return isFitRst;
	}

	void init(){
		ab_sug_Diploid.emplace_back(1, 0);
		ab_sug_Diploid.emplace_back(8, 1);
		ab_sug_Diploid.emplace_back(2, 1);
		ab_sug_Diploid.emplace_back(1, 1);
		ab_sug_Diploid.emplace_back(1, 2);
		ab_sug_Diploid.emplace_back(1, 3);
		ab_sug_Diploid.emplace_back(1, 4);
		ab_sug_Diploid.emplace_back(1, 5);

		ab_sug_Haploid.emplace_back(1, 0);
		ab_sug_Haploid.emplace_back(8, 1);
		ab_sug_Haploid.emplace_back(1, 1);
		ab_sug_Haploid.emplace_back(1, 2);
		ab_sug_Haploid.emplace_back(1, 3);
		ab_sug_Haploid.emplace_back(1, 4);
		ab_sug_Haploid.emplace_back(1, 5);
	}
};

struct PathCluster{
	int bgnUID;//the begin and end UID
	int endUID;
	//the basic UNITIG nodes and count
	//the loop  UNITIG nodes and count
	bool isCPX = false;

	int basicSVLength;
	int stepSVLength;

	int signalNum = 0;
	Path basePath;

	std::vector<NodeAnalysisINFO> nfl;
	std::vector<Path> path_list;

	//IDs
	double BASE_UNIQ_ID_product;
	//final
	bool isFit;
	int suggestFinalSVLength_N;
	int suggestFinalSVLength;
	int confidenceScore;

	int getPathIndex(int suggestFinalSVLength){
		int path_idx = -1;
		int curP = 0;
		for(Path &p :path_list){
			if(p.SV_len == suggestFinalSVLength)
				path_idx = curP;
			curP++;
		}
		return path_idx;
	}

	void showAllPAth(std::vector<DBG_INFO_item *> &id2info){
		for(Path &p :path_list){
			p.show(id2info);
		}
	}

	void getSuggestSVLen(AB_SuggsetHandler &abH, bool usingHeadNode, bool isDoubleHaplotype, bool & isHOM){
		isHOM = isDoubleHaplotype;
		isFit = false;
		if(signalNum == 1){
			suggestFinalSVLength = basePath.SV_len;
			isFit = false;
			fprintf(stderr, "Skip linear path @ suggestFinalSVLength%d\n", suggestFinalSVLength);
		}
		else if(isCPX){
			suggestFinalSVLength = -1;
		}else{
			AB_SuggsetItem * bsetS = abH.getMaxSuggsetM2(nfl, isDoubleHaplotype, usingHeadNode, bgnUID, endUID, 0.3, 0.4);
			isFit = abH.getSVLen(bsetS, basicSVLength, stepSVLength, suggestFinalSVLength, isHOM, false);
		}
		if(isFit)
			fprintf(stderr, "Analysis END, suggestFinalSVLength %5d , basicSVLength %5d , stepSVLength %5d \n", suggestFinalSVLength, basicSVLength, stepSVLength);

	}
	int repeatResultNumber = 0;

	void updateStepCountByNewSignal();
	void generateBasicAndStepSVLen();
	bool isSAME_ID_SET(Path &p){ return basePath.isSAMENodeSet(p); }
	void addSignal(Path &p);

	void show(std::vector<DBG_INFO_item *> &id2info){
		if(isCPX){
			fprintf(stderr, "is CPX, only show base path\n");
			basePath.show(id2info);
		}else{
			fprintf(stderr, "repeatResultNumber %5d\t", repeatResultNumber);
			if(signalNum == 1){
				fprintf(stderr, "is Single Signal\n");
				basePath.show(id2info);
				//fprintf(stderr, "Single Signal END\n");
			}else{
	            fprintf(stderr, "SINFO_WITHRepeat bgnUID %5d  endUID %5d , basicSVLength %5d , "
	                    "stepSVLength  %5d , signalNum %5d isCPX %d \n", bgnUID,  endUID,
	                    basicSVLength, stepSVLength, signalNum, isCPX);
	            for(NodeAnalysisINFO &nf_: nfl)
	            	nf_.show();
			}
		}
	}
};

struct PathFilter{
	struct ReadOffset{
		int readID;
		int readOffset;
		static inline int cmp_by_position(const ReadOffset &a, const ReadOffset &b){
			return a.readOffset < b.readOffset;
		}
	};

	 void showRead(const char *contig_str, ReadOffset & cro, const char *rstr, int rlen, int contigSize){
		for(int i = 0; i < cro.readOffset; i++) fprintf(stderr, " ");
		for(int i = 0; i < rlen; i++){
			int inContigBaseIdx = cro.readOffset + i;
			if(inContigBaseIdx >= 0 && inContigBaseIdx < contigSize){
				if(rstr[i] == contig_str[cro.readOffset + i]){
					fprintf(stderr, "-");
				}
				else{
					fprintf(stderr, "%c", rstr[i]);
				}
			}
		}
		fprintf(stderr, "cro.readID %d\n", cro.readID);
	 }

	 struct ContigBaseCounter{
			 std::vector<std::vector<int>> baseN;
			 int contigSize;
			 void clear(int contigSize){
				 this->contigSize = contigSize;
				baseN.resize(4);
				for(int i = 0; i < 4; i++){
					baseN[i].resize(contigSize);
					memset(&(baseN[i][0]),0,sizeof(int)*contigSize);
				}
			 }

			 void addToContigCounter(const char *contig_str, ReadOffset & cro, const char *rstr, int rlen){
				for(int i = 0; i < rlen; i++){
					int binBase = getBinBase(rstr[i]);
					int inContigBaseIdx = cro.readOffset + i;
					if(inContigBaseIdx >= 0 && inContigBaseIdx < contigSize && cro.readOffset + i >= 0){
						if(binBase <= 3){
							baseN[binBase][cro.readOffset + i]++;
						}
					}
				}
			 }

			 void showContigCounter(const char *contig){
				for(int j = 0; j < contigSize; j++){
					fprintf(stderr, "  %c ", contig[j]);
				}
				fprintf(stderr, "\n");
				for(int i = 0; i < 4; i++){
					for(int j = 0; j < contigSize; j++){
						fprintf(stderr, "%3d ", baseN[i][j]);
					}
					fprintf(stderr, "\n");
				}
			 }
		 };

	 struct AVEKmerNum{
		 int uid;
		 int kmerSize;
		 int depthSUM;
		 int minDepth;
		 AVEKmerNum(int uid, int kmerSize, int depthSUM, int minDepth){
			 this->uid = uid;
			 this->kmerSize = kmerSize;
			 this->depthSUM = depthSUM;
			 this->minDepth = minDepth;
		 }
		 void show(){
			 fprintf(stderr, "old_unitID %5d, N %5d, depthSUM %5d dep %f @ min %5d\n", uid, kmerSize, depthSUM, (float)depthSUM/kmerSize, minDepth);
		 }
	 };

	 struct ContigKmerCounter{
		 std::vector<int> contigKmerCounter;
		 std::vector<int> contigKmerCounterK100;
		 std::unordered_set<int> uniqReadSet;
		 std::vector<int> contigKmerCounterUNIQ;
		 std::vector<int> unitigIDList;
		 int contigSize;
		 int kmerLen;
		 int kmerN;
		 int kmerN_100;
		 int searchBgn = 0;
		 int searchEnd = 0;
		 int KmerK100L;
		 void clear(int contigSize, int kmerLen, int searchBgn, int searchEnd, int KmerK100L){
			 this->KmerK100L = KmerK100L;
			 this->contigSize = contigSize;
			 this->kmerLen = kmerLen;
			 kmerN = contigSize - kmerLen + 1;
			 kmerN_100 = contigSize - KmerK100L + 1;
			 contigKmerCounter.resize(kmerN);
			 contigKmerCounterUNIQ.resize(kmerN);
			 unitigIDList.resize(kmerN);
			 contigKmerCounterK100.resize(kmerN_100);
			 if(kmerN_100 > 0)
				 memset(&(contigKmerCounterK100[0]),0,sizeof(int)*kmerN_100);
			 memset(&(contigKmerCounter[0]),0,sizeof(int)*kmerN);
			 memset(&(contigKmerCounterUNIQ[0]),0,sizeof(int)*kmerN);
			 this->searchBgn = searchBgn;
			 this->searchEnd = searchEnd;
			 uniqReadSet.clear();
		 }

		 bool addSignal(ReadOffset & cro, int rlen, bool readNeedUniq){
			int kmerReadN = rlen - kmerLen + 1;
			for(int i = 0; i < kmerReadN; i++){
				if(cro.readOffset + i >= 0){
					if(cro.readOffset + i >= 0 && cro.readOffset + i < kmerN)
						contigKmerCounter[cro.readOffset + i]++;
				}
			}
			int kmerReadN_100 = rlen - KmerK100L + 1;
			for(int i = 0; i < kmerReadN_100; i++){
				if(cro.readOffset + i >= 0){
					if(cro.readOffset + i >= 0 && cro.readOffset + i < kmerN_100)
						contigKmerCounterK100[cro.readOffset + i]++;
				}
			}

			if(readNeedUniq && uniqReadSet.find(cro.readID) == uniqReadSet.end()){
				for(int i = 0; i < kmerReadN; i++){
					if(cro.readOffset + i >= 0 && cro.readOffset + i < kmerN)
						contigKmerCounterUNIQ[cro.readOffset + i]++;
				}
				uniqReadSet.emplace(cro.readID);
				return true;
			}
			return false;
		 }

		 void addUID(std::vector<UNITIG_INFO_item *> &unitigPL){
			 int global_k_i = 0;
			for(uint curUniIdx = 0; curUniIdx < unitigPL.size(); curUniIdx++){
				UNITIG_INFO_item* up = unitigPL[curUniIdx];
				int ksize = up->kmer_list.size();
				int uid = up->id;
				for(int i = 0; i < ksize; i++){
					if(i + global_k_i >= 0 && i + global_k_i < kmerN){
						//xassert(i + global_k_i >= 0 && i + global_k_i < kmerN, "");
						unitigIDList[i + global_k_i] = uid;
					}
				}
				global_k_i += ksize;
			}
		 }

		 bool Kmer100DepthCheckForUnitig(int minDepth, bool isDEL, int SV_contig_bg, int SV_contig_ed){
			 int kmer_region_bg = SV_contig_bg - KmerK100L/2 - 20;
			 int kmer_region_ed = SV_contig_ed - KmerK100L/2 + 20;

			 int searchEnd100 = searchEnd + kmerLen - KmerK100L - 10;

			 kmer_region_bg = MAX(kmer_region_bg, 0);
			 kmer_region_ed = MIN(kmer_region_ed, searchEnd100);

			 int failNodeNum = 0;
			 for(int j = kmer_region_bg ; j < kmer_region_ed; j++){
				 if(contigKmerCounterK100[j] < minDepth){
					 failNodeNum++;
				 }
			 }
			 fprintf(stderr, "VNTR SV filter: SV_contig_bg %5d, SV_contig_ed %5d, kmer_region_bg %5d , kmer_region_ed %5d, KmerK100L %5d\n",
					 SV_contig_bg, SV_contig_ed, kmer_region_bg , kmer_region_ed, KmerK100L );
			for(int j = kmer_region_bg; j < kmer_region_ed; j++){
				int showBase = MIN(70, contigKmerCounterK100[j]);
				fprintf(stderr, "%c", '0' + showBase);
			}
			 fprintf(stderr, "\n");
			 if(failNodeNum >= 1){
				 fprintf(stderr, "Kmer100DepthCheckForUnitig checking fail, failNodeNum %d",
						 failNodeNum);
				 return false;
			 }
			 if(kmer_region_ed - kmer_region_bg < 10){
				 fprintf(stderr, "Kmer100DepthCheckForUnitig checking fail, contig is too short\n");
				 return false;
			 }
			 return true;
		 }

		 void calculateAveKmerForUnitig(std::vector<AVEKmerNum> & aveKmerNumL){
			aveKmerNumL.clear();
			int old_unitID = -1;
			int depthSUM = 0, N = 0; int minDepth = 999;
			for(int j = searchBgn; j < searchEnd; j++){
				if(unitigIDList[j] != old_unitID || j == searchEnd){
					if(old_unitID != -1)
						aveKmerNumL.emplace_back(old_unitID, N, depthSUM, minDepth);
					old_unitID = unitigIDList[j];
					depthSUM = 0;
					N = 0;
					minDepth = 999;
				}
				N++;
				minDepth = MIN(minDepth,contigKmerCounterUNIQ[j]);
				depthSUM += contigKmerCounterUNIQ[j];
			}
			aveKmerNumL.emplace_back(old_unitID, N, depthSUM, minDepth);
			for(AVEKmerNum & a: aveKmerNumL){
				a.show();
			}
		 }

		 void showContigCounter(const char *contig){
			 fprintf(stderr, "S1 ");
			for(int j = 0; j < contigSize; j++){
				fprintf(stderr, "%c", contig[j]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "S2 ");
			for(int j = 0; j < kmerN; j++){
				int showBase = MIN(70, contigKmerCounter[j]);
				fprintf(stderr, "%c", '0' + showBase);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "S3 ");
			for(int j = 0; j < kmerN; j++){
				int showBase = MIN(70, contigKmerCounterUNIQ[j]);
				fprintf(stderr, "%c", '0' + showBase);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "S4 ");
			for(int j = 0; j < kmerN_100; j++){
				int showBase = MIN(70, contigKmerCounterK100[j]);
				fprintf(stderr, "%c", '0' + showBase);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "S5 ");
			for(int j = 0; j < kmerN; j++){
				int showBase = MIN(70, unitigIDList[j]);
				if(j < searchBgn || j > searchEnd)
					fprintf(stderr, "#");
				else
					fprintf(stderr, "%c", '0' + showBase);
			}
			fprintf(stderr, "\n");
		 }
	 };

	ContigBaseCounter contigBaseCounter;
	ContigKmerCounter contigKmerCounter;
	std::vector<NodeAnalysisINFO> nfl;
	std::vector<AVEKmerNum> aveKmerNumL;
	std::set<int>uniqUntigSet;

	bool runFilterM2(
			std::vector<DBG_INFO_item *> &id2info, std::string & contig,
			std::vector<std::string> &reads, int curKmerLen,
			std::vector<UNITIG_INFO_item *> &unitigPL,
			bool showRst, bool isDEL, int SV_contig_bg, int SV_contig_ed){
		std::vector<int> readID_l;
		fprintf(stderr, "collect read set for unitig!\n");
		int usize = unitigPL.size();
		int contigSize = contig.size();
		//show path
		if(showRst){
			fprintf(stderr, "\nshow path");
			for(int curUniIdx = 0; curUniIdx < usize; curUniIdx++){
				UNITIG_INFO_item* up = unitigPL[curUniIdx];
				fprintf(stderr, "%d(%ld)---", up->id, up->kmer_list.size());
			}
			fprintf(stderr, "\n");
		}
		int global_pos = 0; int globalBgn = 0, globalEnds = 0;
		int MAGIC_NUM = 300;
		int edge_searchBase = 200;
		for(int curUniIdx = 0; curUniIdx < usize; curUniIdx++){
			UNITIG_INFO_item* up = unitigPL[curUniIdx];
			int search_begin = 0; int search_end = up->kmer_list.size();
			if(curUniIdx == 0){//the begin
				search_begin = up->kmer_list.size() - edge_searchBase;
				search_begin = MAX(0,search_begin);
				global_pos = search_begin;
				globalBgn = global_pos;
			}if(curUniIdx == usize - 1){//the end
				search_end = MIN(search_end, edge_searchBase);
				globalEnds = global_pos + search_end;
			}
			for(int i = search_begin; i < search_end; i += 1){
				for(ReadPos &r:id2info[up->kmer_list[i]]->kcp->NGS_rp)
					readID_l.emplace_back((r.readID << 16) + (global_pos - r.readPos + MAGIC_NUM));
				global_pos++;
			}
		}
		fprintf(stderr, "\nreadID_set is: ");
		std::sort(readID_l.begin(), readID_l.end());
		std::vector<ReadOffset> ro;
		for(int i:readID_l){
			int readID =(i>>16);
			int readOffset = i &0xffff;
			if(ro.empty() || ro.back().readID != readID || ro.back().readOffset != readOffset - MAGIC_NUM){
				ro.emplace_back();
				ro.back().readID = readID;
				ro.back().readOffset = readOffset - MAGIC_NUM;
			}
		}
		std::sort(ro.begin(), ro.end(), ReadOffset::cmp_by_position);
		contigBaseCounter.clear(contigSize);

		const char *contig_str = contig.c_str();
		for(ReadOffset & cro:ro)
			contigBaseCounter.addToContigCounter(contig_str, cro,
					reads[cro.readID].c_str(), reads[cro.readID].size());

		if(showRst){
			if(false)
				for(ReadOffset & cro:ro)
					showRead(contig_str, cro, reads[cro.readID].c_str(), reads[cro.readID].size(), contigSize);
			contigBaseCounter.showContigCounter(contig_str);
		}

		int MIN_BRANCH_READ_NUM = 3;
		if(isDEL)
			contigKmerCounter.clear(contigSize, curKmerLen, globalBgn, globalEnds, 110);
		else
			contigKmerCounter.clear(contigSize, curKmerLen, globalBgn, globalEnds, 90);
		for(ReadOffset & cro:ro){
			bool readIsUsed = true;
			const char *rstr = reads[cro.readID].c_str();
			int rlen = reads[cro.readID].size();
			for(int i = 0; i < rlen; i++){
				int binBase = getBinBase(rstr[i]);
				if(cro.readOffset + i >= 0 && cro.readOffset + i < contigSize && rstr[i] != contig_str[cro.readOffset + i]){
					if(binBase <= 3 && contigBaseCounter.baseN[binBase][cro.readOffset + i] >= MIN_BRANCH_READ_NUM){
						readIsUsed = false;
						break;
					}
				}
			}
			if(readIsUsed){
				contigKmerCounter.addSignal(cro, rlen, false);
				if(showRst)	showRead(contig_str, cro, rstr, rlen, contigSize);
			}
		}
	//	contigKmerCounter.addUID(unitigPL);
	//	contigKmerCounter.calculateAveKmerForUnitig(aveKmerNumL);

		if(showRst)
			contigKmerCounter.showContigCounter(contig_str);
		return contigKmerCounter.Kmer100DepthCheckForUnitig(3, isDEL, SV_contig_bg, SV_contig_ed);
	}

	bool runFilter(std::vector<DBG_INFO_item *> &id2info, std::string & contig,
			std::vector<std::string> &reads, int curKmerLen, PathCluster &pc, std::vector<UNITIG_INFO_item *> &unitigPL,
			AB_SuggsetHandler &abH, bool usingHeadNode, bool isDoubleHap,
			int minSVLen, int &new_SV_len,  int old_SV_len, int linearMinDepth,
			float linearMIN_AVE_KMER_DEPTH, bool showRst, bool using_Linear){
		std::vector<int> readID_l;
		bool passFilter = true;
		fprintf(stderr, "collect read set for unitig!\n");
		int usize = unitigPL.size();
		int contigSize = contig.size();
		//show path
		if(showRst){
			fprintf(stderr, "\nshow path");
			for(int curUniIdx = 0; curUniIdx < usize; curUniIdx++){
				UNITIG_INFO_item* up = unitigPL[curUniIdx];
				fprintf(stderr, "%d(%ld)---", up->id, up->kmer_list.size());
			}
			fprintf(stderr, "\n");
		}
		int global_pos = 0; int globalBgn = 0, globalEnds = 0;
		int MAGIC_NUM = 300;
		int edge_searchBase = 200;
		for(int curUniIdx = 0; curUniIdx < usize; curUniIdx++){
			UNITIG_INFO_item* up = unitigPL[curUniIdx];
			int search_begin = 0; int search_end = up->kmer_list.size();
			if(curUniIdx == 0){//the begin
				search_begin = up->kmer_list.size() - edge_searchBase;
				search_begin = MAX(0,search_begin);
				global_pos = search_begin;
				globalBgn = global_pos;
			}if(curUniIdx == usize - 1){//the end
				search_end = MIN(search_end, edge_searchBase);
				globalEnds = global_pos + search_end;
			}
			for(int i = search_begin; i < search_end; i += 1){
				for(ReadPos &r:id2info[up->kmer_list[i]]->kcp->NGS_rp)
					readID_l.emplace_back((r.readID << 16) + (global_pos - r.readPos + MAGIC_NUM));
				global_pos++;
			}
		}
		fprintf(stderr, "\nreadID_set is: ");
		std::sort(readID_l.begin(), readID_l.end());
		std::vector<ReadOffset> ro;
		for(int i:readID_l){
			int readID =(i>>16);
			int readOffset = i &0xffff;
			if(ro.empty() || ro.back().readID != readID || ro.back().readOffset != readOffset - MAGIC_NUM){
				ro.emplace_back();
				ro.back().readID = readID;
				ro.back().readOffset = readOffset - MAGIC_NUM;
			}
		}
		std::sort(ro.begin(), ro.end(), ReadOffset::cmp_by_position);
		contigBaseCounter.clear(contigSize);

		const char *contig_str = contig.c_str();
		for(ReadOffset & cro:ro)
			contigBaseCounter.addToContigCounter(contig_str, cro,
					reads[cro.readID].c_str(), reads[cro.readID].size());

		if(showRst){
			for(ReadOffset & cro:ro)
				showRead(contig_str, cro, reads[cro.readID].c_str(), reads[cro.readID].size(), contigSize);
			contigBaseCounter.showContigCounter(contig_str);
		}

		int MIN_BRANCH_READ_NUM = 3;
		contigKmerCounter.clear(contigSize, curKmerLen, globalBgn, globalEnds, 90);
		for(ReadOffset & cro:ro){
			bool readIsUsed = true;
			const char *rstr = reads[cro.readID].c_str();
			int rlen = reads[cro.readID].size();
			for(int i = 0; i < rlen; i++){
				int binBase = getBinBase(rstr[i]);
				if(cro.readOffset + i >= 0 && cro.readOffset + i < contigSize && rstr[i] != contig_str[cro.readOffset + i]){
					if(binBase <= 3 && contigBaseCounter.baseN[binBase][cro.readOffset + i] >= MIN_BRANCH_READ_NUM){
						readIsUsed = false;
						break;
					}
				}
			}
			if(readIsUsed){
				bool isUniq = contigKmerCounter.addSignal(cro, rlen, true);
				if(showRst && isUniq)	showRead(contig_str, cro, rstr, rlen, contigSize);
			}
		}
		contigKmerCounter.addUID(unitigPL);
		contigKmerCounter.calculateAveKmerForUnitig(aveKmerNumL);

		if(showRst)
			contigKmerCounter.showContigCounter(contig_str);

		//reset depth for each UNITIG
		fprintf(stderr, "pc.signalNum, %5d\n", pc.signalNum);
		if(pc.signalNum > 1){//path is repeat path
			nfl.clear();
			nfl = pc.nfl;//copy
			for(NodeAnalysisINFO &n:nfl)//clear depth
				n.depth = 0;
			//set new depth
			for(AVEKmerNum & a : aveKmerNumL){
				for(NodeAnalysisINFO &n:nfl){
					if(n.up->id == a.uid){
						n.depth += (float)a.depthSUM/a.kmerSize;
						break;
					}
				}
			}
			fprintf(stderr, "New depth is:\n");
			for(NodeAnalysisINFO &n:nfl) n.show();
			AB_SuggsetItem * bsetS = abH.getMaxSuggsetM2(nfl, isDoubleHap, usingHeadNode, pc.bgnUID, pc.endUID, 0.01, 0.3);
			bool isHOM = false;
			bool isFit = abH.getSVLen(bsetS, pc.basicSVLength, pc.stepSVLength, new_SV_len, isHOM, using_Linear);
			if(isFit){
				fprintf(stderr, "New SV len is %d (isHOM %d)\n", new_SV_len, isHOM);
				if(new_SV_len > -minSVLen && new_SV_len < minSVLen)
					passFilter = false;
			}else{
				fprintf(stderr, "No new SVs\n");
				passFilter = false;
			}

		}else{//path is linear path
			uniqUntigSet.clear();
			int bgn_check = 1; int end_check = aveKmerNumL.size() - 1;
			if(aveKmerNumL.size() <= 3){
				bgn_check = 0; end_check = aveKmerNumL.size();
			}
			for(int i = bgn_check ; i < end_check; i++){
				if(uniqUntigSet.find(aveKmerNumL[i].uid) != uniqUntigSet.end()){
					passFilter = false;
					break;
				}else{
					uniqUntigSet.emplace(aveKmerNumL[i].uid);
				}
				if(aveKmerNumL[i].minDepth < 1 ||
						(float)aveKmerNumL[i].depthSUM/aveKmerNumL[i].kmerSize <= linearMIN_AVE_KMER_DEPTH){
					passFilter = false;
					fprintf(stderr, "Linear path fail, at ");
					aveKmerNumL[i].show();
				}
			}
			//fail when too much node is used
			if(uniqUntigSet.size() >= 13){
				passFilter = false;
			}
//			//fail when kmer length is high while is linear deletions
//			if(curKmerLen > 70 && old_SV_len < 0){
//				passFilter = false;
//			}
		}
		return passFilter;
	}
};

struct PathD{
	Path *p1;
	Path *p2;
	int pathID1;
	int pathID2;
	float depthDiffScore;
	int unusedPathN;
	static inline int cmp_by_DiffSum(const PathD &a, const PathD &b){
		//if(a.unusedPathN != b.unusedPathN) return a.unusedPathN < b.unusedPathN;
		return a.depthDiffScore < b.depthDiffScore;
	}

	void show(){
		fprintf(stderr, "pathID1 %d, pathID2 %d depthDiff %.3f unusedPathN %5d \np1->SV_len %d ", pathID1, pathID2, depthDiffScore, unusedPathN, p1->SV_len);
		for(UNITIG_INFO_item* up: p1->unitigPL){
			fprintf(stderr, "%d--",up->id);
		}
		fprintf(stderr, "\np2->SV_len %d ", p2->SV_len);
		for(UNITIG_INFO_item* up: p2->unitigPL){
			fprintf(stderr, "%d--",up->id);
		}
		fprintf(stderr, "\n");
	}
};

#define PathReadFilterMAX_LEVEL 10
#define PathReadFilterMAX_LENGTH 50

//all path cluster begin at specific key node
struct KeyNodePath{
	int beginModeID;
	std::set<int> nodeIDSetOut;
	//std::vector<PathCluster> pc_l;
	std::vector<Path> p_l;
	std::map<uint64_t, int> usedReadPath;
	int allResultNum;
	bool withOtherBranch = false;
	int linearPathN = 0;
	int repeatPathN = 0;

	bool isHOM_ALT = false;
	std::vector<PathD> pathD_l;

	float calFitRateM2(Path *h1, Path *h2,
			std::vector<DBG_INFO_item *> &id2info,  std::vector<UNITIG_INFO_item>& unitig_l, int &unusedPathN, bool showLog){

		int KMER_len = unitig_l[0].unitig.size() - unitig_l[0].kmer_list.size() + 1;
		std::map<int, float> count;
		int totalCount = 0;
//		for(int mode = 0; mode < 2; mode++){
//			Path *ch = (mode == 0)?h1:h2;
//			for(PathNodeCount &c: ch->pnc){
//				std::map<int, float>::iterator it = count.find(c.up->id);
//				totalCount += c.count;
//				if(it == count.end()){
//					count[c.up->id] = c.count;
//				}else{
//					it->second += c.count;
//				}
//			}
//		}

		for(int mode = 0; mode < 2; mode++){
			Path *ch = (mode == 0)?h1:h2;
			for(uint i = 0; i < ch->unitigPL.size(); i++){
				int id = ch->unitigPL[i]->id;
				//if(i >0 && ch->unitigPL[i]->id == ch->unitigPL[i - 1]->id && ch->unitigPL[i]->kmer_list.size() <= 4 )
				//	continue;
				std::map<int, float>::iterator it = count.find(id);
				totalCount += 1;
				if(it == count.end()){
					count[id] = 1;
				}else{
					it->second += 1;
				}
			}
		}

		float totalDiffRate = 0;
		int calNum = 0;
		std::set<int>usedSignal;
		for(int mode = 0; mode < 2; mode++){
			Path *ch = (mode == 0)?h1:h2;
			float preCount = count.find(ch->unitigPL[0]->id)->second;
			for(uint i = 1; i < ch->unitigPL.size(); i++){
				int curID = (ch->unitigPL[i-1]->id << 16) + (ch->unitigPL[i]->id);
				float curCount = count.find(ch->unitigPL[i]->id)->second;
				if(usedSignal.find(curID) == usedSignal.end()){
					usedSignal.emplace(curID);
					float preDepth = ch->unitigPL[i-1]->averageKmerDepth_read_tail;
					float curDepth = ch->unitigPL[i]->averageKmerDepth_read_head;
					calNum += 1;
					float diffCur = ABS_U((preDepth*curCount),(preCount*curDepth))/60;///((preDepth*curCount)+(preCount*curDepth));
					totalDiffRate += diffCur;//log((double)diffCur);
					if(showLog)fprintf(stderr, "[%d(%f: %.4f %ld)--%d(%f: %.4f %ld): %05f]",
							ch->unitigPL[i - 1]->id, preCount, preDepth, ch->unitigPL[i - 1]->kmer_list.size(),
							ch->unitigPL[i]->id, curCount, curDepth, ch->unitigPL[i]->kmer_list.size(), diffCur);
					//}
				}
				preCount = curCount;
			}
			if(showLog)fprintf(stderr, "%d\n",ch->SV_len);
		}

		//L3
		std::set<uint64_t> curUsedReadPath;
		//used read path filter:
		for(int mode = 0; mode < 2; mode++){
			Path *ch = (mode == 0)?h1:h2;
			for(uint i = 1; i < ch->unitigPL.size(); i++){
				if(i >= 1){
					uint64_t link2 = ch->unitigPL[i - 1]->id; link2 <<= 10;
							 link2+= ch->unitigPL[i    ]->id;
					curUsedReadPath.emplace(link2);
				}
				//L3~6
				for(uint linkMode = 3; linkMode <= PathReadFilterMAX_LEVEL; linkMode++){
					if(i >= linkMode - 1 && KMER_len <= 60){
						int middleSUM = 0;
						for(uint j = 0;j < linkMode - 2; j++){
							middleSUM += ch->unitigPL[i - j - 1]->kmer_list.size();
						}
						uint64_t linkNum = 0;
						if(linkMode <= 6){
							for(int j = linkMode - 1;j >= 0 ; j--){
								linkNum <<= 10;
								linkNum += ch->unitigPL[i - j]->id;
							}
						}else{
							for(int j = linkMode - 1;j >= 0 ; j--){
								linkNum ^= (linkNum << 3); linkNum += ch->unitigPL[i - j]->id;
							}
						}
						curUsedReadPath.emplace(linkNum);
					}
				}
			}
		}
		unusedPathN = 0;
		for(std::map<uint64_t, int>::iterator it = usedReadPath.begin(); it != usedReadPath.end(); it++){
			if(curUsedReadPath.find(it->first) == curUsedReadPath.end()){
				if(showLog){
					fprintf(stderr, "%ld ", it->first);
				}
				unusedPathN+= it->second;
			}
		}
		if(showLog)
			fprintf(stderr, "\n");

		//bonus for HOM-SV
		//fprintf(stderr, "unusedPathN %5d\n", unusedPathN);
		float score = totalDiffRate/calNum;
		if(ABS_U(h1->SV_len, h2->SV_len) < 5){//same length
			score = score*0.7;
		}else if(ABS(h1->SV_len) < 5 && ABS(h2->SV_len) < 5){//0/0
			score = score*0.6;
		}else if(ABS(h1->SV_len) < 5 || ABS(h2->SV_len) < 5){//0/1
			score = score*0.85;
		}
		if((h1 == h2 && calNum > 15) || (calNum > 25))
			score = score*4;//penalty for too-complex paths
		score += (0.03*unusedPathN);
		return score;
	}
	bool run_length_recover_ref(int ref_begin, int ref_end, std::vector<uint8_t> &RunLength_count_ref, std::string & reference_str, std::string &ref_str){
		for(int i = ref_begin; i < ref_end; i++){
			char base_char = reference_str[i];
			if(base_char < 'A' || base_char > 't')	return false;
			for(int j = 0; j < RunLength_count_ref[i]; j++){
				ref_str += base_char;
			}
		}
		return true;
	}
	void run_length_recover_contig(bool print_log,int contig_begin, int contig_end, std::vector<ReadPos> &support_read_l,
			std::vector<std::vector<uint8_t>> &RunLength_count_read, std::vector<std::string> &reads,
			 std::string &alt_str ){
		//float read_s = support_read_l.size();
		for(int i = contig_begin; i < contig_end; i++){
			int total_base_count = 0;
			char base_char = 'N';
			int read_s = 0;
			for(ReadPos & r: support_read_l){
				if(r.readPos + i >= (int)RunLength_count_read[r.readID].size()){
					continue;
				}

				if(r.readPos + i < 0){
					total_base_count += 0;
				}else{
					read_s++;
					total_base_count += RunLength_count_read[r.readID][r.readPos + i];
				}

				if(r.readPos + i < 0){
					//do nothing
				}
				else if(base_char == 'N')//blank
					base_char = reads[r.readID][r.readPos + i];
				else if(base_char != reads[r.readID][r.readPos + i]){
					fprintf(stderr, "WARNNING: Contig base ERROR!");
				}
				if(base_char < 'A' || base_char > 'Z'){
					fprintf(stderr, " ");
				}
				if(false) fprintf(stderr, "[%d:%d:%c:%d] ", r.readID, r.readPos + i, reads[r.readID][r.readPos + i], RunLength_count_read[r.readID][r.readPos + i]);
			}
			if(read_s == 0) break;
			float ave_count = (float)total_base_count/read_s;
			int ave_count_int = (int)(ave_count + 0.5);
			if(false) fprintf(stderr, " %f %d\n", ave_count, ave_count_int);
			for(int i = 0; i < ave_count_int; i++){
				alt_str += base_char;
			}
		}
		if(print_log) fprintf(stderr, "ALT: %s\n", alt_str.c_str());
	}

	int main_path_selection(bool print_log, std::vector<DBG_INFO_item *> &id2info, int &total_read_n, int &unknown_read_n){
		int final_path_number = 0;
		if(print_log){
			fprintf(stderr, "\nShow All Path list: path num is %ld, begin at %d, end at ",p_l.size(), beginModeID);
			for(int i : nodeIDSetOut){
				fprintf(stderr, "%d\t", i);
			}
			fprintf(stderr, "\n The path list is below:" );
			for(Path &pp : p_l){
				pp.showSim(id2info);
			}
		}
		if(print_log) fprintf(stderr, "\n The main_path_selection begin\n" );

		std::sort(p_l.begin(), p_l.end(), Path::cmp_by_SV_LEN);
		total_read_n = 0;
		unknown_read_n = 0;

		if( p_l.size() > 0){
			//INIT
			for(uint i = 0; i < p_l.size(); i++){
				if(p_l[i].support_read_l.size() >= 3){
					p_l[i].is_final_path = true;
				}else{
					p_l[i].is_final_path = false;
				}
			}

			int cur_SV_len = p_l[0].SV_len;
			int max_read_size = p_l[0].support_read_l.size();
			int max_path_idx = 0;

			for(uint i = 1; i < p_l.size(); i++){
				Path *ch = &(p_l[i]);
				if(ABS(cur_SV_len - ch->SV_len) < 10){//similar length
					if(max_read_size < (int)p_l[i].support_read_l.size()){
						max_read_size = p_l[i].support_read_l.size();
						max_path_idx = i;
					}
				}else{
					p_l[max_path_idx].is_final_path = true;
					cur_SV_len = ch->SV_len;
					max_read_size = p_l[i].support_read_l.size();
					max_path_idx = i;
				}
			}
			p_l[max_path_idx].is_final_path = true;

			//count the path number the first time:
			for(uint i = 0; i < p_l.size(); i++)
				if(p_l[i].is_final_path)
					final_path_number ++;

			//select the first 2 MAX path when path number is over 2
			if(final_path_number > 2){//re-select when with too
				//count total read number
				fprintf(stderr, "The path number is over size ORI: \n");
				for(uint i = 0; i < p_l.size(); i++){
					Path *ch = &(p_l[i]);
					fprintf(stderr, "%d %ld %d\n", ch->SV_len, p_l[i].support_read_l.size(),p_l[i].is_final_path);
				}
				std::sort(p_l.begin(), p_l.end(), Path::cmp_by_support_read_number);
				for(uint i = 0; i < p_l.size(); i++){
					if(i == 0 || i == 1)	p_l[i].is_final_path = true;
					else					p_l[i].is_final_path = false;
				}
				//count total read number
				fprintf(stderr, "The path number is over size: \n");
				for(uint i = 0; i < p_l.size(); i++){
					Path *ch = &(p_l[i]);
					fprintf(stderr, "%d %ld %d\n", ch->SV_len, p_l[i].support_read_l.size(),p_l[i].is_final_path);
				}
				fprintf(stderr, "The path number is over size END\n");
			}
			final_path_number = 0;
			//count total read number
			for(uint i = 0; i < p_l.size(); i++){
				Path *ch = &(p_l[i]);
				total_read_n += p_l[i].support_read_l.size();
				if(!p_l[i].is_final_path)
					unknown_read_n += p_l[i].support_read_l.size();
				else
					final_path_number ++;
				fprintf(stderr, "%d %ld %d\n", ch->SV_len, p_l[i].support_read_l.size(),p_l[i].is_final_path);
			}

		}


		if(print_log) fprintf(stderr, "\n The main_path_selection END:  total_read_n %d, unknown_read_n %d\n", total_read_n, unknown_read_n);
		return final_path_number;
	}


	int main_path_selection_M1(bool print_log, std::vector<DBG_INFO_item *> &id2info, int &total_read_n, int &unknown_read_n){
		int final_path_number = 0;
		if(print_log){
			fprintf(stderr, "\nShow All Path list: path num is %ld, begin at %d, end at ",p_l.size(), beginModeID);
			for(int i : nodeIDSetOut){
				fprintf(stderr, "%d\t", i);
			}
			fprintf(stderr, "\n The path list is below:" );
			for(Path &pp : p_l){
				pp.showSim(id2info);
			}
		}
		if(print_log) fprintf(stderr, "\n The main_path_selection begin\n" );

		std::sort(p_l.begin(), p_l.end(), Path::cmp_by_SV_LEN);
		total_read_n = 0;
		unknown_read_n = 0;

		if( p_l.size() > 0){
			//INIT
			for(uint i = 0; i < p_l.size(); i++){
				if(p_l[i].support_read_l.size() >= 3){
					p_l[i].is_final_path = true;
				}else{
					p_l[i].is_final_path = false;
				}
			}

			int cur_SV_len = p_l[0].SV_len;
			int max_read_size = p_l[0].support_read_l.size();
			int max_path_idx = 0;

			for(uint i = 1; i < p_l.size(); i++){
				Path *ch = &(p_l[i]);
				if(ABS(cur_SV_len - ch->SV_len) < 10){//similar length
					if(max_read_size < (int)p_l[i].support_read_l.size()){
						max_read_size = p_l[i].support_read_l.size();
						max_path_idx = i;
					}
				}else{
					p_l[max_path_idx].is_final_path = true;
					cur_SV_len = ch->SV_len;
					max_read_size = p_l[i].support_read_l.size();
					max_path_idx = i;
				}
			}
			p_l[max_path_idx].is_final_path = true;
			for(uint i = 0; i < p_l.size(); i++){
				Path *ch = &(p_l[i]);
				total_read_n += p_l[i].support_read_l.size();
				if(!p_l[i].is_final_path)
					unknown_read_n += p_l[i].support_read_l.size();
				else
					final_path_number ++;
				fprintf(stderr, "%d %ld %d\n", ch->SV_len, p_l[i].support_read_l.size(),p_l[i].is_final_path);
			}
		}
		if(print_log) fprintf(stderr, "\n The main_path_selection END:  total_read_n %d, unknown_read_n %d\n", total_read_n, unknown_read_n);
		return final_path_number;
	}

	void analysisAndGenerateSV_TGS(
			 bool print_log,
			 int region_ID_global,
		std::vector<uint8_t> &RunLength_count_ref,
		std::vector<std::vector<uint8_t>> &RunLength_count_read,
		std::vector<NOVA_SV_FINAL_RST_item> &SVE_SVs,
		std::vector<DBG_INFO_item *> &id2info, std::vector<std::string> &reads, AB_SuggsetHandler &abH,
		int curKmerLen,
		std::string & reference_str, Contig_String_aligner *ca, RefRegion &refRegion,
		std::vector<UNITIG_INFO_item>& unitig_l, Genotyping_read_aligner *gra, Bam_file *bam_f,
		bool is_clip, bool output_small_var, int MIN_sv_len
	){
		//print_log = true;
		int total_read_n;//total read number in all path
		int unknown_read_n;//total read number NOT in ANY MAIN path
		//depth analysis
		//generate the result of out node
		//contig selection and remove duplications
		int final_path_number = main_path_selection(print_log, id2info, total_read_n, unknown_read_n);

		//generate the SVs
		int haplotype_ID = -1;
		for(uint i = 0; i < p_l.size(); i++){
			Path *ch = &(p_l[i]);
			if(p_l[i].is_final_path == false)
				continue;
			haplotype_ID ++;
			std::string contig;
			if(print_log){
				fprintf(stderr, "Handle Hap\n");
			}
			if(ABS(ch->SV_len) < 15)
				continue;

			int supp_read_n = ch->support_read_l.size() ;//total read number NOT in ANY MAIN path
			int not_supp_read_n = total_read_n - supp_read_n - unknown_read_n;
			//get the reads position to each path
			{
				std::vector<ReadPos> read_contig_pos;
				ch->getReadPosition_to_contig(id2info, read_contig_pos);
				if(false) {for(ReadPos & r: read_contig_pos) r.show(stderr);}
				for(ReadPos & sr: ch->support_read_l){
					for(ReadPos & rcp: read_contig_pos){
						if(sr.readID == rcp.readID){
							sr.readPos = rcp.readPos;
							break;
						}
					}
				}
				for(ReadPos & r: ch->support_read_l){
					r.show(stderr);
				}
			}

			std::vector<uint8_t> bin_contig;
			std::vector<uint8_t> bin_ref;
			std::string full_contig;
			std::string full_ref_str;

			ch->getContig(id2info, contig);
			if(print_log) ch->showContig(id2info, contig);
			int refPosBgn; int refPosEnd;
			ch->getReferenceRegion(curKmerLen, refPosBgn, refPosEnd);
			if(print_log)
				ch->showReference(refPosBgn,refPosEnd, reference_str);
			//run alignment
			//store contig to binary format
			if(print_log)fprintf(stderr, "Alignment Begin\n");
			run_length_recover_contig(print_log, 0, contig.size(), ch->support_read_l, RunLength_count_read, reads, full_contig);
			if(full_contig.empty())
				continue;
			store_bin_contig(full_contig, bin_contig);
			if(false == run_length_recover_ref(refPosBgn, refPosEnd, RunLength_count_ref, reference_str, full_ref_str))
				fprintf(stderr, "Warning: Ref length not enough!\n");
			store_bin_contig(full_ref_str, bin_ref);
			ca->set_tseq(&(bin_ref[0]));

			uint32_t rl = bin_ref.size();
			uint32_t cl = bin_contig.size();
			bool is_aln_super_length = ((ABS_U(rl, cl) > 1.5*(MIN(rl, cl))) && is_clip);
			int addition_load = ABS_U(rl, cl) * 1.2;
			if(is_aln_super_length){
				ca->align_non_splice_super_long(bin_contig, bin_ref,addition_load);
			}else{//middle
				ca->align_non_splice(&(bin_contig[0]), bin_contig.size(), 0, bin_ref.size(), 0);
			}

			fprintf(stderr, "CIGAR before adjust:\n");
			ca->printCIGAR(stderr);
			ca->show_score();
			int refPosBgn_LOCAL  = 0;
			refPosBgn_LOCAL += ca->adjustCIGAR();
			//ca->printCIGAR(stderr);
			if(print_log)ca->printf_alignment_detail(stderr, refPosBgn_LOCAL, NULL, 0);

			if(is_aln_super_length){
				bool right_aln = ca->align_non_splice_super_long_adjust_cigar(addition_load);
				if(!right_aln)
					continue;
				ca->printCIGAR(stderr);
			}

			//store SVs as result vector, then convert to SVE_SVs
			std::vector<SV_REGION_TGS> sv_r;
			//truSV_num is the number of SVs that not shorter than MIN_sv_len
			int truSV_num = ca->get_canditate_SVs_TGS(false, MIN_sv_len, refPosBgn_LOCAL, sv_r, output_small_var, full_ref_str, full_contig);
			if(truSV_num != 0){
				//generate SVs
				int SV_in_HAP_ID = -1;
				//S1: ref position
				uint true_ref_position = 0;
				for(int i = 0; i < refPosBgn; i++){
					true_ref_position += RunLength_count_ref[i];
				}
				for(SV_REGION_TGS & r: sv_r){
					r.show();
					SV_in_HAP_ID++;
					//recover from the running-length compact:
					//for deletion:
					std::string ref_str;
					std::string alt_str;
					if(r.ref_begin < 1)
						continue;

					alt_str += full_ref_str[r.ref_begin - 1];
					ref_str += full_ref_str[r.ref_begin - 1];

					uint32_t* cigar;
					int cigar_num = ca->get_cigar(&cigar);

					if(r.ref_begin != r.ref_end){//deletion
						//run_length_recover_ref(r.ref_begin, r.ref_end, RunLength_count_ref, reference_str, ref_str);
						ref_str += full_ref_str.substr(r.ref_begin, r.ref_end - r.ref_begin);
						//xassert(r.ref_end <= full_ref_str.size(),"");

						if(!ref_str.empty() && !alt_str.empty()){
							NOVA_SV_FINAL_RST_item::add_to_vector(SVE_SVs, ca->chr_ID, true_ref_position + r.ref_begin + ca->global_ref_pos, "DEL", &(ref_str[0]), &(alt_str[0]),
									alt_str.size() - ref_str.size(), &(bin_contig[0]), bin_contig.size(),
									cigar, cigar_num, r.cigar_idx, r.cigar_idx, refPosBgn_LOCAL, refRegion.st_pos + true_ref_position);
							SVE_SVs.back().setTGS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, final_path_number, SV_in_HAP_ID, sv_r.size());
						}
					}
					else if(r.contig_begin != r.contig_end){//insertion
						alt_str += full_contig.substr(r.contig_begin, r.contig_end - r.contig_begin);
						if(!ref_str.empty() && !alt_str.empty()){
							//run_length_recover_contig(print_log, r.contig_begin, r.contig_end, ch->support_read_l, RunLength_count_read, reads, alt_str);
							NOVA_SV_FINAL_RST_item::add_to_vector(SVE_SVs, ca->chr_ID, true_ref_position + r.ref_begin + ca->global_ref_pos, "INS", &(ref_str[0]), &(alt_str[0]),
									alt_str.size() - ref_str.size(), &(bin_contig[0]), bin_contig.size(),
									cigar, cigar_num, r.cigar_idx, r.cigar_idx, refPosBgn_LOCAL, refRegion.st_pos + true_ref_position);
							SVE_SVs.back().setTGS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, final_path_number, SV_in_HAP_ID, sv_r.size());
						}
					}
					//S3: insertion string
					else{//SNPs
						alt_str.clear(); alt_str += full_contig[r.contig_begin];
						ref_str.clear(); ref_str += full_ref_str[r.ref_begin];
						if(!ref_str.empty() && !alt_str.empty()){
							//run_length_recover_contig(print_log, r.contig_begin, r.contig_end, ch->support_read_l, RunLength_count_read, reads, alt_str);
							NOVA_SV_FINAL_RST_item::add_to_vector(SVE_SVs, ca->chr_ID, true_ref_position + r.ref_begin + ca->global_ref_pos + 1, "SNP", &(ref_str[0]), &(alt_str[0]),
									alt_str.size() - ref_str.size(), NULL, 0,
									NULL, 0, 0, 0, 0, 0);
							SVE_SVs.back().setTGS_INFO(region_ID_global, supp_read_n, not_supp_read_n, unknown_read_n, haplotype_ID, final_path_number, SV_in_HAP_ID, sv_r.size());
						}
					}
				}
			}
		}
		return;
	}

	void analysisAndGenerateSV_M2(
			std::vector<NOVA_SV_FINAL_RST_item> &SVE_SVs,
			std::vector<DBG_INFO_item *> &id2info, std::vector<std::string> &reads, AB_SuggsetHandler &abH,
			int curKmerLen,
			std::string & reference_str, Contig_String_aligner *ca, RefRegion &refRegion, bool showLogs,
			std::vector<UNITIG_INFO_item>& unitig_l, Genotyping_read_aligner *gra, Bam_file *bam_f){
		//depth analysis
		//generate the result of out node
		PathFilter pf;
		int MIN_SV_LEN = 30;
		isHOM_ALT = false;
		fprintf(stderr, "\nShow All Path list: path num is %ld, begin at %d, end at ",p_l.size(), beginModeID);
		for(int i : nodeIDSetOut){
			fprintf(stderr, "%d\t", i);
		}
		fprintf(stderr, "\n The path list is below:" );

		std::vector<Path*> zeroPath_l;
		for(Path &pp : p_l){
			pp.showSim(id2info);
		}

		for(Path &pp : p_l){
			if(ABS(pp.SV_len) <=5)
				zeroPath_l.emplace_back(&pp);
		}

		//combine
		int unusedPathN = 0;
		for(uint i = 0; i < p_l.size(); i++){
			for(uint j = i; j < p_l.size(); j++){
				float totalDiffRateScore = calFitRateM2(&(p_l[i]),&(p_l[j]), id2info,unitig_l, unusedPathN, false);
				pathD_l.emplace_back();
				pathD_l.back().depthDiffScore = totalDiffRateScore;
				pathD_l.back().p1 = &(p_l[i]);
				pathD_l.back().p2 = &(p_l[j]);
				pathD_l.back().pathID1 = i;
				pathD_l.back().pathID2 = j;
				pathD_l.back().unusedPathN = unusedPathN;
			}
		}
		std::sort(pathD_l.begin(), pathD_l.end(), PathD::cmp_by_DiffSum);
		int maxOutput = 0;
		for(PathD & pd: pathD_l){
			pd.show();
			calFitRateM2(pd.p1,pd.p2, id2info,unitig_l, unusedPathN, true);
			if(maxOutput ++ > 3 && ! showLogs){
				break;
			}
		}
		if(pathD_l.empty()){
			return;
		}

		bool resultIsPassFilter = true;

		int SVQUAL = 0;
		if(pathD_l.size() == 1){
			if(pathD_l[0].depthDiffScore > 0.60)
				resultIsPassFilter = false;
			SVQUAL = 30;
		}else if(pathD_l.size() <= 3){
			if(pathD_l[0].depthDiffScore > 0.40)
				resultIsPassFilter = false;
			SVQUAL = 30;
		}else{
			if(pathD_l[0].depthDiffScore > 0.3)
				resultIsPassFilter = false;
			float S1 = pathD_l[0].p1->SV_len; S1 = ABS(S1);
			float S3 = pathD_l[0].p2->SV_len; S3 = ABS(S3);
			if(S1 > S3) std::swap(S1, S3);
			int NOT_same_rst_ID = -1;
			for(uint i = 1;  i < pathD_l.size(); i++){
				float S2 = pathD_l[i].p1->SV_len; S2 = ABS(S2);
				float S4 = pathD_l[i].p2->SV_len; S4 = ABS(S4);
				if(S2 > S4) std::swap(S2, S4);
				bool similarRst1 = ((S1 >= 0.8*S2 && S1 <= 1.2*S2) || (S1 >= S2 - 10 && S1 <= S2 + 10));
				bool similarRst2 = ((S3 >= 0.8*S4 && S3 <= 1.2*S4) || (S3 >= S4 - 10 && S3 <= S4 + 10));
				bool similarRst = (similarRst1 && similarRst2);
				if(!similarRst){
					NOT_same_rst_ID = i;
					break;
				}
			}
			if(NOT_same_rst_ID == -1)//all is similar result
				SVQUAL = 30;
			else{//with other result:
				float scoreDiff = pathD_l[NOT_same_rst_ID].depthDiffScore/pathD_l[0].depthDiffScore;
				if(scoreDiff < 1.1)		 SVQUAL = 0;
				else if(scoreDiff < 1.2) SVQUAL = 3;
				else if(scoreDiff < 1.5) SVQUAL = 5;
				else if(scoreDiff < 2)	 SVQUAL = 15;
				else					 SVQUAL = 30;
			}
		}

		if(SVQUAL == 0){
			fprintf(stderr, "Result is filter out because low QAUL\n");
			resultIsPassFilter = false;
		}

		if(ABS(pathD_l[0].p1->SV_len) < MIN_SV_LEN && ABS(pathD_l[0].p2->SV_len) < MIN_SV_LEN){
			resultIsPassFilter = false;
		}
//		if(pathD_l[0].p1->unitigPL.size() > 20 || pathD_l[0].p2->unitigPL.size() > 20){
//			resultIsPassFilter = false;
//		}

		if(!resultIsPassFilter)
			return;
		bool isHOM = false;
		float LEN_DIFF1 = pathD_l[0].p1->SV_len - pathD_l[0].p2->SV_len;
		LEN_DIFF1 = ABS(LEN_DIFF1);
		float LEN_DIFF2 = ((float)pathD_l[0].p1->SV_len + 0.2)/((float) pathD_l[0].p2->SV_len + 0.1 );
		LEN_DIFF2 = ABS(LEN_DIFF2);
		if(LEN_DIFF2 > 1) LEN_DIFF2 = 1/LEN_DIFF2;

		if(LEN_DIFF1 < 5 || LEN_DIFF2 > 0.9){
			isHOM = true;
		}

		for(int mode = 0; mode < 2; mode++){
			Path *ch = (mode == 0)?pathD_l[0].p1:pathD_l[0].p2;
			std::string contig;
			if(showLogs){
				fprintf(stderr, "Handle Hap\n");
			}
			if(ABS(ch->SV_len) < MIN_SV_LEN)
				continue;
			if(isHOM && mode == 1)
				continue;
			ch->getContig(id2info, contig);
			if(showLogs) ch->showContig(id2info, contig);
			int refPosBgn; int refPosEnd;
			ch->getReferenceRegion(curKmerLen, refPosBgn, refPosEnd);
			if(showLogs) ch->showReference(refPosBgn,refPosEnd, reference_str);
			//run alignment
			//store contig to binary format
			if(showLogs)fprintf(stderr, "Alignment Begin\n");
			std::vector<uint8_t> bin_contig;
			store_bin_contig(contig, bin_contig);
			ca->align_non_splice(&(bin_contig[0]), bin_contig.size(), refPosBgn, refPosEnd, 0);
			refPosBgn += ca->adjustCIGAR();
			ca->printCIGAR(stderr);
			if(showLogs)ca->printf_alignment_detail(stderr, refPosBgn, NULL, 0);
			int match_size;
			int contig_NM = ca->get_contig_NM(refPosBgn, match_size);
			if(showLogs) fprintf(stderr, "Contig filter: contig_NM %d match_size %d \n", contig_NM, match_size);
			if(contig_NM * 20 > match_size){
				fprintf(stderr, "Contig filter: Too mant NM: contig_NM %d error rate: %f%%(over 5%%)\n",
						contig_NM, (float)contig_NM/match_size*100);
				continue;
			}
			std::vector<int> suggest_SV_length;
			suggest_SV_length.emplace_back(ch->SV_len);
			int old_SVE_SVs_size = SVE_SVs.size();
			ca->get_canditate_SVs(false, SVE_SVs, 30, refPosBgn, suggest_SV_length, refRegion.st_pos);

			uint32_t SV_final_size = old_SVE_SVs_size;
			while(SV_final_size < SVE_SVs.size()){
				auto &c_sv = SVE_SVs[SV_final_size];
				c_sv.setGenotype(isHOM);
				c_sv.setVNTR_QUAL(SVQUAL);
				//kmer filter:
				//SV_type_int == 0: INS ; SV_type_int == 1:DEL
				int SV_contig_bg = 0; int SV_contig_ed = 0; int SV_type_int = 0;
				c_sv.get_sv_ST_EN_in_contig(SV_contig_bg, SV_contig_ed, SV_type_int);

				//run additional read filter for find results
				bool passFilter = pf.runFilterM2(id2info, contig, reads, curKmerLen, ch->unitigPL,
						showLogs, SV_type_int == 1, SV_contig_bg, SV_contig_ed);
				if(!passFilter){
					fprintf(stderr, "Result is filter out because low depth\n");
					SVE_SVs[ SV_final_size].passVNTR_depthFilter = false;
					//SVE_SVs.erase(SVE_SVs.begin() + SV_final_size);
				}
				{
					SV_final_size++;
				}
				if(false){
					//this filter is not used
					for(Path * zeroP:zeroPath_l){
						//force calling begin:
						fprintf(stderr, "Force calling Begin\n");
						showLogs = true;
						std::string MateContig;
						zeroP->getContig(id2info, MateContig);
						//if(showLogs) zeroP->showContig(id2info, MateContig);
						//int refPosBgn; int refPosEnd;
						std::vector<uint8_t> bin_contig_mate;
						store_bin_contig(MateContig, bin_contig_mate);
						gra->setRef(&(bin_contig_mate[0]), bin_contig_mate.size(), &(contig[0]), bin_contig.size());
						if(showLogs){
							fprintf(stderr, "RefString is: ");
							for(uint i = 0; i < bin_contig_mate.size(); i++)
								fprintf(stderr, "%c", "ACGT"[bin_contig_mate[i]]);
							fprintf(stderr, "\nContig String is: ");
							for(uint i = 0; i < bin_contig.size(); i++)
								fprintf(stderr, "%c", "ACGT"[bin_contig[i]]);
							fprintf(stderr, "\n");
						}

						bool region_is_overlap = false;
						//load reference in break point 1:
						int edge_len = 148;
						int st_POS = c_sv.st_pos;
						int ed_POS = c_sv.st_pos + (c_sv.SV_length>=0)?1:(-c_sv.SV_length);
						RefRegion main(c_sv.chr_ID, st_POS - 10, st_POS + edge_len);
						RefRegion supp(c_sv.chr_ID, ed_POS - 10, ed_POS + edge_len);
						if(main.region_overlap(supp)){ main.Combine(supp, true); region_is_overlap = true; }
						std::cerr << "Region 1: " << main;
						if(!region_is_overlap) std::cerr << "Region 2: " << supp;
						std::cerr << std::endl;
						uint8_t qseq_buff[1024]; int left_clip = 0;
						fprintf(stderr, "SV breakpoint1 %d SV; breakpoint2 %d\n", st_POS, ed_POS);

						//suggest contig_region_st for BP1/2
						int contig_pos_bp1; int contig_pos_bp2;
						c_sv.get_contig_global_position(&contig_pos_bp1, &contig_pos_bp2);
						int region_support_number[3][3];
						for(int region_ID = 0; region_ID < 3; region_ID++)
							for(int read_genotype_type = 0; read_genotype_type < 3; read_genotype_type++)
								region_support_number[region_ID][read_genotype_type] = 0;

						for(int region_ID = 0; region_ID < 2; region_ID++){ //mode == 0 for region1, mode == 1 for region2
							fprintf(stderr, "region_ID %d\n", region_ID);
							if(region_is_overlap == true && region_ID == 1) continue;
							R_region region;
							region.chr_ID = main.chr_ID;
							region.st_pos = ((region_ID == 0)?main.st_pos:supp.st_pos) + 1;
							region.ed_pos = ((region_ID == 0)?main.ed_pos:supp.ed_pos) + 1;
	//
							resetRegion_ID(bam_f, &region);	//reset region
	//						//reference check:
							while (bam_next(bam_f)) {
								bam1_t *br = &(bam_f->_brec);
								if(bam_is_secondary(br))		continue;
								if(bam_is_supplementary(br))	continue;
								if(bam_is_duplicate(br))       continue;

								int overlap_mode = c_sv.read_overlap_breakpoint(br, region_ID, !region_is_overlap, br->core.l_qseq, &left_clip);
								if(overlap_mode == 0) {
									if(false) fprintf(stderr, "@Read name: %s Overlap_mode %d; SKIP\n", bam_qname(br), overlap_mode);
									continue;
								}

								get_bam_seq_bin(0, br->core.l_qseq, qseq_buff, br);
								//realigned for contig regions
								int read_contig_score = 0;

								int true_read_pos_bg; int true_read_pos_ed;
								c_sv.get_true_read_pos(br, &true_read_pos_bg, &true_read_pos_ed);

								if(false) fprintf(stderr, "\ttrue_read_pos: [%d, %d]\t", true_read_pos_bg, true_read_pos_ed);

								int read_in_contig_st_pos_bp1_rbg = (true_read_pos_bg) - contig_pos_bp1;
								int read_in_contig_st_pos_bp2_rbg = (true_read_pos_bg) - contig_pos_bp2;
								if(showLogs) fprintf(stderr, "\ttrue_read_pos_bg %5d, contig_pos_bp1 %5d, read_in_contig_st_pos_bp1_rbg %5d\t", true_read_pos_bg, contig_pos_bp1, read_in_contig_st_pos_bp1_rbg);
								if(showLogs) fprintf(stderr, "\ttrue_read_pos_bg %5d, contig_pos_bp1 %5d, read_in_contig_st_pos_bp1_rbg %5d\t", true_read_pos_bg, contig_pos_bp2, read_in_contig_st_pos_bp2_rbg);

								int read_in_contig_st_pos_bp1_red = (true_read_pos_ed - br->core.l_qseq) - contig_pos_bp1;
								int read_in_contig_st_pos_bp2_red = (true_read_pos_ed - br->core.l_qseq) - contig_pos_bp2;

								int read_skip_left, read_skip_right;
								int minMismatchQUAL = 12;
								int true_used_read_in_contig = 0;
								if(overlap_mode == 1){
									read_contig_score = c_sv.get_contig_alignment_scoreM2(gra, br, qseq_buff,
											read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red, &read_skip_left,&true_used_read_in_contig, &read_skip_right, minMismatchQUAL);
								}else if(overlap_mode == 2){
									read_contig_score = c_sv.get_contig_alignment_scoreM2(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red,
											&read_skip_left, &true_used_read_in_contig, &read_skip_right, minMismatchQUAL);
								}else if(overlap_mode == 3){
									int read_contig_score1 = c_sv.get_contig_alignment_scoreM2(gra, br, qseq_buff, read_in_contig_st_pos_bp1_rbg, read_in_contig_st_pos_bp1_red,
											&read_skip_left, &true_used_read_in_contig,  &read_skip_right, minMismatchQUAL);
									int read_skip_left_bp2, read_skip_right_bp2; int true_used_read_in_contig2;
									int read_contig_score2 = c_sv.get_contig_alignment_scoreM2(gra, br, qseq_buff, read_in_contig_st_pos_bp2_rbg, read_in_contig_st_pos_bp2_red,
											&read_skip_left_bp2,&true_used_read_in_contig2,  &read_skip_right_bp2, minMismatchQUAL);
									if(read_contig_score1 >= read_contig_score2){
										read_contig_score = read_contig_score1;
									}else{
										read_skip_left = read_skip_left_bp2; read_skip_right = read_skip_right_bp2;
										read_contig_score = read_contig_score2;
										true_used_read_in_contig = true_used_read_in_contig2;
									}
								}
								int read_ori_score = c_sv.get_ori_alignment_score(gra, br, true_read_pos_bg, true_read_pos_bg - true_used_read_in_contig,
										qseq_buff,
										read_skip_left, read_skip_right, minMismatchQUAL, showLogs);

								if(showLogs) {
									fprintf(stderr,"Read_contig_score %d, Read_ori_score is %d,diff is %d ", read_contig_score, read_ori_score, read_contig_score - read_ori_score);
									fprintf(stderr, "@Read name: %s POS: %d Overlap_mode %d;\t", bam_qname(br), br->core.pos, overlap_mode);
									if(true){
										fprintf(stderr, "@Read string: ");
										for(int i = 0; i < br->core.l_qseq; i++)
											fprintf(stderr, "%c", "ACGT"[qseq_buff[i]]);
									}

									fprintf(stderr, "\n");
								}
								int min_score = gra->getGenotypingMinScore(br->core.l_qseq - read_skip_left - read_skip_right);
								if(read_contig_score > read_ori_score + 4 && read_contig_score > min_score)		region_support_number[region_ID][0]++;
								else if(read_contig_score + 4 < read_ori_score && read_ori_score> min_score)	region_support_number[region_ID][1]++;
								else																		region_support_number[region_ID][2]++;
							}
						}
						fprintf(stderr,"Region analysis: for region(1): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[0][0],
								region_support_number[0][1], region_support_number[0][2]);
						fprintf(stderr,"Region analysis: for region(2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[1][0],
								region_support_number[1][1], region_support_number[1][2]);
						fprintf(stderr,"Region analysis: for region(1+2): ALT:REF:UNKNOW [%d, %d, %d]\n", region_support_number[2][0],
								region_support_number[2][1], region_support_number[2][2]);
					}
				}
				//force calling filter:
			}
		}
		return;
	}


	void showBranchNum(std::vector<DBG_INFO_item *> &id2info){
		fprintf(stderr, "Analysis begin: (bgnID %5d)  outNode: ", beginModeID);
		for(int id:nodeIDSetOut){
			fprintf(stderr, "%5d ", id);
		}
		fprintf(stderr, "l_PN %5d, r_PN %5d RstSUM %5d withOtherBranch %s\t", linearPathN , repeatPathN, allResultNum,  withOtherBranch?"YES":"NO ");
	}

	void show(std::vector<DBG_INFO_item *> &id2info){
		showBranchNum(id2info);
	}
};

struct PathReadFilter{
	std::vector<std::map<uint64_t, int>> readPathLink;//link 2/3/4/5/6
	std::vector<std::vector<uint64_t>> readPathLinkBig;

	std::vector<uint64_t> readPath;
	int KMER_len;
	bool isPathFilter(Path &p, std::map<uint64_t, int>& usedReadPath, bool print_log){
		//link 2
		for(uint i = 0; i < p.unitigPL.size(); i++){
			//L2
			if(i >= 1){
				uint64_t link2 = p.unitigPL[i - 1]->id; link2 <<= 10;
						 link2+= p.unitigPL[i    ]->id;
				std::map<uint64_t, int>::iterator it = readPathLink[0].find(link2);
				if(it == readPathLink[0].end() || it->second < 3){
					if(print_log)fprintf(stderr, "link2:i %d\n", i);
					return false;
				}else{
					if(it->second >= 3)
						usedReadPath[link2] = it->second;
				}
			}
			//L3~6
			for(uint linkMode = 3; linkMode <= PathReadFilterMAX_LEVEL; linkMode++){
				std::map<uint64_t, int> &CurLink = readPathLink[linkMode - 2];
				if(i >= linkMode - 1 && KMER_len <= 60){
					int middleSUM = 0;
					for(uint j = 0;j < linkMode - 2; j++){
						middleSUM += p.unitigPL[i - j - 1]->kmer_list.size();
					}
					uint64_t linkNum = 0;
					if(linkMode <= 6){
						for(int j = linkMode - 1;j >= 0 ; j--){
							linkNum <<= 10; linkNum += p.unitigPL[i - j]->id;
						}
					}else{
						for(int j = linkMode - 1;j >= 0 ; j--){
							linkNum ^= (linkNum << 3); linkNum += p.unitigPL[i - j]->id;
						}
					}
					std::map<uint64_t, int>::iterator it = CurLink.find(linkNum);
					if(middleSUM < PathReadFilterMAX_LENGTH && (it == CurLink.end() || it->second < 3)){
						if(print_log)fprintf(stderr, "link%d:i %d\n", linkMode, i);
						return false;
					}
					if(it != CurLink.end() && it->second >= 3)
						usedReadPath[linkNum] = it->second;
				}
			}
		}
		return true;
	}

	void showStringPath(const char * title, std::string &reads, int KMER_len,
			std::unordered_map<std::string, DBG_INFO_item> &kmerIdx,
			std::vector<UNITIG_INFO_item> &unitig_l, bool showLog){

		/// Construct k-mer maps
		/// wordCount: k-mer ==> number of reads containing the k-mer
		/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
		this->KMER_len = KMER_len;
		const unsigned kmerSize = KMER_len;
		// stores the index of a kmer in a read sequence
		const std::string &seq = reads;
		const unsigned seqLen = seq.size();
		// this read is unusable for assembly:
		if (seqLen < kmerSize) return;
		// track all words from the read, including repetitive k
		unsigned kmer_number = seqLen - kmerSize + 1;
		if(showLog)fprintf(stderr, "%s: %s\n", title, reads.c_str());
		for (unsigned j = 0; j < kmer_number; ++j) {
			const std::string word(seq.substr(j, kmerSize));
			//fprintf(stderr, "%s \n",word.c_str());
			std::unordered_map<std::string, DBG_INFO_item>::iterator it = kmerIdx.find(word);
			if(it != kmerIdx.end()){
				fprintf(stderr, "%2d--", it->second.unitigId);
				j += unitig_l[it->second.unitigId].kmer_list.size() - it->second.unitigPos - 1;
			}
			else {
				if(showLog){ fprintf(stderr, "#");}
				//break;
			}
		}
		fprintf(stderr, "\n");
	}

	void generateLink(std::vector<std::string> &reads, int KMER_len,
			std::unordered_map<std::string, DBG_INFO_item> &kmerIdx,
			std::vector<UNITIG_INFO_item> &unitig_l, bool showLog){
		/// Construct k-mer maps
		/// wordCount: k-mer ==> number of reads containing the k-mer
		/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
		if(readPathLink.empty()){
			readPathLink.resize(PathReadFilterMAX_LEVEL - 1);
		}
		//link 2/3/4/5/6
		for(int linkMode = 2; linkMode <= PathReadFilterMAX_LEVEL; linkMode++){
			readPathLink[linkMode - 2].clear();
		}
		readPathLinkBig.clear();

		this->KMER_len = KMER_len;
		const unsigned kmerSize = KMER_len;
		const unsigned strCount = reads.size();
		for (unsigned strIndex = 0; strIndex < strCount; ++strIndex){
			readPath.clear();
			// stores the index of a kmer in a read sequence
			const std::string &seq = reads[strIndex];
			const unsigned seqLen = seq.size();
			// this read is unusable for assembly:
			if (seqLen < kmerSize) continue;
			// track all words from the read, including repetitive k
			unsigned kmer_number = seqLen - kmerSize + 1;
			if(showLog)fprintf(stderr, "read ID %d\t", strIndex);
			for (unsigned j = 0; j < kmer_number; ++j) {
				const std::string word(seq.substr(j, kmerSize));
				//fprintf(stderr, "%s \n",word.c_str());
				std::unordered_map<std::string, DBG_INFO_item>::iterator it = kmerIdx.find(word);
				if(it != kmerIdx.end()){
					//fprintf(stderr, "%2d<%d>--", it->second.unitigId, it->second.unitigPos);
					readPath.emplace_back(it->second.unitigId);
					j += unitig_l[it->second.unitigId].kmer_list.size() - it->second.unitigPos - 1;
				}
				else {
					if(showLog){ fprintf(stderr, "#");}
					break;
				}
			}
			if(showLog){
				for(int uid:readPath)
					fprintf(stderr, "%2d--", uid);
				fprintf(stderr, "\n");
			}
			if(readPath.size() > 1){
				//store link2
				for(uint i = 1; i < readPath.size(); i++){
					uint64_t link2 = (readPath[i - 1] << 10) + readPath[i];
					std::map<uint64_t, int>::iterator it = readPathLink[0].find(link2);
					if(it == readPathLink[0].end())
						readPathLink[0][link2] = 1;
					else
						it->second++;
				}
				//store Link 3~6:
				for(uint linkMode = 3; linkMode <= PathReadFilterMAX_LEVEL; linkMode++){
					std::map<uint64_t, int> &CurLink = readPathLink[linkMode - 2];
					for(uint i = linkMode - 1; i < readPath.size(); i++){
						uint64_t linkNum = 0;
						if(linkMode <= 6){
							for(int j = linkMode - 1;j >= 0 ; j--){
								linkNum <<= 10;
								linkNum += readPath[i - j];
							}
						}else{
							for(int j = linkMode - 1;j >= 0 ; j--){
								linkNum ^= (linkNum << 3); linkNum += readPath[i - j];
							}
						}
						std::map<uint64_t, int>::iterator it = CurLink.find(linkNum);
						if(it == CurLink.end())
							CurLink[linkNum] = 1;
						else
							it->second++;
					}
				}
			}
			if(readPath.size() > 6){
				readPathLinkBig.emplace_back(readPath);
			}
		}
	}

	void show(){
		fprintf(stderr, "Show small links\n");
		for(uint linkMode = 2; linkMode <= 6; linkMode++){
			for(std::map<uint64_t, int>::iterator it = readPathLink[linkMode - 2].begin(); it != readPathLink[linkMode - 2].end(); it++){
				uint64_t linkNum = it->first;
				for(uint i = 0; i< linkMode; i++){
					fprintf(stderr, "%ld--", (linkNum >> (10*(linkMode - i - 1)))&0x3ff);
				}
				fprintf(stderr, "{%5d}\n", it->second);
			}
		}
		fprintf(stderr, "Show big links\n");
		for(std::vector<uint64_t>& rp: readPathLinkBig){
			for(uint64_t r: rp){
				fprintf(stderr, "%ld--", r);
			}
			fprintf(stderr, "\n");
		}
	}
};

struct Kmer_ANA{
	int KMER_N;
	int KMER_ALL;
	int KMER_NEW;
	void clear(){
		KMER_N = 0;
		KMER_ALL = 0;
		KMER_NEW = 0;
	}
	void add_ANA_All(int kmer_number){KMER_ALL += kmer_number;}
	void add_ANA_N(int kmer_number){KMER_N += kmer_number;}
	void add_ANA_New(int kmer_number){KMER_NEW += kmer_number;}
	void show(const char * type){
		if(KMER_ALL == 0)
			return ;
		fprintf(stderr, "Kmer counting %s: (ALL %d : new %d : N %d)\n", type, KMER_ALL, KMER_NEW, KMER_N);
	}
};

int N_count(const std::string &word);

/// Information added to each read in the process of assembly
struct DbgAssemblyReadInformation {
	/// If true, the read was an assembled contig
	int isTGS_read = false;
};

struct MainDBGHandler1 {
	std::vector<std::string> reads_str;
	std::vector<DbgAssemblyReadInformation> read_info;
	//
	//reference handler
	RefRegion refRegion;
	std::vector<std::string> reference;
	void addRef(uint8_t *tseq, int tlen, Contig_String_aligner *ca);
	int curKmerLen;
	Contig_String_aligner *ca;

	KmerCounterTable kmerCountingTable;
	std::unordered_map<std::string, DBG_INFO_item> kmerIdx;

	std::vector<DBG_INFO_item *> id2info;
	std::vector<UNITIG_INFO_item> unitig_l;

	std::vector<KEY_UNITIG> keyPosUnitigL;//unitigs in key position
	std::set<int> keyPosUnitigIDset;

	//BUFFs
	std::vector<int> nextNodeBuff;
	uint64_t random_seed;
	AB_SuggsetHandler abH;

	//path buffs
	std::vector<Path> pathList;
	//std::vector<PathCluster> pathClusterList;
	//
	PathReadFilter prFilter;
	//the SV results
	Genotyping_read_aligner gra;

	/**---------------------------basic function------------------------------------**/
	void init(){ random_seed = get_random_seed_by_time(); abH.init(); gra.init();}
	void clear(){reads_str.clear(); reference.clear(); kmerCountingTable.clear(); }
	void showAllRead(){
		const unsigned strCount = reads_str.size();
		for (unsigned strIndex = 0; strIndex < strCount; ++strIndex){
			fprintf(stderr, "SHOW READs:  %5d %s\n", strIndex, reads_str[strIndex].c_str());;
		}
	}
	/**---------------------------DBG building function------------------------------------**/
	/// Construct k-mer maps
	/// wordCount: k-mer ==> number of reads containing the k-mer
	/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
	/// Construct k-mer maps
	/// wordCount: k-mer ==> number of reads containing the k-mer
	/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
	void kmerCounting(const unsigned kmerSize, bool IsStoreReadInfo) {
		Kmer_ANA k_a[3];
		for (int i = 0; i < 3; i++)
			k_a[i].clear();
		curKmerLen = kmerSize;
		//mode is 0: read; mode is 1: ref; mode is 2: TGS-reads
		for (int mode = 0; mode < 2; mode++) {
			std::vector<std::string> *c_reads = NULL;
			if (mode == 0)
				c_reads = &(reference);
			else
				c_reads = &(reads_str);
			const unsigned strCount = c_reads->size();
			for (unsigned read_idx = 0; read_idx < strCount; ++read_idx) {
				//ref: KMER mode is 0;
				//NGS read: KMER mode is 1;
				//TGS read: KMER mode is 2;
				int kmer_mode = mode;
				if(mode == 1 && read_info[read_idx].isTGS_read)//NGS read
					kmer_mode = 2;
				// stores the index of a kmer in a read sequence
				const std::string &seq = ((*c_reads)[read_idx]);
				const unsigned seqLen = seq.size();
				// this read is unusable for assembly:
				if (seqLen < kmerSize)
					continue;
				// track all words from the read, including repetitive k
				unsigned kmer_number = seqLen - kmerSize + 1;
				k_a[kmer_mode].add_ANA_All(kmer_number);
				for (unsigned kmer_idx = 0; kmer_idx < kmer_number;
						++kmer_idx) {
					const std::string word(seq.substr(kmer_idx, kmerSize));
					//fprintf(stderr, "%s \n",word.c_str());
					// filter words with "N" (either directly from input alignment or marked due to low basecall quality:
					int N_N = N_count(word);
					if (N_N > 0) {
						k_a[kmer_mode].add_ANA_N(1);
						if (false)
							fprintf(stderr, "%s %d \n", word.c_str(), N_N);
						continue;
					}
					KmerCounterTable::iterator it = kmerCountingTable.find(
							word);
					if (it == kmerCountingTable.end()) {
						k_a[kmer_mode].add_ANA_New(1);
						kmerCountingTable[word].add_new_kmer(kmer_mode, kmer_idx, read_idx,
								IsStoreReadInfo);
					} else {
						it->second.add_new_kmer(kmer_mode, kmer_idx, read_idx,
								IsStoreReadInfo);
					}
				}
			}
		}
		k_a[0].show("Ref");
		k_a[0].show("Read NGS");
		k_a[0].show("Read TGS");
	}
	void kmerCountingTest(const unsigned kmerSize);
	/// Construct k-mer maps
	/// wordCount: k-mer ==> number of reads containing the k-mer
	/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
	void errorCorrection(const unsigned kmerSize, int minRight, int maxError);
	//output:
	//	dbg_map_t kmerIdx;//the kmerString to kmer info
	//	id2info: the kmerID to kmer info
	//	unitig_l: the unitig list
	bool buildDBG(int minDepthNGS, int minDepthTGS);
	void showDBG_graph();

	/**---------------------------DBG path search function------------------------------------**/
	std::map<UNITIG_INFO_item *, int > pathNodeCountBuff;
	std::vector<KeyNodePath> pc_l_final;
	void setKmerDepth();
	int randomGetNextNodeID(int curID, bool & withRandomBranch);
	std::vector<UNITIG_INFO_item *> randomSearchBuff;
	bool randomGetPath(int beginNodeID, std::vector<Path> &pathList, int searchPathNum, std::set<int> &stopNodeSet, std::vector<UNITIG_INFO_item *> &randomSearchBuff);
	std::vector<Path> removeRepeatPathBuff;
	void sortAndRemoveDupPath(std::vector<Path> &pathList, bool withRandomPath, int pathNumberPerTime);
	void generatePathCluster(std::vector<Path> &pathList, std::vector<PathCluster> &pathClusterList);

	//deepest-first search the DBG to find all the path in the DBG
	// the counter for how many times the path cross the NODE
	bool D_searchPath(int headID, int curID, std::map<int, int> &node_reach_count, std::set<int> &stopNodeSet,
			std::vector<Path> &pathList, std::vector<UNITIG_INFO_item *> &randomSearchBuff, int MAX_nodePassTime, int &ALL_NODE_USED){
		if(ALL_NODE_USED++ > 30000)
			return false;
		if(pathList.size() > 120)
			return false;
		if(randomSearchBuff.size() > 100)
			return false;
		bool stopAtThisNode = false;
		//stop when reach the BLANK node
		if(curID == -1){
			randomSearchBuff.emplace_back((UNITIG_INFO_item *)NULL);
			stopAtThisNode = true;
		}else
			randomSearchBuff.emplace_back(&(unitig_l[curID]));
		//stop when reach the node without read support
		if(!unitig_l[curID].with_read)
			stopAtThisNode = true;
		//stop when pass the node too many times
		std::map<int, int> ::iterator it = node_reach_count.find(curID);
		if(it == node_reach_count.end()){		node_reach_count[curID] = 1;	}
		else if(it->second++ > MAX_nodePassTime)
			stopAtThisNode = true;
		bool isFullPath = false;
		//stop when reach the node in the END list
		if(curID != headID && stopNodeSet.find(curID) != stopNodeSet.end()){
			stopAtThisNode = true;
			isFullPath = true;
		}
		if(stopAtThisNode){
			if(isFullPath){
				pathList.emplace_back();
				pathList.back().unitigPL = randomSearchBuff;//copy
				pathList.back().endUID = curID;//copy
			}
		}else{
			std::vector<int> &branchUni = unitig_l[curID].outBranchUni;
			for(int outID :branchUni){
				bool notCPX = D_searchPath(headID, outID, node_reach_count, stopNodeSet, pathList, randomSearchBuff, MAX_nodePassTime, ALL_NODE_USED);
				if(!notCPX) return false;
			}

		}
		randomSearchBuff.pop_back();
		node_reach_count[curID]--;
		return true;
	}

	bool getDBGPathM1(bool print_log){
		bool withCPXNode = false;
		pc_l_final.clear();
		std::map<int, int> node_reach_count;
		//for each KEY node
		int keyPosNum = keyPosUnitigL.size();
		for(int i = 0; i < keyPosNum - 1; i++){
			UNITIG_INFO_item *up =  keyPosUnitigL[i].p;
			int beginNodeID = up->id;
			if(unitig_l[beginNodeID].with_read == false){
				continue;
			}
			if(keyPosUnitigIDset.find(beginNodeID) == keyPosUnitigIDset.end())
				continue;
			bool notCPX = true;
			for(int complexMode = 0; complexMode < 4; complexMode++){
				pathList.clear();
				int ALL_NODE_USED = 0;
				node_reach_count.clear();
				randomSearchBuff.clear();
				int MAX_DEEP = 13 - complexMode*4;//13~9~5~1~
				notCPX = D_searchPath(beginNodeID, beginNodeID, node_reach_count,
						keyPosUnitigIDset, pathList, randomSearchBuff, MAX_DEEP, ALL_NODE_USED);
				if(notCPX){
					for(Path & p: pathList){
						p.setSVLen(id2info);
						p.showSim(id2info);
					}
					break;
				}else{
					withCPXNode = true;
					fprintf(stderr, "Path begin at %5d is too complex, using lower deep MAX_DEEP %d\n", beginNodeID, MAX_DEEP - 4);
				}
			}
			if(!notCPX){
				fprintf(stderr, "Path begin at %5d is too complex, SKIP\n", beginNodeID);
				continue;
			}

			//run path filter
			std::vector<Path> pathListTmp;
			std::map<uint64_t, int> usedReadPath;
			std::map<uint64_t, int> usedReadPathTmp;
			for(uint i = 0; i < pathList.size(); i++){
				Path &p = pathList[i];
				//p.show(id2info);
				usedReadPathTmp.clear();
				bool passFilter = prFilter.isPathFilter(p, usedReadPathTmp, print_log);
//				usedReadPathTmp.clear();
//				passFilter = prf.isPathFilter(p, usedReadPathTmp, print_log);
				if(!passFilter){
					if(print_log){
						fprintf(stderr, "Path not pass filter:\n");
						p.show(id2info);
						fprintf(stderr, "Path not pass filter END:\n");
					}
				}
				else{
					//map merging
					for(std::map<uint64_t, int>::iterator it = usedReadPathTmp.begin();it != usedReadPathTmp.end();it++){
						usedReadPath[it->first] = it->second;
					}
					//usedReadPath.merge(usedReadPathTmp);
					pathListTmp.emplace_back();
					std::swap(pathListTmp.back(), p);
				}
			}
			std::swap(pathListTmp, pathList);
			//
			pc_l_final.emplace_back();
			std::swap(pc_l_final.back().usedReadPath, usedReadPath);
			std::swap(pc_l_final.back().p_l, pathList);
			pc_l_final.back().beginModeID = beginNodeID;

			//generate the MAX_out_keyNODE
			int MAX_REF_POS = -1000; int maxRPOS_UID = -1;
			for(Path &pc: pc_l_final.back().p_l){
				if(MAX_REF_POS < unitig_l[pc.endUID].refPos){
					MAX_REF_POS = unitig_l[pc.endUID].refPos;
					maxRPOS_UID = pc.endUID;
				}
				pc_l_final.back().nodeIDSetOut.emplace(pc.endUID);
			}
			if(pc_l_final.back().nodeIDSetOut.size() > 1){
				bool withRemove = false;
				for(int removeID: pc_l_final.back().nodeIDSetOut){
					if(removeID != maxRPOS_UID){
						keyPosUnitigIDset.erase(removeID);
						withRemove = true;
					}
				}
				if(withRemove){
					fprintf(stderr, "All Path Re-calculated begin at %5d\n", beginNodeID);
					pc_l_final.pop_back();
					i --;
				}
			}
		}
		return withCPXNode;
	}

	//deepest-first search the DBG to find all the path in the DBG ONLY using TGS-reads
	void D_searchPath_TGS(int headID, int curID, std::set<int> &stopNodeSet, std::vector<ReadPos> & cur_r_p_l,
			std::vector<Path> &pathList, std::vector<UNITIG_INFO_item *> &randomSearchBuff, int &ALL_NODE_USED){
		bool stopAtThisNode = false;
		std::vector<ReadPos>  ori_r_p_l = cur_r_p_l;//copy

		//stop when reach the BLANK node
		if(curID == -1){
			randomSearchBuff.emplace_back((UNITIG_INFO_item *)NULL);
			stopAtThisNode = true;
			cur_r_p_l.clear();
		}else{
			randomSearchBuff.emplace_back(&(unitig_l[curID]));
			//FOR TGS ONLY: stop when no read support this path
			//TODO::
			std::vector<ReadPos> & TGS_pos = unitig_l[curID].getAllReadTGS(id2info);
			//show all the read TGS
			if(false){for(ReadPos & p: TGS_pos) p.show(stderr);}
			//update the support read list:
			if(curID == headID){//is head
				cur_r_p_l = unitig_l[curID].getAllReadTGSTail(id2info);
			}
			else{
				for(uint i = 0; i < cur_r_p_l.size(); i++){
					ReadPos & rp = cur_r_p_l[i];
					bool truePath = false;
					//todo:: this method is slow
					for(ReadPos & node_rp: TGS_pos){
						//check the node ID
						if(rp.readID == node_rp.readID && rp.readPos == node_rp.readPos){
							truePath = true;
							break;
						}
					}
					if(truePath == false){
						rp.readID = -1;
						cur_r_p_l.erase(cur_r_p_l.begin() + i); i--;
					}
				}
			}
			if(cur_r_p_l.empty())
				stopAtThisNode = true;
			else{
				uint kmer_num = unitig_l[curID].kmer_list.size();
				if(curID == headID){kmer_num = 1;}
				for(ReadPos & node_rp: cur_r_p_l){
					node_rp.readPos += kmer_num;
				}
			}
		}

		bool isFullPath = false;
		//stop when reach the node in the END list
		if(curID != headID && stopNodeSet.find(curID) != stopNodeSet.end() && !cur_r_p_l.empty()){
			stopAtThisNode = true;
			isFullPath = true;
		}

		if(stopAtThisNode){
			if(isFullPath){
				pathList.emplace_back();
				pathList.back().unitigPL = randomSearchBuff;//copy
				pathList.back().endUID = curID;//copy
				pathList.back().support_read_l = cur_r_p_l;//copy
			}
		}else{
			std::vector<int> &branchUni = unitig_l[curID].outBranchUni;
			for(int outID :branchUni){
				//check the branch:
				D_searchPath_TGS(headID, outID, stopNodeSet, cur_r_p_l, pathList, randomSearchBuff, ALL_NODE_USED);
			}
		}
		//restore the read list
		std::swap(ori_r_p_l, cur_r_p_l);//copy
		randomSearchBuff.pop_back();
	}

	void getDBGPathCCS(bool print_log){
		pc_l_final.clear();
		std::map<int, int> node_reach_count;
		//the map to store the current read ID and positions
		std::vector<ReadPos> cur_r_p;
		//for each KEY node
		int keyPosNum = keyPosUnitigL.size();
		for(int i = 0; i < keyPosNum - 1; i++){
			UNITIG_INFO_item *up =  keyPosUnitigL[i].p;
			int beginNodeID = up->id;
			if(unitig_l[beginNodeID].with_read == false){
				continue;
			}
			if(keyPosUnitigIDset.find(beginNodeID) == keyPosUnitigIDset.end())
				continue;
			cur_r_p.clear();
			pathList.clear();
			int ALL_NODE_USED = 0;
			node_reach_count.clear();
			randomSearchBuff.clear();
			//store the begin read list:
			D_searchPath_TGS(beginNodeID, beginNodeID,
					keyPosUnitigIDset, cur_r_p, pathList, randomSearchBuff, ALL_NODE_USED);
			for(Path & p: pathList){
				p.setSVLen(id2info);
				if(print_log) p.showSim(id2info);
			}

			pc_l_final.emplace_back();
			//std::swap(pc_l_final.back().usedReadPath, usedReadPath);
			std::swap(pc_l_final.back().p_l, pathList);
			pc_l_final.back().beginModeID = beginNodeID;
			//generate the MAX_out_keyNODE
			int MAX_REF_POS = -1000; int maxRPOS_UID = -1;
			for(Path &pc: pc_l_final.back().p_l){
				if(MAX_REF_POS < unitig_l[pc.endUID].refPos){
					MAX_REF_POS = unitig_l[pc.endUID].refPos;
					maxRPOS_UID = pc.endUID;
				}
				pc_l_final.back().nodeIDSetOut.emplace(pc.endUID);
			}
			if(pc_l_final.back().nodeIDSetOut.size() > 1){
				bool withRemove = false;
				for(int removeID: pc_l_final.back().nodeIDSetOut){
					if(removeID != maxRPOS_UID){
						keyPosUnitigIDset.erase(removeID);
						withRemove = true;
					}
				}
				if(withRemove){
					if(print_log) fprintf(stderr, "All Path Re-calculated begin at %5d\n", beginNodeID);
					pc_l_final.pop_back();
					i --;
				}
			}



		}
	}


};

#endif /* CPP_LIB_ASSEMBLER_DBGASSEMBLER_HPP_ */
