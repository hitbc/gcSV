/*
 * DBGAssembler.cpp
 *
 *  Created on: 2024-1-26
 *      Author: fenghe
 */

#include "DBGAssembler.hpp"

// int mode
// p: the position of KMER in reference or read
// i: the index of read in read list, ID of read

void KmerCounter::add_new_kmer(int mode, int p, int i, bool IsStoreReadInfo)
{
	if (mode == 0)
	{ // is REF
		refN++;
		if (refN > 1)
		{
			this->ref_pos = -2;
		}
		else
			this->ref_pos = p;
	}
	else if (mode == 1)
	{ // NGS read
		read_Srs_n++;
		if (IsStoreReadInfo)
		{
			NGS_rp.emplace_back(i, p);
		}
	}
	else if (mode == 2)
	{ // NGS read
		read_Lrs_n++;
		if (IsStoreReadInfo)
		{
			LRS_rp.emplace_back(i, p);
		}
	}
}

void DBG_INFO_item::showData()
{
	if (depth_read <= 1 && depth_ref == 0)
	{
		fprintf(stderr, "SKIP-D1");
		return;
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "KMER: [ %s\t%d\t%d\t%d\t] ", kmer.c_str(), kmerId, depth_read, depth_ref);
	fprintf(stderr, "[ %d\t:%d\t] ", IN_NUM, OUT_NUM);
	fprintf(stderr, "[ %d\t:%d\t:%d\t:%d\t]", unitigId, isRepeat, is_UNITIG_begin, is_UNITIG_end);

	fprintf(stderr, "IN: ");
	for (int i = 0; i < IN_NUM; i++)
		fprintf(stderr, "%d ", in[i]);
	fprintf(stderr, "OUT: ");
	for (int i = 0; i < OUT_NUM; i++)
		fprintf(stderr, "%d ", out[i]);
}

void DBG_INFO_item::showDataSimple()
{
	if (depth_read <= 1 && depth_ref == 0)
		return;
	fprintf(stderr, "[%d\t%d\t]\n ", depth_read, depth_ref);
}

void DBG_INFO_item::showDataDetail()
{
	if (depth_read <= 1 && depth_ref == 0)
		return;
	fprintf(stderr, "[%d\t%d\t] readSet is", depth_read, depth_ref);
	for (ReadPos &k : kcp->NGS_rp)
	{
		fprintf(stderr, "(%4d %3d) ", k.readID, k.readPos);
	}
	fprintf(stderr, "\n");
}

void UNITIG_INFO_item::setAverageKmerDepth(std::vector<DBG_INFO_item *> &id2info)
{
	if (totalDepthRead > 0)
	{
		averageKmerDepth_read = ((float)totalDepthRead) / kmer_list.size();
		with_read = (averageKmerDepth_read >= 1); // depth_read = 100*readNum/SIZE
	}
	uint MIN_SIZE = 20;
	if (kmer_list.size() > MIN_SIZE)
	{
		int ksize = kmer_list.size();
		float depth_head = 0;
		float depth_tail = 0;
		for (int i = 0; i < MIN_SIZE; i++)
		{
			depth_head += id2info[kmer_list[i]]->depth_read;
			depth_tail += id2info[kmer_list[ksize - i - 1]]->depth_read;
		}
		averageKmerDepth_read_head = depth_head / MIN_SIZE;
		averageKmerDepth_read_tail = depth_tail / MIN_SIZE;
		;
	}
	else
	{
		averageKmerDepth_read_head = averageKmerDepth_read;
		averageKmerDepth_read_tail = averageKmerDepth_read;
	}

	if (totalDepthRef > 0)
	{
		float ave_refD = ((float)totalDepthRef) / kmer_list.size();
		averageKmerDepth_ref = int(ave_refD + 0.9);
		with_ref = (averageKmerDepth_ref != 0);
	}
}

void UNITIG_INFO_item::setWithRead(int headNodeId, int endNodeId, std::vector<DBG_INFO_item *> &id2info, int MIN_DEPTH)
{
	bool withLowDepthKmerRead = false;
	int searchBg = 0;
	int searchEnd = kmer_list.size();
	if (this->id == headNodeId)
	{
		searchBg = MAX(searchEnd - 20, 0);
	}
	else if (this->id == endNodeId)
	{
		searchEnd = MIN(20, searchEnd);
	}
	for (int i = searchBg; i < searchEnd; i++)
		if (id2info[kmer_list[i]]->depth_read <= MIN_DEPTH)
			withLowDepthKmerRead = true;
	if (withLowDepthKmerRead)
	{
		if (false)
			fprintf(stderr, "unitig (%d) is set unused because with low depth kmer\n", id);
	}
	with_read = !withLowDepthKmerRead;
}

int UNITIG_INFO_item::get_ref_pos(std::vector<DBG_INFO_item *> &id2info)
{
	int pos_ref = id2info[kmer_list[0]]->kmerPOS_ref;
	if (pos_ref == -1)
	{
		int last_pos_ref = id2info[kmer_list.back()]->kmerPOS_ref;
		if (last_pos_ref != -1)
			pos_ref = last_pos_ref - kmer_list.size() + 1;
	}
	return pos_ref;
}

void UNITIG_INFO_item::setString(std::vector<DBG_INFO_item *> &id2info)
{
	unitig = id2info[kmer_list[0]]->kmer;
	unitig.pop_back();
	for (int i : kmer_list)
		unitig.push_back(id2info[i]->kmer.back());
}

void UNITIG_INFO_item::setInOutKmer(std::vector<DBG_INFO_item *> &id2info)
{
	DBG_INFO_item *beginKmer = id2info[kmer_list[0]];
	DBG_INFO_item *endKmer = id2info[kmer_list.back()];

	for (int i = 0; i < beginKmer->IN_NUM; i++)
		inBranchKmer.emplace_back(beginKmer->in[i]);

	for (int i = 0; i < endKmer->OUT_NUM; i++)
		outBranchKmer.emplace_back(endKmer->out[i]);
}

void UNITIG_INFO_item::setBranchUni(std::vector<DBG_INFO_item *> &id2info)
{
	for (int i : inBranchKmer)
		inBranchUni.emplace_back(id2info[i]->unitigId);
	if (inBranchUni.empty())
		inBranchUni.emplace_back(-1);
	for (int i : outBranchKmer)
		outBranchUni.emplace_back(id2info[i]->unitigId);
	if (outBranchUni.empty())
		outBranchUni.emplace_back(-1);
	inBranch_n = inBranchKmer.size();
	outBranch_n = outBranchKmer.size();

	unitig = id2info[kmer_list[0]]->kmer;
	unitig.pop_back();
	for (int i : kmer_list)
		unitig.push_back(id2info[i]->kmer.back());
}

void UNITIG_INFO_item::show(std::vector<DBG_INFO_item *> &id2info)
{
	// fprintf(stderr, "\n");
	if (totalDepthRead == (int)kmer_list.size())
	{
		fprintf(stderr, "SKIP1\t");
		return;
	}
	fprintf(stderr, "[BASE PART][ID %3d total %5ld, DP: %5d %5d (%5.2f ((%5d %5d)(%.2f %.2f)) %5.2f ((%5d %5d)))] [Repeat: %3d %3d]\t",
			id, kmer_list.size(), totalDepthRead, totalDepthRef,
			(float)totalDepthRead / kmer_list.size(),
			id2info[kmer_list[0]]->depth_read,
			id2info[kmer_list.back()]->depth_read,
			averageKmerDepth_read_head, averageKmerDepth_read_tail,

			(float)totalDepthRef / kmer_list.size(),
			id2info[kmer_list[0]]->depth_ref,
			id2info[kmer_list.back()]->depth_ref,

			isRepeat, finalRootId);

	fprintf(stderr, "[refPos:  %5d ]\n", get_ref_pos(id2info));

	fprintf(stderr, "[IN-OUT PART]IN:");
	for (int i : inBranchUni)
		fprintf(stderr, "[%3d->%3d]\t ", i, id);
	fprintf(stderr, "OUT:");
	for (int i : outBranchUni)
		fprintf(stderr, "[%3d->%3d]\t ", id, i);
	fprintf(stderr, "\n");
	fprintf(stderr, "[Depth part]");
	for (int i : kmer_list)
	{
		int dr = MIN(id2info[i]->depth_read, 87);
		fprintf(stderr, "%c", '#' + dr);
	}
	fprintf(stderr, "\n[String part]");
	fprintf(stderr, "\t%s\t:", unitig.c_str());
	fprintf(stderr, "\n");

	fprintf(stderr, "[TGS read position]");
	// the read position
	std::vector<ReadPos> &rpl_head = getAllReadTGSHead(id2info);
	int max_show = 100;
	for (ReadPos &p : rpl_head)
	{
		if (max_show-- < 0)
			break;
		p.show(stderr);
	}

	fprintf(stderr, "[TGS read position END]");
	// the read position
	std::vector<ReadPos> &rpl_tail = getAllReadTGSTail(id2info);
	max_show = 100;
	for (ReadPos &p : rpl_tail)
	{
		if (max_show-- < 0)
			break;
		p.show(stderr);
	}

	fprintf(stderr, "\n");
}

void MainDBGHandler1::addRef(bool print_log, uint8_t *tseq, int tlen, Contig_String_aligner *ca)
{
	std::string ref;
	for (int i = 0; i < tlen; i++)
	{
		ref += ("ACGTNNNNN"[tseq[i]]);
	}
	reference.emplace_back(ref);
	if (print_log)
		fprintf(stderr, "Ref is %s\n", ref.c_str());
	this->ca = ca;
}

bool Path::isSAMEwith(Path &c)
{
	return (bgnUID == c.bgnUID && endUID == c.endUID && SV_len == c.SV_len && unitigPL.size() == c.unitigPL.size() && type == c.type && ID_SUM == c.ID_SUM && pathKey == c.pathKey);
}

void Path::setNodeCount(std::map<UNITIG_INFO_item *, int> &pathNodeCountBuff)
{
	UNIQ_ID_product = 1;
	pathNodeCountBuff.clear();
	for (UNITIG_INFO_item *up : unitigPL)
	{
		std::map<UNITIG_INFO_item *, int>::iterator it = pathNodeCountBuff.find(up);
		if (it != pathNodeCountBuff.end())
		{
			it->second++;
		}
		else
			pathNodeCountBuff[up] = 1;
	}
	for (std::map<UNITIG_INFO_item *, int>::iterator it = pathNodeCountBuff.begin();
		 it != pathNodeCountBuff.end(); it++)
	{
		pnc.emplace_back();
		pnc.back().up = it->first;
		pnc.back().count = it->second;
		if (UNIQ_ID_product < 1e+150)
			UNIQ_ID_product *= (it->first->id + 2);
		UNIQ_ID_SUM += it->first->id;
	}
}

void Path::setSVLen(std::vector<DBG_INFO_item *> &id2info)
{
	int contigSpan = 0;
	ID_SUM = 0;
	for (UNITIG_INFO_item *up : unitigPL)
	{
		contigSpan += up->kmer_list.size();
		ID_SUM += up->id;
	}
	contigSpan -= unitigPL.back()->kmer_list.size();
	int pos_ref_begin = unitigPL[0]->get_ref_pos(id2info);
	int pos_ref_end = unitigPL.back()->get_ref_pos(id2info);
	int refSpan = pos_ref_end - pos_ref_begin;
	SV_len = contigSpan - refSpan;
}

bool Path::isSAMENodeSet(Path &c)
{
	return (bgnUID == c.bgnUID && endUID == c.endUID && UNIQ_ID_product == c.UNIQ_ID_product && UNIQ_ID_SUM == c.UNIQ_ID_SUM);
}

void Path::showSim(std::vector<DBG_INFO_item *> &id2info)
{
	fprintf(stderr, "\nSV_len %5d  ", SV_len);
	if (true)
	{
		for (UNITIG_INFO_item *up : unitigPL)
		{
			if (up == NULL)
				fprintf(stderr, "NA--");
			else
				fprintf(stderr, "%d--", up->id);
		}
	}
	if (true)
		fprintf(stderr, "\t");

	// show support read list:
	if (true)
		fprintf(stderr, "Show all support read for this path ");
	for (ReadPos &rp : support_read_l)
	{
		rp.show(stderr);
	}
	if (true)
		fprintf(stderr, "REMOVE reads ");
	for (ReadPos &rp : removed_read_l)
	{
		rp.show(stderr);
	}

	if (true)
		fprintf(stderr, "\n");
}

void Path::show(std::vector<DBG_INFO_item *> &id2info)
{
	fprintf(stderr, "%s SV signal in loop bgnUID %d, endUID %d node size %5ld \nSV_len %5d , %s, repeatRstNum %5d ID_SUM %5d UNIQSUM[ %f %5d ]\n",
			(type == 1) ? "Repeat" : "Linear", bgnUID, endUID, unitigPL.size(), SV_len, (SV_len > 30 || SV_len < -30) ? "SINFO_WITHSV" : "NOSV", repeatRstNum, ID_SUM, UNIQ_ID_product, UNIQ_ID_SUM);
	int global_pos = 0;
	int offset = 0;
	if (true)
	{
		for (UNITIG_INFO_item *up : unitigPL)
		{
			int pos_ref = up->get_ref_pos(id2info);
			if (pos_ref == -2)
				offset = -999; // repeat ref
			else if (pos_ref == -1)
				offset = -1999; // not include ref
			else
				offset = pos_ref - global_pos;
			fprintf(stderr, "%d--", up->id);
			if (false)
				fprintf(stderr, "[%2d %3ld (refPOS %5d : %5d %5d) (repeat %d %3d)]--->\n",
						up->id, up->kmer_list.size(), pos_ref,
						offset, global_pos, up->isRepeat, up->finalRootId);
			global_pos += up->kmer_list.size();
		}
		//			int mintotalRepeatNumRead = 0;
		//
		//			for(std::map<UNITIG_INFO_item *, int >::iterator it = pathNodeCount.begin();it != pathNodeCount.end(); it++){
		//
		//				int mintotalRepeatNumRead = 0;
		//			}
		for (PathNodeCount &c : pnc)
		{
			// it->first->show(id2info);
			fprintf(stderr, "\nNode %5d * %5d (%.1f %d %5ld)", c.up->id, c.count,
					c.up->averageKmerDepth_read,
					c.up->totalDepthRead,
					c.up->kmer_list.size());
		}
	}
	if (true)
		fprintf(stderr, "\n");
}

int getGreatestCommonDivisor(int m, int n)
{
	if (m == 0 || n == 0)
		return MAX(m, n);
	int gcd;
	while (n != 0)
	{
		gcd = m % n;
		m = n;
		n = gcd;
	}
	return m;
}

void PathCluster::updateStepCountByNewSignal()
{
	int greatestCommonDivisor = 0;
	for (NodeAnalysisINFO &nf_ : nfl)
	{
		if (nf_.currentCountDiff != 0)
		{
			greatestCommonDivisor =
				getGreatestCommonDivisor(nf_.currentCountDiff, greatestCommonDivisor);
		}
	}
	if (greatestCommonDivisor == 0)
		return;
	uint nodeSize = basePath.pnc.size();
	for (uint i = 0; i < nodeSize; i++)
	{
		NodeAnalysisINFO &nf_ = nfl[i];
		int new_step_Count = nf_.currentCountDiff / greatestCommonDivisor;
		if (nf_.step_Count == -1)
		{
			nf_.step_Count = new_step_Count;
		}
		else if (nf_.step_Count != new_step_Count)
		{
			isCPX = true;
		}
	}
	for (NodeAnalysisINFO &nf_ : nfl)
	{
		int new_step_Count = nf_.currentCountDiff / greatestCommonDivisor;
		if (nf_.step_Count == -1)
		{
			nf_.step_Count = new_step_Count;
		}
		else if (nf_.step_Count != new_step_Count)
		{
			isCPX = true;
		}
	}
}

void PathCluster::generateBasicAndStepSVLen()
{
	if (isCPX)
		return;
	int baseXInt = 0;
	if (signalNum == 1)
	{
		stepSVLength = 0;
	}
	else
	{
		uint nodeSize = basePath.pnc.size();
		float baseX = MAX_int32t;
		for (uint i = 0; i < nodeSize; i++)
		{
			NodeAnalysisINFO &nf = nfl[i];
			if (nf.step_Count != 0)
			{
				float curBaseX = basePath.pnc[i].count / nf.step_Count;
				baseX = MIN(baseX, curBaseX);
			}
		}
		baseXInt = baseX;
		stepSVLength = 0;
		BASE_UNIQ_ID_product = 1;
		for (uint i = 0; i < nodeSize; i++)
		{
			NodeAnalysisINFO &nf = nfl[i];
			nf.basic_Count = basePath.pnc[i].count - baseXInt * nf.step_Count;
			stepSVLength += basePath.pnc[i].up->kmer_list.size() * nf.step_Count;
			nf.depth = basePath.pnc[i].up->averageKmerDepth_read;
			if (nf.basic_Count != 0)
			{
				if (BASE_UNIQ_ID_product < 1e+150)
					BASE_UNIQ_ID_product *= (nf.up->id + 2);
			}
		}
	}
	basicSVLength = basePath.SV_len - stepSVLength * baseXInt;
}

void PathCluster::addSignal(Path &p)
{
	signalNum++;
	repeatResultNumber += p.repeatRstNum;
	path_list.emplace_back(p);
	if (signalNum == 1)
	{
		std::swap(basePath, p);
		uint nodeSize = basePath.pnc.size();
		// fprintf(stderr, "BASIC ");
		bgnUID = basePath.bgnUID;
		endUID = basePath.endUID;
		for (uint i = 0; i < nodeSize; i++)
		{
			nfl.emplace_back();
			nfl.back().up = basePath.pnc[i].up;
		}
	}
	else
	{
		// update currentCountDiff
		uint nodeSize = basePath.pnc.size();
		for (uint i = 0; i < nodeSize; i++)
		{
			if (basePath.pnc[i].up->id != p.pnc[i].up->id)
			{
				isCPX = true;
				break;
			}
			nfl[i].currentCountDiff = p.pnc[i].count - basePath.pnc[i].count;
			// fprintf(stderr, "%5d ", nfl[i].currentCountDiff);
		}
		if (!isCPX)
			updateStepCountByNewSignal();
	}
}

/// Returns a suffix (isTail is true) or prefix (otherwise) of the input contig with the specified length.
std::string getContigTail(const std::string &contig, const unsigned length,
						  const bool mode)
{
	xassert(length <= contig.size(), "");
	if (mode)
		return contig.substr(0, length);
	else
		return contig.substr((contig.size() - length), length);
}

std::string addContigBase(const std::string &contig, const char base,
						  const bool mode)
{
	if (mode)
		return base + contig;
	else
		return contig + base;
}

unsigned searchRepeatsUnitig(
	std::vector<UNITIG_INFO_item> &unitig_l,
	UNITIG_INFO_item &cu,
	const int index,
	std::vector<int> &kmerStack)
{

	cu.loopCurr = index;
	cu.loopRoot = index;
	// set the depth index for the current word to the smallest unused index
	unsigned nextIndex = index + 1;
	// push the current KMER into the stack
	kmerStack.push_back(cu.id);
	// deep-first search the DBG
	for (int i : cu.outBranchUni)
	{
		if (i == -1)
		{ /*Do nothing*/
		} // C1: end node
		else if (i == cu.id)
		{ // C2: point to it-self
			cu.isRepeat = true;
			cu.finalRootId = cu.id;
		}
		else if (unitig_l[i].loopCurr == 0)
		{ // the successor word has not been visited
			nextIndex = searchRepeatsUnitig(unitig_l, unitig_l[i], nextIndex, kmerStack);
			cu.loopRoot = MIN(cu.loopRoot, unitig_l[i].loopRoot);
		}
		else if (std::find(kmerStack.begin(), kmerStack.end(), i) != kmerStack.end())
			cu.loopRoot = MIN(cu.loopRoot, unitig_l[i].loopRoot);
		else
		{ /*Do nothing*/
		} // reach other node, but this node is not in the same loop with it
	}

	// if the current UNITIG is a root node:
	// ROOT NOTE: that`s meaning the node is the begin node to step in a cycle OR is a node not form a cycle
	// only the nodes in the middle of a cycle, the CONDITION NOT step in.
	if (cu.loopRoot == index)
	{
		int lastUniId = kmerStack.back();
		// exclude singletons, that is: linear KMERs not form a cycle
		if (lastUniId == cu.id)
			kmerStack.pop_back();
		else
		{ // KMER in a cycle
			int finalRootId = cu.id;
			while (true)
			{
				lastUniId = kmerStack.back();
				kmerStack.pop_back();
				unitig_l[lastUniId].isRepeat = true;
				unitig_l[lastUniId].finalRootId = finalRootId;
				if (lastUniId == cu.id)
					break;
			}
		}
	}
	return nextIndex;
}

bool MainDBGHandler1::buildDBG(int minDepthNGS, int minDepthTGS)
{
	// S1: init the hash table: word-count
	kmerIdx.clear();
	id2info.clear();
	unitig_l.clear();
	// build DBG:
	int kmer_id = 0;
	int kmer_LOW = 0; // used to show the number of low depth KMERs
	// select kmer with enough depth, and get the search directory
	for (KmerCounterTable::value_type &kmerNumItem : kmerCountingTable)
	{
		if (kmerNumItem.second.is_low_depth(minDepthNGS, minDepthTGS))
		{
			kmer_LOW += kmerNumItem.second.read_Srs_n + kmerNumItem.second.read_Lrs_n;
			continue;
		}
		int read_count = 0;
		if (kmerNumItem.second.read_Srs_n >= minDepthNGS)
			read_count = kmerNumItem.second.read_Srs_n;
		kmerIdx[kmerNumItem.first] = DBG_INFO_item();
		kmerIdx[kmerNumItem.first].setNew(kmer_id++, minDepthNGS, minDepthTGS, kmerNumItem.first, &(kmerNumItem.second));
		id2info.emplace_back(&(kmerIdx[kmerNumItem.first]));
	}
	fprintf(stderr, "kmer_LOW %d\n", kmer_LOW);
	const std::string &alphabet = "ACGT";
	// set in and out:
	for (std::unordered_map<std::string, DBG_INFO_item>::value_type &kmerIdx_i : kmerIdx)
	{
		const std::string &kmer_s = kmerIdx_i.first;
		DBG_INFO_item &currKmerInfo = kmerIdx_i.second;
		// deep-first search the DBG
		const std::string tmp(getContigTail(kmer_s, kmer_s.size() - 1, 0));
		for (const char symbol : alphabet)
		{
			// candidate successor of the current word
			const std::string nextKmer(addContigBase(tmp, symbol, 0));
			if (kmerIdx.count(nextKmer) == 0)
				continue; // the successor word does not exist in the reads

			DBG_INFO_item &nextKmerInfo = kmerIdx[nextKmer];
			// nextKmerInfo
			nextKmerInfo.setInOut(true, currKmerInfo.kmerId);
			currKmerInfo.setInOut(false, nextKmerInfo.kmerId);
			//
			if (nextKmerInfo.IN_NUM > 1)
			{
				for (int i = 0; i < nextKmerInfo.IN_NUM; i++)
					id2info[nextKmerInfo.in[i]]->is_UNITIG_end = true;
				nextKmerInfo.is_UNITIG_begin = true;
			}
			if (currKmerInfo.OUT_NUM > 1)
			{
				for (int i = 0; i < currKmerInfo.OUT_NUM; i++)
					id2info[currKmerInfo.out[i]]->is_UNITIG_begin = true;
				currKmerInfo.is_UNITIG_end = true;
			}
		}
	}

	// set the end from the end of read
	for (const std::string &kmer : read_end_kmer)
	{
		auto kmerIdx_i = kmerIdx.find(kmer);
		if (kmerIdx.find(kmer) != kmerIdx.end())
		{
			DBG_INFO_item &currKmerInfo = kmerIdx_i->second;
			for (int i = 0; i < currKmerInfo.OUT_NUM; i++)
				id2info[currKmerInfo.out[i]]->is_UNITIG_begin = true;
			currKmerInfo.is_UNITIG_end = true;
		}
	}

	// set the begin from the begin of read
	for (const std::string &kmer : read_begin_kmer)
	{
		auto kmerIdx_i = kmerIdx.find(kmer);
		if (kmerIdx.find(kmer) != kmerIdx.end())
		{
			DBG_INFO_item &currKmerInfo = kmerIdx_i->second;
			// this kmer is the begin and all kmer point to it is end:
			for (int i = 0; i < currKmerInfo.IN_NUM; i++)
				id2info[currKmerInfo.in[i]]->is_UNITIG_end = true;
			currKmerInfo.is_UNITIG_begin = true;
		}
	}

	// set begin and end
	std::set<int> headKmers;
	std::set<int> tailKmers;
	for (DBG_INFO_item *kmerIdx_i : id2info)
	{
		if (kmerIdx_i->IN_NUM == 0)
			kmerIdx_i->is_UNITIG_begin = true;
		if (kmerIdx_i->OUT_NUM == 0)
			kmerIdx_i->is_UNITIG_end = true;
		if (kmerIdx_i->is_UNITIG_begin)
			headKmers.emplace(kmerIdx_i->kmerId);
		if (kmerIdx_i->is_UNITIG_end)
			tailKmers.emplace(kmerIdx_i->kmerId);
	}

	// build UNITIG:
	int uid = 0;
	for (int kmerId : headKmers)
	{
		unitig_l.emplace_back(uid++);
		UNITIG_INFO_item &cu = unitig_l.back();
		DBG_INFO_item *curKmer = id2info[kmerId];
		for (;; curKmer = id2info[curKmer->out[0]])
		{
			cu.kmer_list.emplace_back(curKmer->kmerId);
			curKmer->unitigId = uid - 1;
			curKmer->unitigPos = cu.kmer_list.size() - 1;
			cu.totalDepthRead += curKmer->depth_read;
			cu.totalDepthRef += curKmer->depth_ref;
			if (curKmer->is_UNITIG_end)
				break;
		}
		cu.setString(id2info);
		cu.setInOutKmer(id2info);
	}
	// build UNITIG-Graph
	for (UNITIG_INFO_item &cu : unitig_l)
	{
		cu.setBranchUni(id2info);
	}
	// search loop in the graph:
	{
		int index = 1;
		std::vector<int> loopStack;
		// for each word, search the index in a single path
		for (UNITIG_INFO_item &cu : unitig_l)
		{
			if (cu.loopCurr == 0)
				index = searchRepeatsUnitig(unitig_l, cu, index, loopStack);
		}
	}
	bool withLoop = false;
	// set repeat info to KMERs
	for (UNITIG_INFO_item &cu : unitig_l)
	{
		bool isRepeat = cu.isRepeat;
		if (isRepeat)
			withLoop = true;
		bool finalRootId = cu.finalRootId;
		for (int i : cu.kmer_list)
		{
			id2info[i]->isRepeat = isRepeat;
			id2info[i]->unitigRootId = finalRootId;
		}
		cu.refPos = cu.get_ref_pos(id2info);
	}

	// set ref pos for all unitig
	// std::sort(unitig_l.begin(), unitig_l.end(), UNITIG_INFO_item::cmp_by_ref_pos);
	// sort UNITIGs by ref-position and repeat STAT
	for (UNITIG_INFO_item &cu : unitig_l)
	{
		bool isRepeat = cu.isRepeat;
		if (isRepeat)
			withLoop = true;
		bool finalRootId = cu.finalRootId;
		for (int i : cu.kmer_list)
		{
			id2info[i]->isRepeat = isRepeat;
			id2info[i]->unitigRootId = finalRootId;
		}
		cu.refPos = cu.get_ref_pos(id2info);
		cu.setAverageKmerDepth(id2info);
	}

	// set key UNITIG vector and set
	{
		keyPosUnitigL.clear();
		for (UNITIG_INFO_item &cu : unitig_l)
		{
			if (cu.isRepeat == false && cu.refPos != -1 && cu.with_read && cu.kmer_list.size() >= 30)
			{
				keyPosUnitigL.emplace_back();
				keyPosUnitigL.back().p = &cu;
				keyPosUnitigL.back().ref_pos = cu.refPos;
			}
		}
		std::sort(keyPosUnitigL.begin(), keyPosUnitigL.end(), KEY_UNITIG::cmp_by_ref_pos);
		keyPosUnitigIDset.clear();
		for (KEY_UNITIG &ku : keyPosUnitigL)
		{
			keyPosUnitigIDset.emplace(ku.p->id);
		}
		if (keyPosUnitigL.empty() == false)
		{
			int headNodeId = keyPosUnitigL[0].p->id;
			int endNodeId = keyPosUnitigL.back().p->id;
			for (UNITIG_INFO_item &cu : unitig_l)
			{
				cu.setWithRead(headNodeId, endNodeId, id2info, 0);
			}
		}
	}

	// debug show all kmers
	if (REPEAT_DEBUG)
	{
		for (DBG_INFO_item *kmerIdx_i : id2info)
			kmerIdx_i->showData();
	}

	// debug show all untigs
	for (UNITIG_INFO_item &cu : unitig_l)
	{
		if (REPEAT_DEBUG)
			cu.show(id2info);
	}

	if (REPEAT_DEBUG)
		printf("END");
	return withLoop;
}

void MainDBGHandler1::showDBG_graph()
{
	if (false)
	{
		fprintf(stderr, "\n\nShow All Kmers\n\n");
		for (KmerCounterTable::iterator it = kmerCountingTable.begin(); it != kmerCountingTable.end(); it++)
		{
			if (it->second.read_Srs_n > 1 || it->second.refN > 0)
				fprintf(stderr, "%s: [ %d | %d ]\n", it->first.c_str(), it->second.read_Srs_n, it->second.refN);
		}
	}
	if (false)
	{
		fprintf(stderr, "\n\nShow All Unitigs\n\n");
		if (true)
			for (UNITIG_INFO_item &cu : unitig_l)
			{
				cu.show(id2info);
			}
		fprintf(stderr, "\n\nShow All Unitigs END \n");
	}
	if (true)
	{
		fprintf(stderr, "\n\nShow All Unitigs LINK\n\n");
		if (true)
			for (UNITIG_INFO_item &cu : unitig_l)
			{
				cu.show_link(id2info);
			}
		fprintf(stderr, "\n\nShow All Unitigs LINK END \n");
	}

	fprintf(stderr, "\n\nShow KEY Unitigs\n\n");
	for (KEY_UNITIG &k : keyPosUnitigL)
	{
		k.p->show(id2info);
	}
	fprintf(stderr, "\n\nEND Show\n\n");
}

bool MainDBGHandler1::randomGetPath(int beginNodeID, std::vector<Path> &pathList, int searchPathNum, std::set<int> &stopNodeSet,
									std::vector<UNITIG_INFO_item *> &randomSearchBuff)
{
	pathList.clear();
	bool withRandomBranch = false;
	for (int i = 0; i < searchPathNum; i++)
	{
		randomSearchBuff.clear();
		int curID = beginNodeID;
		uint64_t pathKey = 0;
		// if(this->curKmerLen == 90 && beginNodeID == 6)
		//	fprintf(stderr, " ");
		randomSearchBuff.emplace_back(&(unitig_l[curID])); // add the begin node
		int totalLength = 0;
		while (totalLength++ < 100)
		{
			pathKey ^= curID;
			pathKey <<= 5;
			pathKey ^= curID;
			pathKey <<= 12;
			pathKey ^= curID;

			curID = randomGetNextNodeID(curID, withRandomBranch);
			if (curID == -1) // reach the end node
			{
				//				pathList.emplace_back();
				//				pathList.back().bgnUID = beginNodeID;
				//				pathList.back().endUID = curID;
				//				pathList.back().type = 1;
				//				std::swap(pathList.back().unitigPL, randomSearchBuff);
				//				pathList.back().setSVLen(id2info);
				break;
			}
			randomSearchBuff.emplace_back(&(unitig_l[curID]));

			// reach one of the key node
			if (stopNodeSet.find(curID) != stopNodeSet.end())
			{
				pathList.emplace_back();
				pathList.back().bgnUID = beginNodeID;
				pathList.back().endUID = curID;
				pathList.back().type = 1;
				pathList.back().pathKey = pathKey;
				std::swap(pathList.back().unitigPL, randomSearchBuff);
				pathList.back().setSVLen(id2info);
				break;
			}
		}
		if (withRandomBranch == false)
		{
			break;
		}
	}
	return withRandomBranch;
}

void MainDBGHandler1::sortAndRemoveDupPath(std::vector<Path> &pathList, bool withRandomPath, int pathNumberPerTime)
{
	std::sort(pathList.begin(), pathList.end(), Path::cmp_by_uid);
	// remove same path
	removeRepeatPathBuff.clear();
	for (uint i = 0; i < pathList.size(); i++)
	{
		if (removeRepeatPathBuff.empty() || !pathList[i].isSAMEwith(removeRepeatPathBuff.back()))
		{
			removeRepeatPathBuff.emplace_back();
			std::swap(removeRepeatPathBuff.back(), pathList[i]);
			removeRepeatPathBuff.back().repeatRstNum = 1;
		}
		else
			removeRepeatPathBuff.back().repeatRstNum++;
	}
	if (!removeRepeatPathBuff.empty() && !withRandomPath)
	{
		removeRepeatPathBuff.back().repeatRstNum = pathNumberPerTime;
	}
	std::swap(removeRepeatPathBuff, pathList);
}

void MainDBGHandler1::generatePathCluster(std::vector<Path> &pathList,
										  std::vector<PathCluster> &pathClusterList)
{
	for (Path &p : pathList)
		p.setNodeCount(pathNodeCountBuff);
	std::sort(pathList.begin(), pathList.end(), Path::cmp_by_idSUM);
	// analysis:
	if (false)
	{
		for (Path &p : pathList)
			p.show(id2info);
	}
	// generate pathCluster
	pathClusterList.clear();
	for (Path &p : pathList)
	{
		if (pathClusterList.empty() || !pathClusterList.back().isSAME_ID_SET(p))
		{
			pathClusterList.emplace_back();
		}
		// todo::
		pathClusterList.back().addSignal(p);
	}
	// update pathCluster
	for (PathCluster &pc : pathClusterList)
		pc.generateBasicAndStepSVLen();
	// absorb linear path using repeat path cluster
	if (pathClusterList.size() <= 30)
	{
		std::vector<PathCluster> pathClusterListTmp;
		// collect all repeat signals
		for (PathCluster &pc : pathClusterList)
		{
			if (pc.signalNum > 1)
			{
				pathClusterListTmp.emplace_back();
				pathClusterListTmp.back().signalNum = 0;
				std::swap(pathClusterListTmp.back(), pc);
			}
		}
		// try to absorb
		for (PathCluster &pc : pathClusterList)
		{
			if (pc.signalNum == 1)
			{
				bool isAbsorbed = false;
				for (PathCluster &pc_r : pathClusterListTmp)
				{
					if (pc.basicSVLength == pc_r.basicSVLength && pc.nfl.size() < pc_r.nfl.size() && pc_r.signalNum > 1)
					{
						// float proDiff = pc_r.basePath.UNIQ_ID_product/pc.basePath.UNIQ_ID_product;
						// test if A is integer
						if (pc_r.BASE_UNIQ_ID_product == pc.basePath.UNIQ_ID_product)
						{
							pc_r.repeatResultNumber += pc.repeatResultNumber;
							isAbsorbed = true;
							pc_r.path_list.emplace_back(pc.basePath);
							// fprintf(stderr, "Node isAbsorbed\n");
							break;
						}
					}
				}
				if (isAbsorbed == false)
				{
					pathClusterListTmp.emplace_back();
					pathClusterListTmp.back().signalNum = 2;
					std::swap(pathClusterListTmp.back(), pc);
				}
			}
		}
		std::swap(pathClusterListTmp, pathClusterList);
	}
	if (false)
		for (PathCluster &pc : pathClusterList)
		{
			pc.show(id2info);
		}
}

int N_count(const std::string &word)
{
	int N = 0;
	const char *s = word.c_str();
	int n = word.size();
	for (int i = 0; i < n; i++)
	{
		if (s[i] == 'N')
			N++;
	}
	return N;
}

int getBinBase(char c)
{
	switch (c)
	{
	case 'A':
		return 0;
		break;
	case 'C':
		return 1;
		break;
	case 'G':
		return 2;
		break;
	case 'T':
		return 3;
		break;
	case 'N':
		return 4;
		break;
	default:
		fprintf(stderr, "c: %c", c);
		xassert(0, "");
		break;
	}
	return 0;
}

/// Construct k-mer maps
/// wordCount: k-mer ==> number of reads containing the k-mer
/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
void MainDBGHandler1::kmerCountingTest(const unsigned kmerSize)
{
	curKmerLen = kmerSize;
	const unsigned strCount = reads_str.size();
	for (unsigned strIndex = 0; strIndex < strCount; ++strIndex)
	{
		// stores the index of a kmer in a read sequence
		const std::string &seq = reads_str[strIndex];
		const unsigned seqLen = seq.size();
		// this read is unusable for assembly:
		if (seqLen < kmerSize)
			continue;
		// track all words from the read, including repetitive k
		unsigned kmer_number = seqLen - kmerSize + 1;
		int usedKmer = 0;
		for (unsigned j = 0; j < kmer_number; ++j)
		{
			const std::string word(seq.substr(j, kmerSize));
			// fprintf(stderr, "%s \n",word.c_str());
			//  filter words with "N" (either directly from input alignment or marked due to low basecall quality:
			KmerCounterTable::iterator it = kmerCountingTable.find(word);
			if (it == kmerCountingTable.end() || it->second.read_Srs_n < 3)
			{
			}
			else
				usedKmer++;
		}
		fprintf(stderr, "strIndex: %d (usedKmer %d)\n", strIndex, usedKmer);
	}
}

/// Construct k-mer maps
/// wordCount: k-mer ==> number of reads containing the k-mer
/// wordSupportReads: k-mer ==> a list of read IDs containg the k-mer
void MainDBGHandler1::errorCorrection(const unsigned kmerSize, int minRight, int maxError)
{
	const unsigned strCount = reads_str.size();
	for (unsigned strIndex = 0; strIndex < strCount; ++strIndex)
	{
		int correctBaseACGT = 0;
		int correctBaseN = 0;
		// stores the index of a kmer in a read sequence
		const std::string seqOri = reads_str[strIndex];
		const std::string &seq = reads_str[strIndex];
		const unsigned seqLen = seq.size();
		// this read is unusable for assembly:
		if (seqLen < kmerSize)
			continue;
		// track all words from the read, including repetitive k
		unsigned kmer_number = seqLen - kmerSize + 1;
		bool withFail = false;
		if (false)
			fprintf(stderr, "Before correction %s \n", seq.c_str());
		for (unsigned j = 0; j < kmer_number; ++j)
		{
			std::string word(seq.substr(j, kmerSize));
			KmerCounterTable::iterator it = kmerCountingTable.find(word);
			if (it == kmerCountingTable.end() || it->second.read_Srs_n <= maxError)
			{ // need correction
				char lastBase = word.back();
				int bestBase = -1;
				int bastCount = 0;
				for (int i = 0; i < 4; i++)
				{
					if (lastBase == "ACCT"[i])
						continue;
					word[word.size() - 1] = "ACGT"[i];
					it = kmerCountingTable.find(word);
					if (it != kmerCountingTable.end() && ((it->second.read_Srs_n >= minRight && it->second.read_Srs_n > bastCount)))
					{ //|| (bestBase != -1 && it->second.refN != 0)
						bestBase = i;
						bastCount = it->second.read_Srs_n;
					}
				}
				if (bestBase != -1)
				{
					if (lastBase != 'N')
						correctBaseACGT++;
					else
						correctBaseN++;

					reads_str[strIndex][j + kmerSize - 1] = "ACGT"[bestBase];
					if (false)
						fprintf(stderr, "Base correct readID %5d pos %3d\n", strIndex, j + kmerSize - 1);
				}
				else
				{
					withFail = true;
					if (false)
						fprintf(stderr, "Base correct FAILED: pos %3d\n", j + kmerSize - 1);
				}
			}
		}
		if (withFail)
		{
			withFail = false;
			for (unsigned j = 0; j < kmer_number; ++j)
			{
				std::string word(seq.substr(seq.size() - j - kmerSize, kmerSize));
				KmerCounterTable::iterator it = kmerCountingTable.find(word);
				if (it == kmerCountingTable.end() || it->second.read_Srs_n <= maxError)
				{ // need correction
					char firstBase = word[0];
					int bestBase = -1;
					int bastCount = 0;
					for (int i = 0; i < 4; i++)
					{
						if (firstBase == "ACCT"[i])
							continue;
						word[0] = "ACGT"[i];
						it = kmerCountingTable.find(word);
						if (it != kmerCountingTable.end() && ((it->second.read_Srs_n >= minRight && it->second.read_Srs_n > bastCount)))
						{ //
							bestBase = i;
							bastCount = it->second.read_Srs_n;
						}
					}
					if (bestBase != -1)
					{
						if (firstBase != 'N')
							correctBaseACGT++;
						else
							correctBaseN++;
						reads_str[strIndex][seq.size() - j - kmerSize] = "ACGT"[bestBase];
						if (false)
							fprintf(stderr, "Base correct (Reverse) readID %5d pos %3ld\n", strIndex, seq.size() - j - kmerSize);
					}
					else
					{
						withFail = true;
						if (false)
							fprintf(stderr, "Base correct (Reverse) FAILED: pos %3ld\n", seq.size() - j - kmerSize);
					}
				}
			}
		}
		// restore original when correct too many bases
		if (correctBaseACGT > 1 || correctBaseN > 1)
		{
			reads_str[strIndex] = seqOri;
			if (false)
				fprintf(stderr, "Base correct readID %5d \t", strIndex);
			if (false)
				fprintf(stderr, "correctBaseN %d correctBaseN %d\n", correctBaseACGT, correctBaseN);
		}
		if (false)
			fprintf(stderr, "correctBaseN %d \n", correctBaseACGT);
		if (false)
			fprintf(stderr, "strIndex %d: After  correction %s \n\n", strIndex, seq.c_str());
	}
}

int MainDBGHandler1::randomGetNextNodeID(int curID, bool &withRandomBranch)
{
	nextNodeBuff.clear();
	std::vector<int> &branchUni = unitig_l[curID].outBranchUni;
	for (int outID : branchUni)
	{
		//		if(read_depth[outID].kmer_list_size == 42 && outID == 11){
		//			fprintf(stderr, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
		//		}
		if (outID != -1 && unitig_l[outID].with_read)
		{
			nextNodeBuff.emplace_back(outID);
		}
	}
	if (nextNodeBuff.empty())
		return -1;
	else if (nextNodeBuff.size() == 1)
	{
		curID = nextNodeBuff[0];
	}
	else
	{
		simple_rand(&random_seed);
		int select = random_seed % nextNodeBuff.size();
		// if(curID == 5 && nextNode.size() == 2 && unitig_l.size() == 9)
		//	fprintf(stderr, "select %d %d mode %d\n", select, nextNode[select],mode);
		curID = nextNodeBuff[select];
		withRandomBranch = true;
	}
	return curID;
}
