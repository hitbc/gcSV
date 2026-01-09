/*
 * assembler.cpp
 *
 *  Created on: 2021-8-9
 *      Author: fenghe
 */

#include "GreedyAssembler.hpp"
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <vector>

// compile with this macro to get verbose output:
//#define DEBUG_ASBL
//#define DEBUG_WALK

/// Adds base @p base to the end (isEnd is true) or start (otherwise) of the contig.
///
/// \return the extended contig.
static std::string addContigBase(const std::string &contig, const char base,
		const bool isEnd) {
	if (isEnd)
		return contig + base;
	else
		return base + contig;
}

/// Returns a suffix (isEnd is true) or prefix (otherwise) of the input contig with the specified length.
static std::string getContigEnd(const std::string &contig, const unsigned length,
		const bool isEnd) {
	const unsigned csize(contig.size());
	assert(length <= csize);

	if (isEnd)
		return contig.substr((csize - length), length);
	else
		return contig.substr(0, length);
}


/// Identify repetitive k-mers
/// i.e. k-mers that form a circular subgraph
///
static unsigned searchRepeatsKmer(const std::string &alphabet, const unsigned index,
		const std::string &word, str_pair_uint_map_t &wordIndices,
		std::vector<std::string> &wordStack,
		std::set<std::string> &repeatWords) {
	// set the depth index for the current word to the smallest unused index
	wordIndices[word] = std::pair<unsigned, unsigned>(index, index);
	unsigned nextIndex = index + 1;
	wordStack.push_back(word);

	const std::string tmp(getContigEnd(word, word.size() - 1, true));
	for (const char symbol : alphabet) {
		// candidate successor of the current word
		const std::string nextWord(addContigBase(tmp, symbol, true));

		// homopolymer
		if (word == nextWord) {
			repeatWords.insert(word);
			continue;
		}

		// the successor word does not exist in the reads
		if (wordIndices.count(nextWord) == 0)
			continue;
		const unsigned nextWordIdx = wordIndices[nextWord].first;
		if (nextWordIdx == 0) {
			// the successor word has not been visited
			// recurse on it
			nextIndex = searchRepeatsKmer(alphabet, nextIndex, nextWord,//
					wordIndices, wordStack, repeatWords);
			// update the current word's lowlink
			const unsigned wordLowLink = wordIndices[word].second;
			const unsigned nextWordLowLink = wordIndices[nextWord].second;
			wordIndices[word].second = std::min(wordLowLink, nextWordLowLink);
		} else {
			const bool isContained(
					std::find(wordStack.begin(), wordStack.end(), nextWord)
							!= wordStack.end());
			if (isContained) {
				// the successor word is in stack and therefore in the current circle of words
				// only update the current word's lowlink
				const unsigned wordLowLink = wordIndices[word].second;
				wordIndices[word].second = std::min(wordLowLink, nextWordIdx);
			}
		}
	}

	// if the current word is a root node,
	if (wordIndices[word].second == index) {
		const std::string &lastWord(wordStack.back());
		// exclude singletons
		bool isSingleton(lastWord == word);
		if (isSingleton) {
			wordStack.pop_back();
		} else {
			// record identified repeat words (i.e. words in the current circle) if the circle is small
			const unsigned lastWordIndex(wordIndices[lastWord].first);
			const bool isSmallCircle((lastWordIndex - index) <= 50);
			while (true) {
				const std::string repeatWd = wordStack.back();
				if (isSmallCircle)
					repeatWords.insert(repeatWd);
				wordStack.pop_back();

				if (repeatWd == word)
					break;
			}
		}
	}

	return nextIndex;
}

#define READ_POS_MAX_SEARCH 5
int search_read_position(std::string &read_seq, std::string &maxWord, int wordLength, int mode, int max_search){
	int read_kmer_num = read_seq.size() - wordLength + 1;
	int max_search_len = MIN(max_search, read_kmer_num);
	int search_len = 0;

	if(mode == 0){//walk to the right
		for(int i = 0; i < max_search_len; i++){
			if(0 == read_seq.compare(i, wordLength, maxWord, 0, wordLength)) break;
			search_len++;
		}
	}else{//walk to the left
		for(int i = max_search_len - 1; i >= 0; i--){
			if(0 == read_seq.compare(i, wordLength, maxWord, 0, wordLength)) break;
			search_len++;
		}
	}
	return (search_len);
}

/// Construct a contig from the seed provided
///
/// The contig is extened in both directions until
/// either the there is no sufficient support evidence
/// or the last k-mer is repeatitive (i.e. part of a bubble in the graph)
///
/// \return True if the contig runs into a repeatitive k-mer when extending in either mode

#define MAX_ALLEL_AS_SNP 2
bool MainAssemblyHandler::walk(const std::string &seed, const unsigned wordLength,
		const str_uint_map_t &wordCount, const str_set_uint_map_t &wordReads,
		const std::set<std::string> &repeatWords,
		std::set<std::string> &unusedWords, AssemblyContig &contig) {

	//load buffs
	std::set<std::string> &wordsInContig = walk_wordsInContig;
	wordsInContig.clear();
	const str_uint_map_t::const_iterator wordCountEnd(wordCount.cend());
	const str_set_uint_map_t::const_iterator wordReadsEnd(wordReads.cend());

	// we start with the seed
	str_set_uint_map_t::const_iterator wordReadsIter = wordReads.find(seed);
	assert(wordReadsIter != wordReadsEnd);
	contig.ass_begin_offset_in_contig = 0;
	contig.supportReads = wordReadsIter->second;
	contig.seedReadCount = contig.supportReads.size();
	contig.seq = seed;
	contig.kmerLength = wordLength;
	contig.actions.clear();
	for (const unsigned read_ID : contig.supportReads) {
		contig.actions.emplace_back(0, read_ID, true);
	}
	contig.ending_reason[0] = 0;
	contig.ending_reason[1] = 0;

	unusedWords.erase(seed);

	if (repeatWords.find(seed) != repeatWords.end()) {
		contig.conservativeRange_bgn = (0);
		contig.conservativeRange_end = (wordLength);
		contig.ending_reason[0] = 1;
		contig.ending_reason[1] = 1;
		return true;
	}

	std::set<unsigned> &maxWordReads = walk_maxWordReads;
	std::set<unsigned> &maxContigWordReads = walk_maxContigWordReads;
	std::set<unsigned> &previousWordReads = walk_previousWordReads;
	std::set<unsigned> &supportReads2Remove = walk_supportReads2Remove;
	std::set<unsigned> &rejectReads2Add = walk_rejectReads2Add;
	std::set<unsigned> &contigWordReads = walk_contigWordReads;
	std::set<unsigned> &sharedReads = walk_sharedReads;
	std::set<unsigned> &toRemove = walk_toRemove;
	std::set<unsigned> &toAdd = walk_toAdd;
	std::set<unsigned> &sharedReads_alleles = walk_sharedReads_alleles;
	std::set<unsigned> &toUpdate = walk_toUpdate;

	// collecting words used to build the contig
	wordsInContig.insert(seed);

	const std::string tmpTrunk = getContigEnd(seed, wordLength - 1, false);
	// collecting rejecting reads for the seed from the unselected branches
	for (const char symbol : o.alphabet) {
		// the seed itself
		if (symbol == seed[wordLength - 1])
			continue;

		// add rejecting reads from an unselected word/branch
		const std::string newKey(addContigBase(tmpTrunk, symbol, true));

		wordReadsIter = wordReads.find(newKey);
		if (wordReadsIter == wordReadsEnd)
			continue;
		const std::set<unsigned> &unselectedReads(wordReadsIter->second);

		contig.rejectReads.insert(unselectedReads.begin(),
				unselectedReads.end());
	}

	bool isRepeatFound(false);
	int &kmer_index = (contig.ass_begin_offset_in_contig);
	// 0 => walk to the right, 1 => walk to the left
	for (unsigned mode(0); mode < 2; ++mode) {
		const bool isEnd(mode == 0);
		unsigned conservativeEndOffset(0);
		kmer_index = 0;
		int kmer_index_step = (mode == 0)?1:-1;

		while (true) {
			const std::string previousWord = getContigEnd(contig.seq, wordLength, isEnd);

			const std::string trunk(getContigEnd(contig.seq, wordLength - 1, isEnd));
			unsigned maxBaseCount(0);
			unsigned maxContigWordReadCount(0);
			char maxBase(o.alphabet[0]);
			std::string maxWord;
			maxWordReads.clear();
			maxContigWordReads.clear();
			previousWordReads.clear();
			supportReads2Remove.clear();
			rejectReads2Add.clear();

			for (const char symbol : o.alphabet) {
				const std::string newKey(addContigBase(trunk, symbol, isEnd));
				const str_uint_map_t::const_iterator wordCountIter(
						wordCount.find(newKey));//
				if (wordCountIter == wordCountEnd)
					continue;
				const unsigned currWordCount(wordCountIter->second);

				wordReadsIter = wordReads.find(newKey);
				if (wordReadsIter == wordReadsEnd)
					continue;
				const std::set<unsigned> &currWordReads(wordReadsIter->second);

				// get the shared supporting reads between the contig and the current word
				contigWordReads.clear();
				std::set_intersection(contig.supportReads.begin(),
						contig.supportReads.end(), currWordReads.begin(),
						currWordReads.end(),
						std::inserter(contigWordReads,
								contigWordReads.begin()));

				// get the shared supporting reads across two alleles
				sharedReads.clear();
				std::set_intersection(maxContigWordReads.begin(),
						maxContigWordReads.end(), currWordReads.begin(),
						currWordReads.end(),
						std::inserter(sharedReads, sharedReads.begin()));

				if (contigWordReads.empty())
					continue;

				const unsigned contigWordReadCount(contigWordReads.size());
				if (contigWordReadCount > maxContigWordReadCount) {
					// the old shared reads support an unselected allele if they don't support the new word
					// remove them from the contig's supporting reads
					if (!maxContigWordReads.empty()) {
						toRemove.clear();
						std::set_difference(maxContigWordReads.begin(),
								maxContigWordReads.end(), sharedReads.begin(),
								sharedReads.end(),
								std::inserter(toRemove, toRemove.begin()));
						if (toRemove.size() >  MAX_ALLEL_AS_SNP)
							supportReads2Remove.insert(toRemove.begin(),
									toRemove.end());
					}

					// the old supporting reads is for an unselected allele if they don't support the new word
					// they become rejecting reads for the currently selected allele
					if (!maxWordReads.empty()) {
						toAdd.clear();
						std::set_difference(maxWordReads.begin(),
								maxWordReads.end(), sharedReads.begin(),
								sharedReads.end(),
								std::inserter(toAdd, toAdd.begin()));
						if (toAdd.size() >  MAX_ALLEL_AS_SNP)
							rejectReads2Add.insert(toAdd.begin(), toAdd.end());
					}
					// new supporting reads for the currently selected allele
					maxWordReads = currWordReads;

					maxContigWordReadCount = contigWordReadCount;
					maxContigWordReads = contigWordReads;
					maxBaseCount = currWordCount;
					maxBase = symbol;
					maxWord = newKey;
				} else {
					toRemove.clear();
					std::set_difference(contigWordReads.begin(),
							contigWordReads.end(), sharedReads.begin(),
							sharedReads.end(),
							std::inserter(toRemove, toRemove.begin()));
					if (toRemove.size() >  MAX_ALLEL_AS_SNP)
						supportReads2Remove.insert(toRemove.begin(),
								toRemove.end());

					toAdd.clear();
					std::set_difference(currWordReads.begin(),
							currWordReads.end(), sharedReads.begin(),
							sharedReads.end(),
							std::inserter(toAdd, toAdd.begin()));
					if (toAdd.size() >  MAX_ALLEL_AS_SNP)
						rejectReads2Add.insert(toAdd.begin(), toAdd.end());
				}
			}

			if (maxBaseCount < o.minCoverage) {
				contig.ending_reason[1 - mode] = 0;
				break;
			}

			// stop walk in the current mode after seeing one repeat word
			if (wordsInContig.find(maxWord) != wordsInContig.end()) {
				isRepeatFound = true;
				contig.ending_reason[1 - mode] = 1;
				break;
			}

			contig.seq = addContigBase(contig.seq, maxBase, isEnd);
			//add kmer count
			kmer_index += kmer_index_step;

			if ((conservativeEndOffset != 0)
					|| (maxBaseCount < o.minConservativeCoverage))
				conservativeEndOffset += 1;

			{
				// walk backwards for one step at a branching point
				if (maxWordReads != previousWordReads) {
					const char tmpSymbol = (isEnd ? previousWord[0] : previousWord[wordLength - 1]);
					for (const char symbol : o.alphabet) {
						// the selected branch: skip the backward word itself
						if (symbol == tmpSymbol)
							continue;

						// add rejecting reads from an unselected branch
						const std::string newKey( addContigBase(trunk, symbol, !isEnd));
						// the selected branch: skip the word just extended
						if (newKey == maxWord)
							continue;

						wordReadsIter = wordReads.find(newKey);
						if (wordReadsIter == wordReadsEnd)
							continue;

						const std::set<unsigned> &backWordReads( wordReadsIter->second);
						// get the shared supporting reads across two alleles
						sharedReads_alleles.clear();
						std::set_intersection(maxContigWordReads.begin(),
								maxContigWordReads.end(), backWordReads.begin(),
								backWordReads.end(),
								std::inserter(sharedReads_alleles,
										sharedReads_alleles.begin()));

						toUpdate.clear();
						std::set_difference(backWordReads.begin(),
								backWordReads.end(),
								sharedReads_alleles.begin(),
								sharedReads_alleles.end(),
								std::inserter(toUpdate, toUpdate.begin()));
						if (toUpdate.size() >  MAX_ALLEL_AS_SNP) {
							rejectReads2Add.insert(toUpdate.begin(),
									toUpdate.end());
							supportReads2Remove.insert(toUpdate.begin(),
									toUpdate.end());
						}
					}
				}
				previousWordReads = maxWordReads;
				// update rejecting reads
				// add reads that support the unselected allele
				for (const unsigned rd : rejectReads2Add) {
					contig.rejectReads.insert(rd);
				}
				// update supporting reads, add reads that support the selected allele
				// normally, the [maxWordReads] is all reads that support a KMER
				for (const unsigned read_ID : maxWordReads) {
					//skip the adding process when the read is already in the support read list
					if(contig.supportReads.find(read_ID) != contig.supportReads.end()) continue;
					bool read_insert_2_support_at_read_begin = (search_read_position(reads[read_ID], maxWord, wordLength, mode, READ_POS_MAX_SEARCH) < READ_POS_MAX_SEARCH);
					//try to add the read into support read set:
					bool read_in_reject_set = (contig.rejectReads.find(read_ID) != contig.rejectReads.end());
					//if read is in reject set, but the kmer is in the begin of the read (or end of the read when reverse extension),
					//it will be removed from the reject set and re-inserted into the support reads list (Unless the read is reject just in this step)
					if(read_in_reject_set){
						bool read_is_not_newly_rejected = (rejectReads2Add.find(read_ID) == rejectReads2Add.end());
						if(read_is_not_newly_rejected){
							//when the new kmer is in the begin of the read
							if(read_insert_2_support_at_read_begin){
								contig.rejectReads.erase(read_ID);
								contig.supportReads.insert(read_ID);
								contig.actions.emplace_back(kmer_index, read_ID, true);
							}
						}
					}
					//if reads is not in rejected set, we check whether the reads is used from the begin (of read), when it is from begin, the reads will be add /
					//into support reads set, otherwise the reads must be undergo read-contig alignment test, when with INDELS or too many wrong bases, it will be reject /
					//when pass the filter, it will be added to support list.
					else{
						bool read_will_be_used = false;
						if(read_insert_2_support_at_read_begin){
							read_will_be_used = true;
						}else{
							int true_read_position = (search_read_position(reads[read_ID], maxWord, wordLength, mode, 1000));
							//match check:
							const char * read_c = reads[read_ID].c_str();
							const char * ref_c  = contig.seq.c_str();
							if(mode == 0){//forward
								int ref_position = contig.seq.size() - wordLength;
								int max_search_len = MIN(true_read_position, ref_position);
								read_c += true_read_position - max_search_len;
								ref_c += ref_position - max_search_len;
								int mismatch_number = 0;
								int max_wrong_base = 6;
								for(int i = 0; i < max_search_len && mismatch_number < max_wrong_base; i++)
									if(read_c[i] != ref_c[i] && read_c[i] != 'N')
										mismatch_number++;
								read_will_be_used = (mismatch_number < max_wrong_base)?true:false;
							}else{//reverse
								read_c += reads[read_ID].size() - true_read_position;
								ref_c += wordLength;
								int max_search_len = contig.seq.size() - wordLength;//max_search_len for reference
								max_search_len = MIN(true_read_position, max_search_len);
								int mismatch_number = 0;
								int max_wrong_base = 6;
								for(int i = 0; i < max_search_len && mismatch_number < max_wrong_base; i++)
									if(read_c[i] != ref_c[i] && read_c[i] != 'N')
										mismatch_number++;
								read_will_be_used = (mismatch_number < max_wrong_base)?true:false;
							}
						}
						if(read_will_be_used){
							if(contig.supportReads.insert(read_ID).second)
								contig.actions.emplace_back(kmer_index, read_ID, true);
						}else{
							contig.actions.emplace_back(kmer_index, read_ID, false);
							contig.rejectReads.insert(read_ID);
						}
					}
				}
				// remove reads that do NOT support the selected allele anymore
				for (const unsigned read_ID : supportReads2Remove) {
					if(contig.supportReads.erase(read_ID) != 0){
						contig.actions.emplace_back(kmer_index, read_ID, false);
					}
				}
			}

			// remove the last word from the unused list, so it cannot be used as the seed in finding the next
			// contig
			unusedWords.erase(maxWord);
			// collect the words used to build the contig
			wordsInContig.insert(maxWord);
		}

		// set conservative coverage range for the contig
		if (mode == 0)
			contig.conservativeRange_end = (conservativeEndOffset);
		else
			contig.conservativeRange_bgn = (conservativeEndOffset);
	}
	if(kmer_index > 0) kmer_index = 0;
	contig.conservativeRange_end = (contig.seq.size()
			- contig.conservativeRange_end);

	return isRepeatFound;
}
/// Construct k-mer maps
/// k-mer ==> number of reads containing the k-mer
/// k-mer ==> a list of read IDs containg the k-mer
void MainAssemblyHandler::getKmerCounts(const unsigned wordLength,
		str_uint_map_t &wordCount, str_set_uint_map_t &wordSupportReads) {
	const unsigned readCount = reads.size();
	std::set<std::string> &readWords = getKmerCounts_readWords_BUFF;
	for (unsigned readIndex = 0; readIndex < readCount; ++readIndex){
		// stores the index of a kmer in a read sequence
		const std::string &seq = reads[readIndex];
		const unsigned readLen = seq.size();
		// this read is unusable for assembly:
		if (readLen < wordLength) continue;
		// track all words from the read, including repetitive words
		readWords.clear();
		unsigned kmer_number = readLen - wordLength + 1;
		for (unsigned j = 0; j < kmer_number; ++j) {
			const std::string word(seq.substr(j, wordLength));
			// filter words with "N" (either directly from input alignment or marked due to low basecall quality:
			if (word.find('N') != std::string::npos) continue;
			readWords.insert(word);
		}

		unsigned wordCountAdd = readInfo[readIndex].base_depth;
		// pseudo reads must have passed coverage check with smaller kmers
		// Assigning minCoverage (instead of 1) to a pseudo read allows the pseudo read to rescue the regions
		// where coverage is as low as minCoverage and reads overlap is small
		if (readInfo[readIndex].isPseudo)
			wordCountAdd = o.minCoverage;
		// total occurrences from this read
		for (const std::string &word : readWords) {
			wordCount[word] += wordCountAdd;
			// record the supporting read
			wordSupportReads[word].insert(readIndex);
		}
	}
}

void MainAssemblyHandler::getRepeatKmers(const str_uint_map_t &wordCount,
		std::set<std::string> &repeatWords) {
	str_pair_uint_map_t &wordIndices = getRepeatKmers_wordIndices;
	wordIndices.clear();
	for (const str_uint_map_t::value_type &wdct : wordCount) {
		wordIndices[wdct.first] = std::pair<unsigned, unsigned>(0, 0);
	}
	unsigned index = 1;
	std::vector<std::string> &wordStack = getRepeatKmers_wordStack;
	wordStack.clear();
	for (const str_pair_uint_map_t::value_type &wdidx : wordIndices) {
		const std::string word = wdidx.first;
		const unsigned wordIdx = wdidx.second.first;
		if (wordIdx == 0)
			index = searchRepeatsKmer(o.alphabet, index, word, wordIndices,
					wordStack, repeatWords);
	}
}

bool MainAssemblyHandler::buildContigs(const unsigned wordLength) {
	tmpContigs.clear();
	bool isAssemblySuccess(true);
	// counts the number of occurrences for each kmer in all reads
	str_uint_map_t &wordCount = buildContigs_wordCount;
	wordCount.clear();
	// records the supporting reads for each kmer
	str_set_uint_map_t &wordSupportReads = buildContigs_wordSupportReads;
	wordSupportReads.clear();
	// identify repeat kmers (i.e. circles from the de bruijn graph)
	std::set<std::string> &repeatWords = buildContigs_repeatWords;
	repeatWords.clear();
	std::set<std::string> &unusedWords = buildContigs_unusedWords;

	// get counts and supporting reads for each kmer
	{
		char log_title[1024]; sprintf(log_title, "[getRepeatKmers]");
		double __cpu_time = cputime(); double __real_time = realtime();

		getKmerCounts(wordLength, wordCount, wordSupportReads);
		getRepeatKmers(wordCount, repeatWords);

		fprintf(stderr, "wordLength %d \n%s: CPU time: %.3f sec; real time: %.3f sec\n", wordLength, log_title, cputime() - __cpu_time, realtime() - __real_time);
	}


	// track kmers can be used as seeds for searching for the next contig
	unusedWords.clear();
	for (const str_uint_map_t::value_type &wdct : wordCount) {
		// filter out kmers with too few coverage
		if (wdct.second >= o.minCoverage)
			unusedWords.insert(wdct.first);
	}
	// limit the number of contigs generated for the seek of speed

	{
		char log_title[1024]; sprintf(log_title, "[walk]");
		double __cpu_time = cputime(); double __real_time = realtime();

		int normal_contig = 0;
		while ((!unusedWords.empty())
				&& (normal_contig < 2 * (int)o.maxAssemblyCount)) {
			std::string maxWord;
			unsigned maxWordCount(0);
			// get the kmers corresponding the highest count
			for (auto &word : unusedWords) {
				assert(wordCount.count(word) > 0);
				const unsigned currWordCount = wordCount.at(word);
				if (currWordCount > maxWordCount) {
					maxWord = word;
					maxWordCount = currWordCount;
				}
			}
			// solve for a best contig in the graph by a heuristic greedy maxflow-ish criteria
			AssemblyContig contig;
			bool isRepeatFound = walk(maxWord, wordLength, wordCount,
					wordSupportReads, repeatWords, unusedWords, contig);
			if (isRepeatFound)
				isAssemblySuccess = false;
			if(contig.seq.size() > wordLength)
				normal_contig++;
			tmpContigs.push_back(contig);
		}
		fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
	}


	// done with this now
	return isAssemblySuccess;
}

void MainAssemblyHandler::selectContigs(const unsigned normalReadCount) {
	unsigned finalContigCount(0);
	// a set of reads that has been used to construct contigs, including pseudo ones.
	std::set<unsigned> &usedReads = selectContigs_usedReads;
	usedReads.clear();
	// a set of pseudo reads that has been used to construct contigs
	std::set<unsigned> &usedPseudoReads = selectContigs_usedPseudoReads;
	usedPseudoReads.clear();
	std::set<unsigned> &contigs2Remove = selectContigs_contigs2Remove;
	std::set<unsigned> &newSupportReads = selectContigs_newSupportReads;

	// contig are selected based on the number of supporting reads that are not pseudo
	while ((!tmpContigs.empty()) && (finalContigCount < o.maxAssemblyCount)) {
		// count unused reads that are not pseudo reads
		const unsigned usedNormalReads = usedReads.size()
				- usedPseudoReads.size();
		const unsigned unusedNormalReads = normalReadCount - usedNormalReads;
		if (unusedNormalReads < o.minUnusedReads)
			return;

		unsigned contigIndex(0);
		contigs2Remove.clear();

		AssemblyContig selectedContig;
		unsigned selectedContigIndex;
		unsigned maxSupport(0);
		unsigned maxLength(0);
		for (const AssemblyContig &contig : tmpContigs) {
			// identify new supporting reads that were not used for the previously identified contigs
			newSupportReads.clear();
			std::set_difference(contig.supportReads.begin(),
					contig.supportReads.end(), usedReads.begin(),
					usedReads.end(),
					std::inserter(newSupportReads, newSupportReads.end()));

			// count the number of new supporting reads that are not pseudo reads
			unsigned newNormalSupport(0);
			for (const unsigned rd : newSupportReads) {
				const AssemblyReadInformation &rinfo(readInfo[rd]);
				if (!rinfo.isPseudo)
					newNormalSupport += rinfo.base_depth;
			}
			if (contigs.size() > 0 && newNormalSupport < o.minSupportReads) {
				contigs2Remove.insert(contigIndex);
				contigIndex++;
				continue;
			}
			//contig.new_support_read = newNormalSupport;
			// either more supporting reads that were not used
			// or the same number of supports but longer contig
			const unsigned currNewSupport = newSupportReads.size();
			const unsigned currContigLen = contig.seq.size();
			const bool isBetterContig(
					(currNewSupport > maxSupport)
							|| ((currNewSupport == maxSupport)
									&& (currContigLen > maxLength)));
			if (isBetterContig) {
				selectedContig = contig;
				selectedContig.new_support_read = newNormalSupport;
				selectedContigIndex = contigIndex;
				maxSupport = currNewSupport;
				maxLength = currContigLen;
			}
			contigIndex++;
		}

		// no more contigs selected, selection is done.
		if (maxSupport == 0)
			break;
		// select one contig
		contigs.push_back(selectedContig);
		contigs2Remove.insert(selectedContigIndex);
		// remove selected & failed contigs
		std::set<unsigned>::reverse_iterator it;
		for (it = contigs2Remove.rbegin(); it != contigs2Remove.rend(); it++) {
			tmpContigs.erase(tmpContigs.begin() + *it);
		}
		// update the info about used reads
		for (const unsigned rd : selectedContig.supportReads) {
			usedReads.insert(rd);
			AssemblyReadInformation &rinfo(readInfo[rd]);
			// read info record the ID of all contigs that the read supports
			//rinfo.contigIds.push_back(finalContigCount);
			if (rinfo.isPseudo)
				usedPseudoReads.insert(rd);
		}
		finalContigCount++;
	}
}


void MainAssemblyHandler::assembley() {
	const unsigned normalReadCount(reads.size());

	assert(o.alphabet.size() > 1);
	//readInfo.clear();
	//readInfo.resize(normalReadCount);
	contigs.clear();

	int wordLength_idx = 0;
	unsigned wordLength = o.word_length_list[0];
	bool AssemblySuccessfully = false;
	reachMaxKmerLength = false;
	for (; wordLength <= o.maxWordLength; wordLength_idx++, wordLength = o.word_length_list[wordLength_idx]) {
		if(wordLength == o.maxWordLength){
			reachMaxKmerLength = true;
		}
		AssemblySuccessfully = buildContigs(wordLength);
		if(AssemblySuccessfully == true)
			break;

		// remove pseudo reads from the previous iteration
		const unsigned readCount(reads.size());
		for (unsigned readIndex(0); readIndex < readCount; ++readIndex) {
			AssemblyReadInformation &rinfo(readInfo[readIndex]);
			if (rinfo.isPseudo) {
				reads.erase(reads.begin() + readIndex, reads.end());
				readInfo.erase(readInfo.begin() + readIndex, readInfo.end());
				break;
			}
		}

		//  Add contigs from the current iteration as pseudo reads
		for (const AssemblyContig &contig : tmpContigs) {
			if (contig.seq.size() > (wordLength + o.word_length_list[wordLength_idx])) {
				reads.push_back(contig.seq);
				readInfo.push_back(AssemblyReadInformation(true));
			}
		}

		//select contig when word length is over 100 bp, in this conditions, reads supported a contig will be not enough and
		//longer word length may not generated good results
		//if(wordLength >= 100) selectContigs(normalReadCount);
	}

	// greedy selection of contigs to be returned
	if(AssemblySuccessfully)
		selectContigs(normalReadCount);

	for (AssemblyContig &contig : contigs) {
		contig.debug_print(stderr);
	}
}
