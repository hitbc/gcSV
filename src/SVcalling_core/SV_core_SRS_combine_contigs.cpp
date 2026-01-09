/*
 * SV_core_SRS_combine_contigs.cpp
 *
 *  Created on: 2025-7-9
 *      Author: fenghe
 */

#include <SVcalling_core/SV_core.hpp>

struct Contig_AT_tail{
	bool with_AT_tail;
	int is_tail_not_head;
	int begin_pos;//the begin position of homopolmer
	int end_pos;//the end position of homopolmer
	int contig_ID;
	char repreat_char;
	void show(){
		if(with_AT_tail)
			fprintf(stderr, "Contig_AT_tail: contig_ID %d [%d~%d] %c repeat char %c\n", contig_ID, begin_pos, end_pos, "HT"[is_tail_not_head], repreat_char);
		else{
			fprintf(stderr, "Contig_AT_tail: with out AT_tail\n");
		}
	}
};

static Contig_AT_tail contig_with_AT_TAIL(bool print_log, int contig_ID, AssemblyContig &contig,
		int ori_read_number, std::vector<std::string> &read_list){
	Contig_AT_tail t;
	t.with_AT_tail = false;
	t.contig_ID = contig_ID;
	//fill in the tail of the contig:
	fprintf(stderr, "contig_ID%d [len %ld]\t", contig_ID, contig.seq.size());
	int contig_len = contig.seq.size();
	int contig_middle = contig_len/2;
	for (auto &ca : contig.actions){
		//filters
		if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
		if(contig.remove_read_set.find(ca.read_ID) != contig.remove_read_set.end()) continue;
		if(ca.wrong_base >= MAX_WRONG_BASE) continue;
		const char * read_str = &(read_list[ca.read_ID][0]);
		int read_len = read_list[ca.read_ID].length();
		int max_repeat_n = 0; int max_repeat_end_pos = -1;	int cur_repeat_n = 0; char repreat_char;

		for(int i = 1; i < read_len; i++){
			if(read_str[i] == read_str[i - 1]){
				cur_repeat_n ++;
				if(cur_repeat_n > max_repeat_n){ max_repeat_n = cur_repeat_n; max_repeat_end_pos = i; repreat_char = read_str[i];}
			}
			else
				cur_repeat_n = 0;
		}

		//ca.print(stderr, true);
		if(max_repeat_n >= 15){
			//ca.print(stderr, false);
			int read_position_in_contig = -contig.ass_begin_offset_in_contig + ca.position_in_contig - ca.position_read;
			int repeat_tail_bg = max_repeat_end_pos - max_repeat_n + read_position_in_contig;
			int repeat_tail_ed = max_repeat_end_pos+ read_position_in_contig;
			t.with_AT_tail = true;
			t.begin_pos = repeat_tail_bg;
			t.end_pos = repeat_tail_ed;
			//xassert(t.end_pos != 40, "");
			t.repreat_char = repreat_char;
			if(repeat_tail_bg < contig_middle){
				t.is_tail_not_head = false;
				fprintf(stderr, "%d \t",  repeat_tail_ed);
			}
			else{
				t.is_tail_not_head = true;
				fprintf(stderr, "%d \t",  repeat_tail_bg);
			}
		}
	}

	//check:
	if(t.is_tail_not_head){// repeat at tail
		if(t.end_pos + 80 < (int)contig.seq.size()){
			t.with_AT_tail = false;
			fprintf(stderr, "FAIL repeat not at edge 1\n");
		}
	}else{// repeat at head
		if(t.begin_pos > 80){
			t.with_AT_tail = false;
			fprintf(stderr, "FAIL repeat not at edge 2\n");
		}
	}
	return t;
}

int get_contig_type(bool print_log, AssemblyContig &contig,int ori_read_number,std::vector<ASS_reads_info> &ass_read_list){

	std::vector<char> read_type;
	for (auto &ca : contig.actions){
		if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
		if(contig.remove_read_set.find(ca.read_ID) != contig.remove_read_set.end()) continue;
		if(ca.wrong_base >= MAX_WRONG_BASE) continue;
		char s_t_c;
		switch (ass_read_list[ca.read_ID].signal_type) {
			case Read_type::SR: s_t_c = 'S'; break;
			case Read_type::DR: s_t_c = 'D'; break;
			case Read_type::UM: s_t_c = 'U'; break;
			case Read_type::TL: s_t_c = 'T'; break;
			default:  s_t_c = 'N'; break;
		}
		read_type.emplace_back(s_t_c);
	}
	read_type.emplace_back(0);
	if(print_log) fprintf(stderr, "BEGIN contig_is_part_ASS_check \t");
	if(print_log) fprintf(stderr, " %s \t", &(read_type[0]));

	struct CONTINUE_REGION_ITEM{char c;	int n;	CONTINUE_REGION_ITEM(char c, int n){ this->c = c; this->n = n; }};
	std::vector<CONTINUE_REGION_ITEM> rl;
	int continue_number = 1;
	for(uint i = 1; i < read_type.size(); i++){
		if(read_type[i]== read_type[i-1])
			continue_number++;
		else{
			//MIN: 2
			if(continue_number >= 2){
				char c = read_type[i-1];
				if(rl.size() > 0 && rl.back().c == c)	rl.back().n += continue_number;
				else									rl.emplace_back(read_type[i-1], continue_number);
			}
			continue_number = 1;
		}
	}
	for(CONTINUE_REGION_ITEM & r:rl)
		fprintf(stderr, "[%c %d]\t",r.c, r.n);
	//get type:
	int contig_type = 0;//
	if(rl.size() == 2){
		if(rl[0].c == 'S' && rl[1].c == 'T'){
			contig_type = 1;//S+T mode left
		}else  if(rl[0].c == 'T' && rl[1].c == 'S'){
			contig_type = 2;//T+S mode: right
		}
	}
	if(print_log) fprintf(stderr, "signal summary: %c \t", "NLR"[contig_type]);
	if(print_log) fprintf(stderr, "\n %s\n", contig.seq.c_str());
	return contig_type;
}


static bool contig_signal_of_tran_based_ins(bool print_log, int contig_ID, AssemblyContig &contig,
		int ori_read_number, std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list,
		RefHandler *ref, std::vector<Tran_ins_suggest_pos> & magic_pos_set, int read_len){

	int contig_type = get_contig_type(print_log, contig, ori_read_number, ass_read_list);
	//process the TL reads list:
	if(contig_type == 1 || contig_type == 2){
		for (auto &ca : contig.actions){
			//filters
			if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
			if(contig.remove_read_set.find(ca.read_ID) != contig.remove_read_set.end()) continue;
			if(ass_read_list[ca.read_ID].signal_type != Read_type::TL)  continue;
			if(ca.wrong_base >= MAX_WRONG_BASE) continue;

			//get source positions
			int contig_POS = -contig.ass_begin_offset_in_contig + ca.position_in_contig;
			ASS_reads_info & r_info = ass_read_list[ca.read_ID];
			bool read_is_revered_in_ori_ref_and_contig = r_info.number_NM;
			int contig_pos_in_ref = 0;//the contig position in the reference
			int ref_chr_ID = r_info.read_st_pos;//var reused
			int ref_POS = r_info.read_ed_pos;
			if(read_is_revered_in_ori_ref_and_contig) contig_pos_in_ref = ref_POS + (contig_POS - ca.position_read) + read_len + 1;
			else								      contig_pos_in_ref = ref_POS - (contig_POS - ca.position_read);
			//
			Tran_ins_suggest_pos sug(ref_chr_ID, contig_pos_in_ref, contig_type, read_is_revered_in_ori_ref_and_contig, contig_ID);
			if(magic_pos_set.empty() || !magic_pos_set.back().is_same(sug)) 	magic_pos_set.emplace_back(sug);
			else																magic_pos_set.back().repeat ++;

			//show logs
			if(false){
				fprintf(stderr, "Print TL reads\t");
				//load reference
				int load_ref_length;
				char * c_reference;
				int contig_len = contig.seq.size() + 500;
				if(read_is_revered_in_ori_ref_and_contig){
					c_reference = ref->load_ref_by_region(ref_chr_ID, contig_pos_in_ref - contig_len, contig_pos_in_ref, &load_ref_length);
					getReverseStr_char(c_reference, load_ref_length);
				}else{
					c_reference = ref->load_ref_by_region(ref_chr_ID, contig_pos_in_ref, contig_pos_in_ref + contig_len, &load_ref_length);
				}
				fprintf(stderr, "%d \t", read_is_revered_in_ori_ref_and_contig);
				fprintf(stderr, "%s\t", c_reference);
				free(c_reference);
				fprintf(stderr, "%d magic %d {%d~%d}, \tPOS CONTIG: %d \t POS read %d \t", read_is_revered_in_ori_ref_and_contig, contig_pos_in_ref, r_info.read_st_pos, r_info.read_ed_pos, contig_POS, ca.position_read);
				fprintf(stderr, "%s\t", read_list[ca.read_ID].c_str());
				ca.print(stderr, true);
			}
		}
	}
	return true;
}

void SV_CALLING_Handler::SRS_combine_repeat_tail_of_contigs(bool print_log, int ori_read_number,
		std::vector<ASS_reads_info> &ass_read_list, std::vector<std::string> &read_list,
		std::vector<AssemblyContig> &contigs){
	int target_tid; int target_ref_pos; int target_ref_len;
	ca.get_ref_info(target_tid, target_ref_pos, target_ref_len);

	//reference check before repeat tail combine
	bool reference_check_fail = false;
	{
		uint8_t * target_ref = ca.get_ref();
		int cur_repeat_n = 0; int repeat_string_count[8] = {0};
		for(int i = 1; i < target_ref_len; i++){
			if(target_ref[i] == target_ref[i - 1]){
				cur_repeat_n ++;
			}
			else{
				if(cur_repeat_n >= 13)
					repeat_string_count[target_ref[i-1]]++;
				cur_repeat_n = 0;
			}
		}
		for(int i = 0; i < 4; i++){
			if(repeat_string_count[i] >= 2)
				reference_check_fail = true;
		}
	}
	if(reference_check_fail) {
		fprintf(stderr, "combine_repeat_tail_of_contigs: reference_check_failm Return\n");
		return;
	}

	//call variant using the combination of multiple Contigs
	//check the AT_tail
	std::vector<Contig_AT_tail> contig_AT_l;
	int contig_ID = -1;
	for (AssemblyContig &contig : contigs) {
		contig_ID++;
		if(contig.contig_is_filtered) continue;
		if(contig.ending_reason[0] != 0 ||  contig.ending_reason[1] != 0) continue;
		contig_AT_l.emplace_back(contig_with_AT_TAIL(print_log, contig_ID, contig, ori_read_number, read_list));
	}
	Contig_AT_tail *first_H = NULL; Contig_AT_tail *first_T = NULL;
	//contig combine using AT_tail
	for(Contig_AT_tail &a: contig_AT_l){
		a.show();
		if(a.with_AT_tail && a.is_tail_not_head && first_T == NULL && !contigs[a.contig_ID].contig_is_filtered){  first_T = &a; }
		if(a.with_AT_tail && !a.is_tail_not_head && first_H == NULL && !contigs[a.contig_ID].contig_is_filtered){ first_H = &a; }
	}
	if(first_T != NULL && first_H != NULL && first_T->repreat_char == first_H->repreat_char && first_H->end_pos >= 0){
		fprintf(stderr, "Try to combine: %d + %d \n", first_T->contig_ID, first_H->contig_ID );
		//combine contig string
		AssemblyContig & contig_L = contigs[first_T->contig_ID];
		AssemblyContig & contig_R = contigs[first_H->contig_ID];

		std::string & contig_seq_L = contig_L.seq;//the first part of combined strings
		std::string & contig_seq_R = contig_R.seq;//the second part of combined strings
		std::string combine_str = contig_seq_L.substr(0, first_T->begin_pos);

		std::vector<char> AT_string; char repeat_char = first_T->repreat_char;
		int repeat_string_size = std::max(first_T->end_pos - first_T->begin_pos, first_H->end_pos - first_H->begin_pos);
		for(int i = 0; i < repeat_string_size; i++){
			AT_string.emplace_back(repeat_char);
		}
		AT_string.emplace_back(0);
		combine_str.append(&(AT_string[0]));
		//combine the read list
		int contig_R_offset = combine_str.size();
		int contig_R_suggset_offset = combine_str.size();
		std::string combine_str_R;
//		if(first_H->end_pos < 0){
//			for(int i = first_H->end_pos; i < 0; i++)
//				combine_str_R.push_back(first_H->repreat_char);
//			first_H->end_pos = 0;
//		}
		//fprintf(stderr, "%d %d\n", first_H->end_pos, contig_seq_R.size());
		if(first_H->end_pos > contig_seq_R.size())
			return;
		//xassert(first_H->end_pos <= contig_seq_R.size(), "");//this assert is used to make sure that: the RIGHT-part SEQ is not empty
		combine_str_R += contig_seq_R.substr(first_H->end_pos, contig_seq_R.size());
		contig_R_offset -= (first_H->end_pos) + (contig_R.ass_begin_offset_in_contig - contig_L.ass_begin_offset_in_contig);
		contig_R_suggset_offset -= (first_H->end_pos);
		combine_str += combine_str_R;
		//
		//show
		if(print_log)fprintf(stderr, "combine_str is %s\n", combine_str.c_str());
		//combine the action list:

		for (AssemblyReadAction &ca : contig_R.actions){
			//filters
			if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
			if(contig_R.remove_read_set.find(ca.read_ID) != contig_R.remove_read_set.end()) continue;
			if(ca.wrong_base >= MAX_WRONG_BASE) continue;
			ca.position_in_contig += contig_R_offset;
			ca.suggest_contig_offset_in_ref -= contig_R_suggset_offset;
			contig_L.actions.emplace_back(0,0,0);
			std::swap(contig_L.actions.back(), ca);
		}
		std::swap(combine_str, contig_L.seq);
		for(auto & s: contig_R.SUGGEST_pos_list)
			s.suggest_pos -= contig_R_suggset_offset;
		contig_L.SUGGEST_pos_list.insert(contig_L.SUGGEST_pos_list.end(), contig_R.SUGGEST_pos_list.begin(), contig_R.SUGGEST_pos_list.end());
		if(contig_L.SUGGEST_pos_list.empty())
			return;

		if(print_log){
			contig_L.debug_print(stderr);
			for (auto &ca : contig_L.actions){
				if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
				{
					int contig_pos = -contig_L.ass_begin_offset_in_contig + ca.position_in_contig;
					if(ca.isAdd){
						for(int i = 0; i < contig_pos - ca.position_read ;i++)fprintf(stderr, " ");
						fprintf(stderr, "%s\t", read_list[ca.read_ID].c_str());
					}
					ass_read_list[ca.read_ID].print(stderr, 0);
					fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
					ca.print(stderr, true);
				}
			}
		}
		//get variants
		store_bin_contig(contig_L.seq, bin_contig);
		SRS_getSuggestSVlength(contig_L);
		SRS_alignment_and_get_var(print_log, 0, contig_ID, contig_L.SUGGEST_pos_list[0].suggest_pos - target_ref_pos, contig_L.seq.size(), target_ref_len);
		//erase contig R
		contigs.erase(contigs.begin() + first_H->contig_ID);
	}
	else{
		fprintf(stderr, "Combine fail\n");
	}
}

void SV_CALLING_Handler::SRS_get_BP_from_tran_based_ins(bool print_log, std::vector<Tran_ins_suggest_pos> &suggest_pos, std::vector<AssemblyContig> &contigs, std::vector<TRANS_INS_info> &bp_l){
	int target_tid; int target_ref_pos; int target_ref_len;
	ca.get_ref_info(target_tid, target_ref_pos, target_ref_len);
	uint8_t * target_ref = ca.get_ref();
	if(false){
		fprintf(stderr, "Target ref basic: [%d:%d-%d] \n", target_tid, target_ref_pos, target_ref_pos + target_ref_len);
		for(int i = 0; i < target_ref_len; i++) fprintf(stderr, "%c", "ACGT"[ ca.getTseq(i)]);
	}
	//summary : magic_POS
	//Tran_ins_suggest_pos * sug_L = NULL;
	//Tran_ins_suggest_pos * sug_R = NULL;

	int sug_L_idx = -1;
	int sug_R_idx = -1;

	for(uint i = 1; i < suggest_pos.size(); i++){
		if(suggest_pos[i].is_paired(suggest_pos[i - 1])){
			sug_L_idx = i - 1;
			sug_R_idx = i;
			if(suggest_pos[sug_L_idx].contig_type  != 1) std::swap(sug_L_idx, sug_R_idx);
//			sug_L = &(suggest_pos[i - 1]);
//			sug_R = &(suggest_pos[i]);
//			if(sug_L->contig_type != 1)	std::swap(sug_L, sug_R);
			fprintf(stderr, "\n\nPosition paired: \n");
			break;
		}
	}

	if(sug_L_idx == -1 || sug_R_idx == -1){//no paired
		//search the position with the max repeat number
		int max_repeat_idx = -1; int max_repeat_n = 0;
		for(uint i = 0; i < suggest_pos.size(); i++){
			if(suggest_pos[i].repeat > max_repeat_n){
				max_repeat_idx = i;
				max_repeat_n = suggest_pos[i].repeat;
			}
		}

		//try the get new paired positions
		//load the reference near the positions

		if(max_repeat_idx != -1){
			fprintf(stderr, "Load MAX Trans regions\t");
			Tran_ins_suggest_pos & max_sug = suggest_pos[max_repeat_idx];
			int contig_base_len = contigs[max_sug.contig_ID].seq.size();
			//load reference
			int source_load_len = 0;
			char * source_reference = NULL;
			int ref_chr_ID = max_sug.ref_chr_ID;
			int suggest_load_length = 10000;
			int begin_load_position = max_sug.SUG_CON_POS_in_ref;
			//contig_type = 1;//S+T mode left ; contig_type = 2;//T+S mode: right
			if(max_sug.contig_type == 1 && !max_sug.contig_is_revered_to_ref){//L + forward
				begin_load_position = max_sug.SUG_CON_POS_in_ref;
			}else if(max_sug.contig_type == 1 && max_sug.contig_is_revered_to_ref){//L + reverse
				begin_load_position = max_sug.SUG_CON_POS_in_ref - suggest_load_length;
			}else if(max_sug.contig_type == 2 && !max_sug.contig_is_revered_to_ref){//R + forward
				begin_load_position = max_sug.SUG_CON_POS_in_ref- suggest_load_length;
			}
			else{ //R + reverse
				begin_load_position = max_sug.SUG_CON_POS_in_ref;
			}
			begin_load_position -= contig_base_len;
			int end_load_position = begin_load_position + suggest_load_length;
			if(begin_load_position < 1)
				begin_load_position = 1;
			source_reference = ref_handler->load_ref_by_region(ref_chr_ID, begin_load_position, end_load_position, &source_load_len);
			if(max_sug.contig_is_revered_to_ref) getReverseStr_char(source_reference, source_load_len);
			std::vector<uint8_t> source_reference_bin;
			if(source_reference == NULL || source_load_len == 0)
				return;
			//fprintf(stderr, "S3, source_reference %ld source_load_len %d\n", source_reference, source_load_len);
			store_bin_from_char(source_reference, source_load_len, source_reference_bin);//1

			//get the other contig:
			for(uint i = 0; i < suggest_pos.size(); i++){
				if(suggest_pos[i].contig_ID != max_sug.contig_ID && suggest_pos[i].contig_type != max_sug.contig_type){
					if(sug_L_idx != -1 && sug_R_idx != -1)
						break;
					//search the first contig that type is not same
					//
					AssemblyContig & mate_contig = contigs[suggest_pos[i].contig_ID];
					//align the contig to the reference ?? how to ??
					//show the reference and the contig
					if(print_log){
						fprintf(stderr, "Contig try to located: %s\n\n", mate_contig.seq.c_str());
						fprintf(stderr, "Load_reference: %s\n\n", source_reference);
						fprintf(stderr, "Contig base: %s\n\n", contigs[max_sug.contig_ID].seq.c_str());
					}
					store_bin_contig(mate_contig.seq, bin_contig);
					//alignment for it
					ca_re_locate.setRef(&(source_reference_bin[0]), source_load_len, ref_chr_ID, begin_load_position);
					//fprintf(stderr, "S5, &(bin_contig[0]) %ld, mate_contig.seq.size()  %ld , 0, source_load_len  %ld\n", &(bin_contig[0]), mate_contig.seq.size(), 0, source_load_len);
					ca_re_locate.align_non_splice_default_ref(&(bin_contig[0]), mate_contig.seq.size(), 0, source_load_len, 400);
					int contig_st_in_ref = ca_re_locate.adjustCIGAR();
					if(print_log) ca_re_locate.printf_alignment_detail(stderr, contig_st_in_ref, NULL, 0);
					if(false){
						int match_size = 0;
						int contig_NM = ca_re_locate.get_contig_NM(contig_st_in_ref, match_size);
						fprintf(stderr, "Contig filter: contig_NM %d match_size %d \n", contig_NM, match_size);
						if(contig_NM * 20 > match_size){
							fprintf(stderr, "Contig filter(try the get paired positions): Too mant NM: contig_NM %d error rate: %f%%(over 5%)\n", contig_NM, (float)contig_NM/match_size*100);
							continue;
						}
					}
					int seq_i = 0;
					uint32_t* bam_cigar;
					int cigar_len = ca_re_locate.get_cigar(&bam_cigar);
					int output_index = contig_st_in_ref;
					for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
						//if(sug_L_idx != -1 && sug_R_idx != -1)
//							break;
						int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
						int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
						switch(type){
						case 0:	case 7:	case 8:	case 3:	case 4:
							if(cigar_len > 50){
								int contig_type = (max_sug.contig_type == 1)?2:1;
								int contig_pos_in_ref = 0;
								if(max_sug.contig_is_revered_to_ref)//reverse
									contig_pos_in_ref = begin_load_position + source_load_len - (output_index-seq_i);
								else//forward
									contig_pos_in_ref = begin_load_position+(output_index-seq_i);
								suggest_pos.emplace_back(ref_chr_ID, contig_pos_in_ref, contig_type, max_sug.contig_is_revered_to_ref, suggest_pos[i].contig_ID);

								sug_L_idx = max_repeat_idx;
								sug_R_idx = suggest_pos.size() - 1;
								//sug_L = &(max_sug);
								//sug_R = &(suggest_pos.back());
								if(suggest_pos[sug_L_idx].contig_type != 1)	std::swap(sug_L_idx, sug_R_idx);

							}
							output_index += cigar_len; seq_i += cigar_len; break;//match
							case 1: seq_i += cigar_len; break;//insertion
							case 2:	output_index += cigar_len; break; //deletion
							default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len); break;
						}
					}
					break;
				}
			}
			free(source_reference);
		}
	}
	bp_l.clear();
	if(sug_L_idx != -1 && sug_R_idx != -1){//successfully paired
		fprintf(stderr, "get BP in paired ends\t");
		for(int mode = 0; mode < 2; mode++){
			int cur_sug_idx = (mode ==0)?sug_L_idx:sug_R_idx;
			Tran_ins_suggest_pos * sug_c = &(suggest_pos[cur_sug_idx]);
			AssemblyContig & contig = contigs[sug_c->contig_ID];
			std::vector<SUGGEST_POS_LIST_ITEM> &local_ref_l = contig.SUGGEST_pos_list;
			if(local_ref_l.empty())
				continue;
			int contig_len = contig.seq.size();
			int source_tid = sug_c->ref_chr_ID;
			int source_pos = sug_c->SUG_CON_POS_in_ref;
			int source_load_len = 0;
			char *source_reference;
			int source_load_bg = 0;
			if(sug_c->contig_is_revered_to_ref){
				source_load_bg = source_pos - contig_len;
				if(source_load_bg < 1)	source_load_bg = 1;
				source_reference = ref_handler->load_ref_by_region(source_tid, source_load_bg, source_load_bg + contig_len, &source_load_len);
				getReverseStr_char(source_reference, source_load_len);
			}else{
				source_load_bg = source_pos;
				if(source_load_bg < 1)	source_load_bg = 1;
				source_reference = ref_handler->load_ref_by_region(source_tid, source_load_bg, source_load_bg + contig_len, &source_load_len);
			}
			std::vector<uint8_t> source_reference_bin;
			if(source_reference == NULL)
				return;
			store_bin_from_char(source_reference, source_load_len, source_reference_bin);
			//target:
			//
			int contig_pos_in_target = local_ref_l[0].suggest_pos - target_ref_pos;
			contig_pos_in_target = MAX(contig_pos_in_target, 0);
			std::vector<uint8_t> target_reference_bin;
			if(target_ref_len > contig_pos_in_target){
				target_reference_bin.resize(target_ref_len - contig_pos_in_target);
				memcpy(&(target_reference_bin[0]), target_ref + contig_pos_in_target, target_ref_len- contig_pos_in_target);
			}
			//combine:
			SRS_REF_STRING_COMBINE_Handler ref_combine;
			int aln_target_pos = target_ref_pos + contig_pos_in_target;
			int aln_source_pos = source_load_bg;
			//LEFT and RIGHT
			if(mode ==0) ref_combine.store_ref(target_reference_bin, true, aln_target_pos,source_reference_bin ,!sug_c->contig_is_revered_to_ref, aln_source_pos, false); //LEFT
			else 		 ref_combine.store_ref(source_reference_bin, !sug_c->contig_is_revered_to_ref, aln_source_pos, target_reference_bin,true, aln_target_pos, false); //RIGHT

			ca_bnd.setRef(&(ref_combine.s[0]), ref_combine.s.size(), 0,0);
			//get bin contig
			store_bin_contig(contig.seq, bin_contig);
			ca_bnd.align_non_splice_default_ref(&(bin_contig[0]), contig.seq.size(), 0, ca_bnd.get_tlen(), 400);
			int contig_st_in_ref = ca_bnd.adjustCIGAR();
			if(print_log && false) ca_bnd.printf_alignment_detail(stderr, contig_st_in_ref, NULL, 0);
			//analysis alignment detail
			int seq_i = 0;
			uint32_t* bam_cigar;
			int cigar_len = ca_bnd.get_cigar(&bam_cigar);
			int output_index = contig_st_in_ref;
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
					case 0:
					case 7:
					case 8:
					case 3:
					case 4:
						output_index += cigar_len; seq_i += cigar_len; break;//match
					case 1: seq_i += cigar_len; break;//insertion
					case 2://deletion:
					{
						int source_BP_pos; int target_BP_pos;
						if(ref_combine.check_and_store_TRAN_INS(output_index, output_index + cigar_len, mode, source_BP_pos, target_BP_pos)){
							bp_l.emplace_back(); //store the trans ins results
							bp_l.back().set(&contig, !sug_c->contig_is_revered_to_ref, mode, seq_i, source_tid, source_BP_pos, target_tid, target_BP_pos);
						}
					}
					output_index += cigar_len;
					break;
					default: fprintf(stderr, "ERROR CIGAR  %d %d ", type, cigar_len); break;
				}
			}
			if(print_log && false){
				fprintf(stderr, "BP %c \t source ref:", "LR"[mode]);
				sug_c->print();
				fprintf(stderr, "\t target ref: ");
				for(SUGGEST_POS_LIST_ITEM &p: local_ref_l){
					p.printf(stderr);
				}
				fprintf(stderr, "\n");
				ref_combine.print_ref();
			}
		}
	}
}

void SV_CALLING_Handler::SRS_handle_tran_based_ins(bool print_log, int ori_read_number, std::vector<ASS_reads_info> &ass_read_list,
		std::vector<std::string> &read_list, std::vector<AssemblyContig> &contigs)
{
	uint32_t SV_size_ori = SV_result_TMP_buff.size();
	int target_tid; int target_ref_pos; int target_ref_len;
	ca.get_ref_info(target_tid, target_ref_pos, target_ref_len);
	int contig_ID = 0;
	//call variant using the combination of multiple Contigs
	std::vector<Tran_ins_suggest_pos> suggest_pos;
	bool with_part_ass_contig = false;
	for (AssemblyContig &contig : contigs) {
		if(contig.ending_reason[0] != 0 ||  contig.ending_reason[1] != 0) continue;
		if(contig.contig_is_filtered) continue;
		with_part_ass_contig |= contig_signal_of_tran_based_ins(print_log, contig_ID, contig, ori_read_number, ass_read_list,read_list, ref_handler, suggest_pos, sig_para.MaxReadLen);
		contig_ID++;
	}
	std::sort(suggest_pos.begin(), suggest_pos.end(), Tran_ins_suggest_pos::cmp_by_ref_position);
	//remove duplication positions
	for(uint i = 1; i < suggest_pos.size(); i++){
		if(suggest_pos[i].is_same(suggest_pos[i-1])){
			suggest_pos[i - 1].repeat += suggest_pos[i].repeat;
			suggest_pos.erase(suggest_pos.begin()+ i);
		}
	}
	std::vector<TRANS_INS_info> bp_l;
	SRS_get_BP_from_tran_based_ins(print_log, suggest_pos, contigs, bp_l);

	//generate variants
	if(with_part_ass_contig){
		//begin the calling process
		if(print_log){
			for(TRANS_INS_info & bp: bp_l)
				bp.print();
		}

		//generate SVs:
		if(bp_l.size() == 2 &&  bp_l[0].is_R_not_L == false && bp_l[1].is_R_not_L == true){
			fprintf(stderr, "Source: [%d:%d ~ %d:%d] len %d \t", bp_l[0].source_tid, bp_l[0].source_BP_pos, bp_l[1].source_tid, bp_l[1].source_BP_pos, bp_l[1].source_BP_pos - bp_l[0].source_BP_pos);
			fprintf(stderr, "Target: [%d:%d ~ %d:%d] len %d \t", bp_l[0].target_tid, bp_l[0].target_BP_pos, bp_l[1].target_tid, bp_l[1].target_BP_pos, bp_l[1].target_BP_pos - bp_l[0].target_BP_pos);
			fprintf(stderr, "DIR: [%c %c]\n", "RF"[bp_l[0].source_is_forward], "RF"[bp_l[1].source_is_forward]);
			//S1: generate combined contigs
			//
			//bp_l[0].source_tid, bp_l[0].source_BP_pos, bp_l[1].source_tid, bp_l[1].source_BP_pos;
			AssemblyContig & contig_L = *(bp_l[0].contig_p);
			AssemblyContig & contig_R = *(bp_l[1].contig_p);

			std::string & contig_seq_L = contig_L.seq;
			std::string & contig_seq_R = contig_R.seq;
			std::string contig_string_middle;
			std::string combine_str = contig_seq_L;

			//get the middle string:
			//from t1 to t2:

			int source_tid = bp_l[0].source_tid;
			int source_load_bg = bp_l[0].source_BP_pos; int source_load_ed = bp_l[1].source_BP_pos; int source_load_len;
			if(source_load_bg > source_load_ed)
				std::swap(source_load_bg, source_load_ed);
			char * source_reference = ref_handler->load_ref_by_region(source_tid, source_load_bg, source_load_ed, &source_load_len);
			if(!bp_l[0].source_is_forward)
				getReverseStr_char(source_reference, source_load_len);
			contig_string_middle.append(source_reference);
			// combine the middle:
			int replace_L_len = bp_l[0].source_contig_ed - bp_l[0].source_contig_bg;
			int replace_R_len = bp_l[1].source_contig_ed - bp_l[1].source_contig_bg;
			int middle_len = contig_string_middle.size();
			if(middle_len - replace_L_len - replace_R_len <= 0){
				int remove_size = -( middle_len - replace_L_len - replace_R_len);
				int remove_pos = combine_str.size() - remove_size;
				if((int)combine_str.size() < remove_size){
					combine_str.clear();
				}
				else
					combine_str.erase(remove_pos,remove_size);
			}
			else{
				//xassert(replace_L_len < contig_string_middle.size(),"");
				combine_str += contig_string_middle.substr(replace_L_len, middle_len - replace_L_len - replace_R_len);
			}
			//combine the string R
			int contig_R_offset = combine_str.size() - (contig_R.ass_begin_offset_in_contig - contig_L.ass_begin_offset_in_contig);//;// - bp_l[0].source_contig_bg
			int contig_R_suggest_offset = combine_str.size();//;// - bp_l[0].source_contig_bg

			//combine the right
			combine_str += contig_seq_R;
			//store the new string to contig_L:
			std::swap(contig_L.seq, combine_str);
			//combine the action list:
			for (AssemblyReadAction &ca : contig_R.actions){
				//filters
				if(ca.read_ID >= ori_read_number || !ca.isAdd)continue;
				if(contig_R.remove_read_set.find(ca.read_ID) != contig_R.remove_read_set.end()) continue;
				if(ca.wrong_base >= MAX_WRONG_BASE) continue;
				ca.position_in_contig += contig_R_offset;
				ca.suggest_contig_offset_in_ref -= contig_R_suggest_offset;
				contig_L.actions.emplace_back(0,0,0);
				std::swap(contig_L.actions.back(), ca);
			}
//
			for(auto & s: contig_R.SUGGEST_pos_list)
				s.suggest_pos -= contig_R_suggest_offset;
			contig_L.SUGGEST_pos_list.insert(contig_L.SUGGEST_pos_list.end(), contig_R.SUGGEST_pos_list.begin(), contig_R.SUGGEST_pos_list.end());

			if(print_log){
				fprintf(stderr, "contig_L.debug_print");
				contig_L.debug_print(stderr);
			}

			if(print_log){
				for (auto &ca : contig_L.actions){
					if(ca.read_ID < ori_read_number)// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
					{
						int contig_pos = -contig_L.ass_begin_offset_in_contig + ca.position_in_contig;
						if(ca.isAdd){
							for(int i = 0; i < contig_pos - ca.position_read ;i++)fprintf(stderr, " ");
							fprintf(stderr, "%s\t", read_list[ca.read_ID].c_str());
						}
						ass_read_list[ca.read_ID].print(stderr, 0);
						fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
						ca.print(stderr, true);
					}
				}
			}
//			//get variants
			store_bin_contig(contig_L.seq, bin_contig);
			SRS_getSuggestSVlength(contig_L);

			if(middle_len > 800){
				if(middle_len < 5000)
					set_region_addition_load_extramly_long();
				else if(middle_len < 50000)
					set_region_addition_load_super_super_long();
			}

			SRS_alignment_and_get_var(print_log, 0, -1,
					contig_L.SUGGEST_pos_list[0].suggest_pos - target_ref_pos, contig_L.seq.size(), target_ref_len);
			//
			//filters
			//remove the contig R
			contigs.erase(contigs.begin() + (&contig_R - &(contigs[0])));
		}else{
			fprintf(stderr, "Source and target un-balance\n");
			//generate:
		}
	}

	//filters
	if(SV_result_TMP_buff.size() != SV_size_ori){
		for(uint i = SV_size_ori; i < SV_result_TMP_buff.size(); i++){
			SV_result_TMP_buff[i].set_SVLEN_is_Imprecise();
		}
	}
}

