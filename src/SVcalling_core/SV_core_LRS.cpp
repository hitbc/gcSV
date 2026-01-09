/*
 * SveHandler_LRS.cpp
 *
 *  Created on: 2025/7/9
 *      Author: fenghe
 */

#include <SVcalling_core/SV_core.hpp>

void SPECIAL_BRANCH_HANDLER::special_branch_check(bool print_log, int minDepth_LRS){
	need_reclassify = false;
	int support_read_num = 1;
	int begin_idx = 0;
	std::sort(long_branch_check_list.begin(), long_branch_check_list.end(), LONG_BRANCH_CHECK::cmp_by_uid);
	if(print_log){
		for(uint i = 0; i < long_branch_check_list.size(); i++)
			long_branch_check_list[i].show();
	}

	for(uint i = 1; i < long_branch_check_list.size() + 1; i++){
		if(i < long_branch_check_list.size() && long_branch_check_list[i].query_is_absorb == true &&
				long_branch_check_list[i].is_similar(long_branch_check_list[begin_idx])){
			//long_branch_check_list[i].show();
			support_read_num ++;
			if(print_log)fprintf(stderr, "SAME: support_read_num %d \n", support_read_num);
		}else{
			if(print_log)fprintf(stderr, "RST: support_read_num %d \n", support_read_num);
			if(support_read_num >= minDepth_LRS){
				LONG_BRANCH_CHECK & c = long_branch_check_list[begin_idx];
				special_branch_check_list.emplace_back();
				special_branch_check_list.back().set(c.query_total_node, c.targt_total_node, c.begin_UID);
				need_reclassify = true;
			}
			begin_idx = i;
			support_read_num = 1;
		}
	}
}

void TARGET_path_and_index::get_the_contig_string(bool print_log, std::vector<UNITIG_INFO_item> &unitig_l,
		std::vector<std::string> & read_string){
	print_log = false;
	xassert(!unitig_l.empty(), "");
	int kmer_size = unitig_l[0].unitig.size()  - unitig_l[0].kmer_list.size() + 1;
	bool is_head = true;
	contig.clear();
	//todo::
	if(print_log)
		fprintf(stderr, "\nGenerate the contig strings: \n");
	int total_kmer_size = 0;
	int old_read_id = -1;
	for(int target_node_id = the_bg_node_id; target_node_id != -1; target_node_id = target_link_path[target_node_id].next_node_id){
		ASS_TARGET_ITEM & a = target_link_path[target_node_id];
		int uid = a.UID;
		if(uid == -1)
			continue;
		std::string add_string;
		//loading the first kmer in the head:

		std::string & unitig_str = unitig_l[uid].unitig;
		if(is_head){
			//kmer: 19 bp:
			add_string += unitig_str.substr(0, kmer_size - 1);
			if(print_log) fprintf(stderr, "0: %s\n", add_string.c_str());
			is_head = false;
		}
		bool FULL_conver = (a.U_pos_bg == 0 && a.U_pos_ed == a.unitig_kmer_number);
		if(FULL_conver){
			add_string += unitig_str.substr(kmer_size - 1, unitig_str.size() - kmer_size + 1);
			if(print_log) fprintf(stderr, "1: %s\n", add_string.c_str());
			if(print_log && false){
				fprintf(stderr, "support read list: \n");
				for(READ_PATH_COMBINE_ITEM * query_p : a.query_node_info){
					fprintf(stderr, "Read ID %d; uid %d, R_pos:%d-%d U_pos:%d-%d\n", query_p->read_id, query_p->UID, query_p->read_pos_bg, query_p->read_pos_ed,  query_p->U_pos_bg,  query_p->U_pos_ed);
				}
			}
		}
		else{
			int bg_pos_unitig = a.U_pos_bg + kmer_size - 1;
			int length_unitig = a.U_pos_ed - a.U_pos_bg;
			bool read_switch = (a.ori_read_id != old_read_id);
			if(read_switch){
				bg_pos_unitig = kmer_size - 1;
				length_unitig = a.U_pos_ed;
			}
			add_string += unitig_str.substr(bg_pos_unitig, length_unitig);
			if(print_log){
				fprintf(stderr, "2: read_switch %d, new_read_id %d, old_read_id %d, %s\n", read_switch, a.ori_read_id, old_read_id, add_string.c_str());
				fprintf(stderr, "support read list: \n");
				for(READ_PATH_COMBINE_ITEM * query_p : a.query_node_info){
					fprintf(stderr, "Read ID %d; uid %d, R_pos:%d-%d U_pos:%d-%d\n", query_p->read_id, query_p->UID, query_p->read_pos_bg, query_p->read_pos_ed,  query_p->U_pos_bg,  query_p->U_pos_ed);
				}
			}
		}
		old_read_id = a.ori_read_id;
		bool with_outside_gap = (a.gap_kmer_num_outside != 0);
		//string out of the unitig:
		if(with_outside_gap){
			std::string & c_r = read_string[a.ori_read_id];
			int pos_outside_gap = a.read_pos_ed + kmer_size - 1;
			int outside_gap_size = a.gap_kmer_num_outside;
			add_string += c_r.substr(pos_outside_gap, outside_gap_size);
			if(print_log) fprintf(stderr, "3: %s\n", add_string.c_str());
		}
		contig += add_string;
		if(print_log) fprintf(stderr, "4: %s\n", contig.c_str());
		if(print_log){
			total_kmer_size += a.kmer_number + a.gap_kmer_num_inside + a.gap_kmer_num_outside;
			bool with_outside_gap = (a.gap_kmer_num_outside != 0);
			fprintf(stderr, "tid: %d con_len %ld kmer_size %d p_diff: %ld u_diff: %d uid %d (%d:%d:%d) U_POS(%d-%d) U_len %d R_POS(%d:%d-%d)",
					target_node_id, contig.size(),
					total_kmer_size,
					contig.size() - total_kmer_size - kmer_size + 1,
					a.unitig_kmer_number - a.kmer_number - a.gap_kmer_num_inside - a.gap_kmer_num_outside,
					uid,
					a.kmer_number, a.gap_kmer_num_inside, a.gap_kmer_num_outside,
					a.U_pos_bg, a.U_pos_ed,a.unitig_kmer_number,
					a.ori_read_id, a.read_pos_bg, a.read_pos_ed
			);
			if(!FULL_conver)
				fprintf(stderr, "NOT_full_conver;");
			if(with_outside_gap)
				fprintf(stderr, "with_outside_gap;");
			//std::string & c_r = read_string[a.ori_read_id];
			//std::string new_str = c_r.substr(a.read_pos_bg, a.read_pos_ed - a.read_pos_bg);
			fprintf(stderr, "ADD:%s:UNI:%s\n", add_string.c_str(), unitig_l[target_link_path[target_node_id].UID].unitig.c_str());
		}
	}
	if(print_log)
		fprintf(stderr, "final contig strings: %s\n",contig.c_str());
}


void TARGET_path_and_index::show_the_final_result(bool print_log, int cur_target_ID, std::vector<UNITIG_INFO_item> &unitig_l){
	int kmer_size = unitig_l[0].unitig.size()  - unitig_l[0].kmer_list.size() + 1;
	fprintf(stderr, "Show_the_target_path [ID%d]: LEN: %ld \n", cur_target_ID, contig.size());
	if(print_log){
		show_total_length(unitig_l);
		for(int target_node_id = the_bg_node_id; target_node_id != -1; target_node_id = target_link_path[target_node_id].next_node_id){
			int uid = target_link_path[target_node_id].UID;
			if(uid == -1)
				continue;
			fprintf(stderr, "%d:%d:%d:%d(%d):%d:%d:%s:%s:%d-->\n",
					target_link_path[target_node_id].node_idx_position_in_path,
					target_link_path[target_node_id].node_base_position_in_path,
					uid,
					target_link_path[target_node_id].kmer_number,
					(int)unitig_l[uid].unitig.size() -target_link_path[target_node_id].kmer_number,
					target_link_path[target_node_id].gap_kmer_num_inside,
					target_link_path[target_node_id].gap_kmer_num_outside,
					unitig_l[uid].unitig.c_str(),
					unitig_l[uid].unitig.c_str() + kmer_size - 1,
					unitig_l[uid].with_ref);
		}
		fprintf(stderr, "\n final contig is %s\n", contig.c_str());
	}
	//show_the_target_path();
	show_read_list();
	fprintf(stderr, "\n");
}

void TARGET_path_and_index::show_read_list(){
	fprintf(stderr, "The support read is: \n");
	std::sort(support_read_l.begin(), support_read_l.end(), SUPPORT_READ_LIST_ITEM::cmp_by_readid);
	for(SUPPORT_READ_LIST_ITEM & r : support_read_l){
		r.show();
	}
	fprintf(stderr, "\n");
}

void TARGET_path_and_index::show_the_target_path(){
	fprintf(stderr, "Show_the_target_path: \n");
	for(int target_node_id = the_bg_node_id; target_node_id != -1; target_node_id = target_link_path[target_node_id].next_node_id){
		fprintf(stderr, "%d:%d:%d:%d:%d:%d-->\n",
				target_link_path[target_node_id].node_idx_position_in_path,
				target_link_path[target_node_id].node_base_position_in_path,
				target_link_path[target_node_id].UID,
				target_link_path[target_node_id].kmer_number,
				target_link_path[target_node_id].gap_kmer_num_inside,
				target_link_path[target_node_id].gap_kmer_num_outside);
	}
	fprintf(stderr, "\n");
}

void TARGET_path_and_index::show_total_length(std::vector<UNITIG_INFO_item> &unitig_l){
	fprintf(stderr, "\nshow_total_length:\n");
	xassert(!unitig_l.empty(), "");
	int contig_total_len = 0;
	int kmer_size = unitig_l[0].unitig.size()  - unitig_l[0].kmer_list.size() + 1;
	bool is_head = true;
	int total_kmer_size = 0;
	for(int target_node_id = the_bg_node_id; target_node_id != -1; target_node_id = target_link_path[target_node_id].next_node_id){
		ASS_TARGET_ITEM & a = target_link_path[target_node_id];
		int uid = a.UID;
		total_kmer_size += a.kmer_number + a.gap_kmer_num_inside + a.gap_kmer_num_outside;
		if(uid == -1)
			continue;
		if(is_head){
			contig_total_len += unitig_l[uid].unitig.size();
			is_head = false;
		}else
			contig_total_len += unitig_l[uid].unitig.size() - kmer_size + 1;
		bool NOT_full_conver = (a.U_pos_bg != 0 || a.U_pos_ed != a.unitig_kmer_number);
		bool with_outside_gap = (a.gap_kmer_num_outside != 0);

		fprintf(stderr, "tid: %d con_len %d kmer_size %d p_diff: %d u_diff: %d uid %d (%d:%d:%d) U_POS(%d-%d) U_len %d R_POS(%d:%d-%d)",
				target_node_id, contig_total_len,
				total_kmer_size,
				contig_total_len - total_kmer_size - kmer_size + 1,
				a.unitig_kmer_number - a.kmer_number - a.gap_kmer_num_inside - a.gap_kmer_num_outside,
				uid,
				a.kmer_number, a.gap_kmer_num_inside, a.gap_kmer_num_outside,
				a.U_pos_bg, a.U_pos_ed,a.unitig_kmer_number,
				a.ori_read_id, a.read_pos_bg, a.read_pos_ed
		);
		if(NOT_full_conver)
			fprintf(stderr, "NOT_full_conver;");
		if(with_outside_gap)
			fprintf(stderr, "with_outside_gap;");
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "contig_total_len: %d\n", contig_total_len);
}

void TARGET_path_and_index::build(std::vector<READ_PATH_COMBINE_ITEM> & target_seq, int read_id, std::vector<UNITIG_INFO_item> & unitig_l,
		int target_id, int clip_target_mode){
	this->target_id = target_id;
	this->clip_target_mode = clip_target_mode;
	support_read_l.clear();
	support_read_l.emplace_back(read_id);
	int unitig_number = unitig_l.size();
	if(unitig_number == 0)
		return;
	target_index.resize(unitig_number);
	for(int i = 0; i < unitig_number; i++)
		target_index[i].clear();
	uint target_node_id = 0;
	target_total_len = 0;
	for(target_node_id = 0; target_node_id < target_seq.size(); target_node_id++){
		READ_PATH_COMBINE_ITEM * tq_node = &(target_seq[target_node_id]);
		target_link_path.emplace_back();
		target_link_path.back().set(read_id, tq_node, target_node_id, target_total_len);
		target_link_path.back().query_node_info.clear();
		target_link_path.back().query_node_info.emplace_back(&(target_seq[target_node_id]));

		target_total_len += target_link_path.back().kmer_number +  target_link_path.back().gap_kmer_num_outside +  target_link_path.back().gap_kmer_num_inside;
		if(target_node_id != 0){
			target_link_path[target_node_id].prev_node_id = target_node_id - 1;
			target_link_path[target_node_id - 1].next_node_id = target_node_id;
		}
		//skip too short unitig node:
		if(tq_node->UID >= 0 && unitig_l[tq_node->UID].kmer_list.size() >= 3 && (tq_node->gap_kmer_num_outside) * 4 < tq_node->match_kmer_num){
//				xassert(tq_node.UID < target_index.size(), "");
			target_index[tq_node->UID].emplace_back(target_node_id);
		}
	}

	the_bg_node_id = -1;
	if(!target_link_path.empty())
		the_bg_node_id = 0;
}

void TARGET_path_and_index::search_match(std::vector<READ_PATH_COMBINE_ITEM> &query, std::vector<SDP_MATCH_ITEM> &SDP_match_list){
	//function: generate the match list
	//S1: find all match of read and the target
	SDP_match_list.clear();
	//for all read node:
	for(uint read_node_idx = 0; read_node_idx < query.size(); read_node_idx++){
		READ_PATH_COMBINE_ITEM & rn = query[read_node_idx];
		//search the uid of read node in the target index:
		if(rn.UID != -1 && !target_index[rn.UID].empty()){
			//SKIP too repeat UNITIG
			if(target_index[rn.UID].size() > 30){
				//fprintf(stderr, "target_index[rn.UID].size(), %ld\n", target_index[rn.UID].size());

			}else{
				for(int target_node_id : target_index[rn.UID]){
					if(target_node_id != -1){
						SDP_match_list.emplace_back();
						SDP_match_list.back().set(target_node_id, &(query[read_node_idx]));
					}
				}
			}
		}
	}
}

void TARGET_path_and_index::SDP(bool print_log, std::vector<READ_PATH_COMBINE_ITEM> &query, bool LRS_full_begin, bool LRS_full_end, int read_length){
	print_log = false;
	//SDP match generating:
	search_match(query, match_list);
	if(print_log){
		fprintf(stderr, "SDP match result: before SDP \n");
		for(SDP_MATCH_ITEM & m : match_list){
			m.show(target_link_path);
		}
	}

	//sort: THE match is sorted by its position in the read, so it is no need to re-sort
	max_score_final = -100000000;
	max_score_idx = -1;
	uint match_number = match_list.size();
	for(uint handle_match_idx = 0; handle_match_idx < match_number; handle_match_idx++){
		SDP_MATCH_ITEM & cur_handle_match = match_list[handle_match_idx];
		int max_score = cur_handle_match.read_node->match_kmer_num;
		int diff_max_score = MAX_int32t;
		//set the init score, NW-like score setting method.
		if(LRS_full_begin){
			int begin_puanty = target_link_path[cur_handle_match.target_node_id].node_base_position_in_path;
			begin_puanty = MAX(begin_puanty, cur_handle_match.read_node->read_pos_bg);
			max_score -= begin_puanty*begin_puanty;
		}
		int max_score_suffix_idx = -1;
		//search reverse to get the max score:
		int search_end = (int)handle_match_idx - 200;
		search_end = MAX(0, search_end);
		for(int suffix_match_idx = handle_match_idx - 1; suffix_match_idx >= search_end; suffix_match_idx--){
			SDP_MATCH_ITEM & pre_try_match = match_list[suffix_match_idx];
			//break conditions
			//condition1:
			//the position in the target and query must be increased
			//target
			xassert(cur_handle_match.target_node_id >= 0 && cur_handle_match.target_node_id < (int)target_link_path.size(), "");
			xassert(pre_try_match.target_node_id >= 0 && pre_try_match.target_node_id < (int)target_link_path.size(), "");

			//target: the position must be increasing
			if(target_link_path[cur_handle_match.target_node_id].node_idx_position_in_path <=
					target_link_path[pre_try_match.target_node_id].node_idx_position_in_path)
				continue;
			//query: the position must be increasing
			if(cur_handle_match.read_node->read_pos_bg <= pre_try_match.read_node->read_pos_bg)
				continue;

			//condition2:
			//the gap between two Match must less the one value
			int tp1 = target_link_path[pre_try_match.target_node_id].node_base_position_in_path;
			int tp2 = target_link_path[cur_handle_match.target_node_id].node_base_position_in_path;
			int qp1 = pre_try_match.read_node->read_pos_bg;
			int qp2 = cur_handle_match.read_node->read_pos_bg;
			int diff = (tp1 - qp1) - (tp2 - qp2);
			//int ABS_diff = ABS(diff);
			//if(ABS_diff > 100)continue;

			//condition3: using the max score
			int score = pre_try_match.max_score + cur_handle_match.read_node->match_kmer_num - diff*diff;
			if(score > max_score){
				diff_max_score = diff;
				max_score = score;
				max_score_suffix_idx = suffix_match_idx;
			}
		}
		//store the SDP results
		cur_handle_match.max_score = max_score;
		cur_handle_match.diff = diff_max_score;
		cur_handle_match.max_previous_node = max_score_suffix_idx;
	}
	//get the max_score node:
	for(uint handle_match_idx = 0; handle_match_idx < match_number; handle_match_idx++){
		if(LRS_full_end){
			ASS_TARGET_ITEM &a = target_link_path[match_list[handle_match_idx].target_node_id];
			int end_punaty_read = read_length - match_list[handle_match_idx].read_node->read_pos_ed;//todo::
			int end_punaty_target = target_total_len - (a.node_base_position_in_path + a.gap_kmer_num_inside + a.kmer_number);
			int end_punaty = MAX(end_punaty_read, end_punaty_target);
			//fprintf(stderr, "LRS_full_end: target_total_len %d end_punaty_read %d\n, end_punaty_target %d\n", target_total_len, end_punaty_read, end_punaty_target);
			match_list[handle_match_idx].max_score -= end_punaty*end_punaty;
		}
		int score = match_list[handle_match_idx].max_score;
		if(score > max_score_final){
			max_score_final = score;
			max_score_idx = handle_match_idx;
		}
	}
	SDP_link_path.clear();
	int node_idx = max_score_idx;
	for(;node_idx != -1;){
		SDP_link_path.emplace_back(node_idx);
		xassert(node_idx >= 0 && node_idx < (int)match_list.size(), "");
		node_idx = match_list[node_idx].max_previous_node;
	}
	std::reverse(SDP_link_path.begin(), SDP_link_path.end());

	if(print_log){
		fprintf(stderr, "SDP match result: after SDP \n");
		for(SDP_MATCH_ITEM & m : match_list){
			m.show(target_link_path);
		}
	}
}

void TARGET_path_and_index::showSDP_result(READ_PATH_COMBINE_ITEM * q_p, int read_id, int cur_target_ID){
	if(max_score_idx == -1)
		fprintf(stderr, "[SDP log:] [r_id:%d;t_id:%d]No SDP results", read_id, cur_target_ID);
	else
		fprintf(stderr, "[SDP log:] [r_id:%d;t_id:%d]The final max score is %d @ UID: %d \n", read_id, cur_target_ID,
				max_score_final, target_link_path[match_list[max_score_idx].target_node_id].UID);
	if(true){
		if(false) fprintf(stderr, "the SDP path is: \n");
		int total_diff = 0;
		int total_ABS_diff = 0;
		for(uint i = 0; i < SDP_link_path.size(); i++){
			int link_id = SDP_link_path[i];
			int the_uid = target_link_path[match_list[link_id].target_node_id].UID;
			int	idx_target = target_link_path[match_list[link_id].target_node_id].node_idx_position_in_path;
			//int	base_position_target = target_link_path[match_list[link_id].target_node_id].node_idx_position_in_path;
			int idx_query = match_list[link_id].read_node - q_p;
			if(true) fprintf(stderr, "%d:%d:%d:%d:%d-->", the_uid, idx_target, idx_query, match_list[link_id].diff, match_list[link_id].max_score);
			if(i != 0){
				total_ABS_diff += ABS(match_list[link_id].diff);
				total_diff += match_list[link_id].diff;
			}
		}
		fprintf(stderr, "total_ABS_diff %d total_diff %d \n", total_ABS_diff, total_diff);
	}
}


bool TARGET_path_and_index::true_branch_check_is_same_run_length(bool print_log, std::string &a, std::string &b ){
	char old_same_char = -1;
	uint a_i = 0, b_i = 0;
	for(;;){
		if(a_i < a.size() && a[a_i] == old_same_char) {
			a_i++; continue;
		}
		if(b_i < b.size() && b[b_i] == old_same_char) {
			b_i++; continue;
		}
		if(a_i < a.size() && b_i < b.size() ){
			if(a[a_i] != b[b_i]){
				return true;
			}
			else{
				old_same_char = a[a_i]; a_i++; b_i++;
			}
		}else
			break;
	}
	if(a_i == a.size() && b_i == b.size()){
		if(print_log)
			fprintf(stderr, "true_branch_check_is_same_run_length: false %s %s \n",a.c_str(), b.c_str());
		return false;
	}
	return true;

}

bool TARGET_path_and_index::node_is_branch(bool print_log,
		bool is_clip_read,
		SDP_MATCH_ITEM * r, int read_id,
		int target_bg_id, int target_ed_id, READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed, std::vector<UNITIG_INFO_item> & unitig_l,
		std::vector<LONG_BRANCH_CHECK> &long_branch_check_list, std::vector <SPECIAL_DETECT_BRANCH> &special_branch_check_list,
		int & new_long_branch_check_node){
	//if(target_ed_id == -1)	return false;
	int diff = ABS_U(r->target_len_total, r->query_len_total);
	int max_len = MAX(r->target_len_total, r->query_len_total);
	int diff1 = r->read_node->read_pos_bg - target_link_path[r->target_node_id].node_base_position_in_path;
	if(print_log){
		fprintf(stderr, "++diff %d, diff1 %d max_len %d r->target_len_total %d , r->query_len_total %d r->target_node_gap %d,	r->query_node_gap %d\n",
				diff, diff1, max_len, r->target_len_total, r->query_len_total, r->target_node_gap,	r->query_node_gap);
	}
	if((diff *5 >= max_len && diff > 20)){
		if(print_log)fprintf(stderr, "[True branch: length diff HIGH, target_len_total %d, query_len_total %d] \n", r->target_len_total, r->query_len_total);
		return true;
	}
	else if(r->target_node_gap == false && r->query_node_gap == false){
		bool is_branch_node_TMP = false;
		//basic checking
		READ_PATH_COMBINE_ITEM * query_node = query_bg;
		for(int target_node_id = target_bg_id ;target_node_id != target_ed_id && query_node  != query_ed ;
				target_node_id = target_link_path[target_node_id].next_node_id, query_node++){
			if( target_link_path[target_node_id].UID != query_node->UID &&
					(target_link_path[target_node_id].kmer_number >= 3 || query_node->match_kmer_num >= 3)){
				is_branch_node_TMP = true;
			}
		}

		//home-polymer-checking
		if(is_branch_node_TMP){
			std::string s_target;
			std::string s_query;
			int kmer_size = unitig_l[0].unitig.size()  - unitig_l[0].kmer_list.size() + 1;
			for(int target_node_id = target_bg_id; target_node_id != target_ed_id; target_node_id = target_link_path[target_node_id].next_node_id){
				std::string & u_str = unitig_l[target_link_path[target_node_id].UID].unitig;
				if(target_node_id == target_bg_id)	s_target = u_str;
				else								s_target += u_str.substr(kmer_size - 1, u_str.size() - kmer_size + 1);
			}
			for(READ_PATH_COMBINE_ITEM * query_node = query_bg;  query_node != query_ed; query_node++){
				std::string & u_str = unitig_l[query_node->UID].unitig;
				if(query_node == query_bg)	s_query = u_str;
				else						s_query += u_str.substr(kmer_size - 1, u_str.size() - kmer_size + 1);
			}

			if(false == true_branch_check_is_same_run_length(print_log, s_target, s_query)){
				is_branch_node_TMP = false;
			}
		}

		if(is_branch_node_TMP){
			//count node number
			int query_total_node = 0;
			int targt_total_node = 0;
			int query_short_node = 0;
			int targt_short_node = 0;

			if(print_log) fprintf(stderr, "New branch: BG: @ read pos: %d, UID %d query_total_node %d , targt_total_node %d \n", query_bg->read_pos_bg, query_bg->UID, query_total_node, targt_total_node);
			if(print_log){
				for(READ_PATH_COMBINE_ITEM * query_node = query_bg; query_node != query_ed;query_node++)
					fprintf(stderr, "QUE: UID %d %d AVE:%f\n", query_node->UID, query_node->match_kmer_num, unitig_l[query_node->UID].averageKmerDepth_read);
				for(int target_node_id = target_bg_id; target_node_id != target_ed_id ; target_node_id = target_link_path[target_node_id].next_node_id)
					fprintf(stderr, "TAR:UID %d %d AVE:%f\n",
							target_link_path[target_node_id].UID, target_link_path[target_node_id].kmer_number,
							unitig_l[target_link_path[target_node_id].UID].averageKmerDepth_read);
				fprintf(stderr, "\n");
			}

			for(READ_PATH_COMBINE_ITEM * q = query_bg;q != query_ed;q++, query_total_node++)						{ if(q->match_kmer_num < 5) query_short_node ++; }
			for(int t = target_bg_id; t != target_ed_id ; t = target_link_path[t].next_node_id, targt_total_node++)	{ if(target_link_path[t].kmer_number < 5) targt_short_node ++; }

			//check if is is_sp_branch:
			if(!special_branch_check_list.empty()){
				for(SPECIAL_DETECT_BRANCH & s :special_branch_check_list){
					if(s.same_as(query_total_node, targt_total_node,query_bg->UID)){
						if(print_log)fprintf(stderr, "[True branch: sp_branch] \n");
						return true;
					}
				}
			}

			if(query_total_node + targt_total_node > 7){
				if(print_log)fprintf(stderr, "[False branch: too many node:ADD SP: [%d %d] ] \t",  query_total_node, targt_total_node);
				long_branch_check_list.emplace_back();
				new_long_branch_check_node ++;
				long_branch_check_list.back().set(read_id, target_id, query_total_node, targt_total_node, false, query_bg->UID);
			}
			else if(query_short_node != 0 || targt_short_node != 0){
				if(print_log)fprintf(stderr, "[False branch: short node:] \t");
			}
			else{
				if(print_log)fprintf(stderr, "[True branch: new] \t");
				return true;
			}
		}
	}
	return false;
}

void TARGET_path_and_index::target_length_count(SDP_MATCH_ITEM * r, int target_bg_id, int target_ed_id){
	r->target_node_gap = false;
	r->target_len_total = 0;
	for(int target_node_id = target_bg_id; target_node_id != target_ed_id; target_node_id = target_link_path[target_node_id].next_node_id){
		xassert(target_node_id < (int)target_link_path.size(),"");
		ASS_TARGET_ITEM & a = target_link_path[target_node_id];
		r->target_len_total += a.kmer_number + a.gap_kmer_num_inside + a.gap_kmer_num_outside;
		if(a.target_node_with_gap())
			r->target_node_gap = true;
	}
}

void TARGET_path_and_index::query_length_count(SDP_MATCH_ITEM * r, READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed){
	r->query_node_gap = false;
	r->query_len_total = 0;
	for(READ_PATH_COMBINE_ITEM * query_node = query_bg;  query_node != query_ed; query_node++){
		r->query_len_total += query_node->match_kmer_num + query_node->gap_kmer_num_inside + query_node->gap_kmer_num_outside;
		if(query_node->query_node_with_gap())
			r->query_node_gap = true;
	}
}

void TARGET_path_and_index::update_one_region(bool print_log,  int read_id,
		int target_bg_id, int target_ed_id, READ_PATH_COMBINE_ITEM * query_bg, READ_PATH_COMBINE_ITEM * query_ed, std::vector<UNITIG_INFO_item> & unitig_l){
	if(print_log){
		//update length check:
		int target_len_update = 0;
		int query_len_update = 0;
		{
			for(READ_PATH_COMBINE_ITEM * query_node = query_bg; query_node != query_ed;query_node++)
				query_len_update += query_node->match_kmer_num + query_node->gap_kmer_num_inside + query_node->gap_kmer_num_outside;
			for(int target_node_id = target_bg_id; target_node_id != target_ed_id ; target_node_id = target_link_path[target_node_id].next_node_id){
				ASS_TARGET_ITEM & a = target_link_path[target_node_id];
				target_len_update += a.kmer_number + a.gap_kmer_num_inside + a.gap_kmer_num_outside;
			}
		}

		fprintf(stderr, "NODE-UPDATA: %d %d target_len_update %d query_len_update %d \n", query_bg->read_pos_bg, query_bg->UID, target_len_update, query_len_update);
		for(READ_PATH_COMBINE_ITEM * query_node = query_bg; query_node != query_ed;query_node++){
			fprintf(stderr, "QUE: UID %d (%d:%d:%d) AVE:%f\n", query_node->UID, query_node->match_kmer_num, query_node->gap_kmer_num_inside,
					query_node->gap_kmer_num_outside, unitig_l[query_node->UID].averageKmerDepth_read);
		}
		for(int target_node_id = target_bg_id; target_node_id != target_ed_id ; target_node_id = target_link_path[target_node_id].next_node_id){
			fprintf(stderr, "TAR:UID %d (%d:%d:%d) AVE:%f\n",
				target_link_path[target_node_id].UID, target_link_path[target_node_id].kmer_number,
				target_link_path[target_node_id].gap_kmer_num_inside,
				target_link_path[target_node_id].gap_kmer_num_outside,
				unitig_l[target_link_path[target_node_id].UID].averageKmerDepth_read);
		}
		fprintf(stderr, "\n");
	}

	//there is no need to remove the NO-use node, just skip it in the link table
	//remove the nodes in the index
	for(int target_node_id = target_bg_id; target_node_id != target_ed_id; target_node_id = target_link_path[target_node_id].next_node_id){

		if(print_log){
			fprintf(stderr, "remove node: %d\n", target_node_id);
			xassert(target_node_id >= 0 && target_node_id < (int)target_link_path.size(), "");
		}
		//search the nodes in the index and set id to -1 to show the remove
		for(uint i = 0; i < target_index[target_link_path[target_node_id].UID].size(); i++){
			xassert(target_link_path[target_node_id].UID < (int)target_index.size(), "");
			xassert(target_node_id < (int)target_link_path.size(),"");
			if(target_index[target_link_path[target_node_id].UID][i] == target_node_id){
				xassert(i < target_index[target_link_path[target_node_id].UID].size(), "");
				target_index[target_link_path[target_node_id].UID][i] = -1;
				break;
			}
		}
	}

	{
		//add new nodes
		int new_node_num = 0;
		for(READ_PATH_COMBINE_ITEM * query_node = query_bg;  query_node != query_ed; query_node++){
			//new a node and store to the end of the path link table
			target_link_path.emplace_back();
			target_link_path.back().set(read_id, query_node, -1, -1);
			new_node_num ++;
			//add link:
			if(new_node_num != 1){
				xassert(target_link_path.size() - 2 >= 0,"");
				target_link_path.back().prev_node_id = target_link_path.size() - 2;
				target_link_path[target_link_path.size() - 2].next_node_id = target_link_path.size() - 1;
			}
			//store node to index:
			if(unitig_l[query_node->UID].kmer_list.size() >= 3 && query_node->gap_kmer_num_outside == 0)
				target_index[query_node->UID].emplace_back(target_link_path.size() - 1);
		}
		//add the link to the bg node
		int first_new_node_id = target_link_path.size() - new_node_num;
		int last__new_node_id = target_link_path.size() - 1;
		//record the link before the begin and after the end

		int node_before_the_first_id = target_link_path[target_bg_id].prev_node_id;
		int node_after__the_last__id = target_ed_id;

		if(print_log){
			fprintf(stderr, "build link: %d<--%d[target_bg_id %d, target_ed_id %d]\n", first_new_node_id, node_before_the_first_id, target_bg_id, target_ed_id);
			fprintf(stderr, "build link: %d-->%d\n", node_before_the_first_id, first_new_node_id);
			fprintf(stderr, "build link: %d-->%d\n", last__new_node_id, node_after__the_last__id);
			if(node_after__the_last__id >= 0)
				fprintf(stderr, "build link: %d<--%d\n", node_after__the_last__id, last__new_node_id);
		}

		target_link_path[first_new_node_id].prev_node_id = node_before_the_first_id;
		if(node_before_the_first_id >= 0)//only set when data exist
			target_link_path[node_before_the_first_id].next_node_id = first_new_node_id;
		//add the link to the ed node
		target_link_path[last__new_node_id].next_node_id = node_after__the_last__id;
		if(node_after__the_last__id >= 0)//only set when data exist
			target_link_path[node_after__the_last__id].prev_node_id = last__new_node_id;

		if(print_log){
			xassert(target_bg_id < (int)target_link_path.size(),"");
			xassert(first_new_node_id < (int)target_link_path.size(),"");
			xassert(node_before_the_first_id < (int)target_link_path.size(),"");
			xassert(last__new_node_id < (int)target_link_path.size(),"");
			xassert(node_after__the_last__id < (int)target_link_path.size(),"");
			xassert(target_bg_id >= 0,"");
			xassert(last__new_node_id >= 0,"");
		}
	}
}

bool TARGET_path_and_index::SDP_update_and_branch_condition_check(
		bool print_log, int MIN_score_final,
		std::vector<READ_PATH_COMBINE_ITEM> &query, int read_id, bool LRS_full_begin, bool LRS_full_end, bool is_clip_read, //read info
		std::vector<UNITIG_INFO_item> & unitig_l,//background
		std::vector<LONG_BRANCH_CHECK> &long_branch_check_list, std::vector <SPECIAL_DETECT_BRANCH> &special_branch_check_list){//
	if(max_score_final < MIN_score_final){
		bool is_absorb_by_target = false;
		return is_absorb_by_target;
	}
	bool is_absorb_by_target = true;
	//fprintf(stderr, "The final max score is %d @ %d, the path (reverse) is: \n", max_score_final, max_score_idx);
	bool with_update_node = false;
	int new_long_branch_check_node = 0;
	//check all node to show if with other branch?
	for(uint i = 0; i < SDP_link_path.size(); i++){
		xassert(SDP_link_path[i] < (int)match_list.size() && SDP_link_path[i]  >= 0, "");
		if(i != SDP_link_path.size() - 1){
			xassert(SDP_link_path[i + 1] < (int)match_list.size() && SDP_link_path[i + 1]  >= 0, "");
		}
		//the begin
		SDP_MATCH_ITEM * r = &(match_list[SDP_link_path[i]]);
		int target_bg_id = r->target_node_id;
		READ_PATH_COMBINE_ITEM * query_bg = r->read_node;
		//the end
		SDP_MATCH_ITEM * r_next = (i == SDP_link_path.size() - 1)?NULL:(&(match_list[SDP_link_path[i + 1]]));
		int target_ed_id = (r_next == NULL)?(-1):r_next->target_node_id;
		READ_PATH_COMBINE_ITEM * query_ed = (r_next == NULL)?(&(query[0]) + query.size()):(r_next->read_node);
		//count the length and gap
		target_length_count(r, target_bg_id, target_ed_id);
		query_length_count(r, query_bg, query_ed);
		//set flag: is branch
		//check if the node is branch node:
		//condition1: clip read, and at the same time, when checking the first or last node, skip
		bool is_branch_node = false;
		if(is_clip_read && (i == 0 || i == SDP_link_path.size() - 1)){
			/*do nothing*/
		}else{
			is_branch_node = node_is_branch(print_log,is_clip_read, r, read_id, target_bg_id, target_ed_id, query_bg, query_ed,
					unitig_l, long_branch_check_list, special_branch_check_list, new_long_branch_check_node);
		}
		if(is_branch_node)
			is_absorb_by_target = false;
		//set flag: if need update
		if(r->target_node_gap == true && r->query_node_gap == false){
			with_update_node = true;
		}
	}

	if(false == is_absorb_by_target)
		return is_absorb_by_target;
	//when is_absorb_by_target is true, set the info to "long_branch_check_list"
	for(int i = 0; i < new_long_branch_check_node; i++)
		long_branch_check_list[long_branch_check_list.size() - i - 1].query_is_absorb = true;

	if(print_log){
		fprintf(stderr, "Read is absorb to the branch.\n");
		show_total_length(unitig_l);
	}

	support_read_l.emplace_back(read_id);

	//store the node info of query to each target node
	for(uint i = 0; i < SDP_link_path.size(); i++){
		xassert(SDP_link_path[i] < (int)match_list.size() && SDP_link_path[i]  >= 0, "");
		SDP_MATCH_ITEM * r = &(match_list[SDP_link_path[i]]);
		ASS_TARGET_ITEM & target_node = target_link_path[r->target_node_id];
		target_node.query_node_info.emplace_back(r->read_node);
	}

	if(with_update_node){
		//update all node:
		for(uint i = 0; i < SDP_link_path.size(); i++){
			if(false == LRS_full_end && i ==  SDP_link_path.size() - 1)//when the read is not reach the end of reference region, and the same time, it try to update the tail of the contig:
				continue;
			if(match_list[SDP_link_path[i]].node_node_update()){
				//the begin
				SDP_MATCH_ITEM * r = &(match_list[SDP_link_path[i]]);
				int target_bg_id = r->target_node_id;
				READ_PATH_COMBINE_ITEM * query_bg = r->read_node;
				//the end
				SDP_MATCH_ITEM * r_next = (i == SDP_link_path.size() - 1)?NULL:(&(match_list[SDP_link_path[i + 1]]));
				int target_ed_id = (r_next == NULL)?(-1):r_next->target_node_id;
				READ_PATH_COMBINE_ITEM * query_ed = (r_next == NULL)?(&(query[0]) + query.size()):(r_next->read_node);

				update_one_region(print_log, read_id, target_bg_id, target_ed_id, query_bg, query_ed, unitig_l);
			}
		}

		//update all node ID when the path is updated
		int target_bg_id = the_bg_node_id;
		int node_position_in_path = 0;
		target_total_len = 0;

		for(int target_node_id = target_bg_id; target_node_id != -1; target_node_id = target_link_path[target_node_id].next_node_id){
			ASS_TARGET_ITEM & a = target_link_path[target_node_id];
			a.node_idx_position_in_path = node_position_in_path;
			a.node_base_position_in_path = target_total_len;
			node_position_in_path ++;
			target_total_len += a.kmer_number + a.gap_kmer_num_outside + a.gap_kmer_num_inside;
		}
		//show new path:
		if(print_log){
			fprintf(stderr, "The target is updated by new read, new_len:%d. ", target_total_len);
			show_total_length(unitig_l);
			show_the_target_path();
		}
	}
	return true;
}

int run_length_count(std::string & a){
	const char *a_p = &(a[0]);
	int len = a.size();
	if(len == 0) return 0;
	int run_len = 1; uint8_t c = a_p[0];
	if(false) fprintf(stderr, "%c", c);
	for(int i = 1; i < len; i++){
		if(a_p[i] != c){
			run_len ++;
			if(false) fprintf(stderr, "%c(%d)", c, i);
			c = a_p[i];
		}
	}
	if(false) fprintf(stderr, "\n");
	return run_len;
}

void CONTIG_SELECTER::contig_basic_count(bool print_log, std::vector<TARGET_path_and_index> &target_path_idx_l, std::vector<READ_MATCH_NUM> & rml){
	int total_read_number = 0;
	for(uint cur_target_ID = 0; cur_target_ID < target_path_idx_l.size(); cur_target_ID++){
		for(SUPPORT_READ_LIST_ITEM & s : target_path_idx_l[cur_target_ID].support_read_l){
			if(rml[s.read_id].is_full_cover == true || rml[s.read_id].is_clip_read == true){
				if(s.is_best_target == true)
					total_read_number ++;
			}
		}
	}
	select_list.clear();
	for(uint cur_target_ID = 0; cur_target_ID < target_path_idx_l.size(); cur_target_ID++){
		select_list.emplace_back();
		int full_clip_supp_read_n = 0;
		int clip_read_n = 0;
		int clip_read_left_n = 0;
		int clip_read_right_n = 0;
		int full_cover_read_n = 0;

		for(SUPPORT_READ_LIST_ITEM & s : target_path_idx_l[cur_target_ID].support_read_l){
			if(s.is_best_target == true){
				if(rml[s.read_id].is_clip_right == true  && rml[s.read_id].is_clip_left == true)
					continue;//skip when both end clip
				if(rml[s.read_id].is_full_cover == true || rml[s.read_id].is_clip_read == true)
					full_clip_supp_read_n ++;
				if(rml[s.read_id].is_full_cover == true)
					full_cover_read_n ++;
				if(rml[s.read_id].is_clip_read){
					clip_read_n++;
					if(rml[s.read_id].is_clip_left) clip_read_left_n++;
					if(rml[s.read_id].is_clip_right) clip_read_right_n++;
				}
			}
		}
		select_list.back().set(
				full_clip_supp_read_n, cur_target_ID,
				clip_read_n, clip_read_left_n, clip_read_right_n, full_cover_read_n,
				total_read_number - full_clip_supp_read_n,
				run_length_count(target_path_idx_l[cur_target_ID].contig),
				target_path_idx_l[cur_target_ID].contig.size());
	}
}

bool CONTIG_SELECTER::select_contig(bool print_log, int preset_platform){
	if(select_list.empty())
		return false;
	//sort by contig length:
	std::sort(select_list.begin(), select_list.end(), CONTIG_SELECT::cmp_by_contig_len);
	//
	int cur_cmp_SV_len = select_list[0].contig_cmp_len;
	//int contig_len = select_list[0].contig_len;
	int max_read_size = select_list[0].support_read_number;
	int max_path_idx = 0;

	int MAX_DIFF = 0;
	if(preset_platform == PRESET_CCS)	MAX_DIFF = 5;
	else								MAX_DIFF = 10;

	for(uint i = 1; i < select_list.size(); i++){
		CONTIG_SELECT *ch = &(select_list[i]);
		if(ABS_U(cur_cmp_SV_len, ch->contig_cmp_len) < MAX_DIFF){//similar length
			if(max_read_size < ch->support_read_number){
				max_read_size = ch->support_read_number;
				max_path_idx = i;
			}
		}else{
			//store
			select_list[max_path_idx].is_final_path = true;
			//reset
			cur_cmp_SV_len = ch->contig_cmp_len;
			//contig_len = ch->contig_len;
			max_read_size = ch->support_read_number;
			max_path_idx = i;
		}
	}
	select_list[max_path_idx].is_final_path = true;

	int final_path_number = 0;
	//count the path number the first time:
	for(uint i = 0; i < select_list.size(); i++)
		if(select_list[i].is_final_path)
			final_path_number ++;
	//select the first 2 MAX path when path number is over 2
	if(final_path_number > 2){//re-select when with too
		std::sort(select_list.begin(), select_list.end(), CONTIG_SELECT::cmp_by_supp_read);
		select_list[0].is_final_path = true; select_list[1].is_final_path = true;
		for(uint i = 2; i < select_list.size(); i++) select_list[i].is_final_path = false;
	}
	if(print_log){
		fprintf(stderr, "Before contig select: \n\n");
		for(auto & s :select_list)	s.show(true);
	}

	std::sort(select_list.begin(), select_list.end(), CONTIG_SELECT::cmp_by_final);
	if(print_log){
		fprintf(stderr, "Before contig select: \n\n");
		for(auto & s :select_list)	s.show(true);
	}

	//count total read number
	final_path_number = 0;
	int total_read_n = 0;
	int unknown_read_n = 0;
	for(uint i = 0; i < select_list.size(); i++){
		total_read_n += select_list[i].support_read_number;
		if(!select_list[i].is_final_path)
			unknown_read_n += select_list[i].support_read_number;
		else{
			final_path_number++;
		}
	}
	use_contig_num = final_path_number;
	return true;
}

void Clip_contig_combine_Handler::init(){
	REF_IDX_KMER_LEN = 20;
	extern uint64_t kmerMask[33];
	MASK = kmerMask[REF_IDX_KMER_LEN];
}

void Clip_contig_combine_Handler::clear(){
	left_contig_n = 0;
	right_contig_n = 0;
	total_contig_n = 0;

	CLIP_AT_right_contigs.clear();
	CLIP_AT_left_contigs.clear();
	full_contigs.clear();

	CLIP_AT_right_contigs_target_ID.clear();
	CLIP_AT_left_contigs_target_ID.clear();
	full_contigs_target_ID.clear();

	clip_contig_with_both_end = false;
	clip_contig_over2 = false;

	contig_match_l.clear();
	match_MEM.clear();
	SDP_link_path.clear();
	SDP_analysis_l.clear();
	kmer_index.clear();;
	bin_contig.clear();;
	CLIP_AT_right_contigs_n = 0;

	combine_contigs.clear();
}


void Clip_contig_combine_Handler::build_idx_for_right_contig(std::vector<std::string > &right_contig){
	kmer_index.clear();
	CLIP_AT_right_contigs_n = right_contig.size();
	for(uint contig_ID = 0; contig_ID < right_contig.size(); contig_ID++){
		store_bin_contig(right_contig[contig_ID], bin_contig);			//char to bin
		uint8_t * buff_bin = &(bin_contig[0]);
		std::unordered_map<int32_t, std::vector<Index_Item>>::iterator it;
		int kmer_number = bin_contig.size() - REF_IDX_KMER_LEN + 1;
		uint64_t kmer = bit2_nextKmer_init(buff_bin, REF_IDX_KMER_LEN);
		for(int contig_position = 0; contig_position < kmer_number; contig_position++){
			kmer = bit2_nextKmerMASK(buff_bin + contig_position, kmer, REF_IDX_KMER_LEN);
			kmer_index[kmer].emplace_back(contig_position, contig_ID);
		}
	}
}

void Clip_contig_combine_Handler::get_contig_kmer_match(std::string &right_contig){
	store_bin_contig(right_contig, bin_contig);			//char to bin
	uint8_t * buff_bin = &(bin_contig[0]);
	int kmer_number = bin_contig.size() - REF_IDX_KMER_LEN + 1;
	contig_match_l.clear();
	contig_match_l.resize(CLIP_AT_right_contigs_n);

	uint64_t kmer = bit2_nextKmer_init(buff_bin, REF_IDX_KMER_LEN);
	for(int c2_pos = 0; c2_pos < kmer_number; c2_pos++){
		kmer = bit2_nextKmerMASK( buff_bin + c2_pos, kmer, REF_IDX_KMER_LEN);
		std::unordered_map<int32_t, std::vector<Index_Item>>::iterator it;
		it = kmer_index.find(kmer);
		if(it!=kmer_index.end()){
			for(Index_Item & ii: it->second){
				//store the match position
				contig_match_l[ii.contig_id].emplace_back();
				contig_match_l[ii.contig_id].back().set(ii.position, c2_pos);
			}
		}
	}
	for(int ci = 0; ci < CLIP_AT_right_contigs_n; ci++){
		//sort all the results:
		std::sort(contig_match_l[ci].begin(), contig_match_l[ci].end(),
				CONTIG_KMER_MATCH_ITEM::cmp_by_position_diff);
		//show all result:
		if(show_log && false){
			for(CONTIG_KMER_MATCH_ITEM & m: contig_match_l[ci]){
				fprintf(stderr, "m.position_diff %d, m.contig_position %d, m.ref_position %d\n", m.position_diff, m.c2_pos, m.c1_pos);
			}
		}
	}
}

void Clip_contig_combine_Handler::load_MEMs(std::vector<CONTIG_KMER_MATCH_ITEM> & cur_c_r_match){
	match_MEM.clear();
	if(cur_c_r_match.empty())
		return;
	int match_num = cur_c_r_match.size();
	int old_diff = -10000000;
	int old_c1 = -10000000;
	int match_count = 0;
	int c1_pos_bg = 0;
	for(int c2_idx = 0; c2_idx < match_num + 1; c2_idx++){
		if(c2_idx == match_num || old_diff != cur_c_r_match[c2_idx].position_diff || old_c1 + 1 != cur_c_r_match[c2_idx].c1_pos){
			//store the result:
			if(match_count > 20){
				match_MEM.emplace_back();
				match_MEM.back().store_final(old_diff, match_count, c1_pos_bg);
				//match_MEM.back().show(0);
			}
			//re-init
			old_diff = cur_c_r_match[c2_idx].position_diff;
			old_c1 = cur_c_r_match[c2_idx].c1_pos;
			c1_pos_bg = old_c1;
			match_count = 1;
		}
		else{
			match_count ++;
			old_c1 = cur_c_r_match[c2_idx].c1_pos;
		}
	}
	if(match_MEM.size() > 30000){
		std::vector<KMER_MATCH_MEM> match_MEM_TMP;
		std::swap(match_MEM, match_MEM_TMP);
		for(auto & m :match_MEM_TMP){
			if(m.count > 100){
				match_MEM.emplace_back();
				std::swap(match_MEM.back(), m);
			}
		}
	}
}
int Clip_contig_combine_Handler::sdp_score(KMER_MATCH_MEM & pre_try_match, KMER_MATCH_MEM & cur_handle_match){
	//condition2:
	//the gap between two Match must less the one value
	int c1_1 = pre_try_match.c_R_pos_bg;
	int c1_2 = cur_handle_match.c_R_pos_bg;
	int c2_1 = pre_try_match.get_c_L_pos();
	int c2_2 = cur_handle_match.get_c_L_pos();
	int diff = (c1_1 - c2_1) - (c1_2 - c2_2);
	int c1_1_ed = c1_1 + pre_try_match.count;
	int dis_c1 = ABS_U(c1_1_ed ,c1_2);
	//int ABS_diff = ABS(diff);
	//if(ABS_diff > 100)continue;
	//condition3: using the max score
	return pre_try_match.max_score + cur_handle_match.count - diff*diff*0.3 - 0.5*dis_c1;
}

void Clip_contig_combine_Handler::SDP_run(std::vector<SDP_ANALYSIS_ITEM> &SDP_analysis_l,
		int CLIP_AT_right_contig_len, int CLIP_AT_left_contig_len,
		int clip_AT_right_idx, int clip_AT_left_idx){
	if(match_MEM.empty())
		return;
	//sort all MEMs
	std::sort(match_MEM.begin(), match_MEM.end(), KMER_MATCH_MEM::cmp_by_contig_position);
	//SDP:
	//sort: THE match is sorted by its position in the read, so it is no need to re-sort
	int max_score_final = -100000000;
	int max_score_idx = -1;
	uint MEM_number = match_MEM.size();
	for(uint handle_match_idx = 0; handle_match_idx < MEM_number; handle_match_idx++){
		KMER_MATCH_MEM & cur_handle_match = match_MEM[handle_match_idx];
		int max_score = cur_handle_match.count;
		//set the init score, NW-like score setting method.
		int max_score_suffix_idx = -1;
		//search reverse to get the max score:
		int search_end = (int)handle_match_idx - 30000;
		search_end = MAX(0, search_end);
		for(int suffix_match_idx = handle_match_idx - 1; suffix_match_idx >= search_end; suffix_match_idx--){
			KMER_MATCH_MEM & pre_try_match = match_MEM[suffix_match_idx];
			//break conditions
			//condition1:
			//the position in the target and query must be increased
			//target
			//xassert(cur_handle_match.min_contig_position >= 0 && cur_handle_match.min_contig_position < (int)target_link_path.size(), "");
			//xassert(pre_try_match.min_contig_position >= 0 && pre_try_match.min_contig_position < (int)target_link_path.size(), "");

			//target: the position must be increasing
			if(cur_handle_match.c_R_pos_bg < pre_try_match.c_R_pos_bg)
				continue;
			//query: the position must be increasing
			if(cur_handle_match.get_c_L_pos() < pre_try_match.get_c_L_pos())
				continue;
			//condition3: using the max score
			int score = sdp_score(pre_try_match, cur_handle_match);
			if(score > max_score){
				max_score = score;
				max_score_suffix_idx = suffix_match_idx;
			}
		}
		//store the SDP results
		cur_handle_match.max_score = max_score;
		cur_handle_match.max_previous_node = max_score_suffix_idx;
	}
	//get the max_score node:
	for(uint handle_match_idx = 0; handle_match_idx < MEM_number; handle_match_idx++){
		int score = match_MEM[handle_match_idx].max_score;
		if(score > max_score_final){
			max_score_final = score;
			max_score_idx = handle_match_idx;
		}
	}
	if(false){
		for(uint handle_match_idx = 0; handle_match_idx < MEM_number; handle_match_idx++)
			match_MEM[handle_match_idx].show(handle_match_idx);
	}
	SDP_link_path.clear();
	int node_idx = max_score_idx;
	for(;node_idx != -1;){
		SDP_link_path.emplace_back(node_idx);
		xassert(node_idx >= 0 && node_idx < (int)match_MEM.size(), "");
		node_idx = match_MEM[node_idx].max_previous_node;
	}
	std::reverse(SDP_link_path.begin(), SDP_link_path.end());
	SDP_analysis_l.emplace_back();
	SDP_analysis_l.back().ana(match_MEM, SDP_link_path, CLIP_AT_right_contig_len, CLIP_AT_left_contig_len, clip_AT_right_idx, clip_AT_left_idx);
	SDP_analysis_l.back().show();
	//xassert(SDP_analysis_l.back().bg_SDP_cL <= CLIP_AT_left_contigs[SDP_analysis_l.back().clip_AT_left_idx].size(), "");
	;
	if(false){
		for(uint i = 0; i < SDP_link_path.size(); i++)
			match_MEM[SDP_link_path[i]].show(SDP_link_path[i]);
	}
}


void Clip_contig_combine_Handler::classify_contig_ref_match(std::vector<std::string > &CLIP_AT_right_contigs,
		std::vector<std::string > &CLIP_AT_left_contigs, int clip_AT_left_idx){
	int CLIP_AT_left_contig_len = CLIP_AT_left_contigs[clip_AT_left_idx].size();
	SDP_analysis_l.clear();
	//for each right contig, RUN:
	for(int clip_AT_right_idx = 0; clip_AT_right_idx < CLIP_AT_right_contigs_n; clip_AT_right_idx++){
		int CLIP_AT_right_contig_len = CLIP_AT_right_contigs[clip_AT_right_idx].size();
		//fprintf(stderr, "\nCLIP@right contig is %d, size %d\n", clip_AT_right_idx, CLIP_AT_right_contig_len);
		std::vector<CONTIG_KMER_MATCH_ITEM> & cur_c_r_match =  contig_match_l[clip_AT_right_idx];
		load_MEMs(cur_c_r_match);
		SDP_run(SDP_analysis_l, CLIP_AT_right_contig_len, CLIP_AT_left_contig_len, clip_AT_right_idx, clip_AT_left_idx);
		//fprintf(stderr, "[END]The CLIP@right ID is %d\n\n", clip_AT_right_idx);
	}
	//select_the_best right contig:
	int max_score = -1; int max_idx = -1;
	for(uint i = 0;  i< SDP_analysis_l.size();i++){
		if(SDP_analysis_l[i].condition1_absorbed == true ||
			SDP_analysis_l[i].condition2_linked == true){
			if(max_score < SDP_analysis_l[i].max_score && SDP_analysis_l[i].max_score > 1000){
				max_score = SDP_analysis_l[i].max_score;
				max_idx = i;
			}

		}
	}
	//combination with the best contig:
	if(max_idx != -1){
		SDP_ANALYSIS_ITEM &s = SDP_analysis_l[max_idx];
		std::string & rs = CLIP_AT_right_contigs[s.clip_AT_right_idx];
		std::string & ls = CLIP_AT_left_contigs[s.clip_AT_left_idx];
		std::string combine_contig;
		combine_contig.insert(combine_contig.begin(), rs.begin(), rs.begin() + s.bg_SDP_cR);//copy the 1st part
		combine_contig.insert(combine_contig.end(),   ls.begin() + s.bg_SDP_cL, ls.end());//copy the 2nd part
		fprintf(stderr, "rs %s, ls %s, combine_contig %s\n", rs.c_str(), ls.c_str(), combine_contig.c_str());
		//store the results
		combine_contigs.emplace_back();
		std::swap(combine_contigs.back().combine_s, combine_contig);
		combine_contigs.back().target_ID_part1 = CLIP_AT_right_contigs_target_ID[s.clip_AT_right_idx];
		combine_contigs.back().target_ID_part2 = CLIP_AT_left_contigs_target_ID[s.clip_AT_left_idx];
	}
}

void Clip_contig_combine_Handler::combine(){
	if(CLIP_AT_right_contigs.empty()){
		std::swap(full_contigs, CLIP_AT_right_contigs);
		std::swap(full_contigs_target_ID, CLIP_AT_right_contigs_target_ID);
	}else if(CLIP_AT_left_contigs.empty()){
		std::swap(full_contigs, CLIP_AT_left_contigs);
		std::swap(full_contigs_target_ID, CLIP_AT_left_contigs_target_ID);
	}
	build_idx_for_right_contig(CLIP_AT_right_contigs);
	//search the right-clip contig in the index
	int clip_AT_left_idx = -1;
	for(std::string &clc : CLIP_AT_left_contigs){
		clip_AT_left_idx ++;
		get_contig_kmer_match(clc); //search sort and show:
		classify_contig_ref_match(CLIP_AT_right_contigs, CLIP_AT_left_contigs, clip_AT_left_idx); //cluster the matches
	}
	//show final combine results:
	fprintf(stderr, "Show the final results\n");
	for(auto & cc: combine_contigs){
		fprintf(stderr, "[%d %d] %s\n",	cc.target_ID_part1,	cc.target_ID_part2,	cc.combine_s.c_str());
	}
}


void Contig_align_region_Handler::get_contig_ref_match(uint8_t * buff_bin, int kmer_number){
	CONTIG_REF_MATCH_list.clear();
	uint64_t kmer = bit2_nextKmer_init(buff_bin, REF_IDX_KMER_LEN);
	for(int i = 0; i < kmer_number; i++){
		kmer = bit2_nextKmerMASK( buff_bin + i, kmer, REF_IDX_KMER_LEN);
		int query_position = ref_handler->search_ref_idx(kmer);
		if(query_position != -1){
			//store the match position
			CONTIG_REF_MATCH_list.emplace_back();
			CONTIG_REF_MATCH_list.back().set(i,query_position);
		}
	}
	//sort all the results:
	std::sort(CONTIG_REF_MATCH_list.begin(), CONTIG_REF_MATCH_list.end(),
			LONG_SV_CONTIG_REF_MATCH::cmp_by_position_diff);
	//show all result:
	if(show_log && false){
		for(LONG_SV_CONTIG_REF_MATCH & m: CONTIG_REF_MATCH_list){
			fprintf(stderr, "m.position_diff %d, m.contig_position %d, m.ref_position %d\n", m.position_diff, m.contig_position, m.ref_position);
		}
	}
}

bool Contig_align_region_Handler::generate_sv_region(int contig_length){
	//init the results
	r_p1.set_empty();
	r_p2.set_empty();
	using_local_reference = false;
	sv_handle_is_needed = false;
	//additional match cluster:
	//for(LONG_SV_MATCH_CLUSTER & c : match_clusters)			c.show();
	bool with_clip_info = false;
	int clip_diff_position = 0;
	if(is_left_clip_contig && contig_length > 500)
	{
		with_clip_info = true;
		int ref_end = local_ref_region->ed_pos;
		int contig_end = contig_length;
		clip_diff_position = ref_end - contig_end;
		match_clusters.emplace_back();
		match_clusters.back().set(500, clip_diff_position, ref_end - 500, ref_end, contig_end - 500, contig_end);
		std::sort(match_clusters.begin(), match_clusters.end(), LONG_SV_MATCH_CLUSTER::cmp_by_ref_position);
		//remove all signal after it
		for(uint i = 0; i < match_clusters.size(); i++){
			if(match_clusters[i].position_diff == clip_diff_position){
				if(i <  match_clusters.size() - 1)
					match_clusters.erase(match_clusters.begin() + i + 1, match_clusters.end());
				break;
			}
		}
		std::sort(match_clusters.begin(), match_clusters.end(), LONG_SV_MATCH_CLUSTER::cmp_by_diff_position);
	}
	if(is_right_clip_contig && contig_length > 500){
		with_clip_info = true;
		int ref_bgn = local_ref_region->st_pos;
		int contig_bgn = 0;
		clip_diff_position = ref_bgn - contig_bgn;
		match_clusters.emplace_back();
		match_clusters.back().set(500, clip_diff_position, ref_bgn, ref_bgn + 500, contig_bgn, contig_bgn + 500);
		std::sort(match_clusters.begin(), match_clusters.end(), LONG_SV_MATCH_CLUSTER::cmp_by_ref_position);
		//remove all signal before it
		for(uint i = 0; i < match_clusters.size(); i++){
			if(match_clusters[i].position_diff == clip_diff_position){
				if(i > 0)
					match_clusters.erase(match_clusters.begin() , match_clusters.begin() + i);
				break;
			}
		}
		std::sort(match_clusters.begin(), match_clusters.end(), LONG_SV_MATCH_CLUSTER::cmp_by_diff_position);
	}
	for(LONG_SV_MATCH_CLUSTER & c : match_clusters)
		c.show();

	//analysis
	if(match_clusters.size() < 2)
		return false;
	int suggest_sv_len = 0;
	bool with_long_SV = false;
	for(uint i = 1; i < match_clusters.size(); i++){
		int ave_ref_pos_1 = match_clusters[i].get_ave_ref_pos();
		int ave_ref_pos_2 = match_clusters[i-1].get_ave_ref_pos();
		int ave_contig_pos_1 = match_clusters[i].get_ave_con_pos();
		int ave_contig_pos_2 = match_clusters[i-1].get_ave_con_pos();
		suggest_sv_len = (ave_contig_pos_2 - ave_contig_pos_1) - (ave_ref_pos_2 - ave_ref_pos_1);
		if(ave_contig_pos_1 > ave_contig_pos_2)
			suggest_sv_len = - suggest_sv_len;
		int abs_suggest_sv_len = ABS(suggest_sv_len);
		if(abs_suggest_sv_len > min_suggest_sv_len && suggest_sv_len < max_suggest_sv_len){
			//check one region is within the current region?
			//check
			RefRegion r1_try(local_ref_region->chr_ID, match_clusters[i-1].min_ref_position, match_clusters[i-1].max_ref_position);
			RefRegion r2_try(local_ref_region->chr_ID, match_clusters[i].min_ref_position, match_clusters[i].max_ref_position);
			bool is_r1_local = false;
			bool is_r2_local = false;

			int r1_overlap_len = r1_try.region_overlap_length(*local_ref_region);
			int r2_overlap_len = r2_try.region_overlap_length(*local_ref_region);

            if(r1_overlap_len > 400) is_r1_local = true;
            if(r2_overlap_len > 400) is_r2_local = true;

			//if(r1_overlap_len > 1000 || r1_try.Within(*local_ref_region)) is_r1_local = true;
			//if(r2_overlap_len > 1000 || r2_try.Within(*local_ref_region)) is_r2_local = true;

			fprintf(stderr, "r1_overlap %d, r2_overlap %d\n", r1_overlap_len, r2_overlap_len);
			//when no region is within the current region, skip SV calling
			if(false == is_r1_local && false == is_r2_local)
				continue;
//				//contig overlap check
//				RefRegion c1_try(0, match_clusters[i-1].min_contig_position, match_clusters[i-1].max_contig_position);
//				RefRegion c2_try(0, match_clusters[i].min_contig_position, match_clusters[i].max_contig_position);
//				int c_overlap_len = c1_try.region_overlap_length(c2_try);
//				//skip when contig region overlap:
//				if(c_overlap_len > 500)
//					continue;

			sv_handle_is_needed = true;
			with_long_SV = true;
			//C0: when both is within the local region, treat as local SV:
			if(is_r1_local && is_r2_local){
				using_local_reference = true;
				r_p1.copy(*local_ref_region);
				r_p2.copy(*local_ref_region);
				break;
			}
			//the condition when one and only one of region is local
			//C1: region-overlap
			if(r1_try.region_overlap(r2_try)){
				r_p1.copy(r1_try);
				r_p1.Combine(r2_try, true);
				r_p2.copy(r_p1);
				using_local_reference = false;
				break;
			}
			//C2: r2 < r1
			if(r2_try.ed_pos <= r1_try.st_pos)
				std::swap(r2_try, r1_try);
			//C3: r1 < r2
			r_p1.copy(r1_try);
			r_p2.copy(r2_try);
			r_p1.ed_pos += additional_load;
			r_p2.st_pos -= additional_load;
			using_local_reference = false;
			if(r_p1.region_overlap(r_p2)){
				r_p1.Combine(r_p2, true);
				r_p2.copy(r_p1);
			}
			break;
		}
	}
	//condition 2: cliping read not right aligned by aligner, partly aligned by aligner?
	//if(9629 == r_p1.getLen()){
		//fprintf(stderr, " ");
	//}
	int ref_mapping_len = 0, contig_mapping_len = 0;
	if(suggest_sv_len < 0){
		ref_mapping_len =  local_ref_region->getLen() + suggest_sv_len;
		contig_mapping_len = contig_length;
	}else{
		ref_mapping_len =  local_ref_region->getLen();
		contig_mapping_len = contig_length - suggest_sv_len;
	}

	if(contig_mapping_len > 5000 && (contig_mapping_len)> 2 * ref_mapping_len){
		std::sort(match_clusters.begin(), match_clusters.end(), LONG_SV_MATCH_CLUSTER::cmp_by_contig_position);
		//remove the match far away:
		if(with_clip_info){
			std::vector<LONG_SV_MATCH_CLUSTER> match_clusters_tmp;
			//remove:
			for(LONG_SV_MATCH_CLUSTER & mc : match_clusters){
				if(ABS_U(mc.position_diff, clip_diff_position) <= ABS(suggest_sv_len) + 5000){
					match_clusters_tmp.emplace_back();
					std::swap(match_clusters_tmp.back(), mc);
				}
			}
			std::swap(match_clusters_tmp, match_clusters);
		}

		for(LONG_SV_MATCH_CLUSTER & c : match_clusters)
			c.show();

		//PASS
		//the ref begin position:
		int ref_pos_bgn = match_clusters[0].min_ref_position - match_clusters[0].min_contig_position;
		int ref_pos_end = match_clusters.back().max_ref_position + (contig_length - match_clusters.back().max_contig_position);
		if(ref_pos_end - ref_pos_bgn < 2 * contig_length && ref_pos_end - ref_pos_bgn > contig_length / 2) {
			r_p1.chr_ID = local_ref_region->chr_ID;
			r_p1.st_pos = ref_pos_bgn;
			r_p1.ed_pos = ref_pos_end;
			r_p2.copy(r_p1);
			sv_handle_is_needed = true;
			using_local_reference = false;
		}
	}
	r_p1.st_pos = MAX(r_p1.st_pos, 0);
	r_p2.st_pos = MAX(r_p2.st_pos, 0);
	if(sv_handle_is_needed){
		//show the final result:
		fprintf(stderr, "SUGGESTed SV size :%d using_local_reference %d, "
				"r_p1 [%d:%d-%d] len %d, r_p2 [%d:%d-%d] len %d\n",
				suggest_sv_len, using_local_reference,
				r_p1.chr_ID, r_p1.st_pos, r_p1.ed_pos, r_p1.getLen(),
				r_p2.chr_ID, r_p2.st_pos, r_p2.ed_pos, r_p2.getLen()
		);
	}else
		fprintf(stderr, "NO SUGGESTed SV\n");
	return true;
}


void Contig_align_region_Handler::init_cluster_position(){
	min_ref_position = MAX_int32t;
	max_ref_position = -1000000;
	min_contig_position = MAX_int32t;
	max_contig_position = -1000000;
}

void Contig_align_region_Handler::set_cluster_position(LONG_SV_CONTIG_REF_MATCH & m){
	min_ref_position = MIN(min_ref_position,m.ref_position);
	max_ref_position = MAX(max_ref_position,m.ref_position);
	min_contig_position = MIN(min_contig_position,m.contig_position);
	max_contig_position = MAX(max_contig_position,m.contig_position);
	cluster_count++;
}

void Contig_align_region_Handler::classify_contig_ref_match(){
	match_clusters.clear();
	int old_diff = -1000000;
	cluster_count = 0;
	int match_num = CONTIG_REF_MATCH_list.size();
	init_cluster_position();

	for(int i = 0; i < match_num + 1; i++){
		if(i == match_num || old_diff + 100 < CONTIG_REF_MATCH_list[i].position_diff){
			//store the result:
			if(cluster_count > 60){
				match_clusters.emplace_back();
				match_clusters.back().set(cluster_count, old_diff, min_ref_position, max_ref_position, min_contig_position, max_contig_position);
			}
			old_diff = CONTIG_REF_MATCH_list[i].position_diff;
			cluster_count = 1;
			init_cluster_position();
		}
		else
			set_cluster_position(CONTIG_REF_MATCH_list[i]);
	}
}

void Contig_align_region_Handler::store_string(uint8_t * s, int l){
	for(int i = 0; i < l; i++)
		load_ref += ("ACGTNN"[s[i]]);
	int ori_size = bin_ref.size();
	bin_ref.resize(ori_size + l);
	memcpy(&(bin_ref[ori_size]), s, l);
}

void Contig_align_region_Handler::load_ref_string(){
	load_ref.clear();
	bin_ref.clear();
	store_string(ref_handler->getRefStr(r_p1.st_pos), r_p1.getLen()); //load part1:
	if(!r_p1.same_as(r_p2)) //load part2:
		store_string(ref_handler->getRefStr(r_p2.st_pos), r_p2.getLen());
	if(0) fprintf(stderr, "Ref is %s\n", load_ref.c_str());
}

bool Contig_align_region_Handler::get_the_long_range_ref_region(std::vector<uint8_t> &bin_contig){
	//S1:search the match of the contig and the ref index;
	uint8_t * buff_bin = &(bin_contig[0]);
	int kmer_number = bin_contig.size() - REF_IDX_KMER_LEN + 1;
	get_contig_ref_match(buff_bin, kmer_number);
	//S1:search the match of the contig and the ref index;
	if(CONTIG_REF_MATCH_list.empty())
		return false;
	//S2:classify the contig-ref-matchs
	classify_contig_ref_match();
	//add the local position when not exist
	//S3:get the true sv alignment region based on the match result
	generate_sv_region(bin_contig.size());
	return true;
}

bool Contig_align_region_Handler::load_alignment_reference(bool is_long_range_contig, bool is_left_clip_contig, bool is_right_clip_contig,//flag
		RefRegion &local_ref_region, //input ori_reference
		std::vector<uint8_t> &ori_bin_reference, std::string & ori_ref_str,//input ori_reference
		std::vector<uint8_t> &bin_contig, //input ori_contig
		RefRegion &rst_r_p1, RefRegion &rst_r_p2,
		std::string & rst_ref_str, std::vector<uint8_t> &rst_bin_ref,
		std::string &full_ref_string//the full size reference string
		){
	bool ref_need_reload = false;
	if(is_long_range_contig){
		//need to combine the reference region????
		//get the long range refernece region:
		//additional match region:
		this->local_ref_region = &local_ref_region;
		this->is_left_clip_contig = is_left_clip_contig;
		this->is_right_clip_contig = is_right_clip_contig;
		get_the_long_range_ref_region(bin_contig);
		if(false == sv_handle_is_needed){
			fprintf(stderr, "No_sv_can_be handler!\n\n");
			return false;
		}
		if(using_local_reference){
			ref_need_reload = false;
		}else{
			ref_need_reload = true;
			//the reference is not the local reference, re-load is needed
			rst_r_p1 = r_p1;
			rst_r_p2 = r_p2;
			//loading the reference
			load_ref_string();
			rst_ref_str = load_ref;//copy
			rst_bin_ref = bin_ref;//copy
			//load full size string:
			uint8_t *full_p = ref_handler->getRefStr(r_p1.st_pos);
			int load_size = r_p2.ed_pos - r_p1.st_pos;
			full_ref_string.clear();
			for(int i = 0; i < load_size; i++)
				full_ref_string += ("ACGTNN"[full_p[i]]);
		}
	}

	if(ref_need_reload == false){
		rst_r_p1 = (local_ref_region);
		rst_r_p2 = (local_ref_region);
		rst_ref_str = ori_ref_str;//copy
		rst_bin_ref = ori_bin_reference;//copy
		full_ref_string = ori_ref_str;//re-copy
	}
	if(rst_ref_str.empty()){
		fprintf(stderr, " ");
	}
	return true;
}

void get_primary_alignment_position(bam1_t *br, faidx_t * c_ref_idx, int &primary_chrID, int &primary_pos){
	//add the supplementary data:
	//show and get the SA
	primary_chrID = -1;
	primary_pos = -1;
	bool bam_is_primary = ((false == bam_is_supplementary(br)) && (false == bam_is_second(br)));
	if(bam_is_primary){
		primary_chrID = br->core.tid;
		primary_pos = br->core.pos;
	}else{
		char* SA_tag_char = bam_get_string_tag(br, "SA");
		if(SA_tag_char != NULL){
			std::vector<std::string> item_tmp;
			split_string(item_tmp, SA_tag_char, ";");
			if(item_tmp.empty() == false){
				std::string primary_sa = item_tmp[0];
				split_string(item_tmp, primary_sa.c_str(), ",");
				//chr5,47154842,+,87624S77370M374I1749543S,1,12794
				primary_chrID =  faidx_get_chrID(c_ref_idx, item_tmp[0].c_str(), NULL, 0);
				primary_pos = (int)atoi(item_tmp[1].c_str());
			}
		}
	}

}

void SV_CALLING_Handler::LRS_SIG_COLLECTION(bool print_log, BAM_handler &cur_read, LRS_SIG_Handler & sig_h,
		float &average_read_depth) {
	std::vector<LRS_SIG> &CIGAR_sig = sig_h.cigar_sig;
	std::vector<LRS_SIG> &CLIP_sig = sig_h.clip_sig;
	std::vector<LRS_SIG> &SMALL_sig = sig_h.small_sig;

	bool using_small_signal = false;
	bool using_SNP_signal = false;
	bool using_SA_signal = false;
	if(preset_platform == PRESET_ONT_Q20 || preset_platform == PRESET_CCS)
		using_small_signal = true;
	if(preset_platform == PRESET_ONT_Q20 || preset_platform == PRESET_CCS)
		using_SNP_signal = false;//todo::not used now:
	bool using_clip_signal = true;
	if(preset_platform == PRESET_ERR)
		using_clip_signal = false;
	if(preset_platform == PRESET_ASM || true == using_clip_signal)
		using_SA_signal = true;

	Bam_file *c_b = &(cur_read.file);
	cur_read.clear();
	R_region region;
	ref_handler->get_cur_region()->toR_region(region);
	resetRegion_ID(c_b, &region);	//reset region
	read_counter = 0;
	int read_local_c = 0;
	uint32_t total_read_num_primary = 0;
	uint32_t total_read_len_primary = 0;
	LRS_SMALL_SIG_HANDLER ssh;
	snp_sh.clear();

	while (bam_next(c_b)) {
		read_local_c++;
		//new sig for this read
		ssh.clear();
		bam1_t *br = &(c_b->_brec);
		if (print_log)
			fprintf(stderr, "Read Begin NAME %s \n", bam_get_qname(br));
		//filter:
		//basic filter
		if (bam_is_secondary(br))
			continue;
		if (br->core.qual < 1)
			continue;
		//supplementary-min length check
		if (bam_is_supplementary(br) && br->core.l_qseq < 1200)
			continue;

		get_bam_seq_bin(0, br->core.l_qseq, cur_read.storeReadBuff, br);
		//not running quality for ASM dataset
		if(preset_platform != PRESET_ASM && false == LRS_read_quality_check(print_log, br, cur_read.storeReadBuff)){
			continue;
		}
		if (false == bam_is_supplementary(br)) {
			total_read_len_primary += br->core.l_qseq;
			total_read_num_primary++;
		}
		//load primary alignment position
		int primary_chrID = -1;
		int primary_pos = -1;
		get_primary_alignment_position(br, ref_handler->getIdx(), primary_chrID, primary_pos);

		cur_read.store_LRS_ReadCore(br->core.l_qseq, cur_read.storeReadBuff,
				bam_get_qual(br), br->core.n_cigar, bam_get_cigar(br), bam_get_qname(br), primary_chrID, primary_pos, br->core.pos, bam_is_fwd_strand(br));
		int cigar_idx = cur_read.read_list.back().cigar_index;
		int cigar_end = cigar_idx + cur_read.read_list.back().cigar_l;
		int seq_i = 0;
		int ref_i = 0;
		int total_read_len = 0;
		if(using_SNP_signal)
			snp_sh.SNP_read_h.add_sig_read(cur_read.read_list.size() - 1, br->core.pos);

		//get the reference sequence:
		uint8_t *tseq = ref_handler->getRefStr(br->core.pos);
		uint8_t *qseq = cur_read.storeReadBuff;
		bool with_long_CLIP = false;
		for (int i = cigar_idx; i < cigar_end; i++) {
			path_segment *p = cur_read.cigar.a + i;
			if (p->type == align_t::CIGAR_SOFT_CLIP)			//soft clip:
					{
				if (p->length > 1000) {					//at lease 1000 bp
					with_long_CLIP = true;
					if (print_log)
						fprintf(stderr,
								"CLIP_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid,
								br->core.pos + ref_i, p->length,
								pathSTR[p->type], bam_get_qname(br));
					if(using_clip_signal){
						CLIP_sig.emplace_back();
						CLIP_sig.back().store(LRS_SIG_CLIP_RIGHT,
								br->core.tid, br->core.pos + ref_i,
								br->core.pos + ref_i);				//END==POS
						if (i == cigar_idx) {
							CLIP_sig.back().type = LRS_SIG_CLIP_LEFT;
							cur_read.read_list.back().soft_left = p->length;
						} else {
							CLIP_sig.back().type = LRS_SIG_CLIP_RIGHT;
							cur_read.read_list.back().soft_right = p->length;
							;
						}
					}
				}
				seq_i += p->length; total_read_len += p->length;
			} else if (p->type == align_t::CIGAR_INSERT) {
				if (p->length >= (unsigned) MIN_sv_len) {
					if (print_log)
						fprintf(stderr,
								"BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid,
								br->core.pos + ref_i, p->length,
								pathSTR[p->type], bam_get_qname(br));
					CIGAR_sig.emplace_back();
					CIGAR_sig.back().store(LRS_SIG_INS, br->core.tid,
							br->core.pos + ref_i, br->core.pos + ref_i);
				} else if(using_small_signal){
					//store small sigs
					ssh.store_sig(print_log, br->core.tid, br->core.pos + ref_i, p->length, SMALL_sig);
				}
				seq_i += p->length; total_read_len += p->length;
			} else if (p->type == align_t::CIGAR_DELETE) {
				if (p->length >= (unsigned) MIN_sv_len) {
					if (print_log)
						fprintf(stderr,
								"BIG_SIG:idx:%d TID: %d POS %d LEN:%d TYPE:%c @ NAME %s \n",
								i - cigar_idx, br->core.tid,
								br->core.pos + ref_i, p->length,
								pathSTR[p->type], bam_get_qname(br));
					CIGAR_sig.emplace_back();
					CIGAR_sig.back().store(LRS_SIG_DEL, br->core.tid,
							br->core.pos + ref_i,
							br->core.pos + ref_i + p->length);
				} else if(using_small_signal) {
					//store small sigs
					ssh.store_sig(print_log, br->core.tid, br->core.pos + ref_i, p->length, SMALL_sig);
				}
				ref_i += p->length;
			} else if (p->type == align_t::CIGAR_MATCH
					|| p->type == align_t::CIGAR_SEQ_MATCH
					|| p->type == align_t::CIGAR_SEQ_MISMATCH) {
				for(uint i = 0; i < p->length; i++, ref_i++, seq_i++){
					if(qseq[seq_i] != tseq[ref_i]){
						if(using_small_signal)
							ssh.store_sig(print_log, br->core.tid, br->core.pos + ref_i, 1, SMALL_sig);
						if(using_SNP_signal)
							snp_sh.SNP_read_h.add_sig_snp(ref_i + br->core.pos, qseq[seq_i]);
					}
				}
				//analysis the mismatch
				//ref_i += p->length;
				//seq_i += p->length;
				total_read_len += p->length;
			} else if (p->type == align_t::CIGAR_HARD_CLIP) {
				if (p->length > 1000) with_long_CLIP = true; //at lease 1000 bp
				total_read_len += p->length;
			} else {
				//DO nothing
			}
		}
		cur_read.store_READ_END(br->core.pos + ref_i);
		if(using_SNP_signal)
			snp_sh.SNP_read_h.add_sig_read_end(ref_i + br->core.pos);
		if(using_SA_signal && with_long_CLIP){
			//add the supplementary data:
			//show and get the SA
			char* SA_tag_char = bam_get_string_tag(br, "SA");
			if(SA_tag_char != NULL){
				//fprintf(stderr, "SA: Read Begin NAME %s \t", bam_get_qname(br));
				bool bam_is_primary = ((false == bam_is_supplementary(br)) && (false == bam_is_second(br)));
				//ANALYSIS SA CIGAR
				sig_h.sa_handler.store(
						cur_read.read_list.size() - 1, br->core.tid, br->core.pos, bam_is_primary,
						cur_read.cigar.a + cigar_idx, cur_read.cigar.a + cigar_end,
						SA_tag_char, total_read_len, ref_handler->getIdx());
			}
		}
	}
	if(using_SNP_signal && false)
		snp_sh.SNP_read_h.show_all_data();

	//todo::
	int region_len = ref_handler->get_cur_region()->getLen();
	xassert(region_len != 0, "");
	float ave_read_len =
			(total_read_num_primary == 0) ?
					(0) :
					((float) total_read_len_primary / total_read_num_primary);
	average_read_depth = (float) total_read_len_primary
			/ (ave_read_len + region_len);
	fprintf(stderr,
			"[Region read stat calculation:] read_number %d; ave_read_len %f; region_len %d; ave_read_dp %f\n ",
			total_read_num_primary, ave_read_len, region_len,
			average_read_depth);
}

void SV_CALLING_Handler::LRS_SV_SIG_COMBINE(bool print_log,
		LRS_SIG_Handler &sig_h,
			float average_read_depth){
	std::vector<LRS_SIG> &CLIP_sig = sig_h.clip_sig;
	std::vector<LRS_SIG> &CIGAR_sig = sig_h.cigar_sig;
	std::vector<LRS_SIG> &SMALL_sig = sig_h.small_sig;

	fprintf(stderr, "LRS variations calling(CLIP) process begin\n" );
	//basic
	int min_support_read_number = (average_read_depth)*0.1;/////？？？？？？ TODO：：：：；
	min_support_read_number = MAX(min_support_read_number,1);//at least 1
	fprintf(stderr, "average_read_depth %f, min_support_read_number %d, \n", average_read_depth, min_support_read_number);
	int signal_nearby_length = 300;
	int region_edge_length = 600;
	int small_signal_nearby_length = 20;
	std::vector<LRS_SIG> cmb_store_tmp;

	{//P0: handle SMALL signals
		std::sort(SMALL_sig.begin(), SMALL_sig.end(), LRS_SIG::cmp_by_ref_pos);
		std::vector<LRS_SIG> &l= SMALL_sig;
		//SMALL signal combine:
		for(uint i = 0; i < l.size();i++){
			if(l[i].isUsed)
				continue;
			l[i].support_read_Num = 1;
			for(uint j = i+1; j < l.size();j++){
				if(l[i].similar_SIG_SMALL(l[j], small_signal_nearby_length)){
					l[i].POS = MIN(l[i].POS, l[j].POS);
					l[i].END = MAX(l[i].END, l[j].END);
					l[i].support_read_Num ++;
					l[j].isUsed = true;
				}
				else
					break;
			}
			if(l[i].support_read_Num >= min_support_read_number && l[i].support_read_Num >= 2)
				cmb_store_tmp.emplace_back(l[i]);
		}
		l.swap(cmb_store_tmp);
	}

	{	//P1： simple combine and sort the clip and cigar:
		std::vector<LRS_SIG> &l= CLIP_sig;
		l.insert(l.end(), CIGAR_sig.begin(), CIGAR_sig.end());
		std::sort(l.begin(), l.end(), LRS_SIG::cmp_by_ref_pos);
		cmb_store_tmp.clear();
		//cigar simple combine:
		for(uint i = 0; i < l.size();i++){
			if(l[i].isUsed)
				continue;
			l[i].support_read_Num = 1;
			for(uint j = i+1; j < l.size();j++){
				if(l[i].similar_SIG_Normal(l[j], signal_nearby_length)){
					l[i].POS = MIN(l[i].POS, l[j].POS);
					l[i].END = MAX(l[i].END, l[j].END);
					l[i].support_read_Num ++;
					l[j].isUsed = true;
				}
				else
					break;
			}
			if(l[i].support_read_Num >= min_support_read_number)
				cmb_store_tmp.emplace_back(l[i]);
		}
		//simple combine:
		l.swap(cmb_store_tmp);
	}

	//P1: //nearby combine of SV region
	{
		std::vector<LRS_SIG> &l= CLIP_sig;
		l.insert(l.end(), SMALL_sig.begin(), SMALL_sig.end());
		std::sort(l.begin(), l.end(), LRS_SIG::cmp_by_ref_pos);
		cmb_store_tmp.clear();
		for(uint i = 0; i < l.size();i++){
			if(l[i].isUsed)
				continue;
			for(uint j = i+1; j < l.size();j++){
				if(l[i].region_overlap(l[j], region_edge_length)){
					l[i].POS = MIN(l[i].POS, l[j].POS);
					l[i].END = MAX(l[i].END, l[j].END);
					l[i].support_read_Num += l[j].support_read_Num;
					l[j].isUsed = true;
				}
			}
			cmb_store_tmp.emplace_back(l[i]);
		}
		//simple combine:
		l.swap(cmb_store_tmp);
	}
	//store clip to final
	std::swap(sig_h.clip_sig, sig_h.FINAL_sig);
	//show
	if(print_log || true){
		for(LRS_SIG ss:sig_h.FINAL_sig)
			ss.show_simple();
	}
}

int calGT_LRS_CORE(int not_supp_read_n, int supp_read_n){
	int final_genotype = -1;
	if(not_supp_read_n > 6*supp_read_n){
		final_genotype = 0;
	}else if(6*not_supp_read_n < supp_read_n){
		final_genotype = 2;
		if(supp_read_n == 1) final_genotype = 1;
	}else{
		final_genotype = 1;
	}
	return final_genotype;
}

void get_SV_quality_score_core(int sup_ref_number, int sup_alt_number, int sup_both_number, bool SVLEN_is_Imprecise, int final_GT, int &QUAL_INT, int &GQ_INT){
	float RR_log_P = -sup_alt_number*5-5.3;
	float AA_log_P = -sup_ref_number*5;
	float RA_log_P = -0.3*(sup_ref_number + sup_alt_number)-5;
	float QUAL = (MAX(RA_log_P, AA_log_P) - RR_log_P) * 10;
	QUAL = MAX(MIN_QUAL, QUAL);
	if((float)sup_both_number/(float)sup_alt_number > 0.1 && QUAL > 100){
		QUAL = QUAL/10;
	}
	QUAL = MIN(MAX_QUAL, QUAL);
	QUAL_INT = QUAL;
	if(SVLEN_is_Imprecise)
		QUAL_INT = MIN_QUAL;
	float GQ = 0;
		 if(final_GT == 0){ GQ = (RR_log_P - (MAX(AA_log_P,RA_log_P)))*10; }
	else if(final_GT == 1){ GQ = (RA_log_P - (MAX(RR_log_P,AA_log_P)))*10; }
	else{                   GQ = (AA_log_P - (MAX(RR_log_P,RA_log_P)))*10; }
	GQ = MAX(MIN_QUAL, GQ);
	GQ = MIN(MAX_QUAL, GQ);
	GQ_INT = GQ;
}

void SV_CALLING_Handler::LRS_SV_CALLING_germline(bool print_log,
		std::vector<LRS_SIG> &FINAL_sig, float average_read_depth, bool output_pure_contig){
	int min_support_read_number = (average_read_depth)*0.1;/////？？？？？？ TODO：：：：；
	std::vector<LRS_SIG> &l= FINAL_sig;
	//P3： calling
	int signal_list_size = l.size();
	RefRegion ref_r = *(ref_handler->get_cur_region());
	faidx_t * c_ref_idx = ref_handler->getIdx();
	int MAX_CONTIG_SIZE = 50000;
	if(preset_platform == PRESET_ASM){
		int MAX_CONTIG_LEN = 3000000;
		MAX_CONTIG_SIZE = MAX_CONTIG_LEN;
	}
	for(int i = 0; i < signal_list_size; i++){
		//set_region_addition_load_super_super_long();
		set_region_addition_load_short();
		l[i].show_simple();
		//skip when out of range
		RefRegion r(l[i].chrID, l[i].POS,  l[i].END);
		if(!(r.region_overlap(ref_r)))
			continue;
		//skip read number is not enough
		if(l[i].support_read_Num < min_support_read_number){
			if(output_pure_contig){
				fprintf(vcf_w, "REF_REGION=%s:%d-%d\t",	faidx_iseq(c_ref_idx, l[i].chrID), (int)l[i].POS,(int)l[i].END);
				fprintf(vcf_w, "STAT=NO_ENOUGH_READ_SUPPORT\n");
			}
			continue;
		}
		RefRegion loadRef;	std::string full_ref_str; std::vector<uint8_t> bin_ref;
		//output region header
		if(LRS_assembly_variations(print_log, &LRS_read, l[i], loadRef, full_ref_str, bin_ref) == false){
			if(output_pure_contig){
				fprintf(vcf_w, "REF_REGION=%s:%d-%d\t",	faidx_iseq(c_ref_idx, loadRef.chr_ID), (int)l[i].POS,(int)l[i].END);
				fprintf(vcf_w, "STAT=ASSEMBLY_FAIL\n");
			}
			continue;
		}
		if(output_pure_contig){
			int true_out_num = 0;
			for(uint haplotype_ID = 0; haplotype_ID < contig_selecter.select_list.size(); haplotype_ID++){
				if(0 == contig_selecter.select_list[haplotype_ID].support_read_number)		continue;
				if(false == contig_selecter.select_list[haplotype_ID].is_final_path)		continue;
				if((int)target_path_idx_l[contig_selecter.select_list[haplotype_ID].target_ID].contig.size() > MAX_CONTIG_SIZE)//MAX contig length is 50K bp
					continue;
				true_out_num ++;
			}
			if(0 == true_out_num){
				fprintf(vcf_w, "REF_REGION=%s:%d-%d\t",	faidx_iseq(c_ref_idx, loadRef.chr_ID), (int)l[i].POS,(int)l[i].END);
				fprintf(vcf_w, "STAT=NO_CONTIG_SELECTED\n");
			}
			int contig_id = -1;
			for(uint haplotype_ID = 0; haplotype_ID < contig_selecter.select_list.size(); haplotype_ID++){
				if(0 == contig_selecter.select_list[haplotype_ID].support_read_number)
					continue;
				if(false == contig_selecter.select_list[haplotype_ID].is_final_path)
					continue;
				CONTIG_SELECT &c = contig_selecter.select_list[haplotype_ID];
				if((int)target_path_idx_l[c.target_ID].contig.size() > MAX_CONTIG_SIZE)//MAX contig length is 50K bp
					continue;
				contig_id++;
				int supp_read_n, not_supp_read_n;
				contig_selecter.get_haplotype_support_read_n(haplotype_ID, supp_read_n, not_supp_read_n);
				int GT = calGT_LRS_CORE(not_supp_read_n, supp_read_n);
				//QUAL score
				int QUAL_INT = 0;
				int GQ_INT = 0;		//simple GT for TGS reads
				get_SV_quality_score_core( not_supp_read_n, supp_read_n, 0, true, GT, QUAL_INT, GQ_INT);
				std::string GT_STR;
				getGT_string(GT, GT_STR);
				fprintf(vcf_w, "REF_REGION=%s:%d-%d\t",	faidx_iseq(c_ref_idx, loadRef.chr_ID), loadRef.st_pos, loadRef.ed_pos);
				fprintf(vcf_w,
						"STAT=PASS\t"
						"GT:QU:GQ:GI:CI:CN=%s:%d:%d:%d:%d:%d\t"
						"SUP_R:UNSUP_R:FULL_COVER_R:CLIP_AT_LEFT:CLIP_AT_RIGHT:IS_FINAL=%d:%d:%d:%d:%d:%d\t"
						"CONTIG_STR=%s\n",
						GT_STR.c_str(), QUAL_INT, GQ_INT, GT,contig_id,true_out_num,
						supp_read_n, not_supp_read_n, c.full_cover_read_n, c.clip_left_read_number, c.clip_right_read_number, c.is_final_path,
						target_path_idx_l[c.target_ID].contig.c_str());
			}
		}else
			LRS_SV_generating_Germline(print_log, l[i], loadRef, full_ref_str, bin_ref, result_LRS_final, target_path_idx_l, contig_selecter);
	}
}

void SV_CALLING_Handler::LRS_remove_duplications(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> &LRS_r_l){
	std::sort(LRS_r_l.begin(), LRS_r_l.end(), NOVA_SV_FINAL_RST_item::cmp_by_position_TGS);
	if(print_log){
		//set basic GT:
		for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: LRS_r_l)
			sv.calGT_LRS();
		for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: LRS_r_l){
			if(sv.writeVCF_final(&vcfBuffStr, cur_header, NULL))
				fprintf(stderr, "\nBefore remove_TGS_duplications for: %s ", vcfBuffStr.s);
		}
	}
	std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
	for(uint i = 0; i < LRS_r_l.size();i++){
		if(LRS_r_l[i].chr_ID == -1)
			continue;
		for(uint j = i+1; j < LRS_r_l.size();j++){
			if(LRS_r_l[i].st_pos != LRS_r_l[j].st_pos)
				break;
			if(LRS_r_l[i].is_same_var(LRS_r_l[j]) && LRS_r_l[i].is_same_supp_read_LRS(LRS_r_l[j])){
				LRS_r_l[j].chr_ID = -1;
			}
		}
		sv_tmp.emplace_back();
		std::swap(sv_tmp.back(), LRS_r_l[i]);
	}
	//simple combine:
	LRS_r_l.swap(sv_tmp);
}

void  SV_CALLING_Handler::LRS_combine_Homozygous_SVs(std::vector<NOVA_SV_FINAL_RST_item> &LRS_r_l){
	std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
	//simple combine:
	for(uint i = 0; i < LRS_r_l.size();i++){
		if(LRS_r_l[i].chr_ID == -1)
			continue;
		for(uint j = i+1; j < LRS_r_l.size();j++){
			if(LRS_r_l[i].st_pos != LRS_r_l[j].st_pos)
				break;
			if(LRS_r_l[i].is_same_var(LRS_r_l[j]) && LRS_r_l[i].is_same_global_region_ID_TGS(LRS_r_l[j])){
				LRS_r_l[i].combine_supp_read_TGS(LRS_r_l[j]);
				LRS_r_l[j].chr_ID = -1;
			}
		}
		sv_tmp.emplace_back();
		std::swap(sv_tmp.back(), LRS_r_l[i]);
	}
	//simple combine:
	LRS_r_l.swap(sv_tmp);
}


bool SV_CALLING_Handler::LRS_read_quality_check(bool print_log, bam1_t *br, uint8_t * seq_bin){
	bool pass_check = true;
	//quality check:
	int seq_i = 0;
	int ref_i = 0;

	uint8_t * ref_bin = ref_handler->getRefStr(br->core.pos);
	uint32_t * bam_cigar = bam_get_cigar(br);

	int snp_count = 0;
	int match_base = 0;

	int total_match_base = 0;
	int total_snp_count = 0;

	for (uint c_i = 0; c_i <  br->core.n_cigar; c_i ++) {
		int length = (bam_cigar[c_i] >> BAM_CIGAR_SHIFT);
		int type = (int)(1 + (bam_cigar[c_i] & BAM_CIGAR_MASK));
		if (type == align_t::CIGAR_SOFT_CLIP)		{seq_i += length;}
		else if (type == align_t::CIGAR_HARD_CLIP)  {/*DO NOTHING*/}
		else if (type == align_t::CIGAR_INSERT)  	{seq_i += length; snp_count++;}
		else if (type == align_t::CIGAR_DELETE) 		{ref_i += length; snp_count++;}
		else if (type == align_t::CIGAR_MATCH || type == align_t::CIGAR_SEQ_MATCH || type == align_t::CIGAR_SEQ_MISMATCH)
		{
			//SNP check:
			for(int i = 0; i < length; i++){
				match_base ++;
				if(ref_bin[ref_i + i] != seq_bin[seq_i + i]){
					snp_count++;
				}
				if(match_base == 1000){
					float err_rate = ((float)snp_count)/match_base;
					if(err_rate > 0.2){
						if(print_log)
							fprintf(stderr, "LOAD reads LRS[High ERROR rate]: NAME_ %s len %d, ref_pos %d read pos %d, error rate %f\n",
								bam_get_qname(br), br->core.l_qseq, seq_i, ref_i, err_rate);
						pass_check = false;
					}
					total_match_base += match_base;
					total_snp_count += snp_count;
					match_base = 0;
					snp_count = 0;
				}
			}
			ref_i += length; seq_i += length;
		}
		else { /*DO nothing*/}
	}

	float total_err_rate = ((float)total_snp_count)/total_match_base;
	//fprintf(stderr, "LOAD reads LRS[FINAL ERROR rate]: NAME_ %s len %d, error rate %f\n", bam_get_qname(br), br->core.l_qseq, total_err_rate);
	float MAX_ALN_ERROR_RATE = 0.04;
	if(preset_platform == PRESET_ERR)
		MAX_ALN_ERROR_RATE = 0.08;
	if(total_err_rate > MAX_ALN_ERROR_RATE){
		fprintf(stderr, "LOAD reads LRS[FINAL ERROR rate]: NAME_ %s len %d, error rate %f\n",
			bam_get_qname(br), br->core.l_qseq, total_err_rate);
		pass_check = false;
	}
	return pass_check;
}

void SV_CALLING_Handler::LRS_SMALL_SIG_HANDLER::clear(){
	std::queue<int> t1; std::swap(t1, position_queue);//clear the queue
	std::queue<int> t2; std::swap(t2, influence_base_queue);//clear the queue
	total_influence_base = 0;
}
void SV_CALLING_Handler::LRS_SMALL_SIG_HANDLER::store_sig(bool print_log, int chrID,  int POS, int influence_base, std::vector<LRS_SIG> &SMALL_sig){
	position_queue.push(POS);
	influence_base_queue.push(influence_base);
	total_influence_base += influence_base;
	while(POS - position_queue.front() > 200){
		position_queue.pop();
		total_influence_base -= influence_base_queue.front();
		influence_base_queue.pop();
	}
	if(false) {
		fprintf(stderr, "NEW SIG-SMALL-POS @%d [%d %d %ld %d]\n", POS, position_queue.front(), position_queue.back(),position_queue.size(), influence_base);
	}
	if(position_queue.size() > 10 || total_influence_base > 35){
		int pos_bg = position_queue.front();
		int pos_ed = position_queue.back();
		SMALL_sig.emplace_back();
		SMALL_sig.back().store(LRS_SIG_X, chrID, pos_bg, pos_ed);
		if(print_log) fprintf(stderr, "NEW SIG-SMALL, position:[%d:%d-%d] N_SIG: %ld total base: %d\n", chrID, pos_bg, pos_ed, position_queue.size(), total_influence_base);
		clear();
	}
}

void SV_CALLING_Handler::Hybrid_debug_code_load_SVs_from_vcf_f(std::vector<NOVA_SV_FINAL_RST_item> & result_sv_l, const char * vcf_fn, RefHandler *ref){
	result_sv_l.clear();
	BCF_FILE vcf_r;//vcf for read
	VCF_open_read(&vcf_r, vcf_fn);//open for read

	char *c_sv_type = (char *)malloc(1000);
	bcf_hdr_t *header = vcf_r.header;
	int *support_read_number = (int*)xcalloc(3,4);
	int32_t SV_LEN = 0;
	RefRegion * r = ref->get_cur_region();
	char *GT[1]; char GT_1[3]; GT[0] = GT_1;
	do{//read one
		bcf1_t *c_r = &( vcf_r.r);
		if(r->chr_ID != c_r->rid || c_r->pos < r->st_pos || c_r->pos > r->ed_pos)
			continue;
		//unpack the vcf data to get the alt string
		bcf_unpack(c_r, BCF_UN_STR);
		vcf_get_sv_GT(vcf_r.header, c_r, GT);
		vcf_get_sv_type(vcf_r.header, c_r, c_sv_type);
		vcf_get_sv_LENGTH(vcf_r.header, c_r, &SV_LEN);
		vcf_get_sv_SR(vcf_r.header, c_r, support_read_number);
		bool is_vntr = vcf_get_sv_flag(header, c_r, "VNTR");
		fprintf(stderr, "GT:%d,%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t\n", GT[0][0], GT[0][1], is_vntr, c_r->rid , c_r->pos, SV_LEN, c_sv_type, c_r->d.allele[0], c_r->d.allele[1]);// chrID+st+len
		if(c_r->d.allele[1][0] == '<' && c_r->d.allele[1][1] == 'D' && c_r->d.allele[1][2] == 'E'){
			std::string s; s += c_r->d.allele[0][0];
			for(int i = 0;i < -SV_LEN; i++)
				s += 'N';
			NOVA_SV_FINAL_RST_item::add_to_vector(result_sv_l, c_r->rid, c_r->pos + 1, c_sv_type, s.c_str(),
					c_r->d.allele[0], 0,  0, 0, NULL, 0, 0, 0, 0, r->st_pos);
		}else{
			NOVA_SV_FINAL_RST_item::add_to_vector(result_sv_l, c_r->rid, c_r->pos + 1, c_sv_type, c_r->d.allele[0], c_r->d.allele[1],
								0,  0, 0, NULL, 0, 0, 0, 0, r->st_pos);
		}

		result_sv_l.back().isVNTRcallerRst = is_vntr;
		if(*c_r->d.flt == 0)//filter: PASS
			result_sv_l.back().QUAL_presetVNTR = 30;
		else
			result_sv_l.back().QUAL_presetVNTR = 0;
		result_sv_l.back().set_supp_read_NGS(support_read_number[0], support_read_number[1], support_read_number[2]);
		int GT_final = 0;
		if(GT[0][0] == 4) GT_final++;
		if(GT[0][1] == 4) GT_final++;
		result_sv_l.back().setGenotype_directly(GT_final);
		//store data into:
	}while(VCF_next(&vcf_r));
	//close
	bcf_close(vcf_r.file);
}

void SV_CALLING_Handler::Hybrid_remove_similar_SVs_in_Target(bool print_log, std::vector<NOVA_SV_FINAL_RST_item> & query, std::vector<NOVA_SV_FINAL_RST_item> & target){
	std::vector<NOVA_SV_FINAL_RST_item> sv_tmp;
	//simple combine:
	for(uint i = 0; i < query.size();i++){
		NOVA_SV_FINAL_RST_item * SRS_r = &(query[i]);
		for(uint j = 0; j < target.size();j++){
			NOVA_SV_FINAL_RST_item * LRS_r = &(target[j]);
			//position is nearby:
			int ABS_POS = ABS_U(SRS_r->st_pos, LRS_r->st_pos);
			if(SRS_r->chr_ID != LRS_r->chr_ID || ABS_POS > 300)
				continue;
			if(print_log){
				fprintf(stderr, "NGS results is removed:\n");
				SRS_r->printSimple_M2(stderr);
				fprintf(stderr, "When comparing with LRS SV:\n");
				LRS_r->printSimple_M2(stderr);
			}
			SRS_r->chr_ID = -1;
			break;
		}
		if(SRS_r->chr_ID != -1){
			sv_tmp.emplace_back();
			std::swap(sv_tmp.back(), query[i]);
		}
	}
	//simple combine:
	query.swap(sv_tmp);
}

void SV_CALLING_Handler::Hybrid_regenotype_NGS_USING_TGS_data(bool print_log){
	//show all the results
	if(print_log){
		fprintf(stderr, "Combine with LRS begin\n\n");
		showSVList(result_LRS_final, "Show LRS result");
		showSVList(result_SRS_INDEL, "Show SRS result");
	}
	bool is_full_begin = true;
	bool is_full_end = true;
	//merging and remove duplication NGS results that similar with LRS results
	Hybrid_remove_similar_SVs_in_Target(print_log, result_SRS_INDEL, result_LRS_final);
	//
	//genotyping the NGS result using LRS data set
	for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_SRS_INDEL ){
		if(print_log) showSV(sv, "GT NGS SVs using LRS reads[Before]:\n"); //show all the SVs
		int LRS_read_number = 0;
		//get the depth of LRS reads
		for(int mode = 0; mode < 2; mode ++){
			int check_position = sv.st_pos;
			if(mode == 1){
				if(sv.SV_length > 0)	continue;
				else					check_position = sv.st_pos - sv.SV_length;
			}
			int check_edge = 500;
			RefRegion r(sv.chr_ID, check_position - check_edge, check_position + check_edge);
			if(print_log) r.show();
			//load all reads in this regions:
			for(READ_record &rr :LRS_read.read_list){
				std::string s;
				LRS_read.load_read_LRS(s, rr, r.st_pos, r.ed_pos, is_full_begin, is_full_end);
				if(!s.empty() && s.size() > check_edge * 1.9)//only full cover LRS reads is used
					LRS_read_number++;
			}
			if(mode == 1)
				LRS_read_number = LRS_read_number/2;
		}
		sv.set_LRS_INFO(-1, 0, LRS_read_number, 0, -1, -1, -1, -1, 0);
		int max_TGS_read_num = 3;
		if(LRS_read_number > max_TGS_read_num)	sv.setGenotype_directly(0);
		else if(LRS_read_number > 0)			sv.setGenotype_directly(1);
		if(print_log) showSV(sv,  "GT NGS SVs using LRS reads[After]:\n");
	}
	if(print_log) fprintf(stderr, "Show NGS result after merging\n\n");
}
