/*
 * SveHandler_LRS.cpp
 *
 *  Created on: 2025/7/9
 *      Author: fenghe
 */

#include <SVcalling_core/SV_core.hpp>

void SNP_SIG_HANDLER::SNP_Position_Handler::store_all_position_SNP(std::vector<int> &used_read_list, SNP_Read_Handler & snp_read_h){
	ana_l.clear();
	all_read_number = used_read_list.size();
	for(int ri = 0; ri < all_read_number; ri++){
		READ_SNP_L & cur_r = snp_read_h.l[used_read_list[ri]];
		for(uint si = 0; si < cur_r.p.size(); si ++){
			int SNP_position = cur_r.p[si];
			ana_l[SNP_position].emplace_back(ri, cur_r.c[si]);
		}
	}
}

void SNP_SIG_HANDLER::SNP_Position_Handler::store_2_link_format(std::vector<int> &used_read_list, Link_Handler &link_handler, SNP_Read_Handler & snp_read_h){
	int min_support_read_n = 5;
	std::unordered_map<int, std::vector<ANA_ITEM>>::iterator ana_l_it = ana_l.begin();
	for(;ana_l_it != ana_l.end(); ana_l_it++){
		if((int)ana_l_it->second.size() >= min_support_read_n){
			//look for over all depth
			read_id_char_map.clear();
			for(int ri = 0; ri < all_read_number; ri++){
				READ_SNP_L & cur_r = snp_read_h.l[used_read_list[ri]];
				if(cur_r.st_pos < ana_l_it->first && cur_r.ed_pos > ana_l_it->first){
					read_id_char_map[ri] = 5;//set to "N"
				}
			}
			int all_skip_read_num = all_read_number - read_id_char_map.size();
			int all_snp_num = ana_l_it->second.size();
			int all_ref_num = read_id_char_map.size() - all_snp_num;

			if(all_ref_num >= min_support_read_n){
				//calculate all the links
				//std::vector<ANA_ITEM> & crl = ana_l_it->second;//copy:
				//store all the reads
				for(ANA_ITEM & A : ana_l_it->second){
					read_id_char_map[A.read_id] = A.c;//set to the true char
				}
				//store the map to vector
				std::map<int, uint8_t>::iterator it_i, it_j, it_begin = read_id_char_map.begin();
				for(it_i = it_begin; it_i != read_id_char_map.end(); it_i++){
					for(it_j = it_begin; it_j != read_id_char_map.end(); it_j++)
						link_handler.store_link(it_i->first, it_j->first, it_i->second, it_j->second);
				}
				//show all the results
				fprintf(stderr, "The position is %d, \t #SNP %d, #REF %d #SKIP %d; the SNP read list is: ",
						ana_l_it->first, all_snp_num, all_ref_num, all_skip_read_num);
				//show all SNP reads
				for(ANA_ITEM & ai : ana_l_it->second)
					fprintf(stderr, "[%d %c ] ", ai.read_id, "ACGTNNN"[ai.c]);
				fprintf(stderr, "\n");
			}
		}
	}
}

void SNP_SIG_HANDLER::Cluster_Handler::init(int read_number){
	read_final_cluster_id.resize(read_number);
	for(uint i = 0; i < read_final_cluster_id.size(); i++)
		read_final_cluster_id[i] = -1;
	clu_l.clear();
	clu_l.emplace_back();
	clu_l.back().nl.emplace_back(0);//store the first node
	read_final_cluster_id[0] = 0;
}

void SNP_SIG_HANDLER::Cluster_Handler::add_cluster_detection(int rid_add, Link_Handler &link_handler){
	add_clusters.clear();
	//int add_clusters_n = 0;
	int reject_clusters_n = 0;
	int unknown_clusters_n = 0;
	for(uint cid = 0; cid < clu_l.size(); cid++){
		int c_score_sum = 0;
		for(int rid_in_c = 0; rid_in_c < (int)clu_l[cid].nl.size(); rid_in_c++){
			int rid_check = clu_l[cid].nl[rid_in_c].read_id;
			int score_i_j = link_handler.get_link_score(rid_add, rid_check);
			c_score_sum += score_i_j;
		}
		fprintf(stderr, "c_score_sum %d cid %d\n", c_score_sum, cid);
		if(c_score_sum > 0){
			add_clusters.emplace_back();
			add_clusters.back().cid = cid;
			add_clusters.back().c_score_sum = c_score_sum;
		}else if(c_score_sum < 0){
			reject_clusters_n ++;
		}else{
			unknown_clusters_n ++;
		}
	}
	fprintf(stderr, "End handle read id %d add_clusters_n %ld, reject_clusters_n %d, unknown_clusters_n %d\n",
			rid_add, add_clusters.size(), reject_clusters_n, unknown_clusters_n);
	if(add_clusters.empty() && reject_clusters_n > 0){
		//new a cluster
		clu_l.emplace_back();
		clu_l.back().nl.emplace_back(rid_add);//store the first node
		read_final_cluster_id[rid_add] = clu_l.size() - 1;
	}
}

void SNP_SIG_HANDLER::Cluster_Handler::add_new_read_to_cluster(int rid_add, Link_Handler &link_handler){
	//add scores
	//detect the max support cluster
	int max_score = 0;
	Read_Cluster_New_Add * max_rc = NULL;
	for(Read_Cluster_New_Add & rc: add_clusters){
		if(rc.c_score_sum > max_score){
			max_score = rc.c_score_sum ;
			max_rc = &rc;
		}
	}
	if(max_rc != NULL){
		//store read to cluster
		clu_l[max_rc->cid].nl.emplace_back(rid_add);
		//reject the cluster to the read
		read_final_cluster_id[rid_add] = max_rc->cid;
		for(int rid_in_c = 0; rid_in_c < (int)clu_l[max_rc->cid].nl.size(); rid_in_c++){
			int rid_check = clu_l[max_rc->cid].nl[rid_in_c].read_id;
			if(rid_add == rid_check)	clu_l[max_rc->cid].nl[rid_in_c].score = max_rc->c_score_sum;
			else						clu_l[max_rc->cid].nl[rid_in_c].score += link_handler.get_link_score(rid_add, rid_check);
		}
	}
}
void SNP_SIG_HANDLER::Cluster_Handler::calculate_sum_score_for_all_cluster(){
	//get the final score for each cluster
	for(uint cid = 0; cid < clu_l.size(); cid++){
		int clu_sum_score = 0 ;
		int clu_low_score_node = 0 ;
		for(int rid_in_c = 0; rid_in_c < (int)clu_l[cid].nl.size(); rid_in_c++){
			clu_sum_score += clu_l[cid].nl[rid_in_c].score;
			if(clu_l[cid].nl[rid_in_c].score < 40)
				clu_low_score_node++;
		}
		clu_l[cid].clu_sum_score = clu_sum_score;
		clu_l[cid].clu_low_score_node = clu_low_score_node;
	}
}
void SNP_SIG_HANDLER::Cluster_Handler::run(Link_Handler &link_handler ){
	//adding other nodes
	for(int rid_add = 1; rid_add < link_handler.get_node_n(); rid_add++){
		fprintf(stderr, "begin handle read id %d\n", rid_add);
		add_cluster_detection(rid_add, link_handler);
		if(!add_clusters.empty())
			add_new_read_to_cluster(rid_add, link_handler);
	}
	calculate_sum_score_for_all_cluster();
}

void SNP_SIG_HANDLER::analysis(std::vector<int> &used_read_list){
	//collect all signals
	int all_read_number = used_read_list.size();
	//clear
	link_handler.init(all_read_number);
	//store the SNP-reads to SNP-position
	SNP_p_h.store_all_position_SNP(used_read_list, SNP_read_h);
	//store the SNP-position to links
	SNP_p_h.store_2_link_format(used_read_list, link_handler, SNP_read_h);
	//show all the links
	link_handler.show_all_link();
	//clusters
	cluster_handler.init(all_read_number);
	cluster_handler.run(link_handler);
	//show all the clusters
	for(uint cid = 0; cid < cluster_handler.clu_l.size(); cid++)
		cluster_handler.clu_l[cid].show();
	for(uint i = 0;i < cluster_handler.read_final_cluster_id.size(); i ++)
		fprintf(stderr, "[read_id:%d cid:%d]", i, cluster_handler.read_final_cluster_id[i]);
	fprintf(stderr, "\n\n");
}

bool SNP_SIG_HANDLER::is_reject_read_target_pair(int read_id, std::vector<SUPPORT_READ_LIST_ITEM> &support_read_l){
	int support_read_n = 0;
	int reject_read_n = 0;
	for(SUPPORT_READ_LIST_ITEM & r : support_read_l){
		if(cluster_handler.read_final_cluster_id[r.read_id] == cluster_handler.read_final_cluster_id[read_id])
			support_read_n ++;
		else
			reject_read_n++;
	}
	bool is_reject = (support_read_n  < reject_read_n);
	fprintf(stderr, "[is_reject_read_target_pair:]read_id %d support_read_n  %d reject_read_n  %d is_reject %d\n",read_id, support_read_n,reject_read_n, is_reject);
	return is_reject;
	//std::vector<SUPPORT_READ_LIST_ITEM> support_read_l;//support read list:
}

void SV_CALLING_Handler::SRS_show_signals(){
    //main type to call INS and DEL
    SRS_print_signal_list(SIG::DR, SV::DEL);
    SRS_print_signal_list(SIG::SH, SV::INS);
    //other type to call INS
    SRS_print_signal_list(SIG::DR, SV::INS);
    //types to call INV and BND
    SRS_print_signal_list(SIG::DR, SV::INV_1);
    SRS_print_signal_list(SIG::DR, SV::INV_2);
    SRS_print_signal_list(SIG::DR, SV::TRA);
    SRS_print_signal_list(SIG::DR, SV::TRA_INV);
}

void SV_CALLING_Handler::SRS_REMOVE_DUPLICATED_INDEL(){
    fprintf(stderr, "True SVs output\n" );
    std::sort(result_SRS_INDEL.begin(), result_SRS_INDEL.end(),
    		NOVA_SV_FINAL_RST_item::cmp_by_position);
    //part4: remove duplication
    int SV_number = result_SRS_INDEL.size();
    for(int c_SV_idx = 0; c_SV_idx < SV_number - 1; c_SV_idx++)
        for(int cmp_sv_idx = c_SV_idx + 1;
                cmp_sv_idx < SV_number && result_SRS_INDEL[c_SV_idx].duplication_SV_filter(result_SRS_INDEL[cmp_sv_idx]);
                cmp_sv_idx++);
}

void SV_CALLING_Handler::SRS_GENOTYPING_INS_DEL(){
    for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_SRS_INDEL ){
        //running genotyping process:
        //if(sv.SV_is_duplicated()) continue;
        char log_title[1024]; sprintf(log_title, "[genotyping]");
        double __cpu_time = cputime(); double __real_time = realtime();
        if(sv.isVNTRcallerRst)
        	fprintf(stderr, "\nGenotyping for(VNTR): \n ");
        else
        	fprintf(stderr, "\nGenotyping for(short INDEL): \n ");
        sv.print(stderr, SRS_read.file._hdr, ref_handler);
        sv.genotyping(sig_para.MaxReadLen,sig_para.read_depth, sig_para.insert_size_max, &SRS_read.file, &ga, ref_handler);
        fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
    }
}

void SV_CALLING_Handler::SRS_INV_signal_combine(std::vector<SVE> &inv_sve_l){
	while(!SRS_sve[SIG::DR][SV::INV_1].empty() || !SRS_sve[SIG::DR][SV::INV_2].empty()){
		bool with_INV[2] = {0};
		SVE c_sve;
		if(!SRS_sve[SIG::DR][SV::INV_1].empty()){
			std::swap(c_sve, SRS_sve[SIG::DR][SV::INV_1].back());
			SRS_sve[SIG::DR][SV::INV_1].pop_back();
			with_INV[0] = true;
		}else{
			std::swap(c_sve, SRS_sve[SIG::DR][SV::INV_2].back());
			SRS_sve[SIG::DR][SV::INV_2].pop_back();
			with_INV[1] = true;
		}
		//try to combine
		for(int mode = 0; mode < 2; mode++){
			for(SVE & combine:SRS_sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2]){
				if(combine.r1.region_overlap(c_sve.r1) && combine.r2.region_overlap(c_sve.r2)){
					c_sve.r1.Combine(combine.r1, true);
					c_sve.r2.Combine(combine.r2, true);
					std::swap(combine, SRS_sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2].back());
					SRS_sve[SIG::DR][(mode==0)?SV::INV_1:SV::INV_2].pop_back();
					with_INV[mode] = true;
				}
			}
		}
		inv_sve_l.emplace_back();
		std::swap(c_sve, inv_sve_l.back());
		fprintf(stderr, "INV LOGS: \n SIGNAL balance %d %d ; ", with_INV[0], with_INV[1]);
	}
}

bool SV_CALLING_Handler::SRS_INV_ASSEMBLY(SVE &c_sve, int &ASS_SIGNAL_NUM){
    if(SRS_SV_region_depth_filter(c_sve) == false)  return false;;
    int main_pos = c_sve.r1.getMiddle();
    int supp_pos = c_sve.r2.getMiddle();

    bool print_log = false;
    int edge_len = 400;
    RefRegion main_region(c_sve.r1.chr_ID, main_pos - edge_len,main_pos + edge_len);
    RefRegion supp_region(c_sve.r2.chr_ID, supp_pos - edge_len,supp_pos + edge_len);

    //length filter for INV
    int inv_length = ABS_U(main_pos, supp_pos);
    if(inv_length > 50000){//todo::
        fprintf(stderr, "Following INV is skiped because the length is over 50K\n");
        return false;
    }

    //assembly
    ASS_SIGNAL_NUM = SRS_assembly_and_get_breakpoints_INV(print_log, bnd_ass_block, am, main_region, supp_region, main_pos, supp_pos);
    c_sve.r1.st_pos = main_pos; c_sve.r1.st_pos = main_pos;
    c_sve.r2.st_pos = supp_pos; c_sve.r2.st_pos = supp_pos;

    fprintf(stderr, "INV: main_pos %d - supp_pos %d\n", main_pos, supp_pos);
    return true;
}

void SV_CALLING_Handler::SRS_INV_GT(SVE &c_sve, int ASS_SIGNAL_NUM){
	//genotyping
	int main_pos = c_sve.r1.st_pos;
	int supp_pos = c_sve.r2.st_pos;

	int GT_result[3][4];
    SRS_INV_genotyping_and_store(c_sve, true, true, true, true, SRS_read.file._hdr, GT_result[0]);
    std::swap(c_sve.r1, c_sve.r2);
    SRS_INV_genotyping_and_store(c_sve, false, false, false, false, SRS_read.file._hdr, GT_result[1]);
    GT_result[2][0] = GT_result[0][0] + GT_result[1][0];
    GT_result[2][1] = GT_result[0][1] + GT_result[1][1];
    GT_result[2][2] = GT_result[0][2] + GT_result[1][2];
    for(int i = 0; i < 3; i++){
        int *region_support_number = GT_result[i];
        if(region_support_number[0] > region_support_number[1] * 4)     region_support_number[3] = 2;//1/1
        else if(region_support_number[0]*4 < region_support_number[1])  region_support_number[3] = 0;//0/0
        else                                    region_support_number[3] = 1;//0/1
        if((region_support_number[0] + region_support_number[1]) * 3 < region_support_number[2])
            region_support_number[3] = 0;//0/0
        if(region_support_number[0] < 2)
            region_support_number[3] = 0;//0/0
    }

    //filter:
    if(ASS_SIGNAL_NUM != 0 && SRS_store_INV(c_sve.r1.chr_ID, main_pos, supp_pos, result_BGS_INV, GT_result[2], GT_result[2][3])){
        fprintf(stderr, "Read support balance %d:%d:%d %d:%d:%d ", GT_result[0][0], GT_result[0][1], GT_result[0][2], GT_result[1][0], GT_result[1][1], GT_result[1][2]);
        fprintf(stderr, "GT balance %d %d ", GT_result[0][3], GT_result[1][3]);
        fprintf(stderr, "ASS BALANCE: ASS_SIGNAL_NUM %d ; ", ASS_SIGNAL_NUM);
        fprintf(stderr, "\n");
        for(NOVA_SV_FINAL_RST_item &indel:result_SRS_INDEL){
            if(indel.SV_overlap(result_BGS_INV.back()) && indel.get_length() >= 50){
                fprintf(stderr, "INV is removed because it overlap with other INS/DELs\n");
                if(result_BGS_INV.back().writeVCF_final(&vcfBuffStr, SRS_read.file._hdr, NULL))
                    fprintf(stderr, "%s", vcfBuffStr.s);
                if(indel.writeVCF_final(&vcfBuffStr, SRS_read.file._hdr, NULL))
                    fprintf(stderr, "%s", vcfBuffStr.s);
                if(!result_BGS_INV.empty())
                	result_BGS_INV.pop_back();
                break;
            }
        }
    }else{
        fprintf(stderr, "GT balance %d:%d:%d %d:%d:%d ", GT_result[0][0], GT_result[0][1], GT_result[0][2], GT_result[1][0], GT_result[1][1], GT_result[1][2]);
        fprintf(stderr, "ASS BALANCE: ASS_SIGNAL_NUM %d ; ", ASS_SIGNAL_NUM);
        fprintf(stderr, "\n");
    }
}

void SV_CALLING_Handler::SRS_set_min_accpet_possibility(){
	float MAX_POSSIBILITY_DR = 0;
	for(float poss:DR_bp_distribution)
		MAX_POSSIBILITY_DR = MAX(MAX_POSSIBILITY_DR, poss);
	MIN_ACCEPT_POSSIBILITY_DR = 2*MAX_POSSIBILITY_DR;

	float MAX_POSSIBILITY_SH = 0;
	for(float poss:SH_bp_distribution)
		MAX_POSSIBILITY_SH = MAX(MAX_POSSIBILITY_SH, poss);
	MIN_ACCEPT_POSSIBILITY_SH = 2*MAX_POSSIBILITY_SH;

	float MAX_POSSIBILITY_UM = 0;
	for(float poss:UM_stPos_distribution)
		MAX_POSSIBILITY_UM = MAX(MAX_POSSIBILITY_UM, poss);
	MIN_ACCEPT_POSSIBILITY_UM = 3*MAX_POSSIBILITY_UM;
}

void SV_CALLING_Handler::SRS_print_signal_list(SIG::T sigT, SV::T svt){
	SVE_L * cur_var_signals = &(SRS_sve[sigT][svt]);
	fprintf(stderr, "Print %s::%s signals, size %ld\n", SIG::STR[sigT].c_str(), SV::STR[svt].c_str(), cur_var_signals->size());
	for(unsigned int i = 0; i < cur_var_signals->size(); i++)
		std::cerr << cur_var_signals[0][i];
}

bool SV_CALLING_Handler::SRS_store_INV(int chr_ID, int position_bg, int position_ed, std::vector<NOVA_SV_FINAL_RST_item> &region_SVs_TMP, int region_support_number[3], int final_genotype){
	std::string ref_inv;
	std::string alt_inv;
	int INV_len = position_ed - position_bg;
	if(final_genotype == 0) { return false;}//length filter
	if(INV_len < 30 || INV_len > 50000){ return false;}//length filter
	if(INV_len <= 2000){
		int load_len;
		char* load_ref = ref_handler->load_ref_by_region(chr_ID, position_bg, position_ed, &load_len);
		//set REF
		ref_inv.append(load_ref);
		//SET ALT
		alt_inv.resize(load_len + 1);
		char *cur_alt_string_p = &(alt_inv[0]);
		for(int i = 0; i < load_len; i++){
			switch( load_ref[load_len - i - 1]){
			case 'A': case 'a': cur_alt_string_p[i] = 'T'; break;
			case 'C': case 'c': cur_alt_string_p[i] = 'G'; break;
			case 'G': case 'g': cur_alt_string_p[i] = 'C'; break;
			case 'T': case 't': cur_alt_string_p[i] = 'A'; break;
			case 'n': case 'N': cur_alt_string_p[i] = 'N'; break;
			}
		}
		alt_inv[load_len] = 0;
		free(load_ref);
	}else{
		ref_inv.append("<INV_REF>");
		alt_inv.append("<INV_ALT>");
	}
	NOVA_SV_FINAL_RST_item::add_to_vector(region_SVs_TMP, chr_ID, position_bg, "INV",
			&(ref_inv[0]), &(alt_inv[0]), INV_len, NULL, 0, NULL, 0, 0, 0, 0, 0);
	region_SVs_TMP.back().INV_SET_GT(region_support_number, final_genotype);
	return true;
}

//void BND and INVS
void SV_CALLING_Handler::SRS_BND_CALLING_AND_GT_PROCESS(){
	SVE_L * cur_var_signals = NULL;
	int signal_list_size;
	am->setRepeatMode();
	bool read_is_before_breakpoint_in_main;
	bool read_in_main_should_be_forward;
	bool read_is_before_breakpoint_in_supp;
	bool read_in_supp_should_be_forward;
	bool supp_ref_is_forward; bool main_ref_is_forward;
	fprintf(stderr, "BND calling process\n" );
	for(int mode = 0; mode < 2; mode++){
		if(mode == 0)
		{
			cur_var_signals = &(SRS_sve[SIG::DR][SV::TRA]);
			read_is_before_breakpoint_in_main = true;
			read_in_main_should_be_forward = true;
			read_is_before_breakpoint_in_supp = false;
			read_in_supp_should_be_forward = false;
			main_ref_is_forward = true;
			supp_ref_is_forward = true;
		}
		else if(mode == 1)
		{
			cur_var_signals = &(SRS_sve[SIG::DR][SV::TRA_INV]);
			read_is_before_breakpoint_in_main = false;
			read_in_main_should_be_forward = false;
			read_is_before_breakpoint_in_supp = true;
			read_in_supp_should_be_forward = true;
			main_ref_is_forward = true;
			supp_ref_is_forward = false;
		}
		if(!cur_var_signals->empty()){
			for(unsigned int i = 0; i < cur_var_signals->size(); i++)
				std::cerr << cur_var_signals[0][i];
			signal_list_size = cur_var_signals->size();
			for(int i = 0; i < signal_list_size; i++){
				if(SRS_SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
				SRS_assembly_and_genetyping_BND(cur_var_signals[0][i],
						read_is_before_breakpoint_in_main, read_in_main_should_be_forward, main_ref_is_forward,
						read_is_before_breakpoint_in_supp, read_in_supp_should_be_forward, supp_ref_is_forward,
						result_NGS_BND, SRS_read.file._hdr);
			}
		}
	}
}

void SV_CALLING_Handler::SRS_BND_WRITE(){
    std::vector<NOVA_SV_FINAL_RST_item> region_SVs_TMP;
    for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_NGS_BND){
        if(sv.is_BND_TRA()){//only output precise results
            if(!sv.is_BND_IMPRECISE()){
                fprintf(stderr, "OUTPUT BND::TRA\n: \n ");
                if(sv.writeVCF_final(&vcfBuffStr, SRS_read.file._hdr, &global_SV_ID))
                    fprintf(vcf_w, "%s", vcfBuffStr.s);
            }
        }else{
            region_SVs_TMP.emplace_back();
            std::swap(sv, region_SVs_TMP.back());
        }
    }
    std::swap(region_SVs_TMP, result_NGS_BND);
}

//LARGE DELETION calling process
void SV_CALLING_Handler::SRS_LARGE_DELETION_CALLING_PROCESS(){
	SVE_L * cur_var_signals = NULL;
	int signal_list_size;
	 //part 1 large deletions: (using DR/DEL signals)
	am->setNormalMode();
	fprintf(stderr, "Deletions calling process begin\n" );
	cur_var_signals = &(SRS_sve[SIG::DR][SV::DEL]);
	if(!cur_var_signals->empty()){
		for(unsigned int i = 0; i < cur_var_signals->size(); i++)
			std::cerr << cur_var_signals[0][i];
		signal_list_size = cur_var_signals->size();
		for(int i = 0; i < signal_list_size; i++){
			if(SRS_SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
//              if(cur_var_signals[0][i].r1.st_pos == 108733125){
//                  std::cerr << "TTTT" <<cur_var_signals[0][i];
//              }
			set_region_addition_load_long();
			int total_region_length = cur_var_signals[0][i].r2.ed_pos - cur_var_signals[0][i].r1.st_pos;
			if(total_region_length > 1800){
				if(total_region_length < 5000)
					set_region_addition_load_extramly_long();
				else if(total_region_length < 50000)
					set_region_addition_load_super_super_long();
				else
					continue;
			}
			//calling the long deletions
			bool withRepeat;
			SRS_assembly_variations_MEI_AND_LONG_SV(cur_var_signals[0][i], false, withRepeat);
		}
	}
	fprintf(stderr, "Deletions calling process Done\n" );
}

void SV_CALLING_Handler::SRS_SMALL_VAR_CALLING_PROCESS(){
	//part 2: small variations: (using SH/INS signals)
	am->setRepeatMode();
	SVE_L * cur_var_signals = NULL;
	int signal_list_size;
	fprintf(stderr, "Small variations calling process begin\n" );
	cur_var_signals = &(SRS_sve[SIG::SH][SV::INS]);
	for(unsigned int i = 0; i < cur_var_signals->size(); i++)
		std::cerr << cur_var_signals[0][i];
	//std::sort(cur_var_signals.begin(), cur_var_signals.end(), SVE::cmp_by_position);//needed to be sorted?
	if(!cur_var_signals->empty()){
		signal_list_size = cur_var_signals->size();
		for(int i = 0; i < signal_list_size; i++){
			fprintf(stderr, "\nCurrent SVE ID: [%d]", i);
			if(SRS_SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
			set_region_addition_load_short();
			int cur_TL_read_num = 0;
			bool withRepeat = false;
			{
				int MIN_TL_read_number = (sig_para.read_depth *0.15); MIN_TL_read_number = MAX(4, MIN_TL_read_number);//at least 4 TL reads
				fprintf(stderr, "MIN_TL_read_number %d, TL_read_number %d\n", MIN_TL_read_number, cur_TL_read_num);
				//tran-ins and ALU calling process:
				if(cur_TL_read_num >= 0){
					char log_title[1024]; sprintf(log_title, "[INS_DEL_LONG_ALU]");
					double __cpu_time = cputime(); double __real_time = realtime();
					SRS_assembly_variations_MEI_AND_LONG_SV(cur_var_signals[0][i], true, withRepeat);
					fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
				}
			}

			{
				RefRegion &r_main = cur_var_signals[0][i].r1;
				rlr->search(r_main.chr_ID, r_main.st_pos, r_main.ed_pos, rlr_result);
				if(withRepeat){//trigger VNTR caller when repeat is found in assembly caller
					fprintf(stderr, "VNTR Handler is used in repeat region!\n");
					for(RefLocalRepeatItemR &rr : rlr_result)
						rr.show();
					//VNTR calling process
					char log_title[1024]; sprintf(log_title, "[INS_DEL_VNTR]");
					double __cpu_time = cputime(); double __real_time = realtime();
					SRS_assembly_variations_INS_DEL_VNTR(cur_var_signals[0][i],true,true);
					cur_TL_read_num = this->TL_read_number;
					fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
				}
				else{
					fprintf(stderr, "Skip VNTR caller\n");
				}
			}
		}
	}
	fprintf(stderr, "Small variations calling process Done\n" );
}

void SV_CALLING_Handler::SRS_VNTR_VAR_CALLING_PROCESS(){
	am->setRepeatMode();
	SVE_L * cur_var_signals = NULL;
	int signal_list_size;
	//part 2: small variations: (using SH/INS signals)
	fprintf(stderr, "VNTR variations calling process\n" );
	cur_var_signals = &(SRS_sveVNTR);
	for(unsigned int i = 0; i < cur_var_signals->size(); i++)
		std::cerr << cur_var_signals[0][i];
	//std::sort(cur_var_signals.begin(), cur_var_signals.end(), SVE::cmp_by_position);//needed to be sorted?
	if(!cur_var_signals->empty()){
		signal_list_size = cur_var_signals->size();
		for(int i = 0; i < signal_list_size; i++){
			fprintf(stderr, "\nCurrent SVE ID: [%d]", i);
			if(SRS_SV_region_depth_filter(cur_var_signals[0][i]) == false) continue;
			set_region_addition_load_short();
			{
				fprintf(stderr, "VNTR Handler is used in repeat region!\n");
				//VNTR calling process
				char log_title[1024]; sprintf(log_title, "[INS_DEL_VNTR]");
				double __cpu_time = cputime(); double __real_time = realtime();
				SRS_assembly_variations_INS_DEL_VNTR(cur_var_signals[0][i],true,true);
				fprintf(stderr, "%s: CPU time: %.3f sec; real time: %.3f sec\n", log_title, cputime() - __cpu_time, realtime() - __real_time);
			}
		}
	}
	fprintf(stderr, "VNTR variations calling Done\n" );
}

void SV_CALLING_Handler::SRS_INV_CALLING_AND_GT_PROCESS(){
	am->setRepeatMode();
	fprintf(stderr, "INV calling process begin\n" );
   	//signals combine and store data:
	std::vector<SVE> inv_sve_l;
	SRS_INV_signal_combine(inv_sve_l);
    for(SVE &c_sve : inv_sve_l){
    	int ASS_SIGNAL_NUM = 0;
    	if(SRS_INV_ASSEMBLY(c_sve, ASS_SIGNAL_NUM) == false)
    		continue;
    	SRS_INV_GT(c_sve, ASS_SIGNAL_NUM);
    }
    fprintf(stderr, "INV processing Done\n");
}

void SV_CALLING_Handler::SRS_INV_WRITE(){
    //POST_PROCESS_INV();
    for(std::vector<NOVA_SV_FINAL_RST_item>::value_type &sv: result_BGS_INV){
        fprintf(stderr, "OUTPUT BND::INV\n: \n ");
        if(sv.writeVCF_final(&vcfBuffStr, SRS_read.file._hdr, &global_SV_ID))
            fprintf(vcf_w, "%s", vcfBuffStr.s);
    }
}


void SRS_REF_STRING_COMBINE_Handler::print_ref(){
	//fprintf(stderr, "region_overlap %d combined_ref is: \n", region_overlap);
	fprintf(stderr, "combined_ref is: \n");
	for(uint i = 0; i < s.size(); i++)
		fprintf(stderr, "%c",  "ACGTNN"[s[i]]);
	fprintf(stderr, "\n");
}

void SRS_REF_STRING_COMBINE_Handler::print_results(){
	if(with_contig_supp)
		fprintf(stderr, "Final main_pos %d supp_pos %d, length %d \n", main_pos, supp_pos, supp_pos- main_pos);
	else
		fprintf(stderr, "NO data \n");
}

void SRS_REF_STRING_COMBINE_Handler::store_ref(
		std::vector<uint8_t> & s1, bool s1_forward, int s1_ref_st,
		std::vector<uint8_t> & s2, bool s2_forward, int s2_ref_st,
		bool is_reverse
		){
	s.clear();
	BP_region_st[0] = s.size();
	s = s1;
	BP_region_ed[0] = s.size();
	s.emplace_back(4);
	s.emplace_back(4);
	s.emplace_back(4);
	BP_region_st[1] = s.size();
	s.insert(s.end(), s2.begin(), s2.end());
	BP_region_ed[1] = s.size();

	this->s1_ref_st = s1_ref_st;
	this->s2_ref_st = s2_ref_st;
	this->s1_forward = s1_forward;
	this->s2_forward = s2_forward;

	this->is_reverse = is_reverse;
	this->with_contig_supp = false;
	if(is_reverse){
		std::vector<uint8_t> s_rev;
		for(int i = 0; i < (int)s.size(); i++){
			uint8_t c = s[s.size() - i - 1];
			if(c > 3)
				s_rev.emplace_back(4);
			else
				s_rev.emplace_back(3- s[s.size() - i - 1]);
		}

		std::swap(s,s_rev);
		this->s1_forward = !this->s2_forward;
		this->s2_forward = !this->s1_forward;
		std::swap(this->s1_ref_st,this->s2_ref_st);
	}
}

void SRS_REF_STRING_COMBINE_Handler::check_and_store(int BP1, int BP2){
	if(BP1 >= BP_region_st[0] && BP1 <= BP_region_ed[0] && BP2 >= BP_region_st[1] && BP2 <= BP_region_ed[1] ){
		BP2 -= BP_region_st[1]; //align the BP2 to base-0
		main_pos = (s1_forward)?(s1_ref_st + BP1):(s1_ref_st + BP_region_ed[0] - BP_region_st[0] - BP1);
		supp_pos = (s2_forward)?(s2_ref_st + BP2):(s2_ref_st + BP_region_ed[1] - BP_region_st[1] - BP2);
		if(main_pos > supp_pos) std::swap(main_pos, supp_pos);
		fprintf(stderr, "\n After is main_pos(M1) %d supp_pos %d, length %d s1_ref_st %d , s2_ref_st %d , s1_forward %d , s2_forward %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, supp_pos- main_pos, s1_ref_st, s2_ref_st, s1_forward, s2_forward, BP1, BP2);
		with_contig_supp = true;
		fprintf(stderr, "check_and_store Success\n");
	}
}

bool SRS_REF_STRING_COMBINE_Handler::check_and_store_TRAN_INS(int BP1, int BP2, bool is_right_part, int &source_BP_pos, int &target_BP_pos){
	if(BP1 >= BP_region_st[0] && BP1 <= BP_region_ed[0] && BP2 >= BP_region_st[1] && BP2 <= BP_region_ed[1] ){
		BP2 -= BP_region_st[1]; //align the BP2 to base-0
		main_pos = (s1_forward)?(s1_ref_st + BP1):(s1_ref_st + BP_region_ed[0] - BP_region_st[0] - BP1);
		supp_pos = (s2_forward)?(s2_ref_st + BP2):(s2_ref_st + BP_region_ed[1] - BP_region_st[1] - BP2);
		if(is_right_part){ source_BP_pos = main_pos; target_BP_pos = supp_pos;}
		else			 { target_BP_pos = main_pos; source_BP_pos = supp_pos;}

		if(false){
			if(is_right_part)	fprintf(stderr, "\n R(S+T): source_pos %d target_pos %d, s1_ref_st %d , s2_ref_st %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, s1_ref_st, s2_ref_st, BP1, BP2);
			else				fprintf(stderr, "\n L(T+S): target_pos %d source_pos  %d, s1_ref_st %d , s2_ref_st %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, s1_ref_st, s2_ref_st, BP1, BP2);
			fprintf(stderr, "check_and_store Success\n");
		}
		return true;
	}
	return false;
}

void SRS_REF_STRING_COMBINE_Handler::chece_and_store_M2(int BP1, int BP2){
	int region_id = 0;
	if(BP1 <= BP_region_ed[0] && BP2 <= BP_region_ed[0]){
		region_id = 0;
	}
	else if(BP1 >= BP_region_st[1] && BP2 >= BP_region_st[1]){
		region_id = 1;
		BP1 -= BP_region_st[1]; //align the BP2 to base-0
		BP2 -= BP_region_st[1]; //align the BP2 to base-0
	}
	else
		return;
	bool is_forward = (region_id == 0)?s1_forward:s2_forward;
	int ref_st = (region_id == 0)?s1_ref_st:s2_ref_st;
	int region_len = (region_id == 0)?(BP_region_ed[0] - BP_region_st[0]):(BP_region_ed[1] - BP_region_st[1]);
	main_pos = (is_forward)?(ref_st + BP1):(ref_st + region_len - BP1);
	supp_pos = (is_forward)?(ref_st + BP2):(ref_st + region_len - BP2);
	if(main_pos > supp_pos) std::swap(main_pos, supp_pos);
	fprintf(stderr, "After is main_pos(M2) %d supp_pos %d, length %d s1_ref_st %d , s2_ref_st %d , s1_forward %d , s2_forward %d , BP1 %d, BP2 %d\n", main_pos, supp_pos, supp_pos- main_pos, s1_ref_st, s2_ref_st, s1_forward, s2_forward, BP1, BP2);
	with_contig_supp = true;
	fprintf(stderr, "chece_and_store Success\n");
}

