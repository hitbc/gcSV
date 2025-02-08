#include "../SVcalling_core/SveHandler.hpp"


bool SveHandler::assembly_and_genetyping_BND(SVE &sv,
		bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward, bool main_ref_is_forward,
		bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward, bool supp_ref_is_forward,
		std::vector<NOVA_SV_FINAL_RST_item> &region_SVs, bam_hdr_t * header){//buff to store contig
	bool print_log = false;
	if(print_log){
		fprintf(stderr, "assembly_variations_BND BG\n");
	}

	int main_tid = sv.r1.chr_ID;
	int main_pos = (sv.r1.st_pos + sv.r1.ed_pos)/2;

	int supp_tid = sv.r2.chr_ID;
	int supp_pos = (sv.r2.st_pos + sv.r2.ed_pos)/2;

	//length filter for INV
	if(sv.info.type_ID == SV::INV_1 || sv.info.type_ID == SV::INV_2){
		int inv_length = ABS_U(main_pos, supp_pos);
		if(inv_length > 50000){//todo::
			std::cerr << sv;
			fprintf(stderr, "Following INV is skiped because the length is over 50K\n");
			return false;
		}
	}

	//region check
	RefRegion block_region = *(ref->get_cur_region());
	if(block_region.pos_within_same_chr_region(main_tid, main_pos)){
		//do nothing
	}else if( block_region.pos_within_same_chr_region(supp_tid, supp_pos)){
		std::swap(main_tid, supp_tid);
		std::swap(main_pos, supp_pos);
		std::swap(read_is_before_breakpoint_in_main, read_is_before_breakpoint_in_supp);
		std::swap(read_in_main_should_be_forward, read_in_supp_should_be_forward);
		std::swap(main_ref_is_forward, supp_ref_is_forward);
	}
	else
		return false;

	Bam_file *bam_f = &NGS_read.file;
	int normal_read_length = sig_para.MaxReadLen;
	int max_isize = sig_para.insert_size_max;

	bool afte_read_before_BP = (read_is_before_breakpoint_in_main)?false:true;

	RefRegion withInMain(main_tid, main_pos, main_pos);			bool main_dir = read_in_main_should_be_forward;
	RefRegion withInSUPP(supp_tid, supp_pos, supp_pos);			bool supp_dir = read_in_supp_should_be_forward;
	RefRegion withInAfterMain(main_tid, main_pos, main_pos);	bool afte_dir = afte_read_before_BP;

	if(read_is_before_breakpoint_in_main){
		withInMain.st_pos -= (max_isize - normal_read_length);
		withInAfterMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
	}else{
		withInMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
		withInAfterMain.st_pos -= (max_isize - normal_read_length);
	}

	withInSUPP.st_pos -= (max_isize - normal_read_length);
	withInSUPP.ed_pos += (max_isize - normal_read_length - normal_read_length);

	fprintf(stderr, "BND searching: withInMain\t");  withInMain.print(stderr);
	fprintf(stderr, "BND searching: withInSUPP\t");  withInSUPP.print(stderr);
	fprintf(stderr, "BND searching: withInAfterMain\t");  withInAfterMain.print(stderr);

	//only consider reads in the first region:
	R_region bam_load_region;
	bam_load_region.chr_ID = main_tid;
	bam_load_region.st_pos = main_pos - max_isize;
	bam_load_region.ed_pos = main_pos + normal_read_length + normal_read_length;

	//Analyzer
	int read_num_support_BND = 0;
	int64_t total_isize_BND = 0;
	int read_num_support_REF = 0;
	int64_t total_isize_REF = 0;

	//
	resetRegion_ID(bam_f, &bam_load_region);
	while (bam_next(bam_f)) {
		bam1_t *br = &(bam_f->_brec);
		if(bam_is_secondary(br))		continue;
	    if(bam_is_duplicate(br))      continue;
		//if(bam_is_supplementary(br))	continue;
		//get iSIZE
		int isize = br->core.isize;
		int tid = br->core.tid;
		int pos_read_ed = br->core.pos + br->core.l_qseq;
		int pos_read_bg = br->core.pos;
		int pos_read = (read_is_before_breakpoint_in_main)?pos_read_ed:pos_read_bg;
		bool read_forward = bam_is_fwd_strand(br);

		int mtid = br->core.mtid;
		int mpos = br->core.mpos;
		int mpos_read = (read_is_before_breakpoint_in_main)?mpos + br->core.l_qseq:mpos;
		bool mate_forward = bam_is_mate_fwd_strand(br);
		int suppor_int = 0;

		if(false && mate_forward == read_forward){
			fprintf(stderr, "XXX");
			fprintf(stderr, " %d %d %d %d %d ", withInSUPP.pos_within_same_chr_region(mtid, mpos_read),
					main_dir == read_forward, supp_dir == mate_forward,withInSUPP.pos_within_same_chr_region(tid, pos_read),
					pos_read > mpos_read);
		}

		if(	withInSUPP.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && supp_dir == mate_forward){
			if(withInSUPP.pos_within_same_chr_region(tid, pos_read) && (pos_read > mpos_read))
			{ /*DO NOTHING*/}else{
				read_num_support_BND++;
				total_isize_BND += isize;
				suppor_int += 1;
			}
		}

		mpos_read = (afte_read_before_BP)?mpos + br->core.l_qseq:mpos;
		if(withInMain.pos_within_same_chr_region(tid, pos_read) && withInAfterMain.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && afte_dir == mate_forward){//for read pairs that support the reference
			read_num_support_REF++;
			total_isize_REF += isize;
			suppor_int += 2;
		}

		//if(suppor_int == 1){
		bam_load_region.st_pos = main_pos - max_isize;
		bam_load_region.ed_pos = main_pos + normal_read_length + normal_read_length;

		if(print_log && false) {
			fprintf(stderr, "[tid %d, pos_read_bg %d pos_read_end %d dir: %d] , [mtid %d, mpos_bg %d mpos_ed %d dir %d ] the read support: \t", tid, pos_read_bg, pos_read_ed , read_forward, mtid, mpos, mpos + br->core.l_qseq, mate_forward);
			switch(suppor_int){
			case 0: fprintf(stderr, "None\n"); break;
			case 1: fprintf(stderr, "BND\n"); break;
			case 2: fprintf(stderr, "REF\n"); break;
			case 3: fprintf(stderr, "BOTH\n"); break;
			default: fprintf(stderr, "UNKNOWN\n"); break;
			}
		}
	}

	fprintf(stderr, "\nRead_num_support_BND %d %f\n",
			read_num_support_BND, (float)total_isize_BND/read_num_support_BND);
	fprintf(stderr, "Read_num_support_reference(IN BND) %d %f\n",
			read_num_support_REF, (float)total_isize_REF/read_num_support_REF);

	int genotype = 0;
	if((read_num_support_BND >= read_num_support_REF * 3))			genotype = 2;
	else if((read_num_support_BND * 6 <= read_num_support_REF))		genotype = 0;
	else														 	genotype = 1;

	if(read_num_support_BND < 3) //at least 3 read support
		genotype = 0;
	fprintf(stderr, "Final BND genotyping result is(0 for 0/0; 1 for 0/1; 2 for 1/1) : %d\n\n", genotype);

	if(genotype > 0){
//		bool BP_is_polishing = assembly_and_get_breakpoints_TRA(print_log, bnd_ass_block, am,
//				main_tid, main_pos,	read_is_before_breakpoint_in_main, main_ref_is_forward,
//				supp_tid, supp_pos, read_is_before_breakpoint_in_supp, supp_ref_is_forward);

		int load_len; char * load_ref;
		char ref_char_main[2]; load_ref = ref->load_ref_by_region(main_tid, main_pos, main_pos+1, &load_len); ref_char_main[0] = load_ref[0]; ref_char_main[1] = 0;free(load_ref);
		char ref_char_supp[2]; load_ref = ref->load_ref_by_region(supp_tid, supp_pos, supp_pos+1, &load_len); ref_char_supp[0] = load_ref[0]; ref_char_supp[1] = 0;	free(load_ref);

		char BND_str[1024];
		for(int mode = 0; mode < 2; mode ++){
			if(mode == 1){
				if( block_region.pos_within_same_chr_region(supp_tid, supp_pos) && sv.info.readNumB >= 3 && sv.info.readNumE >= 3){
					std::swap(main_tid, supp_tid);
					std::swap(main_pos, supp_pos);
					std::swap(read_is_before_breakpoint_in_main, read_is_before_breakpoint_in_supp);
					std::swap(read_in_main_should_be_forward, read_in_supp_should_be_forward);
					std::swap(main_ref_is_forward, supp_ref_is_forward);
					std::swap(ref_char_main[0], ref_char_supp[0]);
				}
				else
					continue;
			}
			if(sv.info.type_ID == SV::TRA){
				if(read_in_main_should_be_forward)
					sprintf(BND_str, "%c[%s:%d[", ref_char_supp[0], header->target_name[supp_tid], supp_pos);
				else
					sprintf(BND_str, "]%s:%d]%c", header->target_name[supp_tid], supp_pos, ref_char_supp[0]);
			}
			else if(sv.info.type_ID == SV::INV_1)
				sprintf(BND_str, "%c]%s:%d]", ref_char_supp[0], header->target_name[supp_tid], supp_pos);
			else if(sv.info.type_ID == SV::INV_2){
				sprintf(BND_str, "[%s:%d[%c", header->target_name[supp_tid], supp_pos,ref_char_supp[0]);
			}else
				continue;
			NOVA_SV_FINAL_RST_item::store_BND(region_SVs, main_tid, main_pos, supp_pos, ref_char_main, BND_str, SV::STR[sv.info.type_ID].c_str(), read_num_support_BND, read_num_support_REF, 0, genotype, true);
			if(print_log){
				fprintf(stderr, "\n\nEND PROCESS BND: %d %d %s \n\n", main_tid, main_pos, BND_str);
			}
		}
	}

	return true;
}


bool SveHandler::INV_genotyping_and_store(SVE &sv,
		bool read_is_before_breakpoint_in_main, bool read_in_main_should_be_forward,
		bool read_is_before_breakpoint_in_supp, bool read_in_supp_should_be_forward,
		bam_hdr_t * header, int *GT_result){//buff to store contig
	bool print_log = false;
	if(print_log){
		fprintf(stderr, "assembly_variations_BND BG\n");
	}

	int main_tid = sv.r1.chr_ID;
	int main_pos = (sv.r1.st_pos + sv.r1.ed_pos)/2;

	int supp_tid = sv.r2.chr_ID;
	int supp_pos = (sv.r2.st_pos + sv.r2.ed_pos)/2;

	//length filter for INV
	if(sv.info.type_ID == SV::INV_1 || sv.info.type_ID == SV::INV_2){
		int inv_length = ABS_U(main_pos, supp_pos);
		if(inv_length > 50000){//todo::
			std::cerr << sv;
			fprintf(stderr, "Following INV is skiped because the length is over 50K\n");
			return false;
		}
	}

	//region check
	RefRegion block_region = *(ref->get_cur_region());
	if(block_region.pos_within_same_chr_region(main_tid, main_pos)){
		//do nothing
	}else if( block_region.pos_within_same_chr_region(supp_tid, supp_pos)){
		std::swap(main_tid, supp_tid);
		std::swap(main_pos, supp_pos);
		std::swap(read_is_before_breakpoint_in_main, read_is_before_breakpoint_in_supp);
		std::swap(read_in_main_should_be_forward, read_in_supp_should_be_forward);
	}
	else
		return false;

	Bam_file *bam_f = &NGS_read.file;
	int normal_read_length = sig_para.MaxReadLen;
	int max_isize = sig_para.insert_size_max;

	bool afte_read_before_BP = (read_is_before_breakpoint_in_main)?false:true;

	RefRegion withInMain(main_tid, main_pos, main_pos);			bool main_dir = read_in_main_should_be_forward;
	RefRegion withInSUPP(supp_tid, supp_pos, supp_pos);			bool supp_dir = read_in_supp_should_be_forward;
	RefRegion withInAfterMain(main_tid, main_pos, main_pos);	bool afte_dir = afte_read_before_BP;

	if(read_is_before_breakpoint_in_main){
		withInMain.st_pos -= (max_isize - normal_read_length);
		withInAfterMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
	}else{
		withInMain.ed_pos += (max_isize - normal_read_length - normal_read_length);
		withInAfterMain.st_pos -= (max_isize - normal_read_length);
	}

	withInSUPP.st_pos -= (max_isize - normal_read_length);
	withInSUPP.ed_pos += (max_isize - normal_read_length - normal_read_length);

	fprintf(stderr, "BND searching: withInMain\t");  withInMain.print(stderr);
	fprintf(stderr, "BND searching: withInSUPP\t");  withInSUPP.print(stderr);
	fprintf(stderr, "BND searching: withInAfterMain\t");  withInAfterMain.print(stderr);

	//only consider reads in the first region:
	R_region bam_load_region;
	bam_load_region.chr_ID = main_tid;
	bam_load_region.st_pos = main_pos - max_isize;
	bam_load_region.ed_pos = main_pos + normal_read_length + normal_read_length;

	//Analyzer
	int read_num_support_BND = 0;
	int64_t total_isize_BND = 0;
	int read_num_support_REF = 0;
	int read_num_support_BOTH = 0;
	int64_t total_isize_REF = 0;

	//
	resetRegion_ID(bam_f, &bam_load_region);
	while (bam_next(bam_f)) {
		bam1_t *br = &(bam_f->_brec);
		if(bam_is_secondary(br))		continue;
	    if(bam_is_duplicate(br))       continue;

		//if(bam_is_supplementary(br))	continue;
		//get iSIZE
		int isize = br->core.isize;
		int tid = br->core.tid;
		int pos_read_ed = br->core.pos + br->core.l_qseq;
		int pos_read_bg = br->core.pos;
		int pos_read = (read_is_before_breakpoint_in_main)?pos_read_ed:pos_read_bg;
		bool read_forward = bam_is_fwd_strand(br);

		int mtid = br->core.mtid;
		int mpos = br->core.mpos;
		int mpos_read = (read_is_before_breakpoint_in_main)?mpos + br->core.l_qseq:mpos;
		bool mate_forward = bam_is_mate_fwd_strand(br);
		int suppor_int = 0;

		if(false && mate_forward == read_forward){
			fprintf(stderr, "XXX");
			fprintf(stderr, " %d %d %d %d %d ", withInSUPP.pos_within_same_chr_region(mtid, mpos_read),
					main_dir == read_forward, supp_dir == mate_forward,withInSUPP.pos_within_same_chr_region(tid, pos_read),
					pos_read > mpos_read);
		}

		if(	withInSUPP.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && supp_dir == mate_forward){
			if(withInSUPP.pos_within_same_chr_region(tid, pos_read) && (pos_read > mpos_read))
			{ /*DO NOTHING*/}else{
				read_num_support_BND++;
				total_isize_BND += isize;
				suppor_int += 1;
			}
		}

		mpos_read = (afte_read_before_BP)?mpos + br->core.l_qseq:mpos;
		if(withInMain.pos_within_same_chr_region(tid, pos_read) && withInAfterMain.pos_within_same_chr_region(mtid, mpos_read) && main_dir == read_forward && afte_dir == mate_forward){//for read pairs that support the reference
			read_num_support_REF++;
			total_isize_REF += isize;
			suppor_int += 2;
		}
		if(suppor_int >=3)
			read_num_support_BOTH++;

		//if(suppor_int == 1){
		bam_load_region.st_pos = main_pos - max_isize;
		bam_load_region.ed_pos = main_pos + normal_read_length + normal_read_length;

		if(print_log && true) {
			fprintf(stderr, "[tid %d, pos_read_bg %d pos_read_end %d dir: %d] , [mtid %d, mpos_bg %d mpos_ed %d dir %d ] the read support: \t", tid, pos_read_bg, pos_read_ed , read_forward, mtid, mpos, mpos + br->core.l_qseq, mate_forward);
			switch(suppor_int){
			case 0: fprintf(stderr, "None\n"); break;
			case 1: fprintf(stderr, "BND\n"); break;
			case 2: fprintf(stderr, "REF\n"); break;
			case 3: fprintf(stderr, "BOTH\n"); break;
			default: fprintf(stderr, "UNKNOWN\n"); break;
			}
		}
	}
	GT_result[0] = read_num_support_BND;
	GT_result[1] = read_num_support_REF;
	GT_result[2] = read_num_support_BOTH;

	return true;
}


//return: true when successful; false when fail
bool SveHandler::assembly_and_get_breakpoints_TRA(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am,
		int main_tid, int &main_pos, bool read_is_before_breakpoint_in_main, bool main_ref_is_forward,
		int supp_tid, int &supp_pos, bool read_is_before_breakpoint_in_supp, bool supp_ref_is_forward){
	BS.clear();
	if(print_log){
		fprintf(stderr, "main_tid %d, main_pos %d, read_is_before_breakpoint_in_main %d, main_ref_is_forward  %d\n", main_tid, main_pos, read_is_before_breakpoint_in_main, main_ref_is_forward);
		fprintf(stderr, "supp_tid %d, supp_pos %d, read_is_before_breakpoint_in_supp %d, supp_ref_is_forward  %d\n", supp_tid, supp_pos, read_is_before_breakpoint_in_supp, supp_ref_is_forward);
	}
	FILE*log_output = stderr;
	//collecting reads MIAN

	int normal_read_length = sig_para.MaxReadLen;
	int edge_len = 300;
	RefRegion bam_load_region_main(main_tid,main_pos - edge_len,main_pos + edge_len);
	RefRegion bam_load_region_supp(supp_tid,supp_pos - edge_len,supp_pos + edge_len);
	bool region_overlap = bam_load_region_main.region_overlap(bam_load_region_supp);
	if(region_overlap)
		bam_load_region_main.Combine(bam_load_region_supp, true);

	for(int is_main = 1; is_main >= 0; is_main--){
		R_region cur_r;
		if(is_main)
			bam_load_region_main.toR_region(cur_r);
		else{
			if(region_overlap) continue;
			bam_load_region_supp.toR_region(cur_r);
		}
		int ASS_bg_pos = cur_r.st_pos - normal_read_length;
		int ASS_ed_pos = cur_r.ed_pos + normal_read_length;
		Bam_file *bam_f = &NGS_read.file;
		resetRegion_ID(bam_f, &cur_r);
		while (bam_next(bam_f)) {
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
			if(bam_is_duplicate(br))       continue;
			if(br->core.qual < 5)			continue;
			if(br->core.pos > ASS_bg_pos && br->core.pos < ASS_ed_pos){
				BS.storeReadCore(br);
				if(BS.reads.size() > 1000){
					BS.clear();
					break;
				}
			}
		}
	}

	BS.run_assembly(am);
	//collect the reference strings:

	std::vector<uint8_t> combined_ref_string;
	int accurate_BP_search_region_size = 250;
	//load main
	int main_load_st;
	int main_load_ed;
	int main_load_len;
	std::vector<uint8_t> main_ref_string;

	//load SUPP
	int supp_load_st;
	int supp_load_ed;
	int supp_load_len;
	std::vector<uint8_t> supp_ref_string;
	for(int is_main = 1; is_main >= 0; is_main--){
		int cur_tid;
		int cur_pos;
		int *cur_load_st;
		int *cur_load_ed;
		int *cur_load_len;
		bool read_is_before_breakpoint_in_cur_region;
		if(is_main){
			cur_tid = main_tid;
			cur_pos = main_pos;
			cur_load_st = &main_load_st;
			cur_load_ed = &main_load_ed;
			cur_load_len = &main_load_len;
			read_is_before_breakpoint_in_cur_region = read_is_before_breakpoint_in_main;
		}else{
			cur_tid = supp_tid;
			cur_pos = supp_pos;
			cur_load_st = &supp_load_st;
			cur_load_ed = &supp_load_ed;
			cur_load_len = &supp_load_len;
			read_is_before_breakpoint_in_cur_region = read_is_before_breakpoint_in_supp;
		}

		if(read_is_before_breakpoint_in_cur_region){
			*cur_load_st = cur_pos - accurate_BP_search_region_size;
			*cur_load_ed = cur_pos + accurate_BP_search_region_size;
		}else{
			*cur_load_st = cur_pos - accurate_BP_search_region_size;
			*cur_load_ed = cur_pos + accurate_BP_search_region_size;
		}

		*cur_load_st = MAX(10, *cur_load_st);

		char * cur_load_ref = ref->load_ref_by_region(cur_tid, *cur_load_st, *cur_load_ed, cur_load_len);
		uint8_t *cur_ref_string_p = NULL;
		bool cur_ref_is_forward;

		if(is_main){
			main_ref_string.resize(*cur_load_len);
			cur_ref_is_forward = main_ref_is_forward;
			cur_ref_string_p = &(main_ref_string[0]);
		}else{
			supp_ref_string.resize(*cur_load_len);
			cur_ref_is_forward = supp_ref_is_forward;
			cur_ref_string_p = &(supp_ref_string[0]);
		}

		if(cur_ref_is_forward){
			for(int i = 0; i < *cur_load_len; i++){
				switch(cur_load_ref[i]){
				case 'A': case 'a': cur_ref_string_p[i] = 0; break;
				case 'C': case 'c': cur_ref_string_p[i] = 1; break;
				case 'G': case 'g': cur_ref_string_p[i] = 2; break;
				case 'T': case 't': cur_ref_string_p[i] = 3; break;
				case 'n': case 'N': cur_ref_string_p[i] = 4; break;
				}
			}
		}else{
			for(int i = 0; i < *cur_load_len; i++){
				switch( cur_load_ref[*cur_load_len - i - 1]){
				case 'A': case 'a': cur_ref_string_p[i] = 3; break;
				case 'C': case 'c': cur_ref_string_p[i] = 2; break;
				case 'G': case 'g': cur_ref_string_p[i] = 1; break;
				case 'T': case 't': cur_ref_string_p[i] = 0; break;
				case 'n': case 'N': cur_ref_string_p[i] = 4; break;
				}
			}
		}
		free(cur_load_ref);
	}

	int BP1_region_st;
	int BP1_region_ed;

	int BP2_region_st;
	int BP2_region_ed;

	//combine the two strings
	//IF read_is_before_breakpoint_in_main is true:
		//the MAIN part is in the HEAD of combined string, and SUPP is in the TAIL
	//ELSE
		//the SUPP part is in the HEAD of combined string, and MAIN is in the TAIL
	bool main_part_is_in_head = read_is_before_breakpoint_in_main;
	bool successfully_adjusted_breakpoints = false;
	{
		combined_ref_string.clear();
		BP1_region_st = combined_ref_string.size();
		if(main_part_is_in_head)	std::swap(main_ref_string, combined_ref_string);
		else									std::swap(supp_ref_string, combined_ref_string);
		BP1_region_ed = combined_ref_string.size();
		combined_ref_string.emplace_back(4);
		BP2_region_st = combined_ref_string.size();
		if(main_part_is_in_head)	combined_ref_string.insert(combined_ref_string.end(), supp_ref_string.begin(), supp_ref_string.end());
		else									combined_ref_string.insert(combined_ref_string.end(), main_ref_string.begin(), main_ref_string.end());
		BP2_region_ed = combined_ref_string.size();
	}

	ca_bnd.setRef(&(combined_ref_string[0]), combined_ref_string.size(), 0,0);

	//show_reference strings
	if(print_log){
		fprintf(stderr, "combined_ref is: \n");
		for(uint i = 0; i < combined_ref_string.size(); i++)
			fprintf(stderr, "%c",  "ACGTNN"[combined_ref_string[i]]);
		fprintf(stderr, "\n");
	}
	if(print_log)
		BS.show_all_reads();
	//running the realignment:
	{
		//the result:
		AssemblyContig * final_contig = NULL; int final_BP1 = 0; int final_BP2 = 0;
		for(AssemblyContig &contig:BS.contigs){
			if(print_log){
				fprintf(stderr, "New contig\n%s\n", contig.seq.c_str());
			}
			if(print_log){
				contig.remove_read_set.clear();
				for (auto &ca : contig.actions)
					if(ca.read_ID < (int)BS.reads.size())
						ca.set_read_pos( BS.reads[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.kmerLength, contig.remove_read_set, BS.read_list[ca.read_ID].read_in_ref_offset);
				std::sort(contig.actions.begin(), contig.actions.end(), AssemblyReadAction::cmp_by_position);
				for (auto &ca : contig.actions){
					if(ca.read_ID < (int)BS.reads.size())// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
					{
						int contig_pos = -contig.ass_begin_offset_in_contig + ca.position_in_contig;
						if(ca.isAdd){
							for(int i = 0; i < contig_pos - ca.position_read ;i++)fprintf(stderr, " ");
							fprintf(stderr, "%s\t", BS.reads[ca.read_ID].c_str());
						}
						BS.read_list[ca.read_ID].print();
						fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
						ca.print(stderr, true);
					}
				}
			}

			//contig to bin char:
			//alignment to reference:contig.seq.c_str()
			store_bin_contig(contig.seq, bin_contig);
			ca_bnd.align_non_splice(&(bin_contig[0]), contig.seq.size(), 0, ca_bnd.get_tlen(), 400);
			int contig_st_in_ref = ca_bnd.adjustCIGAR();
			if(print_log) ca_bnd.printf_alignment_detail(stderr, contig_st_in_ref, NULL, 0);

			uint32_t* bam_cigar;
			int cigar_len = ca_bnd.get_cigar(&bam_cigar);
		    int output_index = contig_st_in_ref;
		    int seq_i = 0;
		    //for(int i = 0; i < contig_st_in_ref; i++) {fprintf(log_output, "-");  output_index++;}
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
					case 0: case 3: case 4: output_index += cigar_len; seq_i += cigar_len; break;//match
					case 1: seq_i += cigar_len; break;//insertion
					case 2://deletion:
					{
						int BP1 = output_index - 1; int BP2 = output_index + cigar_len + 1;
						//condition:
						if(BP1 > BP1_region_st && BP1 < BP1_region_ed && BP2 > BP2_region_st && BP2 < BP2_region_ed ){
							final_contig = &(contig);
							final_BP1 = BP1; final_BP2 = BP2;
						}
					}
					output_index += cigar_len;
					break;
					default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len); break;
				}
			}
		}
		if(true){
			if(final_contig != NULL)	fprintf(log_output, "\n\nassembly_and_get_breakpoints_BND: final_BP1 %d , final_BP2 %d \n", final_BP1, final_BP2);
			else						fprintf(log_output, "\n\nassembly_and_get_breakpoints_BND: BND no adjust results\n");
		}

		if(final_contig != NULL){
			//convert the BPs to true bps
			//IF read_is_before_breakpoint_in_main is true:
				//the MAIN part is in the HEAD of combined string, and SUPP is in the TAIL
			//ELSE
				//the SUPP part is in the HEAD of combined string, and MAIN is in the TAIL
			final_BP2 -= BP2_region_st;//align the BP2 to base-0
			int final_BP_main = (main_part_is_in_head)?final_BP1:final_BP2;
			int final_BP_supp = (main_part_is_in_head)?final_BP2:final_BP1;
			fprintf(log_output, "\n BND position is adjust by assembly, before is main_pos %d supp_pos %d\n", main_pos, supp_pos);
			main_pos = (main_ref_is_forward)?(main_load_st + final_BP_main):(main_load_ed - final_BP_main);
			supp_pos = (supp_ref_is_forward)?(supp_load_st + final_BP_supp):(supp_load_ed - final_BP_supp);
			fprintf(log_output, "\n After is main_pos %d supp_pos %d\n", main_pos, supp_pos);
			successfully_adjusted_breakpoints = true;
		}
		else{

		}
	}
	return successfully_adjusted_breakpoints;
}


//return: true when successful; false when fail
int SveHandler::assembly_and_get_breakpoints_INV(bool print_log, BND_ASS_Block &BS, MainAssemblyHandler *am,
		RefRegion &main_region, RefRegion &supp_region, int &main_BP, int &supp_BP){
	BS.clear();
	if(print_log){
		fprintf(stderr, "main region is ");	main_region.print(stderr);
		fprintf(stderr, "supp region is ");	supp_region.print(stderr);
		fprintf(stderr, "\n");
	}
	FILE*log_output = stderr;
	//collecting reads MIAN
	int normal_read_length = sig_para.MaxReadLen;
	bool region_overlap = main_region.region_overlap(supp_region);
	if(region_overlap)
		main_region.Combine(supp_region, true);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LOAD READs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	for(int is_main = 1; is_main >= 0; is_main--){
		R_region cur_r;
		if(is_main)									main_region.toR_region(cur_r);
		else{			if(region_overlap) continue;supp_region.toR_region(cur_r); }
		int ASS_bg_pos = cur_r.st_pos - normal_read_length;
		int ASS_ed_pos = cur_r.ed_pos + normal_read_length;
		Bam_file *bam_f = &NGS_read.file;
		resetRegion_ID(bam_f, &cur_r);
		while (bam_next(bam_f)) {
			bam1_t *br = &(bam_f->_brec);
			if(bam_is_secondary(br))		continue;
			if(bam_is_supplementary(br))	continue;
		    if(bam_is_duplicate(br))       continue;
			if(br->core.pos > ASS_bg_pos && br->core.pos < ASS_ed_pos){
				BS.storeReadCore(br);
				if(BS.reads.size() > 1000)	break;
			}
		}
	}
	if(false && print_log)
		BS.show_all_reads();
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ASS READs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	BS.run_assembly(am);
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GET REFERENCE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//collect the reference strings:
	REF_COMBINE ref_combine[4];
	int accurate_BP_search_region_size = 400;
	//load main
	int ref_load_st[2];
	std::vector<uint8_t> main_ref_string[2];
	for(int mode = 0; mode < 2; mode++){
		int cur_tid = (mode == 0)?main_region.chr_ID:supp_region.chr_ID;
		int cur_pos = (mode == 0)?main_region.getMiddle():supp_region.getMiddle();
		int cur_load_st = cur_pos - accurate_BP_search_region_size; cur_load_st = MAX(10, cur_load_st);
		int cur_load_ed = cur_pos + accurate_BP_search_region_size;
		int cur_load_len;
		char * cur_load_ref = ref->load_ref_by_region(cur_tid, cur_load_st, cur_load_ed, &cur_load_len);
		ref_load_char_2_bin((mode == 0)?(true):(false), cur_load_len, main_ref_string[mode], cur_load_ref);
		ref_load_st[mode] = cur_load_st;
		free(cur_load_ref);
	}
	//combine the two strings
	//IF read_is_before_breakpoint_in_main is true:
		//the MAIN part is in the HEAD of combined string, and SUPP is in the TAIL
	//ELSE
		//the SUPP part is in the HEAD of combined string, and MAIN is in the TAIL
	{
		ref_combine[0].store_ref(main_ref_string[0], true, ref_load_st[0], main_ref_string[1],false, ref_load_st[1], false);
		ref_combine[1].store_ref(main_ref_string[0], true, ref_load_st[0], main_ref_string[1],false, ref_load_st[1], true);
		ref_combine[2].store_ref(main_ref_string[1], false, ref_load_st[1], main_ref_string[0],true, ref_load_st[0], false);
		ref_combine[3].store_ref(main_ref_string[1], false, ref_load_st[1], main_ref_string[0],true, ref_load_st[0], true);
	}

	//show_reference strings
	if(print_log){
		ref_combine[0].print_ref();
		ref_combine[1].print_ref();
		ref_combine[2].print_ref();
		ref_combine[3].print_ref();
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Contig re-alignmet~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//running the realignment:
	//the result:
	for(AssemblyContig &contig:BS.contigs){
		//log data
		if(false && print_log){
			fprintf(stderr, "New contig\n%s\n", contig.seq.c_str());
			contig.remove_read_set.clear();
			for (auto &ca : contig.actions)
				if(ca.read_ID < (int)BS.reads.size())
					ca.set_read_pos( BS.reads[ca.read_ID], contig.seq, contig.ass_begin_offset_in_contig, contig.kmerLength, contig.remove_read_set, BS.read_list[ca.read_ID].read_in_ref_offset);
			std::sort(contig.actions.begin(), contig.actions.end(), AssemblyReadAction::cmp_by_position);
			for (auto &ca : contig.actions){
				if(ca.read_ID < (int)BS.reads.size())// && (remove_read_set.find(ca.read_ID) == remove_read_set.end()))
				{
					int contig_pos = -contig.ass_begin_offset_in_contig + ca.position_in_contig;
					if(ca.isAdd){
						for(int i = 0; i < contig_pos - ca.position_read ;i++)fprintf(stderr, " ");
						fprintf(stderr, "%s\t", BS.reads[ca.read_ID].c_str());
					}
					BS.read_list[ca.read_ID].print();
					fprintf(stderr, "POS CONTIG: %d \t POS read %d \t", contig_pos, ca.position_read);
					ca.print(stderr, true);
				}
			}
		}
		//contig to bin char:
		//alignment to reference:contig.seq.c_str()
		store_bin_contig(contig.seq, bin_contig);
		for(int i = 0; i < 4; i++){
			ca_bnd.setRef(&(ref_combine[i].s[0]), ref_combine[i].s.size(), 0,0);
			ca_bnd.align_non_splice(&(bin_contig[0]), contig.seq.size(), 0, ca_bnd.get_tlen(), 400);
			int contig_st_in_ref = ca_bnd.adjustCIGAR();
			if(print_log) ca_bnd.printf_alignment_detail(stderr, contig_st_in_ref, NULL, 0);
			uint32_t* bam_cigar;
			int cigar_len = ca_bnd.get_cigar(&bam_cigar);
			int output_index = contig_st_in_ref;
			int seq_i = 0;
			//for(int i = 0; i < contig_st_in_ref; i++) {fprintf(log_output, "-");  output_index++;}
			for(int cigar_ID = 0;cigar_ID < cigar_len; cigar_ID++){
				int cigar_len = bam_cigar[cigar_ID] >> BAM_CIGAR_SHIFT;
				int type = bam_cigar[cigar_ID] & BAM_CIGAR_MASK;
				switch(type){
					case 0: case 3: case 4: output_index += cigar_len; seq_i += cigar_len; break;//match
					case 1: seq_i += cigar_len; break;//insertion
					case 2://deletion:
					ref_combine[i].check_and_store(output_index, output_index + cigar_len);
					if(region_overlap && cigar_ID != 0){
						int pre_cigar_len = bam_cigar[cigar_ID - 1] >> BAM_CIGAR_SHIFT;
						int pre_type = bam_cigar[cigar_ID - 1] & BAM_CIGAR_MASK;
						if(pre_type == 1 && pre_cigar_len == cigar_len){
							ref_combine[i].chece_and_store_M2(output_index, output_index + cigar_len);
						}
					}
					output_index += cigar_len;
					break;
					default: fprintf(log_output, "ERROR CIGAR  %d %d ", type, cigar_len); break;
				}
			}
		}
	}
	if(true){
		//SUMMARY
		for(int i = 0; i < 4; i++){
			fprintf(stderr, "I: %d \t", i);
			ref_combine[i].print_results();
		}
	}

	//check the results:
	int signal_number = 0;
	main_BP = 0;
	supp_BP = 0;
	int length_inv1 = MAX(ref_combine[0].get_len(), ref_combine[1].get_len());
	int length_inv2 = MAX(ref_combine[2].get_len(), ref_combine[3].get_len());

	if(length_inv1 == 0 || length_inv2 == 0 || length_inv1 > 1.5*length_inv2 || length_inv2 > 1.5*length_inv1){
		fprintf(log_output, "SKIP for UNBALANCE LEN(M1) \n");
		return 0;
	}
	int length_diff = ABS_U(length_inv1, length_inv2);
	if(length_diff > 150){
		fprintf(log_output, "SKIP for UNBALANCE LEN(M2) \n");
		return 0;
	}

	//SUMMARY
	for(int i = 0; i < 4; i++){
		if(ref_combine[i].with_contig_supp){
			main_BP += ref_combine[i].main_pos;
			supp_BP += ref_combine[i].supp_pos;
			signal_number ++;
		}
	}
	if(signal_number != 0){
		main_BP/=signal_number;
		supp_BP/=signal_number;
	}

	return signal_number;
}

