/*
 * sve.cpp
 *
 *  Created on: 2020年4月30日
 *      Author: fenghe
 */

#include "../SVcalling_core/sve.hpp"

#include <algorithm>

//combined signals to begin/end SV type
void sve_combine_STEP1_SH(SVE_L & l, int min_score, bool is_union){
	std::sort(l.begin(), l.end(), SVE::cmp_by_position);
	int store_index = 0;
	auto sve_ed = l.end();//for each sve
	for(auto sve = l.begin(); sve < sve_ed; sve++)	{
		if(sve->info.is_solid == RST::UNKNOWN) continue;//already combined sve
		for(auto sve_try = sve + 1; sve_try < sve_ed && sve->r1.region_overlap(sve_try->r1); sve_try++){//use sve as main, than try to combine
			if(sve->info.is_solid != sve_try->info.is_solid) continue;//condition 1
			if(!sve->r2.region_overlap(sve_try->r2)) 		 continue;//condition 3
			sve->r1.Combine(sve_try->r1, is_union);
			sve->r2.Combine(sve_try->r2, is_union);
			sve->info.combine_score(sve_try->info);
			sve_try->info.is_solid = RST::UNKNOWN;
		}
		if(sve->info.getScore() < min_score && sve->info.getReadNum() < 2) continue;
		//store main evidence to list
		l[store_index++] = sve[0];
	}
	l.resize(store_index);
}
