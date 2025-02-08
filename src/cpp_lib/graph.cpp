/*
 * graph.cpp
 *
 *  Created on: 2021年3月13日
 *      Author: fenghe
 */
#include<stdlib.h>
#include<stdint.h>
#include "graph.hpp"
extern "C"{
#include "../clib/utils.h"
}

int UNI_SEED::cmp(const void *a , const void *b)
{
	UNI_SEED* seed1 = (UNI_SEED *)a;
	UNI_SEED* seed2 = (UNI_SEED *)b;

	if (seed1->ref_end > seed2->ref_end)
		return 1;
	else if (seed1->ref_end < seed2->ref_end)
		return -1;
	else
	{
		if (seed1->ref_begin > seed2->ref_begin)
			return 1;
		else if (seed1->ref_begin < seed2->ref_begin)
			return -1;
		else
			return 0;
	}
}

void  UNI_SEED::show(std::vector<UNI_SEED> &uniseed){
	for (uint32_t i = 0; i < uniseed.size(); ++i)
		fprintf(stderr,"id = %u, read_begin = %u, read_end = %u, cov = %d, seed_id = %u,ref_begin  = %u, ref_end = %u\n",i, uniseed[i].read_begin, uniseed[i].read_end, uniseed[i].cov, uniseed[i].seed_id, \
			uniseed[i].ref_begin, uniseed[i].ref_end);
}


#define PARAM 2
#define INTRON_PENALTY 2

#define MAX_REF_DIS 50
#define MAX_READ_DIS 50
#define MAX_REF_DIS_STR 400
#define MAX_READ_DIS_STR 400
#define MAX_SEARCH_STEP 40 // = max read length / search step length  = 400 / 5 = 80
#define MAX_SEARCH_STEP_STR 80 // = max read length / search step length  = 400 / 5 = 80
#define MAX_ABS_GAP 50
#define MAX_ABS_GAP_STR 20

void Graph_handler::process(std::vector<UNI_SEED>& vertexArr){

	phase_ID = 0;
	vertexArr_size = vertexArr.size();
	vertexArr_ = &vertexArr;
	if(vertexArr_size == 0)
		return;
	vertexArr_size = vertexArr.size();
	qsort(&(vertexArr[0]), vertexArr_size, sizeof(UNI_SEED), UNI_SEED::cmp);

	int max_ref_dis = (readIsSTR == true)?MAX_REF_DIS_STR:MAX_REF_DIS;
	int max_read_dis = (readIsSTR == true)?MAX_READ_DIS_STR:MAX_READ_DIS;
	uint max_search_step = (readIsSTR == true)?MAX_SEARCH_STEP_STR:MAX_SEARCH_STEP;
	uint max_gap = (readIsSTR == true)?MAX_ABS_GAP_STR:MAX_ABS_GAP;
	int search_step = MIN(vertexArr_size, max_search_step);
	bool non_ioslated_point = true;
	//show_uniseed(vertexArr, vertexNum);

	if(graph.size() < vertexArr_size){//malloc for graph
		graph.resize(vertexArr_size);
		dist_path.resize(vertexArr_size);
	}

	// re-inital
	for (uint32_t i = 0; i < vertexArr_size; ++i){
        graph[i].new_node(i, vertexArr[i].cov);
        dist_path[i].new_node(vertexArr[i].cov, -1);
	}

	for (uint32_t target_ID = 0; target_ID < vertexArr_size - 1; ++target_ID){
		uint32_t read_end = vertexArr[target_ID].read_end;
		uint32_t ref_end = vertexArr[target_ID].ref_end;
		uint32_t seed_id = vertexArr[target_ID].seed_id;
		uint32_t search_end = MIN(vertexArr_size, (target_ID + search_step));

		for (uint32_t try_ID = target_ID + 1; try_ID < search_end; ++try_ID)
		{
			if (vertexArr[try_ID].seed_id == seed_id){ continue; }//two mem generated from the same seed will not be connected.
			if (vertexArr[try_ID].ref_end == ref_end){ continue; }
			int32_t dis_ref = (int32_t)(vertexArr[try_ID].ref_begin - ref_end);
			if (dis_ref > max_ref_dis)  	break;
			int32_t dis_read = (int32_t)(vertexArr[try_ID].read_begin - read_end);
			if (dis_read > max_read_dis) 	continue;
			uint32_t abs_gap = ABS_U(dis_read, dis_ref);
			if(abs_gap > max_gap) continue;
			//int32_t gap = (int32_t)(dis_read - dis_ref);
			float penalty = (abs_gap == 0)?0: ((ABS(abs_gap) >> 3) + 3);// ABS(gap) * 0.125 + 3; open:4; ext:0.125
			uint32_t weight = 0;
			//get weight and penalty

			if (dis_read == dis_ref) //without INDEL
				weight = vertexArr[try_ID].cov - MAX(1 - dis_read, 0);
			 else if ((dis_read > 0 && dis_ref > 0))
				weight = vertexArr[try_ID].cov;
			else if((dis_read >= -5 && dis_read <= 0  && dis_ref >= -5))
				weight = vertexArr[try_ID].cov + MIN(dis_read , dis_ref);
			else
				continue;

			graph[try_ID].pre_edge.emplace_back(target_ID, weight, penalty);

			non_ioslated_point = false;

		}
	}

	if (!non_ioslated_point)
		dynamic_programming_path();

	return;
}

void Graph_handler::dynamic_programming_path(){
	// for all For each v∈V\S in Topological Sorting order do dilg(v)=max(u,v)∈E{dilg(u)+w(u, v)}
	for (uint32_t target = 0; target < vertexArr_size; ++target){
		if (graph[target].pre_edge.empty()) 	continue;
		float current_dist = 0;
		int32_t pre_node = -1;
		//travel the processer of vertex target
		for (ANode &arcnode : graph[target].pre_edge){
			int32_t weight = arcnode.weight;
			float penalty = arcnode.penalty;
			uint32_t adjvex = arcnode.adjvex;
			float temp = dist_path[adjvex].dist + weight - penalty;
			if (current_dist <= temp){
				current_dist = temp;
				pre_node = adjvex;
			}
		}
		dist_path[target].dist = current_dist;   // change the distance because there are overlaps between mems
		dist_path[target].pre_node = pre_node; //the front node

//		if (max_distance < dist_path[target].dist){
//			max_distance = dist_path[target].dist;
//			max_index = target;
//		}
	}
}
