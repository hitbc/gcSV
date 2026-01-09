/*
 * graph.hpp
 *
 *  Created on: 2021-3-12
 *      Author: fenghe
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include<vector>
#include<stdlib.h>
#include<string.h>
extern "C"{
#include "../clib/utils.h"
}

///*************************************************************************
//	> File Name: graph.c
//	> Author:
//	> Mail:
// ************************************************************************/
//
//#include <stdio.h>
//#include <stdlib.h>

#include <time.h>

//#include <assert.h>
////#include <malloc.h>
//#include <string.h>
//#include <math.h>
//#include "binarys_qsort.h"
//
//
//#define REFLEN 0Xffffffff
//
//double DAG_time = 0;
//double findpath_time = 0;
//double qsort_time = 0;

struct UNI_SEED
{
	uint32_t read_begin;
	uint32_t read_end;
	uint32_t seed_id; //record the first seed_id of those mems which can be merged
	uint32_t ref_begin;
	uint32_t ref_end;
	uint32_t cov;

	void print(){ fprintf(stderr, "UNI_SEED: read_begin:[%d], read_end:[%d], seed_id:[%d], ref_begin:[%d], ref_end:[%d], cov:[%d], diff:[%d]\n",
				read_begin, read_end, seed_id, ref_begin, ref_end, cov, ref_begin - read_begin);}
	static int cmp(const void *a , const void *b);
	static int cmp_diff(const void *a , const void *b){
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

	static void show(std::vector<UNI_SEED> &uniseed);
};


struct PATH_t
{
	void new_node(float dist_, int32_t pre_node_){
		dist = dist_;
		pre_node = pre_node_;
		already_used = false;
	}
	float dist;
	int32_t pre_node;
	uint8_t already_used;

	static int cmp(const void *a , const void *b){
		PATH_t* seed1 = (PATH_t *)a;
		PATH_t* seed2 = (PATH_t *)b;

		if(seed1->dist == seed2->dist)
			return 0;
		if (seed1->dist < seed2->dist)
			return 1;
		return -1;
	}
    void print(char endl){ fprintf(stderr, "ANode: dist:[%f], pre_node:[%d]%c", dist, pre_node, endl); }

};

struct Graph_handler{
	struct ANode{
		ANode(uint32_t adjvex_, int weight_, float penalty_):
			adjvex(adjvex_), weight(weight_),penalty(penalty_){}
		uint32_t adjvex;
	    int weight;
	    float penalty;

	    void print(char endl){
	    	fprintf(stderr, "ANode: adjvex:[%d] weight:[%d], penalty:[%f]%c", adjvex, weight, penalty, endl);
	    }
	};

	struct VNode{
		void new_node(uint32_t head_vertex_, int weight_){
			head_vertex = head_vertex_;
	        weight = weight_;
			pre_edge.clear();
		}
		uint32_t head_vertex;
	    int weight;
		std::vector<ANode> pre_edge;

	    void print(char endl){
	    	fprintf(stderr, "VNode: head_vertex:[%d], weight:[%d]\t", head_vertex, weight);
	    	for(auto & a :pre_edge)
	    		a.print('\t');
	    	fprintf(stderr, "%c", endl);
	    }
	};

	//information from reads
	uint64_t vertexArr_size;
	//information kept in graph
	std::vector<VNode> graph;//graph node
	std::vector<PATH_t> dist_path;//graph node
	//std::vector<uint8_t> out_degree;//graph node
	std::vector<UNI_SEED>* vertexArr_;
	float max_distance = 0; // store the best path
	uint32_t max_index = 0;
	int phase_ID = 0;
	bool readIsSTR = false;

	//same top max distance count
	int same_top_max_distance_count = 0;
	std::vector<int> same_top_max_distance_id_list;


    void print(){
    	fprintf(stderr, "Graph_handler\n");
    	for(uint32_t i = 0; i < vertexArr_size; i++)	graph[i].print('\n');
    	fprintf(stderr, "\n\n");
    	for(uint32_t i = 0; i < vertexArr_size; i++)	dist_path[i].print('\n');
    	fprintf(stderr, "\n");
    }

    void sort_print(){
    	//output the best
    	if(vertexArr_size == 0){
    		fprintf(stderr, "Graph_handle phase: [%d] :max_distance[ NULL ]@max_index[ NULL ] \n", phase_ID++);
    		return;
    	}
    	max_index = 0; max_distance = 0;
       	for(int i = vertexArr_size - 1; i >= 0; i--){
       		if(dist_path[i].already_used)
       			continue;
    		if (max_distance < dist_path[i].dist){
    			max_distance = dist_path[i].dist;
    			max_index = i;
    		}
       	}

       int read_begin = 0;
       int ref_begin = 0;
       	for(int i = max_index; i != -1; i = dist_path[i].pre_node){
			//fprintf(stderr, "ID: [%d]\t", i);
			//vertexArr_[0][i].print();
			read_begin = vertexArr_[0][i].read_begin;
			ref_begin = vertexArr_[0][i].ref_begin;
			dist_path[i].already_used = true;
		}

       	fprintf(stderr, "Graph_handle phase: [%d]:max_distance[%f]@max_index[%d], read_begin[%d], ref_begin[%d] \n", phase_ID++, max_distance, max_index, read_begin, ref_begin);

    	//for(uint32_t i = 0; i < vertexArr_size; i++)	graph[i].print('\n');



    	//qsort(&(dist_path[0]), vertexArr_size, sizeof(PATH_t), PATH_t::cmp);
    	//for(uint32_t i = 0; i < MIN(vertexArr_size, 5); i++)	dist_path[i].print('\n');

    }

	void process(std::vector<UNI_SEED>& vertexArr);
	void dynamic_programming_path();

};


#endif /* GRAPH_HPP_ */
