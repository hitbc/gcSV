/*
 * Genotyping_handler.hpp
 *
 *  Created on: 2022-10-31
 *      Author: fenghe
 *      Function: handle all calculation of genotyping
 */

#ifndef SRC_SVCALLING_CORE_GT_CAL_HPP_
#define SRC_SVCALLING_CORE_GT_CAL_HPP_

#include<vector>
extern "C"
{
#include "../clib/utils.h"
}

#define GT_CAL_GERMLINE_MODE 0
#define GT_CAL_SOMATIC_MODE 1

//var mode
#define GT_CAL_SHORT_VAR_MODE 0
#define GT_CAL_STRUCT_VAR_MODE 1

#define P_IMPOSSIBLE (-10000)
//each read is a signal
struct GT_signal{
	int haplotype_num;//total hap number
	int P_d_a[10];// the possibility when a read happen when a specific haplotype happen; take log value based 10, then *10

	GT_signal(int haplotype_num_){
		haplotype_num = haplotype_num_;
		for(int i = 0; i < haplotype_num; i++)
			P_d_a[i] = P_IMPOSSIBLE;
	}

	void set_poss_log(int hap_ID, int log_value){
		P_d_a[hap_ID] = log_value;
	}
};

class GT_handler{
	//basic
	std::vector<GT_signal> signal_list;
	//GT mode
	int GT_MODE;
	int hap_number;
	int P_a_G[10][10];//the possible of a haplotype when a specific Genotype happen
	int var_MODE;//"SNP/INDEL" or "SV"
	//P_G: the possible that a Genotype happening
	int genotype_number;
	int P_G[10];
	//GT string
	std::vector<std::string> GT_str;

	//P_D_G: the results
	int P_D_G[10];
	//final results
	int final_GT;
	int final_score;

	//Tips: when P1 and P2 were added,
	//log10(P1 + P2)*10 = MAX(log10(P1)*10, log10(P2)*10) + add_diff[|log10(P1)*10 - log10(P2)*10|]
	//let f(P)= log10(P)*10, the following is always true:
	//f(P1 + P2) = MAX(f(P1), f(P2)) + add_diff[|f(P1) - f(P2)|]
	int add_diff[20];

	void init_str(){
		GT_str.emplace_back("1/1"); //AA
		GT_str.emplace_back("0/1");	//AR
		GT_str.emplace_back("0/0");	//RR

		GT_str.emplace_back("0/0");	//N
		GT_str.emplace_back("1/1");	//S
	}

	void init_add_diff(){
		for(int i = 0; i < 20; i++){
			int a = 0;
			double P1 = pow(10, (double)a/10);
			int b = -i;
			double P2 = pow(10, (double)b/10);

			double P3 = P1+P2;
			double c = log10(P3)*10;
			int cur_add_diff = c + 0.5 - a;
			if(false)
				fprintf(stderr, "current add_diff at %d is %d\n",i, cur_add_diff);
			add_diff[i] = cur_add_diff;
		}
	}

	int get_P_d_G(int GT, GT_signal &cur_sig, int MIN_value){
		int P_d_G = P_IMPOSSIBLE;
		for(int i = 0; i < cur_sig.haplotype_num; i++){
			int P_d_a_mutiply_by_P_a_G = cur_sig.P_d_a[i] + P_a_G[i][GT];
			P_d_G = log_value_add(P_d_G, P_d_a_mutiply_by_P_a_G);
		}
		P_d_G = MAX(P_d_G, MIN_value);
		return P_d_G;
	}

	int get_P_D_G(int GT){
		int P_D_G = 0;
		for(GT_signal &cur_sig:signal_list){
			P_D_G += get_P_d_G(GT, cur_sig, -80);
		}
		return P_D_G;
	}

	//calculate possible value ADD
	//f(P1 + P2) = MAX(f(P1), f(P2)) + add_diff[|f(P1) - f(P2)|]
	int log_value_add(int log_v1, int log_v2){
		if(log_v1 < log_v2)
			std::swap(log_v1, log_v2);
		int cur_add_diff = (log_v1 - log_v2 >= 20)?0:add_diff[log_v1 - log_v2];
		return log_v1 + cur_add_diff;
	}


public:
	//init
	void init(int GT_MODE_, int var_MODE_){
		init_add_diff();
		init_str();
		GT_MODE = GT_MODE_;
		var_MODE = var_MODE_;

		//set hap_number
		if(GT_MODE == GT_CAL_GERMLINE_MODE){
			hap_number = 2; // A or R
			genotype_number = 3;
			if(var_MODE == GT_CAL_SHORT_VAR_MODE){
				P_G[0] = -33;//AA
				P_G[1] = -30;  //AR
				P_G[2] = 0;   //RR
			}
			else if(var_MODE == GT_CAL_STRUCT_VAR_MODE){
				P_G[0] = -53;//XX
				P_G[1] = -50;  //XR
				P_G[2] = 0;   //RR
			}else
				xassert(0, "Unknown VAR mode");
		}
		else if(GT_MODE == GT_CAL_SOMATIC_MODE){
			hap_number = 4; // A / R / N / S
			genotype_number = 5;

			if(var_MODE == GT_CAL_SHORT_VAR_MODE){
				P_G[0] = -33;//XX
				P_G[1] = -30;  //XR
				P_G[2] = 0;   //RR
				P_G[3] = -50;   //N
				P_G[4] = -70;   //S
			}
			else if(var_MODE == GT_CAL_STRUCT_VAR_MODE){
				P_G[0] = -53;//XX
				P_G[1] = -50;  //XR
				P_G[2] = 0;   //RR
				P_G[3] = -50;   //N
				P_G[4] = -70;   //S
			}else
				xassert(0, "Unknown VAR mode");
		}
		else
			xassert(0, "Unknown GT mode");

		//set a_G
		P_a_G[0][0] = 0;				//a=A && GT=AA
		P_a_G[0][1] = -3;				//a=A && GT=AR
		P_a_G[0][2] = P_IMPOSSIBLE;	//a=A && GT=RR
		P_a_G[0][3] = P_IMPOSSIBLE;	//a=A && GT=N
		P_a_G[0][4] = P_IMPOSSIBLE;	//a=A && GT=S

		P_a_G[1][0] = P_IMPOSSIBLE;	//a=R && GT=AA
		P_a_G[1][1] = -3;				//a=R && GT=AR
		P_a_G[1][2] = 0;				//a=R && GT=RR
		P_a_G[1][3] = P_IMPOSSIBLE;	//a=R && GT=N
		P_a_G[1][4] = P_IMPOSSIBLE;	//a=R && GT=S

		P_a_G[2][0] = P_IMPOSSIBLE;			//a=N && GT=AA
		P_a_G[2][1] = P_IMPOSSIBLE;			//a=N && GT=AR
		P_a_G[2][2] = P_IMPOSSIBLE;			//a=N && GT=RR
		P_a_G[2][3] = P_IMPOSSIBLE;			//a=N && GT=N
		P_a_G[2][4] = P_IMPOSSIBLE;//a=N && GT=S

		P_a_G[3][0] = P_IMPOSSIBLE;			//a=S && GT=AA
		P_a_G[3][1] = P_IMPOSSIBLE;			//a=S && GT=AR
		P_a_G[3][2] = P_IMPOSSIBLE;			//a=S && GT=RR
		P_a_G[3][3] = P_IMPOSSIBLE;			//a=S && GT=N
		P_a_G[3][4] = -3;			//a=S && GT=S

	}

	void clear(){
		signal_list.clear();
	}

	//basic
	void new_signal(){
		signal_list.emplace_back(hap_number);
	}

	//P_value = log(P)*10;
	void set_P_log_back(int hap_ID, int P_value_log){
		signal_list.back().P_d_a[hap_ID] = P_value_log;
	}

	void set_P_log_rid(int read_ID, int hap_ID, int P_value){
		signal_list[read_ID].P_d_a[hap_ID] = P_value;
	}

	void calculateGT(){
		int total_P_D_G = P_IMPOSSIBLE;
		int max_GT = 0; int max_P = P_IMPOSSIBLE;
		for(int i = 0; i < genotype_number; i++){
			P_D_G[i] = get_P_D_G(i);
			total_P_D_G = log_value_add(total_P_D_G, P_D_G[i]);
			if(max_P < P_D_G[i]){
				max_P = P_D_G[i];
				max_GT = i;
			}
		}

		int other_p = P_IMPOSSIBLE;
		for(int i = 0; i < genotype_number; i++){
			P_D_G[i] -= total_P_D_G;
			if(i != max_GT){
				other_p = log_value_add(other_p, P_D_G[i]);
			}
		}
		final_GT = max_GT;
		final_score = - other_p;
	}

	int get_score(){ return final_score; }
	int get_GT(){ return final_GT; }

	std::string &get_GT_str(){ return GT_str[final_GT]; }

};

struct GT_CAL_DEMO_10_10_10_10{
	int GT_MODE; int var_MODE;
	void test(){
		//test:
		GT_handler h;
		h.init(GT_CAL_GERMLINE_MODE, GT_CAL_SHORT_VAR_MODE);
		//clear
		h.clear();

		//r1-10, support REF
		for(int i = 0; i < 10; i++){
			h.new_signal(); h.set_P_log_back(0, -40); h.set_P_log_back(1, 0);
		}

		//r11-20, support ALT
		for(int i = 0; i < 2; i++){
			h.new_signal(); h.set_P_log_back(0, 0); h.set_P_log_back(1, -40);
		}

		//r21-30, support Both
		for(int i = 0; i < 0; i++){
			h.new_signal(); h.set_P_log_back(0, 0); h.set_P_log_back(1, 0);
		}

		//r31-40, support neither
		for(int i = 0; i < 100; i++){
			h.new_signal(); h.set_P_log_back(0, -40); h.set_P_log_back(1, -40);
		}

		h.calculateGT();

		fprintf(stderr, "h.get_GT() %s, h.get_score() %d\n", h.get_GT_str().c_str(), h.get_score());

	}
};

#endif /* SRC_SVCALLING_CORE_GT_CAL_HPP_ */
