#include "Solver.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

bool Solver::use_backbone(){
	int total = 0;
	for(auto mus: muses){
		total += mus.int_mus.size();
	}
	float avg = total / float(muses.size());
	cout << "average size: " << avg << endl << endl << endl;
	return avg < 3000;
}

void Solver::marco_base(){
	int crit_candidate = 0;
	bool use = true;
	Formula top = explorer->get_unexplored(1, false);
	while(!top.empty()){
		Formula original_top = top;
		if(is_valid(top, true, true)){
			block_down(top);
			unex_sat++;
			if(model_rotation){
				Formula model = satSolver->get_model();
				int rotation_limit = 100;
				int sat_rotated = explorer->sat_rotation(original_top, model, satSolver->clauses, rotation_limit);
	//			cout << "sat rotated: " << sat_rotated << endl;
			}
		}else{
			MUS mus = shrink_formula(top);
			mark_MUS(mus);
			unex_unsat++;			
			if(variant == 31){
				chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
				if(duration > 30 || muses.size() > 100){
					int rot = recursive_rotation_delta(mus, original_top, 0);
					cout << "recursively rotated: " << rot << endl;
				}
			}
			if(variant == 32){
//				if(muses.size() == 20)
//					use = use_backbone();
				if(use){
					int rot = backbone_mus_rotation(mus, original_top);
					cout << "recursively rotated: " << rot << endl;
				}
			}
			if(variant == 33){
				int rot = backbone_dual_mus(mus, original_top);
				cout << "recursively rotated: " << rot << endl;
			}
			if(variant == 34){
				chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
				if(duration > 30 || muses.size() > 100){
					int rot = backbone_mus_rotation_all(mus, original_top);
					cout << "recursively rotated: " << rot << endl;
				}
			}
		}
		bool crit_attempt = false;
//		while(crit_candidate < dimension){
//			if(explorer->mus_intersection[crit_candidate++]){
//				top = Formula(dimension, true); top[crit_candidate-1] = false; crit_attempt = true; break;
//			}
//		}
		if(!crit_attempt)
			top = explorer->get_unexplored(1, false);
	}
}

