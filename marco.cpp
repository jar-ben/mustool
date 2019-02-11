#include "Master.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

void Master::marco_base(){
	int crit_candidate = 0;
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
				int rot = backbone_mus_rotation(mus, original_top);
				cout << "recursively rotated: " << rot << endl;
			}
		}
		top = explorer->get_unexplored(1, false);
	}
}

