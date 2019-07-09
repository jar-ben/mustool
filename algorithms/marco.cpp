#include "core/Master.h"
#include "core/misc.h"
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
		}else{
			MUS mus = shrink_formula(top);
			mark_MUS(mus);
			unex_unsat++;			
			if(useMatchmaker){
				//we first build a database of at least some MUSes so we do some statistics, 
				//and them use them inside Matchmaker
				chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
				if(duration > 30 || muses.size() > 100){
					recursive_rotation_delta(mus, original_top, 0);
				}
			}
			else if(useBackbone){
				backbone_mus_rotation(mus, original_top);
			}
			else if(useMixedHeuristic){
				mixed_mus_rotation(mus, original_top);
			}
		}
		top = explorer->get_unexplored(1, false);
	}
}

