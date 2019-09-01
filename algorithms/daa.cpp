#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

void Master::daa_base(){
	Formula bot = explorer->get_unexplored(0, false);
	Formula original_top;
	while(!bot.empty()){
		if(!is_valid(bot, false, true)){
			if(useMatchmaker || useBackbone || useMixedHeuristic){
				original_top = explorer->get_top_unexplored(bot);
			}
			MUS mus(bot, -1, muses.size(), count_ones(bot));
			muses.push_back(mus);
			mark_MUS(mus);
			unex_unsat++;
			if(useMatchmaker){
				//we first build a database of at least some MUSes so we do some statistics, 
				//and them use them inside Matchmaker
				chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
				if(duration > 30){// || muses.size() > 100){
					recursive_rotation_delta(mus, original_top, 0);
				}
			}
			else if(useBackbone){
				backbone_mus_rotation(mus, original_top);
			}
			else if(useMixedHeuristic){
				mixed_mus_rotation(mus, original_top);
			}
		}else{
			Formula mss = grow_formula(bot);
			block_down(mss);
			unex_sat++;			
		}
		bot = explorer->get_unexplored(0, false);
	}
}

