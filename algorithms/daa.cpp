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
			MUS mus(bot, -1, muses.size(), count_ones(bot));
			muses.push_back(mus);
			mark_MUS(mus);
			unex_unsat++;
		}else{
			Formula mss = grow_formula(bot);
			block_down(mss);
			unex_sat++;			
		}
		bot = explorer->get_unexplored(0, false);
	}
}

