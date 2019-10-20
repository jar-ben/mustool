#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

// Binary search
pair<Formula, Formula> Master::local_mus(Formula bot, Formula top, int diff){
	if(diff == 1)
		return make_pair(bot, top);
	int ones = 0;
	Formula pivot = bot;
	for(int i = 0; i < dimension; i++)
		if(!pivot[i] && top[i]){
			pivot[i] = true;
			if(++ones >= diff/2)
				break;
		}
	if(is_valid(pivot))
		return local_mus(pivot, top, diff - ones);	
	else
		return local_mus(bot, pivot, ones);
}

// The core of the TOME algorithm
void Master::find_all_muses_tome(){
	Formula bot, top, mss, fmus, original_top;
	vector<Formula> path;
	int diff;
	while(true){
		bot = explorer->get_unexplored(0, false);
		if(bot.empty()) break;
		top = explorer->get_top_unexplored(bot);
		original_top = top;
		if(is_valid(top)){
			block_down(top);
			continue;
		}
		if(!is_valid(bot)){
			MUS mus = MUS(bot, -1, muses.size()); //-1 duration means skipped shrink
			muses.push_back(mus);
			mark_MUS(mus);
			continue;
		}	
		diff =  std::count(top.begin(), top.end(), true) - std::count(bot.begin(), bot.end(), true);
		auto locals = local_mus(bot, top, diff);
		fmus = locals.second;
		mss = locals.first;
		if(domain == "sat" || domain == "smt")
			is_valid(fmus, true, true); //get unsat core of fmus before shrinking
		MUS mus = shrink_formula(fmus); 		
		mark_MUS(mus);
		block_down(mss);
	}
}

