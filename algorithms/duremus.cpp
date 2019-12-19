#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


void Master::duremus_booster(int limit){
	int i = 0;
	while(i < limit){
		Formula top = explorer->get_unexplored(1, false);
		if(top.empty()) return;
		if(is_valid(top, true, false)){
			msses.push_back(MSS(top, -1, msses.size(), count_ones(top)));
			mark_MSS(msses.back());
			guessed++;
		}else{
			block_up(top);
			i++;
		}
	}
}

//remus
//if subset[i] = True then i is fixed to be presented in the current search space
void Master::duremus(Formula subset, Formula conflicts, int depth){
	//cout << "banananas: " << count_ones(conflicts) << " " << depth << " " << count_ones(subset) << endl;
	int subsetSize = count_ones(subset);
	current_depth++;	
	std::vector<int> assumptions;
	for(int i = 0; i < dimension; i++){
		if(subset[i]){
			assumptions.push_back(i + 1);			
		}
	}

	Formula bot;	
	int streak = 0;	
	int iteration = 0;	
	while(true){
		bit++;
		if(bit % 20 == 0){
			duremus_booster(10);
		}
		iteration++;
		if(bit % 10 == 0)
			bot = explorer->get_top_unexplored(assumptions);
		else
			bot = explorer->get_bot_unexplored(assumptions);
		if(bot.empty())	{ 
			current_depth--; 
			return;
		}
		Formula original_bot = bot;		
		if(is_valid(bot, true, true)){
			streak = 0;
			vector<MSS> grown;
			if(satSolver->grow_alg == "mcsls"){
				grown = grow_formulas(bot, conflicts);
			}else{
				grown.push_back(grow_formula(bot, conflicts));
			}
			for(auto m: grown)
				mark_MSS(m);
			if(subsetSize > 0.05 * depthMUS * dimension) continue;	
			//if(depth > depthMUS) continue;
			Formula m = msses.back().bool_mss;
			reduce_mss(original_bot, m);
			duremus(m, conflicts, depth + 1);
		}
		else{	//top is necessarily an MUS of subset
			streak++;
			//prevent getting stuck in a small subset
			block_up(bot);
			if(streak > 10 && depth > 0){ current_depth--; return; }
			//continue;
			vector<int> conf_all;
			for(int i = 0; i < dimension; i++)
				if(!subset[i] && original_bot[i]){
					conf_all.push_back(i);
				}
			if(conf_all.size() == 1)
				continue;
			for(auto conf: conf_all){
				Formula rec_subset = original_bot;
				rec_subset[conf] = false;
				duremus(rec_subset, conflicts, depth + 1);
			}
		}
	}
	current_depth--;
}
	

// helper funcion for the ReMUS algorithm
// modifies mus into a set S such that mus /subseteq S /subseteq top
void Master::reduce_mss(Formula &bot, Formula &mss){
	int ones = count_ones(mss);
	int bot_size = count_ones(bot);
	//int target_size = (1 + (1 - dim_reduction)) * bot_size;
	int target_size = bot_size + 0.05 * dimension;
	if(target_size - bot_size < 10)
		target_size += 10;

	for(int i = 0; i < dimension; i++){
		if(!bot[i] && mss[i]){
			mss[i] = false;
			ones--;
			if(ones <= target_size) 
				break;
		}			
	}
}
