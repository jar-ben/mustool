#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


//remus
//if subset[i] = True then i is fixed to be presented in the current search space
void Master::duremus(Formula subset, Formula conflicts, int depth){
	//cout << "banananas: " << count_ones(conflicts) << " " << depth << " " << count_ones(subset) << endl;
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
	bool bot_found = false;
	while(true){
		iteration++;
		if(!bot_found) bot = explorer->get_bot_unexplored(assumptions);
		bot_found = false;
		if(bot.empty())	{ 
			current_depth--; 
			return;
		}
		Formula original_bot = bot;		
		if(is_valid(bot, true, true)){
			streak = 0;
			MSS mss = grow_formula(bot, conflicts);
			mark_MSS(mss);
			if(depth > depthMUS) continue;
			Formula m = mss.bool_mss;
			reduce_mss(original_bot, m);
			duremus(m, conflicts, depth + 1);
		}
		else{	//top is necessarily an MUS of subset
			streak++;
			//prevent getting stuck in a small subset
			if(streak > 10 && depth > 0){ current_depth--; return; }
			block_up(bot);

			vector<int> conf_all;
			for(int i = 0; i < dimension; i++)
				if(!subset[i] && original_bot[i]){
					conf_all.push_back(i);
				}
			if(conf_all.size() == 1){
				conflicts[conf_all[0]] = true;
				continue;
			}	
			//if(depth > depthMUS) continue;		
			for(auto conf: conf_all){
				Formula rec_subset = original_bot;
				rec_subset[conf] = false;
				Formula rec_conflicts = conflicts;
				rec_conflicts[conf] = true;
				duremus(rec_subset, rec_conflicts, depth + 1);
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
	int target_size = (1 + (1 - dim_reduction)) * bot_size;
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
