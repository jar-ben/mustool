#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


//remus
void Master::find_all_muses_duality_based_remus(Formula subset, Formula crits, int depth){
	current_depth++;	
	std::vector<int> assumptions;
	for(int i = 0; i < dimension; i++){
		if(!subset[i]){
			assumptions.push_back(-1 * i - 1);			
		}
	}

	Formula top;	
	int streak = 0;	
	int iteration = 0;	
	bool top_found = false;
	while(true){
		iteration++;
		if(!top_found) top = explorer->get_top_unexplored(assumptions);
		top_found = false;
		Formula original_top = top;
		if(top.empty())	{ 
			current_depth--; 
			return;
		}
		Formula origin_top = top;		
		if(!is_valid(top, true, true)){
			streak = 0;
			MUS mus = shrink_formula(top, crits);
			mark_MUS(mus);
			if(depth > depthMUS) continue;
			Formula m = mus.bool_mus;
			extend_mus(origin_top, m);
			find_all_muses_duality_based_remus(m, crits, depth + 1);
		}
		else{	//top is necessarily an MSS of subset
			streak++;
			//prevent getting stuck in a small subset
			if(streak > 10 && depth > 0){ current_depth--; return; }
			vector<bool> model;
			//mark_MSS(top);
			block_down(top);	
			vector<int> crit_all;
			for(int i = 0; i < dimension; i++)
				if(subset[i] && !origin_top[i]){
					crit_all.push_back(i);
				}
			if(crit_all.size() == 1){
				crits[crit_all[0]] = true;
				continue;
			}	
			if(depth > depthMUS) continue;		
			for(auto crit: crit_all){
				Formula rec_subset = origin_top;
				rec_subset[crit] = true;
				Formula rec_crits = crits;
				rec_crits[crit] = true;
				find_all_muses_duality_based_remus(rec_subset, rec_crits, depth + 1);
			}
		}
	}
	current_depth--;
}
	

// helper funcion for the ReMUS algorithm
// modifies mus into a set S such that mus /subseteq S /subseteq top
void Master::extend_mus(Formula &top, Formula &mus){
	int origin_top_size = count_ones(top);
	int ones = count_ones(mus);
	for(int i = 0; i < dimension; i++){
		if(top[i] && !mus[i]){
			mus[i] = true;
			ones++;
			if(ones >= origin_top_size * dim_reduction) 
				break;
		}			
	}
}
