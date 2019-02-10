#include "Solver.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


//remus
void Solver::find_all_muses_duality_based_remus(Formula subset, Formula crits, int depth){
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
		current_dimension = count_ones(subset);
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
			Formula rotation_top;
			if(explorer->test_rotation_unex){ rotation_top = explorer->get_top_unexplored(mus.bool_mus); }
			mark_MUS(mus);

			if(variant == 11){
				int rot = recursive_rotation_delta(mus, original_top, 0);
				cout << "recursively rotated: " << rot << endl;
			}
			if(variant == 12){
				int rot = backbone_mus_rotation(mus, original_top);
				cout << "recursively rotated: " << rot << endl;
 			}

			if(depth > depthMUS) continue;
			Formula m = mus.bool_mus;
			extend_mus(origin_top, m);
			find_all_muses_duality_based_remus(m, crits, depth + 1);
		}
		else{	//top is necessarily an MSS of subset
			streak++;
			vector<bool> model;
			if(model_rotation){
				model = satSolver->get_model();
			}
			block_down(top);	
			vector<int> crit_all;
			for(int i = 0; i < dimension; i++)
				if(subset[i] && !origin_top[i]){
					crit_all.push_back(i);
				}
			if(crit_all.size() == 1){
				crits[crit_all[0]] = true;
				if(model_rotation){
					vector<vector<bool>> model_extensions;
					int rotated = satSolver->model_rotation(crits, crit_all[0], subset, model, model_extensions);
					for(auto &extension: model_extensions){
						if(false && !explorer->checkValuation(extension)) explored_extensions++;
						else block_down(extension);
					}					
				}
				continue;
			}	
			if(depth > depthMUS) continue;		
			for(auto crit: crit_all){
				Formula rec_subset = origin_top;
				rec_subset[crit] = true;
				Formula rec_crits = crits;
				rec_crits[crit] = true;
				if(model_rotation){			
					vector<vector<bool>> model_extensions;
					for(auto &extension: model_extensions){
						if(false && !explorer->checkValuation(extension)) explored_extensions++;
						else block_down(extension);
					}
				}
				find_all_muses_duality_based_remus(rec_subset, rec_crits, depth + 1);
			}
		}
	}
	current_depth--;
}
	

// helper funcion for the ReMUS algorithm (experimental, not fully integrated yet)
// modifies mus into a set S such that mus /subseteq S /subseteq top
void Solver::extend_mus(Formula &top, Formula &mus, int dMUS){
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
