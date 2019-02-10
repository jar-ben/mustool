#include "Solver.h"
#include "SimpleExplorer.h"
#include "misc.h"
#include "SimpleSAT.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


void backbone_marco(){
	vector<int> L;
	
}

Formula mapit(Formula f, vector<int> &mapa, int dim){
	Formula result(dim, false);
	for(int i = 0; i < f.size(); i++){
		if(f[i]){
			if(mapa[i] >= 0)
				result[mapa[i]] = true;
		}
	}
	return result;
}


int reduce_vars(vector<vector<int>> &cls, int vars){
	vector<int> mapa(vars + 1,-1);
	int current = 1;

	for(auto &cl: cls){
		for(auto l: cl){
			int var = (l > 0)? l : -l;
			if(mapa[var] == -1){
				mapa[var] = current;
				current++;
			}
		}
	}
	vector<vector<int>> reduced;
	for(auto &cl: cls){
		vector<int> rcl;
		for(auto l: cl){
			bool sign = l > 0;
			int var = (sign > 0)? l : -l;
			rcl.push_back( mapa[var] * (sign)? 1 : -1 );
		}
		reduced.push_back(rcl);
	}
	cls = reduced;
	return current - 1;
}

void Solver::build_mapping(Formula &top, Formula &implied, Formula &values, Formula &core, vector<vector<int>> &reduced,
		vector<int> &map_o_to_n, vector<int> &map_n_to_o, int &max_var, vector<int> &implied_clauses){
	Formula ommit(dimension, false);

	int original_lits = 0;
	int final_lits = 0;

	for(int i = 0; i < dimension; i++){
		if(!top[i]){
			ommit[i] = true;
			continue;
		}
		if(core[i]){
			vector<int> rcl = satSolver->clauses[i];
			rcl.pop_back();
			reduced.push_back(rcl);
			continue;
		}
		auto &cl = satSolver->clauses[i];
		bool skip = false;
		vector<int> rcl;
		int lits_to_add = 0;
		for(int lid = 0; lid < cl.size() - 1; lid++){
			original_lits++;
			bool pos = cl[lid] > 0;
			int var = (pos)? cl[lid] - 1: -cl[lid] - 1;
			if(implied[var] && (pos == values[var]))
				skip = true;
			if(!implied[var]){
				rcl.push_back(cl[lid]);
				lits_to_add++;
			}
		}
		if(!skip){
			reduced.push_back(rcl);
			final_lits += lits_to_add;
		}
		else{
			ommit[i] = true;
			implied_clauses.push_back(i);
		}
	}
	map_o_to_n.resize(dimension, -1);
	int reduced_dimension = reduced.size();
	map_n_to_o.resize(reduced_dimension);

	int current = 0;
	for(int i = 0; i < dimension; i++){
		if(!ommit[i]){
			map_o_to_n[i] = current;
			map_n_to_o[current] = i;
			current++;
		}
	}

	max_var = satSolver->vars;
}

int Solver::backbone_mus_rotation_all(MUS &m1, Formula &top){
	if(!explorer->easy_to_shrink(m1)) return 0;

	int here_found = 0;
	int here_sat = 0;
	int here_limit = 10;
	int here_sat_limit = 30;
	SimpleExplorer exp (dimension);
	exp.block_up(m1.bool_mus);
			

	int vars = satSolver->vars;	
	vector<int> violated;
	vector<int> variables_map;
	vector<int> variables_map_inv;
	backbone_build_literal_map(m1.bool_mus, variables_map, variables_map_inv);
	
	vector<Formula> seeds;

	int main_iteration = 0;
	for(auto c: m1.int_mus){
		if(explorer->is_critical(c,top)) continue;

		main_iteration++;
		vector<int> assumptions;
		for(int i = 0; i < dimension; i++)
			if(m1.bool_mus[i] && i != c)
				assumptions.push_back(i + 1);
		//extend top if possible
		if(explorer->mus_intersection[c]){
			assumptions.push_back(-1 * c - 1);
		}else{
			for(int i = 0; i < dimension; i++)
				if(!top[i] || i == c)
					assumptions.push_back(-1 * i - 1);
		}		

		Formula t = exp.get_top_unexplored(assumptions);
		if(t.empty()) continue;

		Formula core_ = m1.bool_mus; core_[c] = false;
		Formula implied_ (vars, false);
		Formula values_ (vars, false);
		vector<int> lits_left_ (dimension,0);
		vector<bool> satisfied_ (dimension, false);
		vector<int> singletons_;
		vector<int> violated  = backbone_init(core_, implied_, values_, lits_left_, satisfied_, singletons_, t, variables_map, variables_map_inv, c);
		
		if(!violated.empty()){
			for(auto c2: violated){
				Formula seed = core_; 
				seed[c2] = true;
				MUS mus = shrink_formula(seed);
				rotated_muses++;
				mark_MUS(mus);
				exp.block_up(mus.bool_mus);
			}
			t = exp.get_top_unexplored(assumptions);
		}	
		int iteration = 0;	
	
		vector<vector<int>> reduced;
		vector<int> map_o_to_n;
		vector<int> map_n_to_o;
		int max_var;
		vector<int> implied_clauses;
		build_mapping(t, implied_, values_, core_, reduced, map_o_to_n, map_n_to_o, max_var, implied_clauses);
		Formula over_t = t;
		SimpleSAT ssolver(reduced, max_var);
		for(int i = 0; i < implied_.size(); i++)
			if(implied_[i])
				ssolver.add_unit( (values_[i])? i + 1 : -(i + 1) );
		SimpleExplorer sexplorer(reduced.size());
		Formula core_r = mapit(core_, map_o_to_n, reduced.size());

		while(!t.empty()){
			Formula core = core_;
			Formula implied = implied_;
			Formula values = values_;
			vector<int> lits_left = lits_left_;
			Formula satisfied = satisfied_;
			vector<int> singletons = singletons_;

			vector<int> violated = backbone_find_seed_beta( core, implied, values, satisfied, singletons, lits_left, t );
			if(violated.empty()){
				Formula mapped = mapit(t, map_o_to_n, reduced.size());
				if(ssolver.solve(mapped, true, true)){
					for(auto ic: implied_clauses)
						t[ic] = true;
					block_down(t); //blocked after the while cycle via over_t	
					sexplorer.block_down(mapped);
					exp.block_down(t);
					here_sat++;
					if(here_sat > here_sat_limit) return here_found;
				}else{
					Formula unsat_core = mapit(mapped, map_n_to_o, dimension);
					for(int i = 0; i < dimension; i++)
						if(core_[i])
							unsat_core[i] = true;
					MUS mus = shrink_formula(unsat_core);
					mark_MUS(mus);
					rotated_muses++;
					here_found++;
					exp.block_up(mus.bool_mus);
					sexplorer.block_up( mapit(mus.bool_mus, map_o_to_n, reduced.size()) );
					if(here_found > here_limit) return here_found;
				}
			}else{
				int c2 = violated[0];
				Formula seed = core;
				if(c2 >= 0) 
					seed[c2] = true;
				MUS mus = shrink_formula(seed);
				rotated_muses++;
				mark_MUS(mus);
				exp.block_up(mus.bool_mus);
				sexplorer.block_up( mapit(mus.bool_mus, map_o_to_n, reduced.size()) );
			}
			t = sexplorer.get_top_unexplored(core_r);
			if(t.empty()) break;
			t = mapit(t, map_n_to_o, dimension);
		}
		exp.block_down(over_t);
		block_down(over_t);//save to perform provided that all MUSes of over_t have been found
	}
	return here_found;
}

