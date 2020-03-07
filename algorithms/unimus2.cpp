#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include <numeric>

void Master::unimusRec_add_block(MUS &m1, int mid, vector<vector<int>> &blocks, vector<vector<int>> &blocks_hitmap){
	vector<int> block(1, mid);
	for(auto c: muses[mid].int_mus){
		if(!m1.bool_mus[c]){
			block.push_back(c);
			blocks_hitmap[c].push_back(blocks.size());
		}
	}
	if(block.size() > 1){
		blocks.push_back(block);
	}
}

bool Master::unimusRec_isAvailable(int c, Formula &subset, vector<vector<int>> &blocks, vector<vector<int>> &blocks_hitmap){
	for(auto bid: blocks_hitmap[c]){
                bool covered = false;
                for(int i = 1; i < blocks[bid].size(); i++){ //the first element of a block is the mid of the MUS
			auto &c2 = blocks[bid][i];
                        if(c2 != c && !subset[c2]){
                                covered = true;
                                break;
                        }
                }
                if(!covered) return false;
        }
        return true;
}

Formula Master::unimusRec_propagateRefine(Formula &conflict, Formula &base, vector<pair<int,int>> &implied){
	Formula toKeep(dimension, true);
	BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);        
	int finalOne = implied.back().first;
	int removed = 0;
	for(int i = implied.size() - 2; i >= 0; i--){
		bool needed = false;
		int c1 = implied[i].first;
		int lit = implied[i].second;
		for(int j = i + 1; j < implied.size(); j++){
			int c2 = implied[j].first;
			if(!toKeep[c2]) continue;
			vector<int> &cl = bSolver->clauses[c2];
			if(find(cl.begin(), cl.end(), -lit) != cl.end()){
				needed = true;
				break;
			}		
		}
		if(!needed && !base[c1] && c1 != finalOne){
			toKeep[c1] = false;
			removed++;
		}
	}
	for(int i = 0; i < dimension; i++){
		if(!toKeep[i])
			conflict[i] = false;
	}
	if(DBG){
		if(is_valid(conflict)) print_err("satisfiable seed after unimusRec_propagateRefine()");
	}

	return conflict;
}

Formula Master::unimusRec_propagateToUnsat(Formula base, Formula cover, vector<pair<int,int>> implied, vector<vector<int>> &blocks, vector<vector<int>> &blocks_hitmap){
	BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);        
        queue<int> toPropagate;
        for(auto p: implied)
                toPropagate.push(p.second);
        vector<bool> conflict = base;
        vector<int> lits_left (dimension,0);
        for(int i = 0 ; i < dimension; i++){
                lits_left[i] = bSolver->clauses[i].size() - 1; // -1 because there is the activation literal
        }
        vector<bool> satisfied(dimension, false);
        vector<bool> setup(bSolver->vars + 1, false);
        vector<bool> value(bSolver->vars + 1, false);
        int added = 0;
	while(!toPropagate.empty()){
		if(added > 20){
		//	return Formula(); //Optimalization. The seed is usually find either soon or not at all. 
		}
		auto l = toPropagate.front();
                toPropagate.pop();
                int var = (l > 0)? l : -l;
                if(setup[var])
                        continue;
                setup[var] = true;
                value[var] = l < 0;
                if(l < 0){
                        for(auto c1: bSolver->hitmap_neg[var - 1]){ //clauses satisfied by implied negative value of i-th literal
                                satisfied[c1] = true;
                        }
                        for(auto c1: bSolver->hitmap_pos[var - 1]){ //trim clauses that contain positive i-th literal
                                if(satisfied[c1] || !cover[c1]) continue;
                                lits_left[c1]--;
                                if(lits_left[c1] == 1 && unimusRec_isAvailable(c1, conflict, blocks, blocks_hitmap)){
                                        for(auto lit: bSolver->clauses[c1]){
                                                if(lit >= bSolver->vars) break; //the activation literal
                                                int var2 = (lit > 0)? lit : -lit;
                                                if(!setup[var2]){
                                                        toPropagate.push(lit);
                                                        implied.push_back(make_pair(c1,lit));
                                                        break;
                                                }
                                        }
                                        if(!conflict[c1]) added++;
                                        conflict[c1] = true;
                                }
                                if(lits_left[c1] == 0 && unimusRec_isAvailable(c1, conflict, blocks, blocks_hitmap)){
					conflict[c1] = true;
					implied.push_back(make_pair(c1, 0));
					return unimusRec_propagateRefine(conflict, base, implied);
                                }
                        }
                }else{
                        for(auto c1: bSolver->hitmap_pos[var - 1]){ //clauses satisfied by implied negative value of i-th literal
                                satisfied[c1] = true;
                        }
                        for(auto c1: bSolver->hitmap_neg[var - 1]){ //trim clauses that contain positive i-th literal
                                if(satisfied[c1] || !cover[c1]) continue;
                                lits_left[c1]--;
                                if(lits_left[c1] == 1 && unimusRec_isAvailable(c1, conflict, blocks, blocks_hitmap)){
                                        for(auto lit: bSolver->clauses[c1]){
                                                if(lit >= bSolver->vars) break; //the activation literal
                                                int var2 = (lit > 0)? lit : -lit;
                                                if(!setup[var2]){
                                                        toPropagate.push(lit);
                                                        implied.push_back(make_pair(c1,lit));
                                                        break;
                                                }
                                        }
                                        if(!conflict[c1]) added++;
                                        conflict[c1] = true;
                                }
                                if(lits_left[c1] == 0 && unimusRec_isAvailable(c1, conflict, blocks, blocks_hitmap)){
                                        conflict[c1] = true;
					implied.push_back(make_pair(c1, 0));
					return unimusRec_propagateRefine(conflict, base, implied);
                                }
                        }
                }
        }
        return Formula();
}


void Master::unimusRec_rotate_mus(int mid, Formula cover, Formula subset, vector<int> &localMuses){
        MUS m1 = muses[mid];
        vector<vector<int>> blocks;
	vector<vector<int>> blocks_hitmap(dimension);
        for(auto mi: localMuses)
		unimusRec_add_block(m1, mi, blocks, blocks_hitmap);
	
	BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);        
	int iter = 0;
	int limit = 500;
	if(m1.int_mus.size() > limit)
		shuffle_ints(m1.int_mus);
	for(auto c1: m1.int_mus){			
		if(explorer->is_critical(c1, cover)) continue;		
		iter++;
		if(iter > limit) break;
		Formula base = m1.bool_mus;
		base[c1] = false;
		cover[c1] = false;
		vector<pair<int,int>> implied;
		for(int k = 0; k < bSolver->clauses[c1].size() - 1; k++)
			implied.push_back(make_pair(c1, -1 * bSolver->clauses[c1][k]));
		Formula seed = unimusRec_propagateToUnsat(base, cover, implied, blocks, blocks_hitmap);
		cover[c1] = true;
		unimus_attempts++;	
		if(seed.empty()) continue;
		if(DBG){
			if(!is_subset(base, seed)) print_err("the base is not a subset of the seed");
			if(is_valid(seed)) print_err("the rotated seed is valid A");
			if(!is_subset(seed, cover)) print_err("the seed is not a subset of the cover");
			for(auto mi: localMuses){
				if(is_subset(seed, muses[mi].bool_mus)) print_err("the seed is a subset of a local MUS " + to_string(mi));
			}
			for(int mi = 0; mi < muses.size(); mi++){
				if(is_subset(seed, muses[mi].bool_mus)) print_err("the seed is a subset of a MUS " + to_string(mi));
			}
		}		

		if(DBG){
			if(is_valid(seed)) print_err("the rotated seed is valid");
		}
		unimus_rotated++;
		MUS mus = shrink_formula(seed);
		mark_MUS(mus);
		localMuses.push_back(muses.size() - 1);
		unimusRec_add_block(m1, muses.size() - 1, blocks, blocks_hitmap);
		unimus_rotation_stack.push(muses.size() - 1);
	}
}


void Master::unimusRec_mark_mus(MUS &mus, Formula cover, Formula subset){
	vector<int> localMuses {muses.size() - 1};
	unimus_rotation_stack.push(muses.size() - 1);
	int streak = 0;
	while(!unimus_rotation_stack.empty()){
		int mid = unimus_rotation_stack.top();
		unimus_rotation_stack.pop();
		int before = muses.size();
		unimusRec_rotate_mus(mid, cover, subset, localMuses);
		if(muses.size() == before)
			streak++;
		else
			streak = 0;
		if(streak == 5){
			unimus_rotation_stack = stack<int>();
			break;
		}
	}
}

void Master::unimusRec(Formula subset, Formula crits, int depth){
	Formula processed(dimension, false);
	currentSubset = subset;
	unimusRecDepth = depth;
	std::vector<int> assumptions;
        for(int i = 0; i < dimension; i++){
                if(!subset[i]){
                        assumptions.push_back(-1 * i - 1);
                }
        }
	
        Formula top;
        int streak = 0;
        int iteration = 0;
	bool unsatFound = false;
        while(true){
		
		iteration++;
		top = explorer->get_top_unexplored(assumptions);
                if(top.empty()) {
                        return;
                }
                Formula origin_top = top;
                if(!is_valid(top, true, true)){
                        unsatFound = true;
			streak = 0;
                        MUS mus = shrink_formula(top, crits);
                        //MUS mus = shrink_formula(top);
			mark_MUS(mus);
                        //if(depth > 0)
				unimusRec_mark_mus(mus, origin_top, subset);
                }
                else{   //top is necessarily an MSS of subset
                        block_down(top);
                        streak++;
                        //prevent getting stuck in a small subset
                        if(streak > 5 || iteration == 1){ return; }
			vector<int> crit_all;
                        for(int i = 0; i < dimension; i++)
                                if(subset[i] && !origin_top[i]){
                                        crit_all.push_back(i);
                                }
                        if(crit_all.size() == 1){
                                crits[crit_all[0]] = true;
                                continue;
                        }
			if(!unsatFound) continue; // optimization step
                        for(auto crit: crit_all){
				if(processed[crit])
					continue;
				processed[crit] = true;

				BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);
				if(!bSolver->lit_occurences(subset, crit)){
					continue; // this means that the lits in crits are not contained in subset - {crit}
			       			  // and, thus, there is a low chance that we find rotation pairs in subset - {crit}
						  // that are based on {crit}
						  // this is a heuristic backtrack 	
				}

				Formula rec_subset = subset;
                                //rec_subset[crit] = false;

				for(auto crit2: crit_all)
					rec_subset[crit2] = false;
				rec_subset[crit] = true;

				//pick a random subset of the rotationMuses
				unimusRec(rec_subset, crits, depth + 1);
				currentSubset = subset;
				unimusRecDepth = depth;
                        }
                }
        }
}

//return false if Unexplored is empty
bool Master::unimusRecRefine(){
	int foundMUS = 0;
	int foundMSS = 0;
        unimus_refines++;
	Formula seed = explorer->get_unexplored(1, false);
        while(foundMUS < 1 && foundMSS < 50 && !seed.empty()){
                Formula origin_top = seed;
		if(!is_valid(seed, true, false)){
                        MUS mus = shrink_formula(seed);
                        //unimus_mark_mus(mus);
                        mark_MUS(mus);
				unimusRec_mark_mus(mus, origin_top, origin_top);
                        foundMUS++;
                }else{
                        mark_MSS(seed);
			foundMSS++;
                }
                seed = explorer->get_unexplored(1, false);
        }
	return !seed.empty();
}


void Master::unimusRecMain(){
	while(unimusRecRefine()){                
		unimusRec(uni,Formula(dimension, false), 0);
	}
}
