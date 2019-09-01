#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include <unordered_set>


struct OverlapHash{
	public:
		std::size_t operator()(std::vector<int> const& over) const{
			std::size_t ret = over.size();
			for(auto &c: over){
				ret ^= c + 0x9e3779b9 + (ret << 6) + (ret >> 2);
			}
			return ret;
		}
};

int Master::mixed_mus_rotation(MUS &m1, Formula &top){
	int rotated = 0;
	int mid_limit = muses.size() - scope_limit;
	if(mid_limit < 0) mid_limit = 0;

	vector<int> top_int;
	for(int i = 0; i < dimension; i++)
		if(!top[i]) top_int.push_back(i);

	vector<int> indicator(muses.size(), 0);
	for(int mid = muses.size() - 1; mid > mid_limit; mid--){
		auto &m = muses[mid];
		for(auto c2: top_int)
			if(m.bool_mus[c2])
				indicator[mid]++;
	}

	vector<vector<int>> local_muses_overlaps;
	vector<int> local_muses_ids;
	unordered_set<vector<int>, OverlapHash> overlaps;	
	vector<bool> already_used_mus(muses.size(), false);

	vector<int> variables_map;
        vector<int> variables_map_inv;
	backbone_build_literal_map(m1.bool_mus, variables_map, variables_map_inv);

	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	int vars = msSolver->vars;
	int round = 0;
	int mainRound = 0;
	for(auto c1: m1.without_crits){		
		if(explorer->mus_intersection[c1]) continue;
		if(++mainRound > 10) break;
		vector<int> pairs;
		Formula implied(vars, false);
                Formula values(vars, false);
		Formula seed = m1.bool_mus; seed[c1] = false;
		backbone_simplify(seed, c1, implied, values);
		get_backbones_bones(seed, implied, values, variables_map, variables_map_inv);	
		for(int i = 0; i < vars; i++){
			if(implied[i]){
				vector<int> &hits = (values[i])? msSolver->hitmap_pos[i] : msSolver->hitmap_neg[i];
				for(auto c2: hits)
					if(!top[c2]) pairs.push_back(c2);
			}
		}

		for(auto c2: pairs){			
			for(auto mid: explorer->parent_muses[c2]){
				if(indicator[mid] != 1) continue;// || already_used_mus[mid]) continue;
				auto &m2 = muses[mid];
				if(m2.bool_mus[c1]) continue;
				round++;
				vector<int> m1_overlap;
				for(auto o: m2.without_crits)
					if(o != c2 && !m1.bool_mus[o])
						m1_overlap.push_back(o);
				overlaps.insert(m1_overlap);
				bool unexplored = true;					
				for(int i = 0; i < local_muses_overlaps.size(); i++){
					if(muses[local_muses_ids[i]].bool_mus[c1]) continue; //hitten by c1
					auto &overlap = local_muses_overlaps[i];
					bool hit = false;
					for(auto l: overlap){
						if(!m2.bool_mus[l]){
							hit = true; break; 
						}
					}
					if(!hit){ unexplored = false; break; }
				}
				if(!unexplored) continue;
				Formula seed = m1.bool_mus; seed[c1] = false;
				for(auto o: m1_overlap) seed[o] = true;

				MUS m = shrink_formula(seed);
				rotated_muses++;
				mark_MUS(m);
				rotated++;
				already_used_mus[mid] = true;
				//local mus overlap
				vector<int> overlap;
				for(auto l: m.without_crits)
					if(!m1.bool_mus[l])
						overlap.push_back(l);
				local_muses_overlaps.push_back(overlap);
				local_muses_ids.push_back(muses.size() - 1);
			}
		}		
	}
	return rotated;
}

