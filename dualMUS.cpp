#include "Master.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include <unordered_set>

vector<int> Master::minimal_hitting_set(vector<int> local_muses){
	vector<int> hs;
	for(auto mid: local_muses){
		auto &mus = muses[mid];
		bool hitten = false;
		for(auto c: hs){ if(mus.bool_mus[c]){ hitten = true; break; } }
		if(hitten) continue;	
		hs.push_back(mus.without_crits[0]);
	}	
	return hs;
}

int Master::recursive_rotation_delta(MUS m1, Formula &top, int depth){
	if(!explorer->easy_to_shrink(m1)) return 0;
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

	int round = 0;
	for(auto c1: m1.without_crits){
		if(explorer->mus_intersection[c1]) continue;
		vector<int> pairs;	
		for(auto c: top_int){
			if(explorer->is_hitting_pair(c, c1)){ pairs.push_back(c); }
		}		
		for(auto c2: pairs){			
			for(auto mid: explorer->parent_muses[c2]){
				if(indicator[mid] != 1) continue;// || already_used_mus[mid]) continue;
				auto &m2 = muses[mid];
				if(m2.bool_mus[c1] || !explorer->easy_to_shrink(m2)) continue;
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

