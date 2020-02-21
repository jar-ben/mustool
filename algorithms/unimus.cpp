#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


bool Master::unimus_hitting_pair(int mid1, int mid2, int c1, int c2){
	BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);
	vector<int> implied1;
	vector<int> implied2;
	auto search1 = unimus_map[mid1].find(c1);
	if(search1 != unimus_map[mid1].end()){
		implied1 = search1->second;
	}else{
		implied1 = bSolver->get_implied(muses[mid1].bool_mus, c1);
		unimus_map[mid1][c1] = implied1;
	}
	auto search2 = unimus_map[mid2].find(c2);
	if(search2 != unimus_map[mid2].end()){
		implied2 = search2->second;
	}else{
		implied2 = bSolver->get_implied(muses[mid2].bool_mus, c2);
		unimus_map[mid2][c2] = implied2;
	}
	for(auto l: implied1)
		if(find(implied2.begin(), implied2.end(), -l) != implied2.end())
			return true;
	return false;
}

void Master::unimus_add_blocks(MUS &m1, int from, int to, vector<vector<int>> &blocks){
	for(int i = from; i <= to; i++){
		vector<int> block(1, i);
		for(auto c: muses[i].int_mus){
			if(!m1.bool_mus[c]){
				block.push_back(c);
			}
		}
		if(block.size() > 1)
			blocks.push_back(block);
	}
}

void reducePossibilities(Formula &m3, vector<pair<int,int>> &pairs){
	for(int i = pairs.size() - 1; i > -1; i--){		
		if(!m3[pairs[i].first] && !m3[pairs[i].second]){
			swap(pairs[i], pairs[pairs.size() - 1]);
			pairs.pop_back();
		}
	}
}

void Master::unimus_rotate_mus(int mid, int limit){
	MUS m1 = muses[mid];
	vector<vector<int>> blocks;
	unimus_add_blocks(m1, 0, muses.size() - 1, blocks);
	for(int i = mid - 40; i < mid; i++){
	//for(int i = 0; i < mid; i++){
		if(i < 0) continue;
		bool stay = abs(int(m1.int_mus.size() - muses[i].int_mus.size())) < 200;
		stay = stay || (abs(int(m1.int_mus.size() - muses[i].int_mus.size()))/float(dimension)) < 0.1;
		if(!stay) continue;
	

		vector<pair<int,int>> pairs;
		for(auto c1: m1.int_mus){
			for(auto c2: muses[i].int_mus){
				if(!muses[i].bool_mus[c1] && !m1.bool_mus[c2]){
					pairs.push_back(make_pair(c1,c2));
				}
			}
		}
		
		int initSize = pairs.size();
		int iter = 0;
		while(!pairs.empty()){
			iter++;
			pair<int,int> candidate = pairs[pairs.size() - 1];
			pairs.pop_back();
			unimus_attempts++;
			if(!unimus_hitting_pair(mid, i, candidate.first, candidate.second))
				continue;
			//check if the seed is unexplored
			bool ok = true;
			for(auto &block: blocks){
				int musId = block[0];
				ok = muses[musId].bool_mus[candidate.first] || muses[musId].bool_mus[candidate.second];
				if(!ok){
					for(int k = 1; k < block.size(); k++){
						int c = block[k];	
						ok = !muses[i].bool_mus[block[k]];
						if(ok) break;
					}
				}
				if(!ok){
					reducePossibilities(muses[musId].bool_mus, pairs);
					break;
				}
			}
			if(!ok) continue;
			Formula seed = union_sets(m1.bool_mus, muses[i].int_mus);
			seed[candidate.first] = seed[candidate.second] = false;
			unimus_rotated++;
			MUS mus = shrink_formula(seed);					
			mark_MUS(mus);
			unimus_add_blocks(m1, muses.size() - 1, muses.size() - 1, blocks);
			unimus_rotation_queue.push(muses.size() - 1);
		}
//		cout << "malina " << iter << " " << initSize << " " << (float(iter)/initSize) << endl;
	}
}


void Master::unimus_mark_mus(MUS &mus){
	mark_MUS(mus);
	unimus_rotation_queue.push(muses.size() - 1);
	while(!unimus_rotation_queue.empty()){
		int mid = unimus_rotation_queue.front();
		unimus_rotation_queue.pop();
		unimus_rotate_mus(mid, 1);
	}
}

void Master::unimus(){

	Formula top = explorer->get_unexplored(1, false);
        int streak = 0;
	bit = 0;
	while(!top.empty()){
                bit++;
		Formula original_top = top;
                if(is_valid(top, true, true)){
                       	streak++;
		       	block_down(top);
			//MSS mss = grow_formula(top);
			//mark_MSS(mss, true);
                }else{
                        MUS mus = shrink_formula(top);
			unimus_mark_mus(mus);
		}
                if(streak < 10)
			top = explorer->get_top_unexplored_inside(uni);	
		if(top.empty() || streak > 9)
			top = explorer->get_unexplored(1, false);
        }
}


void Master::unimus2(){
	Formula bot = explorer->get_unexplored(0, false);
        int streak = 0;
	bit = 0;
	while(!bot.empty()){
		bot = union_sets(bot,couni);
                bit++;
                if(is_valid(bot, true, true)){
                       	streak++;
		       	//block_down(bot);
			MSS mss = grow_formula(bot);
			mark_MSS(mss, true);
                }else{
                        MUS mus = shrink_formula(bot);
                        unimus_mark_mus(mus);
		}
                if(streak < 10)
			bot = explorer->get_bot_unexplored_inside(uni);	
		if(bot.empty() || streak > 9)
			bot = explorer->get_unexplored(0, false);
        }
}

