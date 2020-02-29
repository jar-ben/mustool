#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


//TODO: return a pointer
vector<int> Master::unimus_get_implied(int mid, int c){
	BooleanSolver *bSolver = static_cast<BooleanSolver*>(satSolver);
	auto search = unimus_map[mid].find(c);
	if(search != unimus_map[mid].end())
		return search->second;
	vector<int> implied = bSolver->get_implied(muses[mid].bool_mus, c);
	unimus_map[mid][c] = implied;
	return implied;
}

bool Master::unimus_hitting_pair(int mid1, int mid2, int c1, int c2){
	vector<int> implied1 = unimus_get_implied(mid1, c1);
	vector<int> implied2 = unimus_get_implied(mid2, c2);
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
	for(int i = max(0,mid - 200); i < mid; i++){
		bool stay = abs(int(m1.int_mus.size() - muses[i].int_mus.size())) < 20;
		int cointersection = 0;
		for(int c = 0; c < dimension; c++)
			if((m1.bool_mus[c] && !muses[i].bool_mus[c]) || (!m1.bool_mus[c] && muses[i].bool_mus[c]))
				cointersection++;
		stay = stay && cointersection < 10;
		//		stay = stay || (abs(int(m1.int_mus.size() - muses[i].int_mus.size()))/float(dimension)) < 0.1;
		if(!stay) continue;
	
		vector<pair<int,int>> pairs;
		unordered_map<int, vector<int>> litHitmap;
		for(auto c1: m1.int_mus){
			if(!muses[i].bool_mus[c1]){
				vector<int> implied = unimus_get_implied(mid, c1);
				for(auto lit: implied){
					if(litHitmap.find(lit) == litHitmap.end())
						litHitmap[lit] = vector<int>(1, c1);
					else
						litHitmap[lit].push_back(c1);
				}
			}
		}
		for(auto c2: muses[i].int_mus){
			if(!m1.bool_mus[c2]){
				vector<int> implied = unimus_get_implied(i, c2);
				for(auto lit: implied){
					auto search = litHitmap.find(-1 * lit);
					if(search != litHitmap.end()){
						for(auto c1: search->second){
							pairs.push_back(make_pair(c1,c2));
						}
					}
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
			//check if the seed is unexplored
			bool ok = true;
			for(auto &block: blocks){
				int musId = block[0];
				ok = muses[musId].bool_mus[candidate.first] || muses[musId].bool_mus[candidate.second];
				if(!ok){
					for(int k = 1; k < block.size(); k++){
						ok = !muses[i].bool_mus[block[k]];
						if(ok) break;
					}
				}
					//reducePossibilities(muses[musId].bool_mus, pairs);
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
			if(unimus_use_stack){
				unimus_rotation_stack.push(muses.size() - 1);
			}
			else{
				unimus_rotation_queue.push(muses.size() - 1);
			}
		}
	}
}


void Master::unimus_mark_mus(MUS &mus){
	mark_MUS(mus);
	//unimus_use_stack = false;
	int streak = 0;
	if(unimus_use_stack){
		unimus_rotation_stack.push(muses.size() - 1);
		while(!unimus_rotation_stack.empty()){
			int mid = unimus_rotation_stack.top();
			unimus_rotation_stack.pop();
			int before = muses.size();
			unimus_rotate_mus(mid, 1);
			if(muses.size() == before)
				streak++;
			else
				streak = 0;
			if(streak == 5){
				unimus_rotation_stack = stack<int>();
				break;
			}
		}
	} else {
		unimus_rotation_queue.push(muses.size() - 1);
		while(!unimus_rotation_queue.empty()){
			int mid = unimus_rotation_queue.front();
			unimus_rotation_queue.pop();
			unimus_rotate_mus(mid, 1);
		}
	}
}

void Master::unimus_refine(){
	int found = 0;
	unimus_refines++;
	while(found < 1){
		Formula seed = explorer->get_unexplored(1, false);
		if(!is_valid(seed, true, false)){
                        MUS mus = shrink_formula(seed);
			unimus_mark_mus(mus);
			//mark_MUS(mus);
			found++;	
		}else{		
			//block_down(seed);	
			mark_MSS(seed);
			//block_down(mss.bool_mss);
			//found++;
		}
	}
	cout << "end of refine" << endl;	
}

/*void Master::unimus_refine(){
	int found = 0;
	while(found < 1){
		Formula bot = explorer->get_unexplored(0, false);
		if(!is_valid(bot, false, true)){
			muses.push_back(MUS(bot, -1, muses.size(), -1));
			MUS mus = muses.back();
			mark_MUS(mus);
		}else{			
			MSS mss = grow_formula(bot);
			mark_MSS(mss);
			//block_down(mss.bool_mss);
			found++;
		}
	}
	cout << "end of refine" << endl;	
}*/

//this is the alternation based version
/*void Master::unimus(){
        int streak = 0;
	bit = 0;
	bool bot = false;
	int bot_muses = 0;
	while(true){
		bit++;
		if(streak == 5){
			unimus_refine();
			streak = 0;
		}
		if(bot){
			Formula seed = explorer->get_bot_unexplored_inside(uni);
			if(seed.empty()){			
				unimus_refine();		
				seed = explorer->get_unexplored(0, false);
				if(seed.empty()) break; //unexplore = \emptyset
			}
			if(is_valid(seed, false, true)){
				streak++;
				block_down(seed);
				bot = false;
			}else{
				streak = 0;
				muses.push_back(MUS(seed, -1, muses.size(), -1));
				MUS mus = muses.back();
				mark_MUS(mus);
				bot_muses++;
				cout << "bot_muses " << bot_muses << endl << endl << endl;
			}
		}else{
			Formula seed = explorer->get_top_unexplored_inside(uni);
			if(seed.empty()){			
				unimus_refine();		
				seed = explorer->get_unexplored(1, false);
				if(seed.empty()) break; //unexplore = \emptyset
			}
			if(is_valid(seed, true, true)){
				streak++;
				block_down(seed);
			}else{
				streak = 0;
				MUS mus = shrink_formula(seed);
				unimus_mark_mus(mus);
			}
			bot = true;
		}
	}
}*/

void Master::unimus(){
	Formula top = explorer->get_unexplored(1, false);
        int streak = 0;
	bit = 0;
	while(!top.empty()){
                bit++;
                if(is_valid(top, true, true)){
                       	streak++;
		       	//block_down(top);
			//cout << "growing" << endl;
			MSS mss = grow_formula(top);
			mark_MSS(mss, true);
                }else{
			streak = 0;
                        MUS mus = shrink_formula(top);
			unimus_mark_mus(mus);
		}
		top = explorer->get_top_unexplored_inside(uni);	
		if(top.empty() || streak > 10){			
			cout << "empty top" << endl;
			unimus_refine();		
			top = explorer->get_unexplored(1, false);
			streak = 0;
		}
	}
}
