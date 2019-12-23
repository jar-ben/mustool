#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

vector<int> Master::extendCS(MSS &mss, Formula &uni, Formula &couni){
	vector<int> added;
	Formula extended = union_sets(couni, mss.int_mss);
	Formula copy = extended;
	int i = 0;
	while(!is_valid(copy, true, false)){
		i++;
		if(count_ones(intersection(copy, couni)) > 10)
			copy = shrink_formula(copy).bool_mus;
		int diff = count_ones(intersection(copy, couni));
		int over = -1;
		for(int c = 0; c < dimension; c++){
			if(copy[c] && couni[c]){
				over = c;
				if(diff < 2) added.push_back(c);
			}
		}
		if(over == -1) print_err("over -1");
		extended[over] = false;
		copy = extended;
	}
	if(i > 0){
		Formula conflicts = subtract(extended, uni);
		vector<MSS> grown = grow_formulas(extended, conflicts, 1);
		for(auto g: grown){
			mark_MSS(g);
		}
	}else{
		mark_MSS(MSS(extended, mss.duration, mss.id, mss.seed_dimension));
	}
	return added;
}


void Master::unimus(){
	Formula whole(dimension, true);
	Formula uni = unibase_union(2);
        cout << "initial uni size: " << count_ones(uni) << ", " << (count_ones(uni)/float(dimension)) << endl;
	Formula couni(dimension, false);
	for(int i = 0; i < dimension; i++)
		if(!uni[i]) couni[i] = true;

	Formula seed = explorer->get_bot_unexplored_inside(uni);	
	while(!seed.empty()){
		bit++;
		if(is_valid(seed, true, true)){
			vector<int> missing = subtract_int(seed, whole);
			if(missing.size() < 5){			
				if(missing.size() > 1)
					grow_combined(seed);
				mark_MSS(MSS(seed, -1, msses.size(), count_ones(seed)));
			}else{
				MSS mss = grow_formula(seed);
				mark_MSS(mss);
			}	
		}else{
			block_up(seed);
		}
		seed = explorer->get_bot_unexplored_inside(uni);
		//the current uni is whole searched
		if(seed.empty()){
			block_up(uni);
			block_down(uni);
			Formula uni2 = unibase_union(2);
			if(count_ones(uni2) == 0) break; //all MSSes found
			uni = union_sets(uni, uni2);
			couni = subtract(uni, whole); 
			seed = explorer->get_bot_unexplored_inside(uni);
		}	
	}	
}
	
