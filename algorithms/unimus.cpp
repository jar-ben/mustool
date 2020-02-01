#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


void Master::unimus(){
	Formula whole(dimension, true);
	unibase_union(1);
        cout << "initial uni size: " << count_ones(uni) << ", " << (count_ones(uni)/float(dimension)) << endl;
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
			unibase_union(1);
			seed = explorer->get_bot_unexplored_inside(uni);
		}	
	}	
}
	
