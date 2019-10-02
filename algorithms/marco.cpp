#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

void Master::marco_base(){
	int crit_candidate = 0;
	Formula top = explorer->get_unexplored(1, false);
	while(!top.empty()){
		Formula original_top = top;
		if(is_valid(top, true, true)){
			block_down(top);
			unex_sat++;
		}else{
			MUS mus = shrink_formula(top);
			mark_MUS(mus);
			unex_unsat++;			
		}
		top = explorer->get_unexplored(1, false);
	}
}

