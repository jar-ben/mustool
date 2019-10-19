#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>

void Master::marco_base(){
	Formula top = explorer->get_unexplored(1, false);
	while(!top.empty()){
		Formula original_top = top;
		if(is_valid(top, true, true)){
			block_down(top);
		}else{
			MUS mus = shrink_formula(top);
			mark_MUS(mus);
		}
		top = explorer->get_unexplored(1, false);
	}
}

