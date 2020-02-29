#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>



void Master::unimusRec(Formula subset, Formula crits, int depth){
        std::vector<int> assumptions;
        for(int i = 0; i < dimension; i++){
                if(!subset[i]){
                        assumptions.push_back(-1 * i - 1);
                }
        }

        Formula top;
        int streak = 0;
        int iteration = 0;
        bool top_found = false;
        while(true){
                iteration++;
		if(!top_found) top = explorer->get_top_unexplored(assumptions);
                top_found = false;
                if(top.empty()) {
                        return;
                }
                Formula origin_top = top;
                if(!is_valid(top, true, true)){
                        streak = 0;
                        MUS mus = shrink_formula(top, crits);
                        mark_MUS(mus);
                }
                else{   //top is necessarily an MSS of subset
                        streak++;
                        //prevent getting stuck in a small subset
                        if(streak > 10){ current_depth--; return; }
                        //mark_MSS(top);
                        block_down(top);
                        vector<int> crit_all;
                        for(int i = 0; i < dimension; i++)
                                if(subset[i] && !origin_top[i]){
                                        crit_all.push_back(i);
                                }
                        if(crit_all.size() == 1){
                                crits[crit_all[0]] = true;
                                continue;
                        }
                        for(auto crit: crit_all){
                                Formula rec_subset = origin_top;
                                rec_subset[crit] = true;
                                Formula rec_crits = crits;
                                rec_crits[crit] = true;
                                unimusRec(rec_subset, rec_crits, depth + 1);
                        }
                }
        }
}

//return false if Unexplored is empty
bool Master::unimusRecRefine(){
	int found = 0;
        unimus_refines++;
	Formula seed = explorer->get_unexplored(1, false);
        while(found < 1 && !seed.empty()){
                if(!is_valid(seed, true, false)){
                        MUS mus = shrink_formula(seed);
                        //unimus_mark_mus(mus);
                        mark_MUS(mus);
                        found++;
                }else{
                        mark_MSS(seed);
                        //found++;
                }
                seed = explorer->get_unexplored(1, false);
        }
	return !seed.empty();
}


void Master::unimusRecMain(){
	while(unimusRecRefine()){                
		cout << "start of bit " << ++bit << endl;		
		unimusRec(uni,Formula(dimension, false),0);
		cout << "end of bit " << bit << endl;
	}
}
