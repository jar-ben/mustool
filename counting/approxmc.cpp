#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include <cmath>

int Master::bsat(float tresh){
	Formula top = explorer->get_unexplored(1, false);
	int i = 0;
	int j = 0;
	while(i < tresh && !top.empty()){
		if(is_valid(top, true, true)){
			block_down(top);
			vector<CMSat::Lit> cl;
			for(int i = 0; i < dimension; i++)
				if(!top[i])
					cl.push_back(CMSat::Lit(i, false));
			blocksDown.push_back(cl);
		}else{
			Formula original = top;
			MUS mus = shrink_formula(top);
			mark_MUS(mus, true, false);
			vector<CMSat::Lit> cl;
			for(int i = 0; i < dimension; i++)
				if(mus.bool_mus[i])
					cl.push_back(CMSat::Lit(i, true));
			blocksUp.push_back(cl);
			i++;
			if(muses.size() % 100 == 0)
				std::cout << "bsat muses size: " << muses.size() << std::endl;
		}
		top = explorer->get_unexplored(1, false);
	}
	return i;
}

Formula Master::xor_shrink(Formula &f, XorExplorer &xe, int &limit){
	Formula top = xe.solve_subset(f);
	while(!top.empty() && limit > 0){
		Formula original = top;
		if(is_valid(top,true,true)){
			block_down(top);
			blocksDown.push_back(xe.block_down(top));
		}else{
			limit--;
			std::cout << count_ones(original) << " -> ";
			MUS mus = shrink_formula(top);
			mark_MUS(mus, true, false);
			if(xe.check(mus.bool_mus)){
				blocksUp.push_back(xe.block_up(mus.bool_mus));
				return mus.bool_mus;
			}
			blocksUp.push_back(xe.block_up(mus.bool_mus));
			Formula rec = xor_shrink(original, xe, limit);
			if(!rec.empty()) return rec;		
		}
	}
	std::cout << std::endl;
	return Formula();
}

int Master::bsat_xor(float tresh, XorExplorer &xe){
	int i = 0;
	std::cout << "muses size: " << muses.size() << std::endl;	
	//first check already identifeis MUSes
	for(auto &mus: muses){
		if(xe.check(mus.bool_mus)){
			xe.block_up(mus.bool_mus);
			i++;
		}
		if(i >= tresh){
			std::cout << "exceeds treshold via inhritance" << std::endl;       
			return i;
		}
	}
	std::cout << "inherited MUSes: " << i << std::endl;

	xe.solver->set_default_polarity(true);
	Formula f = xe.solve();
	int valid = 0;
	int invalid = 0;
	int least_distance = 1000000;
	while(i < tresh && !f.empty()){
		assert(xe.check(f));
		xe.block(f);
		Formula original = f;
		if(is_valid(f, true, true)){
			block_down(f);
			blocksDown.push_back(xe.block_down(f));
			valid++;
		}else{
			invalid++;
			MUS mus = shrink_formula(f);
			if(xe.check(mus.bool_mus) || mus.bool_mus == original){
				i++;
			}
			mark_MUS(mus, true, false);
			blocksUp.push_back(xe.block_up(mus.bool_mus));
		}
		f = xe.solve();
		if((valid + invalid) % 200 == 0 || xe.xorSize > 100){
		      	std::cout << "iter: " << (valid + invalid) << ", muses: " << i << ", valid: " << valid << ", invalid: " << invalid;
	       		cout << "least diff: " << least_distance << std::endl;
			cout << "f size: " << count_ones(f) << std::endl;
		}
	}
	std::cout << "bananas: " << (invalid + valid) << ", muses: " << i << ", valid: " << valid << ", invalid: " << invalid << std::endl;
	return i;
}

vector<vector<int>> generate_As(int dimension){
        int m = dimension - 1;
        vector<vector<int>> a(m, vector<int>());
	int min = 1000000;
        for(int i = 0; i < m; i++){
                for(int j = 0; j < dimension; j++){
                        if(random_bool())
                                a[i].push_back(j);
                }
		min = (min < a[i].size())? min : a[i].size();
		a[i].push_back(random_bool()); //this is the initial a0
        }
	std::cout << "min xor size: " << min << std::endl;
        return a;
}

int Master::logSATSearch(vector<vector<int>> &As, int tresh, int mPrev){
	std::cout << "start of binary search" << std::endl;
	int low = 0; 
	int high = dimension - 1;
	int m = mPrev;
	while((low + 1) != high){
		std::cout << "m: " << m << ", l: " << low << ", h: " << high << std::endl;
		XorExplorer xe(dimension, blocksDown, blocksUp);
		xe.add_xor(m, As);
		int y = bsat_xor(tresh, xe);
		std::cout << "log y: " << y << std::endl;
		if(y >= tresh){
			low = m;
		}
		else{
			high = m;
		}
		m = ceil((low + high)/2);					
	}	
	std::cout << "final m: " << m << std::endl;
	return m;
}

int Master::approxMC2Core(float tresh, int &nCells){
	//choose h at random
	vector<vector<int>> As = generate_As(dimension);
/*	XorExplorer xe(dimension, blocksDown, blocksUp);	
	xe.add_xor(dimension - 1, As);
	int y = bsat_xor(tresh, xe);
	std::cout << "precheck y: " << y << std::endl;
	if(y >= tresh){
		nCells = -1;
		return -1;
	}
*/	int mPrev = log2(nCells);
	int m = logSATSearch(As, tresh, mPrev);
	XorExplorer xe2(dimension, blocksDown, blocksUp);	
	xe2.add_xor(m, As);
	nCells = pow(2,m);
	return bsat_xor(tresh, xe2);
}

void Master::initialBooster(){
	Formula whole(dimension, true);
	for(int i = 0; i < dimension; i++){
		whole[i] = false;
		Formula copy = whole;
		if(is_valid(copy,true,true)){
			block_down(copy);
			vector<CMSat::Lit> cl;
			for(int i = 0; i < dimension; i++)
				if(!copy[i])
					cl.push_back(CMSat::Lit(i, false));
			blocksDown.push_back(cl);
		}else{
			MUS mus = shrink_formula(copy);
			mark_MUS(mus, true, false);
			vector<CMSat::Lit> cl;
			for(int i = 0; i < dimension; i++)
				if(mus.bool_mus[i])
					cl.push_back(CMSat::Lit(i, true));
			blocksUp.push_back(cl);
		}
		whole[i] = true;
	}
	std::cout << "unit crits: " << explorer->criticals << std:: endl;
	bsat(200);
	std::cout << "known MUSes after booster: " << muses.size() << std::endl;
}

int Master::counting(float e, float d){
	//	initialBooster();
	float tresh = 1 + 9.84 * (1 + (e / (1 + e)))*(1 + 1/e)*(1 + 1/e);
	cout << "tresh: " << tresh << endl;
	int y = bsat(tresh);
	cout << "y: " << y << endl;
	if(y < tresh){
		cout << "y < tresh" << endl;
		return y;
	}
	int t = ceil(17 * log2(3 / d));
	cout << "t: " << t << endl;
	vector<int> counts;
	int nCells = 2;
	t = 10; //TOOD: just for testing purposes set to 3
	for(int iter = 0; iter < t; iter++){
		cout << "iteration: " << iter << endl;
		counts.push_back(approxMC2Core(tresh, nCells));
		cout << "muses in the last iteration: " << counts.back() << endl;
	}
	return 0 ;
}

