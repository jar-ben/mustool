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
	while(i < tresh && !top.empty()){
		if(is_valid(top, true, true)){
			block_down(top);
			xe->block_down(top);
		}else{
			Formula original = top;
			MUS mus = shrink_formula(top);
			mark_MUS(mus, true, false);
			if(original == mus.bool_mus)
				i++;
		}
		top = explorer->get_unexplored(1, false);
	}
	return i;
}

int Master::bsat_xor(float tresh){
	Formula f = xe->solve();
	int i = 0;
	while(i < tresh && !f.empty()){
		if(is_valid(f, true, true)){
			block_down(f);
			xe->block_down(f);
		}else{
			Formula original = f;
			MUS mus = shrink_formula(f);
			mark_MUS(mus, true, false);
			if(original == mus.bool_mus)
				i++;
			xe->temporary_block_up(original);
		}
		f = xe->solve();
	}
	std::cout << "iters: " << i << std::endl;
	return i;
}

vector<vector<int>> generate_As(int dimension){
        int m = dimension - 1;
        vector<vector<int>> a(m, vector<int>());
        for(int i = 0; i < m; i++){
                for(int j = 0; j < dimension + 1; j++){
                        if(random_bool())
                                a[i].push_back(j);
                }
        }
        return a;
}

int Master::approxMC2Core(float tresh, int &nCells){
	vector<vector<int>> As = generate_As(dimension);
	vector<bool> alfa(dimension - 1, false);
       	for(int i = 0; i < dimension -1; i++){
		alfa[i] = random_bool();
	}	
	xe->add_xor(dimension - 1, As, alfa);
	int y = bsat_xor(tresh);
	std::cout << "y: " << y << std::endl;	
	if(y >= tresh){
		nCells = -1;
		return -1;
	}

	return 0;
}

int Master::counting(float e, float d){
	float tresh = 1 + 9.84 * (1 + (e / (1 + e)))*(1 + 1/e)*(1 + 1/e);
	cout << "tresh: " << tresh << endl;
//	int y = bsat(tresh);
//	cout << "y: " << y << endl;
//	if(y < tresh){
//		cout << "y < tresh" << endl;
//		return y;
//	}
	int t = ceil(17 * log2(3 / d));
	cout << "t: " << t << endl;
	vector<int> counts;
	int nCells = 2;
	for(int iter = 0; iter < t; iter++){
		cout << "iteration: " << iter << endl;
		counts.push_back(approxMC2Core(tresh, nCells));
		cout << "muses in the last iteration: " << counts.back() << endl;
	}
	return 0 ;
}

