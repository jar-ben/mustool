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
		}
		top = explorer->get_unexplored(1, false);
	}
	return i;
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
	while(i <  tresh){
		Formula f = counting_get_mus(xe);
		if(f.empty()) 
			return i;
		i++;
		validate_mus(f);
		muses.push_back(MUS(f, -1, muses.size(), count_ones(f)));
		mark_MUS(muses.back());
		cout << i << " validated" << endl;
	}
	return i;
}
	
Formula Master::counting_get_mus(XorExplorer &xe){
	Formula whole(dimension, true);
	string constraints_file = "constraints.cnf";
	ofstream cfile;
	cfile.open(constraints_file, ios::out);
	cfile << satSolver->toString(whole);
	cfile.close();

	string unexXor_file = "unex.cnf";
	ofstream file;
	file.open(unexXor_file, ios::out);
	file << "p cnf 0 0\n";
	file << explorer->toCnf();
	file << xe.export_xors();
	file.close();

	stringstream cmd;
	cmd << "python qbf.py " << constraints_file << " " << unexXor_file;
	string result = exec(cmd.str().c_str());

	Formula f(dimension, false);	
	vector<string> lines = getLines(result);
	bool reading = false;
	for(auto &line: lines){
		if(reading){
			vector<string> nums = split(line);
			for(auto n: nums){
				int i = stoi(n);
				if(i == 0) return f;
				if(i < 0) f[(-1 * i) - 1] = false;
				else f[i - 1] = true;
			}
		}
		if(line.find("SOLUTION") != std::string::npos){
			reading = true;
		}		
	}
	return Formula();
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

int Master::counting(float e, float d){
	float tresh = 1 + 9.84 * (1 + (e / (1 + e)))*(1 + 1/e)*(1 + 1/e);
	tresh = 10;
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
	t = 10; //TOOD: just for testing purposes set to 3
	for(int iter = 0; iter < t; iter++){
		cout << "iteration: " << iter << endl;
		counts.push_back(approxMC2Core(tresh, nCells));
		cout << "muses in the last iteration: " << counts.back() << endl;
	}
	return 0 ;
}

