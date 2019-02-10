#include "SatSolver.h"

using namespace std;

//basic, domain agnostic, implementation of the shrink procedure
vector<bool> SatSolver::shrink(std::vector<bool> &f, std::vector<bool> crits){
	if(crits.empty())
		crits = std::vector<bool> (dimension, false);
	vector<bool> s = f;
	for(int i = 0; i < dimension; i++){
		if(s[i] && !crits[i]){
			s[i] = false;
			if(solve(s, true))
				s[i] = true;
		}
	}
	return s;
}

//basic, domain agnostic, implementation of the grow procedure
vector<bool> SatSolver::grow(std::vector<bool> &f){
	vector<bool> s = f;
	for(int i = 0; i < dimension; i++){
		if(!s[i]){
			s[i] = true;
			if(!solve(s))
				s[i] = false;
		}
	}
	return s;
}

