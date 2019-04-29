#include "SatSolver.h"
#include <fstream>

using namespace std;

//basic, domain agnostic, implementation of the shrink procedure
vector<bool> SatSolver::shrink(std::vector<bool> &f, std::vector<bool> crits){
	shrinks++;
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
	grows++;
	vector<bool> s = f;
	for(int i = 0; i < dimension; i++){
		if(!s[i]){
			s[i] = true;
			if(!solve(s, false, true))
				s[i] = false;
		}
	}
	return s;
}

void SatSolver::exportMUS(std::vector<bool> mus, std::string outputFile){
	exported_muses++;
	ofstream file;
	file.open(outputFile, std::ios_base::app);
	file << "MUS #" << exported_muses << "\n";
	file << toString(mus) << "\n";
	file.close();
}
