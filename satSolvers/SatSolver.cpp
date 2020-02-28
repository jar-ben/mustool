#include "SatSolver.h"
#include <fstream>

using namespace std;


bool SatSolver::solve(std::vector<bool> &f, std::vector<int> conflicts, bool core, bool grow){
	return solve(f, core, grow);
}

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

vector<bool> SatSolver::shrink(std::vector<bool> &f, Explorer *e, std::vector<bool> crits){
	return shrink(f, crits);
}


//basic, domain agnostic, implementation of the grow procedure
vector<bool> SatSolver::grow2(std::vector<bool> &f, std::vector<std::vector<bool>> &unsats, std::vector<bool> conflicts){
	grows++;
	if(conflicts.empty())
		conflicts = std::vector<bool> (dimension, false);
	vector<bool> s = f;
	for(int i = 0; i < dimension; i++){
		if(!s[i] && !conflicts[i]){
			s[i] = true;
			std::vector<bool> copy = s;
			if(!solve(copy, true, true)){
				unsats.push_back(s);
				s[i] = false;
			}else{
				s = copy;
			}
		}
	}
	return s;
}


//basic, domain agnostic, implementation of the grow procedure
vector<bool> SatSolver::grow(std::vector<bool> &f, std::vector<bool> conflicts){
	grows++;
	if(conflicts.empty())
		conflicts = std::vector<bool> (dimension, false);
	vector<bool> s = f;
	for(int i = 0; i < dimension; i++){
		if(!s[i] && !conflicts[i]){
			s[i] = true;
			if(!solve(s, false, true))
				s[i] = false;
		}
	}
	return s;
}

vector<vector<bool>> SatSolver::growMultiple(std::vector<bool> &f, std::vector<bool> conflicts, int limit){
	vector<vector<bool>> mcses;
	mcses.push_back(grow(f, conflicts));
	return mcses;
}

void SatSolver::exportMUS(std::vector<bool> mus, std::string outputFile){
	exported_muses++;
	ofstream file;
	file.open(outputFile, std::ios_base::app);
	file << "MUS #" << exported_muses << "\n";
	file << toString(mus) << "\n";
	file.close();
}

void SatSolver::exportMSS(std::vector<bool> mss, std::string outputFile){
	exported_msses++;
	ofstream file;
	file.open(outputFile, std::ios_base::app);
	file << "MSS #" << exported_msses << "\n";
	file << toString(mss) << "\n";
	file.close();
}
