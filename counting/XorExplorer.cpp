#include "XorExplorer.h"
#include "core/misc.h"
#include <random>

using namespace CMSat;
using std::vector;

XorExplorer::XorExplorer(int v, bool verb){
	vars = v;
	dimension = v;
	solver = new SATSolver();
	solver->new_vars(v);
	iter = 0;
}

XorExplorer::~XorExplorer(){
	delete solver;
}
bool XorExplorer::block_up(MUS& mus){
	return true;
}

bool XorExplorer::block_down(Formula f){
	std::vector<Lit> cl;
	for(int i = 0; i < dimension; i++){
		if(!f[i])
			cl.push_back(Lit(i,false));
	}
	solver->add_clause(cl);
	return true;
}

bool XorExplorer::temporary_block_up(Formula f){
	vector<Lit> ban_solution;
	for(int i = 0; i < dimension; i++){
		if(f[i])
			ban_solution.push_back(Lit(i, true));
		else
			ban_solution.push_back(Lit(i, false));
	}
	ban_solution.push_back(Lit(dimension + iter - 1, false));
	solver->add_clause(ban_solution);
}

void XorExplorer::add_xor(int m, std::vector<std::vector<int>> As, vector<bool> alfa){
	std::cout << "adding xors" << std::endl;
	iter++;
	solver->new_vars(1);
	for(int i = 0; i < m; i++){
		vector<uint32_t> new_xor;
		for(auto c: As[i]){
			if(c < m)
				new_xor.push_back(c);
		}
		bool flip = As[i].back() == 1;
		flip ^=alfa[i];
		new_xor.push_back(dimension + iter - 1); //the control variable
		solver->add_xor_clause(new_xor, !flip);
		if(i % 50 == 0)
			std::cout << "i: " << i << std::endl;
	}
	std::cout << "added xors" << std::endl;	
}

Formula XorExplorer::solve(){
	vector<CMSat::Lit> assumptions;
	assumptions.push_back(Lit(dimension + iter - 1, false)); //force the last xor to hold
	auto ret = solver->solve(&assumptions);

	if(ret == l_False)
		return Formula();
	auto m = solver->get_model();
	Formula f(dimension, false);
	for(int i = 0; i < dimension; i++)
		if(m[i] == l_True)
			f[i] = true;
	return f;
}
