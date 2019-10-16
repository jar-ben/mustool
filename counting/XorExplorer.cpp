#include "XorExplorer.h"
#include "core/misc.h"
#include <random>

using namespace CMSat;
using std::vector;

XorExplorer::XorExplorer(int v, vector<vector<Lit>> &blocksDown, vector<vector<Lit>> &blocksUp){
	vars = v;
	dimension = v;
	solver = new SATSolver();
	solver->new_vars(v);
	xorSize = -1;
	for(auto &b: blocksDown)
		solver->add_clause(b);
	for(auto &b: blocksUp)
		solver->add_clause(b);
}

XorExplorer::~XorExplorer(){
	delete solver;
}

std::vector<Lit> XorExplorer::block_up(Formula &f){
	std::vector<Lit> cl;
	for(int i = 0; i < dimension; i++){
		if(f[i])
			cl.push_back(Lit(i,true));
	}
	solver->add_clause(cl);
	return cl;
}

std::vector<Lit> XorExplorer::block_down(Formula &f){
	std::vector<Lit> cl;
	for(int i = 0; i < dimension; i++){
		if(!f[i])
			cl.push_back(Lit(i,false));
	}
	solver->add_clause(cl);
	return cl;
}

bool XorExplorer::block(Formula f){
	vector<Lit> ban_solution;
	for(int i = 0; i < dimension; i++){
		if(f[i])
			ban_solution.push_back(Lit(i, true));
		else
			ban_solution.push_back(Lit(i, false));
	}
	std::cout << "dim: " << dimension << ", nVars: " << solver->nVars() << std::endl;
	solver->add_clause(ban_solution);
}

void XorExplorer::add_xor(int m, std::vector<std::vector<int>> &As){
	for(int i = 0; i < m; i++){
		vector<uint32_t> new_xor;
		for(int j = 0; j < As[i].size() - 1; j++){
			auto c = As[i][j];
			if(c < dimension)
				new_xor.push_back(c);
		}
		bool rhs = As[i].back() == 1;
		solver->add_xor_clause(new_xor, rhs);
	}
	assert(xorSize == -1);
	xorSize = m;
}

Formula randomF(int dimension){
	Formula f(dimension, false);
	for(int i = 0; i < dimension; i++)
		f[i] = random_bool();
	return f;
}

void XorExplorer::checkUniformity(){
	std::cout << "checking uniformity, m = " << xorSize << std::endl;
	int in = 0;
	int rounds = 3000;
	for(int i = 0; i < rounds; i++){
		Formula f = randomF(dimension);
		if(check(f)) in++;
	}
	float exp = (float)1/pow(2,xorSize);
	std::cout << in << " out of " << rounds << ", " << (((float)in)/rounds) << ", expected " << exp << std::endl;
}

bool XorExplorer::check(Formula &f){
	vector<CMSat::Lit> assumptions;
	for(int i = 0; i < dimension; i++){
		if(f[i])
			assumptions.push_back(Lit(i, false));
		else
			assumptions.push_back(Lit(i, true));
	}
	auto ret = solver->solve(&assumptions);
	assert(ret == l_True || ret == l_False);
	return ret == l_True;
}

Formula XorExplorer::solve_subset(Formula &f){
	vector<CMSat::Lit> assumptions;
	for(int i = 0; i < dimension; i++)
		if(!f[i])
			assumptions.push_back(Lit(i,true));
	return solve(assumptions);
}

Formula XorExplorer::solve(vector<CMSat::Lit> assumptions){
	auto ret = solver->solve(&assumptions);

	if(ret == l_False){
		return Formula();
	}
	auto m = solver->get_model();
	Formula f(dimension, false);
	for(int i = 0; i < dimension; i++)
		if(m[i] == l_True)
			f[i] = true;
	return f;
}
