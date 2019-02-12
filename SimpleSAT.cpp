#include "SimpleSAT.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <random>
#include <stdlib.h>   
#include <time.h>   


using namespace CustomMinisat;
using namespace std;

Lit itoLitS(int i){
	bool sign = i < 0;
	int var = (sign)? -i-1 : i-1;
	return (sign) ? ~mkLit(var) : mkLit(var);
}

int LittoiS(Lit l){
	return (var(l)+1) * (sign(l) ? -1 : 1);
}

SimpleSAT::SimpleSAT(vector<vector<int>> cls, int v){
	solver = new Solver();
	dimension = cls.size();
	vars = v;

	for(int i = 0; i < vars; i++)
		solver->newVar();
	for(int i = 0; i < dimension; i++)
		solver->newVar(lbool(uint8_t(1)), true);	// control variables

	for(int i = 0; i < dimension; i++){
		auto &cl = cls[i];
		cl.push_back(vars + i + 1);
		add_clause(cl);
	}
}

SimpleSAT::~SimpleSAT(){
	delete solver;
}

bool SimpleSAT::add_clause(vector<int> cl){
	vec<Lit> msClause;
	for(auto &lit: cl)
		msClause.push(itoLitS(lit));		
	return solver->addClause(msClause);
}

bool SimpleSAT::add_unit(int lit){
	return solver->addClause(itoLitS(lit));
}

void SimpleSAT::get_unsat_core(std::vector<bool> &core){
	for (int i = 0 ; i < solver->conflict.size() ; i++){ 
		int index = var(solver->conflict[i]) - vars;
		if(index >= 0)
			core[var(solver->conflict[i]) - vars] = true;
	}	
};


// check formula for satisfiability using miniSAT
// the core and grow variables controls whether to return an unsat core or model extension, respectively
bool SimpleSAT::solve(vector<bool>& controls, bool unsat_improve, bool sat_improve){
	checks++;
	vec<Lit> lits;
	for(unsigned int i = 0; i < controls.size(); i++){
		if(controls[i])
			lits.push(itoLitS((i + vars + 1) * (-1)));
		else
			lits.push(itoLitS(i + vars + 1 ));
	}
	bool sat = solver->solve(lits);
	if(!sat && unsat_improve){ // extract unsat core
		controls.resize(dimension, false);
	        for (int i = 0 ; i < solver->conflict.size() ; i++){ 
			if(var(solver->conflict[i]) - vars <= dimension && var(solver->conflict[i]) - vars >= 0)
				controls[var(solver->conflict[i]) - vars] = true;
		}
	}				
	return sat;
}



