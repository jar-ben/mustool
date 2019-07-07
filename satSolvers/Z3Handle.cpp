#include "Z3Handle.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>
#include "core/misc.h"
using namespace z3;

void Z3Handle::addExpr(expr &e){
	std::string control_name = "ceasar" + std::to_string(clauses.size()); //name and create control variables
	expr control = ctx.bool_const(control_name.c_str());
	expr impl = implies(control, e);
	controls.push_back(control);
	implications.push_back(impl);
	clauses.push_back(e);
}


void Z3Handle::parseMultiple(Z3_ast_vector &a_vec, int size){
	for(int i = 0; i < size; i++){
		expr e(ctx, Z3_ast_vector_get(ctx, a_vec, i));
		addExpr(e);
	}
}

void Z3Handle::parseAndFormula(Z3_ast_vector &a_vec){
	Z3_ast a = Z3_ast_vector_get(ctx, a_vec, 0);
	expr e(ctx, a);
	for(int i = 1; i < e.num_args(); i++){
		expr ei = e.arg(i);
		addExpr(ei);
	}
}

Z3Handle::Z3Handle(std::string filename): SatSolver(filename){
	Z3_ast_vector a_vec = ctx.parse_file(filename.c_str());
	int size = Z3_ast_vector_size(ctx, a_vec);
	if(size == 1)
		parseAndFormula(a_vec);
	else if(size > 1)
		parseMultiple(a_vec, size);
	else{
		std::cout << "cannot parse the input file" << std::endl;
		exit(1);
	}
	dimension = clauses.size();
	s = new solver(ctx);
	for(auto &impl: implications)
		s->add(impl);
}


std::string Z3Handle::toString(std::vector<bool> &f){	
	solver sol(ctx);
	for(int i = 0; i < f.size(); i++){
		if(f[i]){
			sol.add(clauses[i]);
		}
	}
	return sol.to_smt2();
}

// implements the shrinking procedure
// exploits ability of z3 solver to return unsat cores
std::vector<bool> Z3Handle::shrink(std::vector<bool> &formula, std::vector<bool> crits){
	if(crits.empty())
		crits = std::vector<bool> (dimension, false);	
	std::vector<bool> shrinked = formula;
	int cores = 3;
	for(int i = 0; i < formula.size(); i++){
		if(shrinked[i] && !crits[i]){
			shrinked[i] = false;
			checks++;
			shrinked[i] = solve(shrinked, cores > 0, false);
			cores--;			
		}
	}
	return shrinked;
}

// check formula for satisfiability
// the core and grow variables controls whether to return an unsat core or model extension, respectively
bool Z3Handle::solve(std::vector<bool> &formula, bool core, bool grow){
	checks++;
	bool result;
	std::vector<expr> assumptions;
	for(int i = 0; i < formula.size(); i++)
		if(formula[i])
			assumptions.push_back(controls[i]);
		else
			assumptions.push_back(!controls[i]); 
	switch (s->check(assumptions.size(), &assumptions[0])) {
		case unsat:   result = false; break;
		case sat:     result = true; break;
		case unknown: result = false; break;
	}	

	if(!result && core){ // extract unsat core
		std::fill(formula.begin(), formula.end(), false);
		auto unsat_core = s->unsat_core();
		for(int i = 0; i < unsat_core.size(); i++){
        		std::stringstream fname; fname << unsat_core[i];
			std::string name = fname.str();
			name.erase(0,6);
			int id = std::stoi(name);
			formula[id] = true;
		}		
	}
	if(result && grow){ // extract model extension
		int count = 0;
		model m = s->get_model();	
		for(int i = 0; i < clauses.size(); i++){
			if(!formula[i]){
       				std::stringstream res_s; res_s << m.eval(clauses[i]);
				std::string res = res_s.str();
				if(res == "true"){
					formula[i] = true;
					count++;
				}
			}
		}
	}
	return result;
}

Z3Handle::~Z3Handle(){
	delete s;
}

