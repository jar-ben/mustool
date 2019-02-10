#include "Z3Handle.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <regex>
#include <random>
using namespace z3;

//from z3 import *


void print_bit(std::vector<bool> &f){
	for(int i = 0; i < f.size(); i++)
		std::cout << ((f[i])? "1" : "0");
	std::cout << std::endl;
}

void trim(std::string &f){
	f.erase( std::remove(f.begin(), f.end(), '\n'), f.end() );
	std::regex pattern("  ");
	f = std::regex_replace(f, pattern, "");	
}


Z3Handle::Z3Handle(std::string filename): SatSolver(filename){
	Z3_ast a = Z3_parse_smtlib2_file(ctx, filename.c_str(), 0, 0, 0, 0, 0, 0);    
	expr e(ctx, a);

	if(e.num_args() == 0){
		std::cout << "z3 init error" << std::endl;
		exit(1);
	}

	expr control = ctx.bool_const("ceasar0");
	expr impl = implies(control, e.arg(0));
	controls.push_back(control);
	implications.push_back(impl);
	clauses.push_back(e.arg(0));
	
	s = new solver(ctx);
	s->add(impl);

	for(int i = 1; i < e.num_args(); i++){
		std::string control_name = "ceasar" + std::to_string(i); //name and create control variables
		control = ctx.bool_const(control_name.c_str());		
		impl = implies(control, e.arg(i)); //the need to satisfy individual formulas is controled by the control variables
		s->add(impl);
		controls.push_back(control);
		implications.push_back(impl);
		clauses.push_back(e.arg(i));
	}
	dimension = clauses.size();

	//store individual formulas for possible output of MUSes to a file
	std::stringstream cl;
	std::string pom;
	for(int i = 0; i < e.num_args(); i++){
		cl.str( std::string() );
		cl.clear();
		cl << e.arg(i);
		pom = cl.str();
		trim(pom);
		clauses_string.push_back(pom);
		clauses_string_map[pom] = clauses_string.size() - 1;	
	}	
}

// implements the shrinking procedure
// exploits ability of z3 solver to return unsat cores
std::vector<bool> Z3Handle::shrink(std::vector<bool> &formula, std::vector<bool> crits){
	if(crits.empty())
		crits = std::vector<bool> (dimension, false);	
	std::vector<bool> shrinked = formula;
	for(int i = 0; i < formula.size(); i++){
		if(shrinked[i] && !crits[i]){
			shrinked[i] = false;
			checks++;
			shrinked[i] = solve(shrinked, true, false);				
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
		formula.resize(formula.size(), false);
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

