#include "MSSSolver.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <algorithm>

using namespace Minisat;
using namespace std;

void print_formula3(std::vector<bool> f){
	for(int i = 0; i < f.size(); i++)
		std::cout << int(f[i]);
	std::cout << std::endl;
}

Lit itoLit3(int i){
	bool sign = i < 0;
	int var = (sign)? -i-1 : i-1;
	return (sign) ? ~mkLit(var) : mkLit(var);
}

int Littoi3(Lit l){
	return (var(l)+1) * (sign(l) ? -1 : 1);
}

vector<int>  convert_clause3(string clause){
	vector<int> ret;
	istringstream is(clause);
	int n;
	while(is >> n)
		ret.push_back(n);
	ret.pop_back();
	return ret;
}

MSSSolver::MSSSolver(){
	solver = new Solver();
        solver->phase_saving = 2;
	vars = 0;
}

MSSSolver::~MSSSolver(){
	delete solver;
}

bool MSSSolver::block_down(vector<bool> cl){
        vec<Lit> msClause;
        for(int i = 0; i < cl.size(); i++)
                if(!cl[i])
                        msClause.push(mkLit(vars + i));
        return solver->addClause(msClause);
}

bool MSSSolver::block_up(vector<bool> cl){
        vec<Lit> msClause;
        for(int i = 0; i < cl.size(); i++)
                if(cl[i])
                        msClause.push(~mkLit(vars + i));
        return solver->addClause(msClause);
}

bool MSSSolver::add_clause(vector<int> cl){
	clauses.push_back(cl);
	vec<Lit> msClause;
	for(vector<int>::iterator it = cl.begin(); it != cl.end(); it++)
		msClause.push(itoLit3(*it));
	return solver->addClause(msClause);
}

bool MSSSolver::add_unit(int lit){
	return solver->addClause(itoLit3(lit));
}

bool MSSSolver::parse_dimacs(string path){
        ifstream infile(path, ifstream::in);
        if (!infile.is_open())
               std:: cout << endl << "err parse_dimacs" << endl;

        vector<string> clauses_str;
        string line;
	vector<int> clause;
	string pom;
        while (getline(infile, line))
        {
                if (line[0] == 'p'){
			istringstream is(line);
			is >> pom;	// p
			is >> pom;	// cnf
			is >> vars;	//number of variables	
		}
                else if(line[0] == 'c')
			continue;
		else
			clauses_str.push_back(line);
        }
	for(int i = 0; i < vars; i++){
		solver->newVar(lbool(uint8_t(0)), true);	// clause variables
	}
	for(int i = 0; i < clauses_str.size(); i++)
		solver->newVar(lbool(uint8_t(0)), true);	// control variables
	for(size_t i = 0; i < clauses_str.size(); i++){
		clause = convert_clause3(clauses_str[i]);
		clause.push_back(-1 * (vars + i + 1));	//control variable
		add_clause(clause);
	}
	return true;
}

vector<bool> MSSSolver::solve(vector<bool> &f){

        solver->rnd_pol = false; //default value is not randomly chosen
        for(int i = 0; i < vars + f.size(); i++)
                solver->setPolarity(i, lbool((uint8_t) 1)); //default variable value is true (1)

        vec<Lit> lits;
	for(unsigned int i = 0; i < f.size(); i++)
		if(f[i])
			lits.push(itoLit3(vars + i + 1 ));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> mss(f.size());
        for (int i = 0 ; i < f.size() ; i++) {
             mss[i] = solver->modelValue(vars + i) != l_False;
        }	

        return mss;
}

