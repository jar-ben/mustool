#include "core/misc.h"
#include "satSolvers/GlucoseHandle.h"
#include "satSolvers/mcsmus_handle.h"
#include <sstream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <random>
#include <stdlib.h>   
#include <time.h>   
#include <cstdio>
#include <assert.h>

using namespace CustomGlucose;
using namespace std;

Lit GLitoLit(int i){
	bool sign = i < 0;
	int var = (sign)? -i-1 : i-1;
	return (sign) ? ~mkLit(var) : mkLit(var);
}

int GLLittoi(Lit l){
	return (var(l)+1) * (sign(l) ? -1 : 1);
}

vector<int>  GLconvert_clause(string clause){
	vector<int> ret;
	istringstream is(clause);
	int n;
	while(is >> n)
		ret.push_back(n);
	ret.pop_back();
	return ret;
}

GlucoseHandle::GlucoseHandle(string filename):BooleanSolver(filename){
	//solver = new Solver();
	solver = new SimpSolver();
	vars = 0;
	parse(filename);
	dimension = clauses.size();
	srand (time(NULL));
	rotated_crits = 0;
	flip_edges_computed.resize(dimension, false);
	flip_edges.resize(dimension);
	flip_edges_flatten.resize(dimension);
}

GlucoseHandle::~GlucoseHandle(){
	delete solver;
}


bool GlucoseHandle::add_clause(vector<int> cl){
	std::sort(cl.begin(), cl.end());
	vector<int> copy = cl; 
	copy.pop_back(); //get rid of control variable
	clauses.push_back(cl);
	clauses_map[copy] = clauses.size() - 1; //used for manipulation with single MUS extractions (muser2, dmuser)
	vec<Lit> msClause;
	for(auto &lit: cl)
		msClause.push(GLitoLit(lit));		

	for(auto &lit: copy){
		if(lit > 0)
			hitmap_pos[lit - 1].push_back(clauses.size() - 1);
		else
			hitmap_neg[(-1 * lit) - 1].push_back(clauses.size() - 1);
	}
	return solver->addClause(msClause);
}

bool GlucoseHandle::add_unit(int lit){
	return solver->addClause(GLitoLit(lit));
}

bool GlucoseHandle::parse(string path){
        ifstream infile(path, ifstream::in);
        if (!infile.is_open())
		print_err("wrong input file");

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
		else if(clauses_unique_map.find(line) != clauses_unique_map.end()){
			cout << "a duplicate clause found in the input formula" << endl;
			continue;		
		}
		else{
			clauses_str.push_back(line);
			clauses_unique_map[line] = clauses_str.size() - 1;
		}
        }
	hitmap_pos.resize(vars);
	hitmap_neg.resize(vars);
	for(int i = 0; i < vars; i++){
		solver->newVar();	// clause variables
	}
	for(int i = 0; i < clauses_str.size(); i++)
		solver->newVar(true, true);	// control variables
	for(size_t i = 0; i < clauses_str.size(); i++){
		clause = GLconvert_clause(clauses_str[i]);
		clause.push_back(vars + i + 1);	//control variable
		add_clause(clause); //add clause to the solver
	}
	return true;
}

vector<bool> GlucoseHandle::model_extension(vector<bool> subset, vector<bool> model){
	int flipped = 0;
	vector<bool> extension = subset;
	for(int i = 0; i < extension.size(); i++){
		if(!extension[i]){
			for(int j = 0; j < clauses[i].size() - 1; j++){
				int lit = clauses[i][j];
				if(lit > 0 && model[lit - 1]){
					extension[i] = true;
					flipped++;
					break;
				}
				else if(lit < 0 && !model[(-1 * lit) - 1]){
					extension[i] = true;
					flipped++;
					break;
				}				
			}
		}
	}
	return extension;
}


// check formula for satisfiability using miniSAT
// the core and grow variables controls whether to return an unsat core or model extension, respectively
bool GlucoseHandle::solve(vector<bool>& controls, bool unsat_improve, bool sat_improve){
	checks++;
	vec<Lit> lits;
	for(unsigned int i = 0; i < controls.size(); i++){
		if(controls[i])
			lits.push(GLitoLit((i + vars + 1) * (-1)));
		else
			lits.push(GLitoLit(i + vars + 1 ));
	}
	bool sat = solver->solve(lits);
	if(sat && sat_improve){ // extract model extension		
		for(int f = 0; f < controls.size(); f++){
			if(!controls[f]){
				for(int l = 0; l < clauses[f].size() - 1; l++){
					if(clauses[f][l] > 0){
						if(solver->model[clauses[f][l] - 1] == l_True){
							controls[f] = true;
							break;
						}
					}
					else{
						if(solver->model[-1 * clauses[f][l] - 1] == l_False){
							controls[f] = true;
							break;
						}
					}
				}
			}
		}
	}			
	else if(!sat && unsat_improve){ // extract unsat core
		vector<bool> core = vector<bool> (dimension, false);		
	        for (int i = 0 ; i < solver->conflict.size() ; i++) 
			core[var(solver->conflict[i]) - vars] = true;
		controls = core;		

	}				
	return sat;
}

vector<bool> GlucoseHandle::get_model(){
	vector<bool> model(vars, false);
	for(int i = 0; i < vars; i++){
		if(solver->model[i] == l_True)
			model[i] = true;
		else if(solver->model[i] == l_False)
			model[i] = false;
		else
			cout << "literal not evaluated " << i << endl;
	}
	return model;
}

string GlucoseHandle::toString(vector<bool> &f){
	int formulas = std::count(f.begin(), f.end(), true);
	stringstream result;
	result << "p cnf " << vars << " " << formulas << "\n";
	for(int i = 0; i < f.size(); i++)
		if(f[i]){
			result << clauses_str[i] << "\n";
		}
	return result.str();
}

void GlucoseHandle::export_formula_crits(vector<bool> f, string filename, vector<bool> crits){
	int formulas = std::count(f.begin(), f.end(), true);


	FILE *file;
	file = fopen(filename.c_str(), "w");
	if(file == NULL) cout << "failed to open " << filename << ", err: " << strerror(errno) << endl;
	
	fprintf(file, "p gcnf %d %d %d\n", vars, formulas, formulas);
	for(int i = 0; i < f.size(); i++)
		if(f[i] && crits[i]){
			fprintf(file, "{0} %s\n", clauses_str[i].c_str());
		}
	int group = 1;
	for(int i = 0; i < f.size(); i++)
		if(f[i] && !crits[i]){
			fprintf(file, "{%d} %s\n", group++, clauses_str[i].c_str());
		}
	if (fclose(file) == EOF) {
		cout << "error closing file: " << strerror(errno) << endl;
	}	
}

// implementation of the shrinking procedure
// based on the value of basic_shrink employes either muser2 or dmuser
vector<bool> GlucoseHandle::shrink(std::vector<bool> &f, std::vector<bool> crits){
	shrinks++;
	if(shrink_alg == "custom"){
		return SatSolver::shrink(f, crits); //shrink with unsat cores
	}
	return shrink_mcsmus(f, clauses, crits);
}

int GLmuser_output(std::string filename){
	ifstream file(filename, ifstream::in);
	std::string line;
	while (getline(file, line))
	{
		if (line.find("c Calls") != std::string::npos){
			std::string calls = line.substr(line.find(":") + 2, line.size()); // token is "scott"
			return atoi(calls.c_str());
		}
	}	
	return 0;
}
