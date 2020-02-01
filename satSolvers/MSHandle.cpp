#include "satSolvers/Dimacs.h"
#include "satSolvers/mcsmus_handle.h"
#include "core/misc.h"
#include "satSolvers/MSHandle.h"
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

using namespace CustomMinisat;
using namespace std;

Lit itoLit(int i){
	bool sign = i < 0;
	int var = (sign)? -i-1 : i-1;
	return (sign) ? ~mkLit(var) : mkLit(var);
}

int Littoi(Lit l){
	return (var(l)+1) * (sign(l) ? -1 : 1);
}

vector<int>  convert_clause(string clause){
	vector<int> ret;
	istringstream is(clause);
	int n;
	while(is >> n)
		ret.push_back(n);
	ret.pop_back();
	return ret;
}

MSHandle::MSHandle(string filename):BooleanSolver(filename){
	solver = new Solver();
	vars = 0;
	parse(filename);
	dimension = clauses.size();
	srand (time(NULL));
	rotated_crits = 0;
	flip_edges_computed.resize(dimension, false);
	flip_edges.resize(dimension);
	flip_edges_flatten.resize(dimension);
}

MSHandle::~MSHandle(){
	delete solver;
}


bool MSHandle::add_clause(vector<int> cl){
	std::sort(cl.begin(), cl.end());
	vector<int> copy = cl; 
	copy.pop_back(); //get rid of control variable
	clauses.push_back(cl);
	clauses_map[copy] = clauses.size() - 1; //used for manipulation with single MUS extractions (muser2, dmuser)
	vec<Lit> msClause;
	for(auto &lit: cl)
		msClause.push(itoLit(lit));		

	for(auto &lit: copy){
		if(lit > 0)
			hitmap_pos[lit - 1].push_back(clauses.size() - 1);
		else
			hitmap_neg[(-1 * lit) - 1].push_back(clauses.size() - 1);
	}
	return solver->addClause(msClause);
}

bool MSHandle::add_unit(int lit){
	return solver->addClause(itoLit(lit));
}

bool MSHandle::parse(string path){
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
		solver->newVar(lbool(uint8_t(1)), true);	// control variables
	for(size_t i = 0; i < clauses_str.size(); i++){
		clause = convert_clause(clauses_str[i]);
		clause.push_back(vars + i + 1);	//control variable
		add_clause(clause); //add clause to the solver
	}
	return true;
}

vector<bool> MSHandle::model_extension(vector<bool> subset, vector<bool> model){
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

bool MSHandle::solve(std::vector<bool> &controls, std::vector<int> conflicts, bool unsat_improve, bool sat_improve){
	checks++;
	vec<Lit> lits;
	for(unsigned int i = 0; i < controls.size(); i++){
		if(controls[i])
			lits.push(itoLit((i + vars + 1) * (-1)));
		else{
			lits.push(itoLit(i + vars + 1 ));
		}
	}

	//for(auto c: conflicts){
	//	for(auto l: clauses[c]){
	//		lits.push(itoLit(l)); 
	//	}
	//}

	bool sat = solver->solve(lits);
	if(sat && sat_improve){ // extract model extension		
		int initSize = count_ones(controls);
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
		int finalSize = count_ones(controls);
	}			
	else if(!sat && unsat_improve){ // extract unsat core
		vector<bool> core = vector<bool> (dimension, false);		
	        for (int i = 0 ; i < solver->conflict.size() ; i++) 
			//if(var(solver->conflict[i]) - vars > 0 && var(solver->conflict[i]) - vars < dimension)
				core[var(solver->conflict[i]) - vars] = true;
		controls = core;		

	}				
	return sat;

}

// check formula for satisfiability using miniSAT
// the core and grow variables controls whether to return an unsat core or model extension, respectively
bool MSHandle::solve(vector<bool>& controls, bool unsat_improve, bool sat_improve){
	return solve(controls, vector<int>(), unsat_improve, sat_improve);
}

vector<bool> MSHandle::get_model(){
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

string MSHandle::toString(vector<bool> &f){
	int formulas = std::count(f.begin(), f.end(), true);
	stringstream result;
	result << "p cnf " << vars << " " << formulas << "\n";
	for(int i = 0; i < f.size(); i++)
		if(f[i]){
			result << clauses_str[i] << "\n";
		}
	return result.str();
}

void MSHandle::export_formula_crits(vector<bool> f, string filename, vector<bool> crits){
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

vector<bool> MSHandle::import_formula_crits(string filename){
	vector<bool> f(dimension, false);
	vector<vector<int>> cls;
	ReMUS::parse_DIMACS(filename, cls);
	for(auto cl: cls){
		sort(cl.begin(), cl.end());
		if(clauses_map.count(cl))
			f[clauses_map[cl]] = true;
		else { assert(false); }		
	}	
	return f;
}

// implementation of the shrinking procedure
// based on the value of basic_shrink employes either muser2 or dmuser
vector<bool> MSHandle::shrink(std::vector<bool> &f, std::vector<bool> crits){
	shrinks++;
	if(shrink_alg == "custom"){
		return SatSolver::shrink(f, crits); //shrink with unsat cores
	}
#ifdef UMCSMUS
	if(shrink_alg == "default"){
		return shrink_mcsmus(f, clauses, crits);
	}
#endif
	stringstream exp;			
	exp << "./tmp/f_" << hash << ".cnf";			
	export_formula_crits(f, exp.str(), crits);	

	return shrink_muser(exp.str(), hash);
}

int muser_output(std::string filename){
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

vector<bool> MSHandle::shrink_muser(string input, int hash2){
	stringstream cmd, muser_out, imp;
	muser_out << "./tmp/f_" << hash << "_output";
	imp << "./tmp/f_" << hash << "_mus";
	cmd << "./muser2-para -grp -wf " << imp.str() << " " << input << " > " << muser_out.str();// */ " > /dev/null";
	int status = system(cmd.str().c_str());
	if(status < 0){
		std::cout << "Invalid muser return code" << std::endl; exit(0);
	}
	imp << ".gcnf";
	vector<bool> mus = import_formula_crits(imp.str());
	int sat_calls = muser_output(muser_out.str());	
//	checks += sat_calls;
	remove(imp.str().c_str());
	remove(muser_out.str().c_str());
	remove(input.c_str());
	return mus;
}

std::vector<bool> MSHandle::grow(std::vector<bool> &f, std::vector<bool> conflicts){
	grows++;
	if(grow_alg == "default"){
		return SatSolver::grow(f, conflicts);
	}else if(grow_alg == "uwr"){
		return grow_uwrmaxsat(f, conflicts);
	}else if(grow_alg == "mcsls"){
		std::vector<std::vector<bool>> ms = growMultiple(f, conflicts, 1);
		if(!ms.empty()) // there is a chance that mcsls fail; in such a case, we proceed with grow_cmp
			return ms.back();
	}
	return grow_cmp(f, conflicts);
}

std::vector<bool> MSHandle::grow_cmp(std::vector<bool> &f, std::vector<bool> &conflicts){
	cout << "growing cmp" << endl;
	stringstream ex;			
	ex << "./tmp/f_" << hash << ".wcnf";			
	vector<int> mapa = export_formula_wcnf(f, conflicts, ex.str());
	stringstream cmd;
	cmd << "./cmp_linux -strategyCoMss=" << growStrategy << " " << ex.str();
	string result = exec(cmd.str().c_str());
	
	std::string line;    
	bool reading = false;
	std::vector<int> res_mcs;
	std::istringstream res(result);
	while (std::getline(res, line)) {
		if (line.find("UNSATISFIABLE") != std::string::npos)
			reading = true;			
		else if(reading && line[0] == 'v'){
			line.erase(0, 2);
			res_mcs = convert_clause(line);
		}
	}
	cout << "mcs size: " << res_mcs.size() << endl;
	if(!reading)
		print_err("cmp failed to find a mcs");
	vector<bool> mss(dimension, true);
	for(int i = 0; i < dimension; i++)
		if(conflicts[i])
			mss[i] = false;
	for(auto l:res_mcs){
		mss[mapa[l - 1]] = false;
	}
	return mss;
}

std::vector<bool> MSHandle::satisfied(std::vector<int> &valuation){
	std::vector<bool> result(dimension, false);
	for(auto l: valuation){
		if(l > 0){
			for(auto c: hitmap_pos[l - 1]){
				result[c] = true;
			}
		}else{
			for(auto c: hitmap_neg[(-1 * l) - 1]){
				result[c] = true;
			}
		}	
	}
	return result;
}

std::vector<bool> MSHandle::grow_uwrmaxsat(std::vector<bool> &f, std::vector<bool> &conflicts){
	stringstream ex;			
	ex << "./tmp/uwr_" << hash << ".wcnf";			
	vector<int> mapa = export_formula_wcnf(f, conflicts, ex.str());
	stringstream cmd;
	cmd << "./uwrmaxsat -no-msu -m -v0 " << ex.str();
	string result = exec(cmd.str().c_str());
	
	std::string line;    
	std::vector<int> res_vars;
	std::istringstream res(result);
	while (std::getline(res, line)) {
		if(line[0] == 'v' && line[1] == ' '){
			line.erase(0, 2);
			istringstream is(line);
			int n;
			while(is >> n)
				res_vars.push_back(n);
		}
	}
	return satisfied(res_vars);	
}

std::vector<std::vector<bool>> MSHandle::growMultiple(std::vector<bool> &f, std::vector<bool> conflicts, int limit){
	vector<bool> original = f;
	stringstream ex;			
	ex << "./tmp/mcsls_" << hash << ".wcnf";			
	vector<int> mapa = export_formula_wcnf(f, conflicts, ex.str());
	stringstream cmd;
	cmd << "./mcsls -num " << limit << " " << mcslsArgs << " " << ex.str();
	cout << cmd.str() << endl;
	string result = exec(cmd.str().c_str());
	
	std::string line;    
	std::vector<std::vector<bool>> res_msses;
	std::istringstream res(result);
	while (std::getline(res, line)) {
		if(starts_with(line, "c MCS")){
			vector<bool> mss(dimension, true);
			for(int  i = 0; i < dimension; i++)
				if(conflicts[i]) mss[i] = false;
			line.erase(0, 6);
			istringstream is(line);
			int n;
			while(is >> n)				
				mss[mapa[n - 1]] = false;
		res_msses.push_back(mss);	
		}
	}
	grows += res_msses.size();
	return res_msses;
}

vector<int> MSHandle::export_formula_wcnf(std::vector<bool> f, std::vector<bool> &conflicts, std::string filename){
	int formulas = dimension - std::count(conflicts.begin(), conflicts.end(), true);
	int hard = std::count(f.begin(), f.end(), true);
	int hh = dimension * dimension;

	FILE *file;
	file = fopen(filename.c_str(), "w");
	if(file == NULL) cout << "failed to open " << filename << ", err: " << strerror(errno) << endl;
	vector<int> mapa;	
	fprintf(file, "p wcnf %d %d %d\n", vars, formulas, hh);
	for(int i = 0; i < f.size(); i++){
		if(f[i]){
			fprintf(file, "%d %s\n", hh, clauses_str[i].c_str());
			mapa.push_back(i);
		}else if(!f[i] && !conflicts[i]){
			fprintf(file, "1 %s\n", clauses_str[i].c_str());
			mapa.push_back(i);
		}
	}
	if (fclose(file) == EOF) {
		cout << "error closing file: " << strerror(errno) << endl;
	}
	return mapa;	
}

