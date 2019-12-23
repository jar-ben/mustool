#include "CadicalHandle.h"
#include "core/misc.h"
#include <sstream>
#include <iostream>

using namespace std;


CadicalHandle::CadicalHandle(string filename):BooleanSolver(filename){
        solver = new CaDiCaL::Solver;
        vars = 0;
        parse(filename);
        dimension = clauses_str.size();
}

bool CadicalHandle::parse(string path){
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
                        is >> pom;      // p
                        is >> pom;      // cnf
                        is >> vars;     //number of variables
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
        for(size_t i = 0; i < clauses_str.size(); i++){
                clause = convert_clause(clauses_str[i]);
                clauses.push_back(clause);
		clause.push_back(vars + i + 1); //control variable
                for(auto l: clause){
			solver->add(l);
		}
		solver->add(0);
        }
	solver->write_dimacs("dim.cnf");
        return true;
}

vector<int> CadicalHandle::convert_clause(string clause){
        vector<int> ret;
        istringstream is(clause);
        int n;
        while(is >> n)
                ret.push_back(n);
	ret.pop_back();
        return ret;
}

CadicalHandle::~CadicalHandle(){
}

string CadicalHandle::toString(vector<bool> &mus){
	stringstream result;
	for(int i = 0; i < dimension; i++){
		if(mus[i]){
			result << clauses_str[i] << "\n";
		}
	}
	return result.str();
}

bool CadicalHandle::solve(std::vector<bool> &f, bool shrink, bool grow){
	checks++;
	for(int i = 0; i < dimension; i++){
		if(!f[i]){
//			solver->assume(vars + i + 1);
		}
		else{
			solver->assume((vars + i + 1) * -1);
		}
	}
	bool sat = solver->solve() == 10;
	if(sat && grow){
		int initSize = count_ones(f);
		for(int i = 0; i < f.size(); i++){
                        if(!f[i]){
                                for(int l = 0; l < clauses[i].size() - 1; l++){
	  				if(clauses[i][l] > 0){
						if(solver->val(clauses[i][l]) > 0){                  				
		  					f[i] = true;
                                                        break;
                                                }
                                        }
                                        else{
						if(solver->val(-1 * clauses[i][l]) < 0){
							f[i] = true;
                                                        break;
                                                }
                                        }
                                }
                        }
                }
		int finalSize = count_ones(f);
	}
	if(!sat && shrink){
	}
	return sat;
}
