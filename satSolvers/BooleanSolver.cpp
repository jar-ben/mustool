#include "core/misc.h"
#include "satSolvers/BooleanSolver.h"
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
#include <queue>

using namespace std;


BooleanSolver::BooleanSolver(string filename):SatSolver(filename){
}

BooleanSolver::~BooleanSolver(){
}

void print_clause(vector<int> cl){
	cout << "clause: ";
	for(auto l: cl)
		cout << l << " ";
	cout << endl;
}

void print_model(vector<int> model){
	for(auto l: model)
		cout << l << " ";
	cout << endl;
}


bool BooleanSolver::lit_occurences(vector<bool> subset, int c2){
	for(auto l: clauses[c2]){
		int var = (l > 0)? l : -l;
		if(var > vars) break; //the activation literal
		int count = 0;
		if(l > 0){
			for(auto c1: hitmap_pos[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				if(subset[c1] && c1 != c2){ count++; }
			}
		} else {
			for(auto c1: hitmap_neg[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				if(subset[c1] && c1 != c2){ count++; }
			}
		}
		if(count == 0) return false;
	}
	return true;
}


vector<bool> BooleanSolver::propagateToUnsat(vector<bool> base, vector<bool> cover, vector<int> implied){
	queue<int> toPropagate;
	for(auto l: implied)
		toPropagate.push(l);
	vector<bool> conflict = base;
	vector<int> lits_left (dimension,0);
	for(int i = 0 ; i < dimension; i++){
		lits_left[i] = clauses[i].size() - 1; // -1 because there is the activation literal
	}
	vector<bool> satisfied(dimension, false);
	vector<bool> setup(vars + 1, false);
	vector<bool> value(vars + 1, false);
	while(!toPropagate.empty()){
		auto l = toPropagate.front();
		toPropagate.pop();
		int var = (l > 0)? l : -l;
		if(setup[var])
			continue;
		setup[var] = true;
		value[var] = l < 0;
		if(l < 0){
			for(auto c1: hitmap_neg[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				satisfied[c1] = true;
			}
			for(auto c1: hitmap_pos[var - 1]){ //trim clauses that contain positive i-th literal
				if(satisfied[c1] || !cover[c1]) continue;
				lits_left[c1]--;
				if(lits_left[c1] == 1){
					for(auto lit: clauses[c1]){
						if(lit >= vars) break; //the activation literal
						int var2 = (lit > 0)? lit : -lit;
						if(!setup[var2]){
							toPropagate.push(lit);
							implied.push_back(lit);
							break;
						}
					}
					conflict[c1] = true;
				}
				if(lits_left[c1] == 0){
					conflict[c1] = true;
					return conflict;
				}
			}
		}else{
			for(auto c1: hitmap_pos[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				satisfied[c1] = true;
			}
			for(auto c1: hitmap_neg[var - 1]){ //trim clauses that contain positive i-th literal
				if(satisfied[c1] || !cover[c1]) continue;
				lits_left[c1]--;
				if(lits_left[c1] == 1){
					for(auto lit: clauses[c1]){
						if(lit >= vars) break; //the activation literal
						int var2 = (lit > 0)? lit : -lit;
						if(!setup[var2]){
							toPropagate.push(lit);
							implied.push_back(lit);
							break;
						}
					}
					conflict[c1] = true;
				}
				if(lits_left[c1] == 0){
					conflict[c1] = true;
					return conflict;
				}
			}
		}
	}
	return vector<bool>();
}
vector<int> BooleanSolver::get_implied(vector<bool> mus, int c){
	if(!mus[c]) print_err("c is not in mus");
	mus[c] = false;
	vector<int> lits_left (dimension,0);
	for(int i = 0 ; i < dimension; i++){
		lits_left[i] = clauses[i].size() - 1; // -1 because there is the activation literal
	}
	vector<bool> satisfied(dimension, false);
	queue<int> to_propagate;
	vector<bool> setup(vars + 1, false);
	vector<bool> value(vars + 1, false);
	vector<int> model;
	for(auto l: clauses[c]){
		int var = (l > 0)? l : -l;
		if(var <= vars){
			to_propagate.push(-l);
			model.push_back(-l);
		}
	}
	while(!to_propagate.empty()){
		auto l = to_propagate.front();
		to_propagate.pop();
		int var = (l > 0)? l : -l;
		if(setup[var])
			continue;
		setup[var] = true;
		value[var] = l < 0;
		model.push_back(l);
		if(l < 0){
			for(auto c1: hitmap_neg[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				satisfied[c1] = true;
			}
			for(auto c1: hitmap_pos[var - 1]){ //trim clauses that contain positive i-th literal
				if(!mus[c1] || satisfied[c1]) continue;
				lits_left[c1]--;
				if(lits_left[c1] == 1){
					for(auto lit: clauses[c1]){
						int var2 = (lit > 0)? lit : -lit;
						if(!setup[var2]){
							to_propagate.push(lit);
							break;
						}
					}
				}
				if(lits_left[c1] == 0){
					print_err("conflict during propagation");
				}
			}
		}else{
			for(auto c1: hitmap_pos[var - 1]){ //clauses satisfied by implied negative value of i-th literal
				satisfied[c1] = true;
			}
			for(auto c1: hitmap_neg[var - 1]){ //trim clauses that contain positive i-th literal
				if(!mus[c1] || satisfied[c1]) continue;
				lits_left[c1]--;
				if(lits_left[c1] == 1){
					for(auto lit: clauses[c1]){
						int var2 = (lit > 0)? lit : -lit;
						if(!setup[var2]){
							to_propagate.push(lit);
							break;
						}
					}
				}
				if(lits_left[c1] == 0){
					print_err("conflict during propagation");
				}
			}
		}
	}
	return model;
}

void BooleanSolver::compute_flip_edges(int c){
	if(flip_edges_computed.empty()){
		flip_edges_computed.resize(dimension, false);
		flip_edges.resize(dimension);
		flip_edges_flatten.resize(dimension);
	}
	if(flip_edges_computed[c]) return;
	flip_edges_computed[c] = true;

	vector<bool> flatten(dimension, false);
	for(int l = 0; l < clauses[c].size() - 1; l++){
		auto lit = clauses[c][l];
		vector<int> edges;
		if(lit > 0){
			for(auto &h: hitmap_neg[lit - 1]){
				if(h != c){
					edges.push_back(h);
					flatten[h] = true;
				}
			}
		}
		else{
			for(auto &h: hitmap_pos[(-1 * lit) - 1]){
				if(h != c){
					edges.push_back(h);			
					flatten[h] = true;
				}
			}
		}
		flip_edges[c].push_back(edges);
	}
	for(int i = 0; i < dimension; i++)
		if(flatten[i])
			flip_edges_flatten[c].push_back(i);
}

void BooleanSolver::criticals_rotation(vector<bool>& criticals, vector<bool> subset){
	vector<int> criticals_int;
	for(int i = 0; i < dimension; i++)
		if(criticals[i] && subset[i]) criticals_int.push_back(i);

	for(int i = 0; i < criticals_int.size(); i++){
		int c = criticals_int[i];
		compute_flip_edges(c); //TODO: better encansulape
		for(auto &lit_group: flip_edges[c]){
			int count = 0;			
			int flip_c;	
			for(auto c2: lit_group){
				if(subset[c2]){ count++; flip_c = c2; }
			}
			if(count == 1 && !criticals[flip_c]){
				criticals_int.push_back(flip_c);
				criticals[flip_c] = true;
			}
		}
	}
}

vector<int> BooleanSolver::critical_extension_clause(vector<bool> &f, vector<bool> &crits, int c){
	vector<int> extensions;
	for(auto l: clauses[c]){
		vector<int> hits;
		int var = (l > 0)? l - 1 : (-l) - 1;
		if(var >= vars) break; //these are the control variables
		if(l > 0){
			for(auto c2: hitmap_neg[var]){
				if(f[c2])
					hits.push_back(c2);
			}
		}else{
			for(auto c2: hitmap_pos[var]){
				if(f[c2])
					hits.push_back(c2);
			}
		}
		if(hits.size() == 1 && !crits[hits[0]]){
			extensions.push_back(hits[0]);
			crits[hits[0]] = true;
		}
	}
	return extensions;
}

int BooleanSolver::critical_extension(vector<bool> &f, vector<bool> &crits){
	int initSize = count_ones(crits);
        queue<int> S;
        for(int i = 0; i < dimension; i++){
                if(crits[i]){
                        S.push(i);
                }
        }

        while(!S.empty()){
                int c = S.front();
                S.pop();
		for(auto e: critical_extension_clause(f, crits, c))
			S.push(e);
        }
	return count_ones(crits) - initSize;
}

int BooleanSolver::critical_propagation(vector<bool> &f, vector<bool> &crits, int cl){
        int extensions = 0;
	queue<int> S;
	S.push(cl);
        while(!S.empty()){
		extensions++;
                int c = S.front();
                S.pop();
		for(auto e: critical_extension_clause(f, crits, c))
			S.push(e);
        }
	return extensions - 1; //exclude the initical clause cl
}

//basic, domain agnostic, implementation of the shrink procedure
vector<bool> BooleanSolver::shrink(std::vector<bool> &f, std::vector<bool> crits){
        shrinks++;
        if(crits.empty())
                crits = std::vector<bool> (dimension, false);
        vector<bool> s = f;
	int extensions = 0;
        for(int i = 0; i < dimension; i++){
                if(s[i] && !crits[i]){
                        s[i] = false;
                        if(solve(s, true, false)){
                                s[i] = true;
				crits[i] = true;
				extensions += critical_propagation(f, crits, i);
			}
                }
        }
	cout << "crit extensions: " << extensions << endl;
        return s;
}

vector<bool> BooleanSolver::shrink(std::vector<bool> &f, Explorer *e, std::vector<bool> crits){
        shrinks++;
        if(crits.empty())
                crits = std::vector<bool> (dimension, false);
        vector<bool> s = f;
	int extensions = 0;
        for(int i = 0; i < dimension; i++){
                if(s[i] && !crits[i]){
                        s[i] = false;
                        if(solve(s, true, false)){
                                s[i] = true;
				crits[i] = true;
				extensions += critical_propagation(f, crits, i);
			}
                }
        }
	cout << "crit extensions: " << extensions << endl;
        return s;
}

