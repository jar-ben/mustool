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

using namespace std;


BooleanSolver::BooleanSolver(string filename):SatSolver(filename){
}

BooleanSolver::~BooleanSolver(){
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

