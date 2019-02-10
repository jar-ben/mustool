#include "RotationSolver.h"
#include <iostream>
#include <fstream>

using namespace Minisat;
using namespace std;

void RotationSolver::block_up(std::vector<bool> &f){	
        vec<Lit> msClause;
        for(int i = 0; i < dimension; i++)
                if(f[i])
                        msClause.push(~mkLit(i));
        solver->addClause(msClause);
}

void RotationSolver::add_rotation_seed(std::vector<bool> &f){
	controls++;
	controls_type.push_back(true);
	solver->newVar(lbool(uint8_t(2)), true);
	
	//first part clause_k implies p_k
        vec<Lit> msClause;
        for(int i = 0; i < dimension; i++)
                if(f[i])
                        msClause.push(~mkLit(i));
		else
			msClause.push(mkLit(i));
	msClause.push(mkLit((dimension + controls) -1));
        solver->addClause(msClause);

	//second part, p_k implies clause_k
	for(int i = 0; i < dimension; i++){
		vec<Lit> msClause;
		msClause.push(~mkLit((dimension + controls) -1));
		if(f[i]) 
			msClause.push(mkLit(i));
		else
			msClause.push(~mkLit(i));
		solver->addClause(msClause);
	}	
}

std::vector<bool> union_sets3(std::vector<bool> a, std::vector<bool> b){
	int dimension = a.size();
	std::vector<bool> res(dimension, false);
	for(int i = 0; i < dimension; i++)
		if(a[i] || b[i])	res[i] = true;
	return res;
}

bool is_rotation_pair(std::vector<std::vector<std::vector<int>>> &flip_edges, int c, int c2){
	for(auto &edges: flip_edges[c]){
		for(auto &edge: edges){
			if(edge == c2) return true;
		}
	}			
	return false;
}

void RotationSolver::add_mus(
	std::vector<bool> &mus, 
	std::vector<std::vector<std::vector<int>>> &flip_edges, 
	std::vector<std::vector<int>> &parent_muses, 
	std::vector<std::vector<bool>> &muses,
	std::vector<std::vector<int>> muses_int,
	std::vector<bool> &critical){
	
	std::cout << "adding mus to rotation solver" << std::endl;
	block_up(mus);
	int counter = 0;
	for(int i = 0; i < dimension; i++){
		if(mus[i] && !critical[i]){
			//for(int m_id = 0; m_id < muses.size(); m_id++){
			for(int m_id = muses.size() - 1; m_id >= 0 && m_id >= muses.size() - 11; m_id--){
				auto &m_bool = muses[m_id];
				if(m_bool[i]) continue;
				auto &m_int = muses_int[m_id];
				for(auto &c: m_int){
					counter++;
					if(!mus[c] && is_rotation_pair(flip_edges, c, i)){ //(c,i) is a rotation pair		
						std::vector<bool> united = union_sets3(mus, m_bool);		
						united[i] = united[c] = false;
						add_rotation_seed(united);	
						break;
					}
				}
			}
		}
	}


	std::cout << "counter: " << counter << std::endl;
}

std::vector<bool> RotationSolver::get_seed(){	
        solver->rnd_pol = true; //default value is randomly chosen

	//control variable for the control clause
	controls++;
	controls_type.push_back(false);
	solver->newVar(lbool(uint8_t(2)), true);

	//add the p_1 || p_2 || ... || p_k control clause	
        vec<Lit> msClause;
        for(int i = 0; i < controls_type.size(); i++)
                if(controls_type[i])
                        msClause.push(mkLit(dimension + i));
	msClause.push(mkLit((dimension + controls) - 1));
        solver->addClause(msClause);
	

	vec<Lit> assumption; assumption.push(~mkLit((dimension + controls) - 1)); //the current control clause cannot be satisfied using the kvazi-control lit 
        if(!solver->solve(assumption))
                return vector<bool>();

        vector<bool> unexplored(dimension);
        for (int i = 0 ; i < dimension ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}

RotationSolver::RotationSolver(int dim): dimension(dim){
	solver = new Solver();
        for(int i = 0; i < dimension; i++){
                solver->newVar(lbool(uint8_t(2)), true);
	}
	controls = 0;
}

RotationSolver::~RotationSolver(){
	delete solver;
}
