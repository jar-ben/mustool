#include "SimpleExplorer.h"
#include "misc.h"

using namespace Minisat;
using namespace std;

Lit itoLit30(int i){
        bool sign = i < 0;
        int var = (sign)? -i-1 : i-1;
        return (sign) ? ~mkLit(var) : mkLit(var);
}

inline int Littoi30(Lit l) {
    return (var(l)+1) * (sign(l) ? -1 : 1);
}

SimpleExplorer::SimpleExplorer(int v, bool verb){
	vars = v;
	dimension = v;
	solver = new Solver();
        for(int i = 0; i < vars; i++){
                solver->newVar(lbool(uint8_t(2)), true);
	}
	added_clauses = 0;
}

SimpleExplorer::~SimpleExplorer(){
	delete solver;
}

bool SimpleExplorer::block_down(Formula cl){

	//add to SAT unex solver
	added_clauses++;
        vec<Lit> msClause;
        for(int i = 0; i < cl.size(); i++)
                if(!cl[i])
                        msClause.push(mkLit(i));
        return solver->addClause(msClause);
}

bool SimpleExplorer::block_up(Formula mus){
	added_clauses++;
        vec<Lit> msClause;
	for(int i = 0; i < mus.size(); i++)
		if(mus[i])
			msClause.push(~mkLit(i));
        return solver->addClause(msClause);
}


vector<bool> SimpleExplorer::get_bot_unexplored(vector<bool> top){
        solver->rnd_pol = false;
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool(uint8_t(0)));

        vec<Lit> lits;
        for(int i = 0; i < vars; i++)
                if(!top[i])
                        lits.push(~mkLit(i));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
		if(top[i])
             		unexplored[i] = solver->modelValue(i) != l_False;
		else
			unexplored[i] = false;
        }
        return unexplored;
}


vector<bool> SimpleExplorer::get_top_unexplored(vector<bool> bot){
        solver->rnd_pol = false;
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool(uint8_t(1)));

	vec<Lit> lits;
	for(int i = 0; i < vars; i++)
		if(bot[i])
			lits.push(mkLit(i));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
		if(bot[i])
			unexplored[i] = true;
		else
	        	unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}

bool SimpleExplorer::check(vector<bool> f){
	vec<Lit> lits;
	for(int i = 0; i < vars; i++)
		if(f[i])
			lits.push(mkLit(i));
		else
			lits.push(~mkLit(i));

        return solver->solve(lits);
}

std::vector<bool> SimpleExplorer::get_unexplored(vector<int> assumptions){
        solver->rnd_pol = true; //default value is randomly chosen

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit30(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}

std::vector<bool> SimpleExplorer::get_top_unexplored(vector<int> assumptions){
        solver->rnd_pol = false; //default value is not randomly chosen
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool((uint8_t) 1)); //default variable value is true (1)

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit30(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}


std::vector<bool> SimpleExplorer::get_bot_unexplored(vector<int> assumptions){
        solver->rnd_pol = false; //default value is not randomly chosen
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool((uint8_t) 0)); //default variable value is false (0)

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit30(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}


std::vector<bool> SimpleExplorer::get_unexplored(uint8_t polarity, bool rnd_pol){

        solver->rnd_pol = rnd_pol;
	if(!rnd_pol)
	        for(int i = 0; i < vars; i++)
        	        solver->setPolarity(i, lbool(polarity));

        if(!solver->solve())
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}

