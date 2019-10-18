#include "Explorer.h"
#include "misc.h"
#include <sstream>

using namespace CustomMinisat;
using namespace std;

Lit itoLit2(int i){
        bool sign = i < 0;
        int var = (sign)? -i-1 : i-1;
        return (sign) ? ~mkLit(var) : mkLit(var);
}

inline int Littoi2(Lit l) {
    return (var(l)+1) * (sign(l) ? -1 : 1);
}

Explorer::Explorer(int v, bool verb){
	vars = v;
	dimension = v;
	solver = new Solver();
	botSolver = new Solver();
	topSolver = new Solver();
        for(int i = 0; i < vars; i++){
                solver->newVar(lbool(uint8_t(2)), true);
		botSolver->newVar(lbool(uint8_t(0)), true);
		topSolver->newVar(lbool(uint8_t(1)), true);	
	}

	//init local critical structures
	critical.resize(dimension, false);
	maybe_critical.resize(dimension, false);
	parent_mcses.resize(dimension, std::vector<int>());
	criticals = 0;

	//init local available structrures
	maybe_available.resize(dimension, false);
	parent_muses.resize(dimension, std::vector<int>());
	mus_intersection.resize(dimension, true);
	mus_union.resize(dimension, false);
}

Explorer::~Explorer(){
	delete solver;
	delete botSolver;
	delete topSolver;
}

string Explorer::toCnf(){
	stringstream result;
	for(auto &mcs: mcses){
		for(auto l: mcs)	
			result << (l + 1) << " ";
		result << "0\n";
	}
	for(auto &mus: muses){
		for(auto l: mus.int_mus)
			result << (-1 * (l + 1)) << " ";
		result << "0\n";
	}
	return result.str();
}

bool Explorer::is_critical(int c, std::vector<bool> &subset){
	if(critical[c]) return true;
	if(maybe_critical[c]){		
		for(auto &mcs: parent_mcses[c]){ //"mcs" is an id of a mcs
			bool crit = true;
			for(auto &c2: mcses[mcs]){
				if(c2 != c && subset[c2]){
					crit = false;
					break;
				}
			}
			if(crit) return true;
		}
	}	
	return false;	
}

bool Explorer::is_available(int c, std::vector<bool> &subset){
	if(!maybe_available[c]) return true;	

	for(auto &mus: parent_muses[c]){ //"mus" is an id of a mus
		bool covered = false;
		for(auto &c2: muses[mus].int_mus){
			if(c2 != c && !subset[c2]){
				covered = true;
				break;
			}
		}
		if(!covered) return false;
	}	
	return true;	
}

std::vector<bool> Explorer::get_unexplored(std::vector<bool> top, std::vector<bool> &mus){
	std::vector<int> top_int;
	for(int i = 0; i < dimension; i++)
		if(!top[i]) top_int.push_back(i);

	bool found = false;
	for(int i = 0; i < dimension; i++){
		if(mus[i] && !is_critical(i, top)){
			top[i] = false;
			found = true;	
			break;
		}
	}

	if(!found)
		return std::vector<bool>();
	return top;
}

bool seeds_sort (const std::pair<std::vector<bool>, int> &s1, const std::pair<std::vector<bool>, int> &s2){ 
	return std::count(s1.first.begin(), s1.first.end(), true) < std::count(s2.first.begin(), s2.first.end(), true);
} 

bool Explorer::block_down(Formula cl){
	//add to local critical structures
	std::vector<int> mcs; 
        for(int i = 0; i < cl.size(); i++)
                if(!cl[i])
			mcs.push_back(i);
	int mcs_id = mcses.size();
	mcses.push_back(mcs);
	if(mcs.size() == 1){
		critical[mcs[0]] = true;
		criticals++;
	}
	for(auto &c: mcs){
		maybe_critical[c] = true;
		parent_mcses[c].push_back(mcs_id);
	}
	

	//add to SAT unex solver
        vec<Lit> msClause;
        for(int i = 0; i < cl.size(); i++)
                if(!cl[i])
                        msClause.push(mkLit(i));
        return solver->addClause(msClause) && botSolver->addClause(msClause);
}

bool Explorer::block_up(MUS& mus){
	//add to local available structures
	std::vector<int> mus_no_crit;
	for(auto c: mus.int_mus)
		if(!critical[c])
			mus_no_crit.push_back(c);

	mus.without_crits = mus_no_crit;
	muses.push_back(mus);
	for(auto &c: mus.int_mus){
		if(critical[c]) continue;
		maybe_available[c] = true;
		parent_muses[c].push_back(mus.id);
	}

	for(int i = 0; i < dimension; i++){
		if(!mus.bool_mus[i]) mus_intersection[i] = false;
		else mus_union[i] = true;
	}

        vec<Lit> msClause;
	for(auto c: mus.without_crits)
		msClause.push(~mkLit(c));
        return solver->addClause(msClause) && topSolver->addClause(msClause);
}


vector<bool> Explorer::get_bot_unexplored(vector<bool> top){
        botSolver->rnd_pol = false;
        for(int i = 0; i < vars; i++)
                botSolver->setPolarity(i, lbool(uint8_t(0)));

        vec<Lit> lits;
        for(int i = 0; i < vars; i++)
                if(!top[i])
                        lits.push(~mkLit(i));

        if(!botSolver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
		if(top[i])
             		unexplored[i] = botSolver->modelValue(i) != l_False;
		else
			unexplored[i] = false;
        }
        return unexplored;
}


vector<bool> Explorer::get_top_unexplored(vector<bool> bot){
        topSolver->rnd_pol = false;
        for(int i = 0; i < vars; i++)
                topSolver->setPolarity(i, lbool(uint8_t(1)));

	vec<Lit> lits;
	for(int i = 0; i < vars; i++)
		if(bot[i])
			lits.push(mkLit(i));

        if(!topSolver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
		if(bot[i])
			unexplored[i] = true;
		else
	        	unexplored[i] = topSolver->modelValue(i) != l_False;
        }
        return unexplored;
}

std::vector<bool> Explorer::get_unexplored(vector<int> assumptions){
        solver->rnd_pol = true; //default value is randomly chosen

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit2(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}

std::vector<bool> Explorer::get_top_unexplored(vector<int> assumptions){
        solver->rnd_pol = false; //default value is not randomly chosen
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool((uint8_t) 1)); //default variable value is true (1)

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit2(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}


std::vector<bool> Explorer::get_bot_unexplored(vector<int> assumptions){
        solver->rnd_pol = false; //default value is not randomly chosen
        for(int i = 0; i < vars; i++)
                solver->setPolarity(i, lbool((uint8_t) 0)); //default variable value is false (0)

        vec<Lit> lits;
	for(auto &assumption: assumptions)
		lits.push(itoLit2(assumption));

        if(!solver->solve(lits))
                return vector<bool>();

        vector<bool> unexplored(vars);
        for (int i = 0 ; i < vars ; i++) {
             unexplored[i] = solver->modelValue(i) != l_False;
        }
        return unexplored;
}


std::vector<bool> Explorer::get_unexplored(uint8_t polarity, bool rnd_pol){

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

bool Explorer::isUnexplored(vector<bool> valuation){
        vec<Lit> lits;
        for(unsigned int i = 0; i < dimension; i++){
                if(valuation[i])
                        lits.push(itoLit2(i + 1 ));
                else
                        lits.push(itoLit2((i + 1) * (-1)));

        }
        return solver->solve(lits);
}

bool Explorer::isUnexploredSat(vector<bool> f){
	for(auto &mcs: mcses){
		if(is_disjoint(f, mcs))
			return false;
	}
	return true;
}

int Explorer::getImplied(std::vector<bool>& implied, std::vector<bool>& f){
        for(int i = 0; i < dimension; i++){
		if(f[i] && is_critical(i, f))
			implied[i] = true;
	}
	return 0;
}
