#include "Explorer.h"
#include "misc.h"

using namespace Minisat;
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
	verbose = false;
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
	bot_calls = top_calls = 0;
	use_implications = true;
	rotation_edge_uses = 0;
	rotated_unsat_seeds = 0;
	added_clauses = 0;
	reduced_unsat = 0;
	mus_based_round = 0;
	rotation_pairs = 0;
 	mus_size_limit = 10000;

	//init local critical structures
	critical.resize(dimension, false);
	maybe_critical.resize(dimension, false);
	parent_mcses.resize(dimension, std::vector<int>());
	found_conflicting = 0;
	criticals = 0;

	//init local available structrures
	maybe_available.resize(dimension, false);
	parent_muses.resize(dimension, std::vector<int>());
	sat_rotated = 0;
	seed_priority_queue = false;

	minimum_duration = 1000;
	maximum_duration = 0;
	total_duration = 0;
	durations = 0;
	total_iters = 0;
	digits = 10;
	int pom = dimension;
	while(pom > 10){
		digits *= 10;
		pom /= 10;
	}
	mus_intersection.resize(dimension, true);
	mus_union.resize(dimension, false);
}

Explorer::~Explorer(){
	delete solver;
	delete botSolver;
	delete topSolver;
}

/*
std::vector<bool> Explorer::get_unex_impl_sat(std::vector<bool> top, std::vector<bool> &searchspace){
	//get new unexplored subset
	int covering_lit = -1;
	Implications copy;
	vector<bool> new_top;
	int i = (impl.positive_clauses.size() + impl.negative_clauses.size()) % dimension;
	for(int j = 0; j < dimension; j++){
		i = (i + 1) % dimension;
		if(!top[i] && !impl.implied[i] && searchspace[i]){ //da se pridat, t.j. neni v top a neni -crit
			new_top = top; new_top[i] = true;
			copy = impl;

			bool covered = true; //are all clauses covered?
			for(auto &cl: copy.negative_in[i]){
				bool already = false;
				for(auto &lit: copy.negative_clauses[cl]){
					if(!new_top[lit]){ already = true; break; }
				}
				if(already) continue;
				bool cov = false;
				for(auto &lit: copy.negative_clauses[cl]){
					if(lit != i && !copy.isCrit(lit, new_top)){
						new_top[lit] = false;
						cov = true; break;
					}
				}
				if(!cov){ covered = false; break; } // this clause is not covered, i.e. the i-th literal cannot be added to top
			}
			if(covered){ covering_lit = i; break; }
		}
	}	

	if(covering_lit == -1){ //covering literal not found
		cout << "covering literal not found" << endl;	
		return vector<bool>();
	}
	//maximize_impl_unex(new_top, copy, searchspace);
	return new_top;
}
*/

bool Explorer::is_hitting_pair(int c1, int c2){
	for(auto l1: satSolver->clauses[c1]){
		for(auto l2: satSolver->clauses[c2]){
			if(l1 == (-1 * l2)) return true;
		}
	}
	return false;
}

bool Explorer::is_backbone_pair(int c2, vector<bool> &implied, vector<bool>& values){
	for(auto l: satSolver->clauses[c2]){
		if(l == 0 || l >= implied.size()) continue;
		int var = (l > 0)? l - 1 : ((l + 1) * -1);
		bool value = l > 0;
		if(implied[var] && (values[var] == value)){
			return true;
		}
	}
	return false;
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

//post extension
int Explorer::sat_rotation(vector<bool>& subset, vector<bool>& model, std::vector<std::vector<int>>& clauses, int &limit){
	int rotated = 0;
	if(limit < 0) return rotated;
	int rotated_maximal = 0;
	std::vector<int> mcs;
	for(int i = 0; i < dimension; i++)
		if(!subset[i]) mcs.push_back(i);
	for(auto &c: mcs){
		satSolver->compute_flip_edges(c); //TODO: better encansulape
		for(int i = 0; i < satSolver->flip_edges[c].size(); i++){ //iterate over literals in c
			int literal = clauses[c][i];
			int lit_index = (literal > 0)? (literal - 1) : ((-1 * literal) - 1);				
			model[lit_index] = !model[lit_index]; //rotate model
			vector<int> &literal_edges = satSolver->flip_edges[c][i]; //edges grouped by literals
			vector<int> unsat_clauses;
			for(auto &c2: literal_edges){ //individual edges
				bool sat = false;
				for(int j = 0; j < clauses[c2].size() - 1; j++){ // be carefull of the "- 1", there is a control literal at the end of each clause
					int lit = clauses[c2][j];		
					if(lit != literal * -1 && ((lit > 0 && model[lit - 1]) || (lit < 0 && !model[(-1 * lit) - 1]))){
						sat = true;
						break;
					}
				}
				if(!sat){ unsat_clauses.push_back(c2); }
			}
			std::vector<int> temporal_c2;
			for(auto &c2: mcs){ //the rotation might satisfy more constraints than just c
				if(c2 == c) continue;
				for(int j = 0; j < clauses[c2].size() - 1; j++){ // be carefull of the "- 1", there is a control literal at the end of each clause
					int lit = clauses[c2][j];		
					if(lit == literal){ subset[c2] = true; temporal_c2.push_back(c2); break; }
				}
			}
			for(auto &c2: unsat_clauses) subset[c2] = false;
			bool covered = true;
			for(auto &c2: unsat_clauses){ // check if the subset is unexplored
				for(auto &mcs_id: parent_mcses[c2]){
					covered = false;
					for(auto &c3: mcses[mcs_id]){
						if(subset[c3] || c3 == c){ covered = true; break; }
					}
					if(!covered) break;
				}
				if(!covered) break;
			}
			if(covered && !unsat_clauses.empty()){
//				cout << "mcs: " << mcs.size() << ", rotated: " << std::count(subset.begin(), subset.end(), false) << endl;
//				cout << "HU" << endl;
				subset[c] = true;				
				if(check_maximality(subset)){
					block_down(subset);
//					cout << "HA" << endl;
					limit--;
					sat_rotated++;
					rotated_maximal++;	
					rotated++;
					rotated += sat_rotation(subset, model, clauses, limit);
					if(limit < 0) return rotated;
				}
				subset[c] = false;
			}
			//undo the rotation
			for(auto &c2: unsat_clauses) subset[c2] = true;
			model[lit_index] = !model[lit_index];
			for(auto c2: temporal_c2) subset[c2] = false;
		}	
	}
	return rotated;
}

bool Explorer::check_maximality(std::vector<bool> &unex){
	std::vector<int> avail;
	for(int i = 0; i < dimension; i++)
		if(!unex[i])
			avail.push_back(i);
	for(auto &c: avail)
		if(is_available(c, unex))
			return false;
	return true;
}

bool Explorer::maximize_unex(std::vector<bool> &unex, int c1, int c2){
	std::vector<int> avail;
	for(int i = 0; i < dimension; i++)
		if(!unex[i] && i != c1 && i != c2)
			avail.push_back(i);
	for(auto &c: avail)
		if(is_available(c, unex))
			unex[c] = true;

	if(c1 >= 0  && is_available(c1, unex)){
		cout << "conflict not preserved" << endl;
		unex[c1] = true;
		return false;
	}
	if(c2 >= 0  && is_available(c2, unex)){
		cout << "conflict not preserved" << endl;
		unex[c2] = true;
		return false;
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
//	maximize_unex(top);
	return top;
}

bool seeds_sort (const std::pair<std::vector<bool>, int> &s1, const std::pair<std::vector<bool>, int> &s2){ 
	return std::count(s1.first.begin(), s1.first.end(), true) < std::count(s2.first.begin(), s2.first.end(), true);
} 

vector<bool> Explorer::init_rotation_single(Formula &top, MUS &m1){
	std::vector<int> top_int;
	for(int i = 0; i < dimension; i++)
		if(!top[i]) top_int.push_back(i);

	int c1, c2;
	for(auto &mus: muses){
		if(!easy_to_shrink(mus)) continue;
		int counter = 0;
		for(auto &c: top_int){
			if(mus.bool_mus[c]){
				if(++counter > 1) break;
				c2 = c;
			}
		}
		if(counter == 1){
			for(auto &c1: m1.int_mus){
				if(!mus.bool_mus[c1] && is_hitting_pair(c2, c1)){
					vector<bool> seed = union_sets(mus.bool_mus, m1.int_mus); seed[c1] = seed[c2] = false;
					top[c1] = false; top[c2] = false;
					return seed;
				}
			}
		}		
	}
	return vector<bool>();
}

//returns a vector of pairs seed, c2
std::vector<std::pair<std::vector<bool>, int>> Explorer::init_rotation(Formula &top, MUS &m2){
	std::vector<std::pair<std::vector<bool>, int>> seeds;
	std::vector<int> top_int;
	for(int i = 0; i < dimension; i++)
		if(!top[i]) top_int.push_back(i);

	std::vector<int> m2_int;
	for(auto c: m2.int_mus)
		if(!critical[c]) m2_int.push_back(c);	

	std::vector<bool> m1;
	int c1;
	bool found = false;
	int options = 0;
	for(auto &mus: muses){
		if(!easy_to_shrink(mus)) continue;
		int counter = 0;
		for(auto &c: top_int){
			if(mus.bool_mus[c]){
				if(++counter > 1) break;
				c1 = c;
			}
		}
		if(counter == 1){
			int current_rotation_pairs = 0;
			for(auto &c: m2_int){
				if(!mus.bool_mus[c] && is_hitting_pair(c1, c)){
					std::vector<bool> seed = union_sets(mus.bool_mus, m2_int); seed[c1] = seed[c] = false;
					auto seed_pair = std::make_pair(seed,c);
					seeds.push_back(seed_pair);
					//return seeds;
					found = true;
					if(!use_implications) break;
					current_rotation_pairs++;
				}
			}
			if(current_rotation_pairs > rotation_pairs) rotation_pairs = current_rotation_pairs;
		}
		if(!use_implications && found) break;
	}
	std::sort(seeds.begin(), seeds.end(), seeds_sort);
	return seeds;
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
		//flip_edges[mcs[0]].clear(); //since mcs[0] is critical, it cannot be flipped (ommited) during mus rotation
	}
	for(auto &c: mcs){
		maybe_critical[c] = true;
		parent_mcses[c].push_back(mcs_id);
	}
	

	//add to SAT unex solver
	added_clauses++;
        vec<Lit> msClause;
        for(int i = 0; i < cl.size(); i++)
                if(!cl[i])
                        msClause.push(mkLit(i));
        return solver->addClause(msClause) && botSolver->addClause(msClause);
}

bool Explorer::block_up_unex(std::vector<int> &mus){
	//add to SAT unex solver
	added_clauses++;
        vec<Lit> msClause;
	for(auto c: mus)
		msClause.push(~mkLit(c));
        return solver->addClause(msClause) && topSolver->addClause(msClause);
}

bool Explorer::block_up(MUS& mus, bool block_unex){
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

	if(mus.duration > 0 && mus.duration < minimum_duration)
		minimum_duration = mus.duration;
	if(mus.duration > maximum_duration)
		maximum_duration = mus.duration;
	if(mus.duration > 0){
		durations++;
		total_duration += mus.duration;
	}

	for(int i = 0; i < dimension; i++){
		if(!mus.bool_mus[i]) mus_intersection[i] = false;
		else mus_union[i] = true;
	}

	rotation_start.push_back(muses.size() - 1);
	if(seed_priority_queue)
		add_seeds_to_queue(mus);

	if(block_unex)
		return block_up_unex(mus_no_crit);
	return true;
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

std::vector<bool> Explorer::get_co_unexplored(){
//								top = union_sets(muses[m1_id].bool_mus, muses[m2_id].int_mus);
	if(muses.size() < 2) return std::vector<bool>();
	auto &m1 = muses[muses.size() - 1];
	auto &m2 = muses[muses.size() - 2];
	vector<int> assumptions;
	bool found = false;
	for(auto c1: m1.int_mus){
		if(!m2.bool_mus[c1]){
			for(auto c2: m2.int_mus){
				if(!m1.bool_mus[c2]){
					if(is_hitting_pair(c1,c2)){
						assumptions.push_back(c1);
						assumptions.push_back(c2);
						found = true; break;
					}
			
				}
			}
			if(found) break;
		}
	}
	if(found){
		cout << "found A" << endl;
		std::vector<bool> top = get_top_unexplored(assumptions);
		if(is_available(assumptions[0], top) || is_available(assumptions[1], top))
			return std::vector<bool>();
		cout << "found B" << endl;
		if(rotate_co(top))
			return top;
		else
			return std::vector<bool>();
	}
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

bool Explorer::checkValuation(vector<bool> valuation){
        vec<Lit> lits;
        for(unsigned int i = 0; i < dimension; i++){
                if(valuation[i])
                        lits.push(itoLit2(i + 1 ));
                else
                        lits.push(itoLit2((i + 1) * (-1)));

        }
        return solver->solve(lits);
}

bool Explorer::checkValuationTopSolver(vector<bool> valuation){
        vec<Lit> lits;
        for(unsigned int i = 0; i < dimension; i++){
                if(valuation[i])
                        lits.push(itoLit2(i + 1 ));
                else
                        lits.push(itoLit2((i + 1) * (-1)));

        }
        return topSolver->solve(lits);
}

bool Explorer::checkValuationBotSolver(vector<bool> valuation){
        vec<Lit> lits;
        for(unsigned int i = 0; i < dimension; i++){
                if(valuation[i])
                        lits.push(itoLit2(i + 1 ));
                else
                        lits.push(itoLit2((i + 1) * (-1)));

        }
        return botSolver->solve(lits);
}

int Explorer::get_implies(std::vector<bool>& implied, std::vector<bool>& f){
        for(int i = 0; i < dimension; i++){
		if(f[i] && is_critical(i, f))
			implied[i] = true;
	}
	return 0;

	vec<Lit> assumptions;
        for(int i = 0; i < dimension; i++){
                if(!f[i])
                        assumptions.push(itoLit2((i + 1) * (-1)));
	}
        vec<Lit> outvec;
        solver->implies(assumptions, outvec, true);
        int impl = 0;
	int new_crits = 0;
        for (int i = 0 ; i < outvec.size(); i++) {
		if(Littoi2(outvec[i]) > 0){
			if(implied[Littoi2(outvec[i]) - 1] != true)
				new_crits++;
			implied[Littoi2(outvec[i]) - 1] = true;
			impl++;
		}
        }
	int f_size =  std::count(f.begin(), f.end(), true);
	//std::cout << (f_size - impl) << ", size: " << f_size << ", implied: " << impl << ", new crits: " << new_crits << std::endl;
        return impl;
}

bool Explorer::is_unexplored(std::vector<bool> seed){
	int seed_size = std::count(seed.begin(), seed.end(), true);
	for(auto &mus: muses){
		if(mus.dimension > seed_size) continue;
		bool unexplored = false;
		for(auto &c: mus.int_mus){
			if(!seed[c]){ unexplored = true; break; }
		}
		if(!unexplored) return false;
	}
	return true;
}

bool Explorer::is_unexplored(std::vector<bool> &m1, std::vector<bool> &m2, int c1, int c2){
	for(auto &mus: muses){
		bool unexplored = false;
		for(auto &c: mus.without_crits){
			if((!m1[c] && !m2[c]) || c == c1 || c == c2 ){ unexplored = true; break; }
		}
		if(!unexplored) return false;
	}
	return true;
}

bool Explorer::rotate_co(std::vector<bool>& top){
	//try to find rotation pair
	vector<int> co;
	for(int i = 0; i < dimension; i++)
		if(!top[i])
			co.push_back(i);
	for(int i = 0; i < co.size(); i++){
		if(maybe_critical[i]) continue;
		for(int j  = i + 1; j < co.size(); j++){
			if(maybe_critical[j]) continue;
			int c1 = co[i]; int c2 = co[j];
			if(is_hitting_pair(c1, c2)){
				//found two corresponding MUSes to reduce the top
				for(auto &m1_id: parent_muses[c1]){
					if(muses[m1_id].bool_mus[c2]) continue;				
					for(auto &m2_id: parent_muses[c2]){
						if(muses[m2_id].bool_mus[c1]) continue;					
						bool ok = true;
						for(auto c: co){
							if(c != c1 && c != c1 && (muses[m1_id].bool_mus[c] || muses[m2_id].bool_mus[c]))
								ok = false;
						}
						if(ok){
							cout << "ok " << m1_id << " " << m2_id << endl;		
							top = union_sets(muses[m1_id].bool_mus, muses[m2_id].int_mus);
							top[c1] = top[c2] = false;
							if(is_hitting_pair(c1,c2)) cout << "it is a hitting pair..." << endl;
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

std::vector<bool> Explorer::get_seed_explicit(){
	int dequeued = 0;
	while(!seeds_queue.empty()){
		dequeued++;
		auto seed = seeds_queue.top();
		seeds_queue.pop();
		auto m1_id = std::get<0>(seed);
		auto m2_id = std::get<1>(seed);
		auto c1 = std::get<2>(seed);
		auto c2 = std::get<3>(seed);
		if(is_unexplored(muses[m1_id].bool_mus,muses[m2_id].bool_mus, c1, c2)){
			auto united = union_sets(muses[m1_id].bool_mus, muses[m2_id].int_mus);		
			united[c1] = united[c2] = false;
//			int size = std::count(united.begin(), united.end(), true);
//			cout << "diff: " << (size - muses[m1_id].dimension) << " " << (size - muses[m2_id].dimension) << endl;			
//			cout << "dequeued: " << dequeued << endl;
			return united;
		}
	}
	return std::vector<bool>();
}

bool Explorer::easy_to_shrink(MUS &mus){
	float limit = (total_duration / durations) * 3;
	return mus.duration < limit;
}

/*
void Explorer::add_seeds_to_queue(MUS &m1){
	if(!easy_to_shrink(m1)) return;
	
	int new_seeds = 0;
	for(auto &m2: muses){
		if(!easy_to_shrink(m2)) continue;
		vector<int> candidates;
		for(auto c2: m2.without_crits){
			if(!m1.bool_mus[c2])
				candidates.push_back(c2);
		}
		for(auto c1: m1.without_crits){
			if(!m2.bool_mus[c1]){
				for(auto c2: candidates){
					if(is_hitting_pair(c1, c2)){
						//int size = m1.dimension + m2.dimension; //cheap to compute approximation
						int size = 0;
						for(int i = 0; i < dimension; i++)
							if(m1.bool_mus[i] || m2.bool_mus[i])
								size++;
						auto seed = std::make_tuple(m1.id, m2.id, c1, c2, size);
						seeds_queue.push(seed);
						if(++new_seeds > 100) return;
					}
				}
			}
		}
	}
}*/

void Explorer::add_seeds_to_queue(MUS &m1){
	//first check whether we should even add the seed
	if(!easy_to_shrink(m1)) return;
	
	int new_seeds = 0;
	for(auto c1: m1.without_crits){
		satSolver->compute_flip_edges(c1); //TODO: better encansulape
		for(auto c2: satSolver->flip_edges_flatten[c1]){
			bool found = false;
			if(m1.bool_mus[c2]) continue;
			for(auto &m2_id: parent_muses[c2]){
				auto &m2 = muses[m2_id];
				if(!easy_to_shrink(m2)) continue;
				if(m2.bool_mus[c1]) continue;
				//int size = m1.dimension + m2.dimension; //cheap to compute approximation
				int size = 0;
				for(int i = 0; i < dimension; i++)
					if(m1.bool_mus[i] || m2.bool_mus[i])
						size++;
				auto seed = std::make_tuple(m1.id, m2.id, c1, c2, size);
				seeds_queue.push(seed);
				found = true;
				if(++new_seeds > 100) return;					
			}
			if(found) break;
		}
	}
}

/*void Explorer::add_seeds_to_queue(MUS &m1){
	//first check whether we should even add the seed
	if(!easy_to_shrink(m1)) return;
	
	int limit = muses.size() - 2;
	int new_seeds = 0;
	for(int m2_id = limit; m2_id > 0; m2_id--){
		MUS& m2 = muses[m2_id];
		if(!easy_to_shrink(m2)) continue;
		for(auto &l: m1.int_mus){
			bool found = false;
			if(m2.bool_mus[l]) continue;
			for(auto &c: m2.without_crits){
				if(!m1.bool_mus[c] && is_hitting_pair(c, l)){ //(c,l) is a rotation pair							
					int size = m1.dimension + m2.dimension;
					auto seed = std::make_tuple(m1.id, m2.id, c, l, size);
					seeds_queue.push(seed);
					new_seeds++;
					found = true;
					break;
				}
			}
			if(found) break;
		}
		if(new_seeds > 100) break;
	}	
	cout << "queue size: " << seeds_queue.size() << endl;
}*/
