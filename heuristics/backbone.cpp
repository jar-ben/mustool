#include "core/Master.h"
#include "core/misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>
#include "Bones.h"

//return ids of violated clause if there is some
vector<int> Master::propagate_backbone(int i, bool positive, vector<int> &lits_left, vector<bool> &satisfied, vector<bool> &top, vector<int> &singletons, Formula &seed, int c, 	vector<int> &influenced, bool mark_influenced){
	vector<int> violated;
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	if(positive){//implied positive value of i-th literal
		for(auto c1: msSolver->hitmap_pos[i]){ //clauses satisfied by implied positive value of i-th literal
			satisfied[c1] = true;
		}
		for(auto c1: msSolver->hitmap_neg[i]){ //trim clauses that contain negated i-th literal
			if(!top[c1] || satisfied[c1]) continue;
			if(mark_influenced) influenced.push_back(c1);
			lits_left[c1]--;
			if(lits_left[c1] == 1){				
				singletons.push_back(c1);
			}
			if(lits_left[c1] == 0){
				if(seed[c1]) return vector<int> (1, -1 * (c1 + 1) ); 
				else violated.push_back(c1);				
			}
		}
	}else{//implied negative value of i-th literal
		for(auto c1: msSolver->hitmap_neg[i]){ //clauses satisfied by implied negative value of i-th literal
			satisfied[c1] = true;
		}
		for(auto c1: msSolver->hitmap_pos[i]){ //trim clauses that contain positive i-th literal
			if(!top[c1] || satisfied[c1]) continue;
			lits_left[c1]--;
			if(mark_influenced) influenced.push_back(c1);
			if(lits_left[c1] == 1){
				singletons.push_back(c1);
			}
			if(lits_left[c1] == 0){
				if(seed[c1]) return vector<int> (1, -1 * (c1 + 1) ); 
				else violated.push_back(c1);				
			}
		}
	}
	return violated;
}

vector<int> Master::backbone_init(Formula &seed, Formula &implied, Formula &values, vector<int> &lits_left, Formula &satisfied, vector<int> &singletons, Formula &top,
			vector<int> &variables_map, vector<int> &variables_map_inv, int c){
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	if(c >= 0)
		backbone_simplify(seed, c, implied, values);	
	get_backbones_bones(seed, implied, values, variables_map, variables_map_inv);

	//get_backbones(seed, implied, values, variables_map, variables_map_inv);
	int vars = msSolver->vars;
	vector<int> violated;
	for(int i = 0 ; i < dimension; i++){
		lits_left[i] = msSolver->clauses[i].size() - 1; // -1 because there is the activation literal
	}		
	//propagate the backbone literals		
	vector<int> influenced;
	for(int i = 0; i < vars; i++){
		if(implied[i]){			
			vector<int> new_violated = propagate_backbone(i, values[i], lits_left, satisfied, top, singletons, seed, 0, influenced);
			violated.insert(violated.end(), new_violated.begin(), new_violated.end());
		}
	}
	return violated;
}

vector<int> Master::backbone_find_seed(Formula &seed, Formula &implied, Formula &values, Formula satisfied, vector<int> singletons, vector<int> lits_left, Formula &top){
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	int iters = 0;
	int vars = msSolver->vars;
	vector<int> violated;
	Formula seedc = seed;
	Formula impliedc = implied;
	Formula valuesc = values;
	Formula satisfiedc = satisfied;
	vector<int> singletonsc = singletons;
	vector<int> lits_leftc = lits_left;

	//slicing structure
	vector<int> unit_lits(dimension, 0);
	vector<int> added_clauses;

	vector<int> influenced;
	while(!singletons.empty()){
		iters++;
		int singleton = singletons.back();
		singletons.pop_back();
		if(!top[singleton]) continue;
		//find the new implied literal
		int lit = 0;
		for(int lid = 0; lid < msSolver->clauses[singleton].size() - 1; lid++){
			auto l = msSolver->clauses[singleton][lid];
			if( l > 0 && !implied[l - 1] ){	lit = l; break;	}
			else if( l < 0 && !implied[(-1 * l) - 1] ){ lit = l; break; }
		}
		if(lit == 0) continue; // this singleton is already satisfied, i.e. there were two same singletons
		seed[singleton] = true;
		unit_lits[singleton] = lit;
		added_clauses.push_back(singleton);

		int var = (lit > 0)? lit - 1 : (-1 * lit) - 1;
		int positive = lit > 0;
		implied[var] = true;
		values[var] = positive;
		violated = propagate_backbone(var, positive, lits_left, satisfied, top, singletons, seed, singleton, influenced);
		if(!violated.empty()) break;
	}
	vector<bool> extra(added_clauses.size(), false);
	for(int i = added_clauses.size() - 1; i >= 0; i--){
		auto c1 = added_clauses[i];

		bool influence = false;
		//first check whether the clause influence one of the violated clauses
		for(auto c2: violated){			
			if(c2 < 0) c2 = (c2 + 1) * -1;
			for(int lid = 0; lid < msSolver->clauses[c2].size() - 1; lid++){
				auto l = msSolver->clauses[c2][lid];
				if( -1 * l == unit_lits[c1] ){ influence = true; break; }
			}
			if(influence) break;
		}
		if(influence) continue;

		//49 50 51
		//check the other clauses
		for(int j = i + 1; j < added_clauses.size(); j++){
			auto c2 = added_clauses[j];	
			if(!seed[c2]) continue;
			for(int lid = 0; lid < msSolver->clauses[c2].size() - 1; lid++){
				auto l = msSolver->clauses[c2][lid];
				if( -1 * l == unit_lits[c1] ){ influence = true; break; }
			}
			if(influence) break;
		}
		if(!influence){ seed[c1] = false; extra[i] = true; }
	}

	seed = seedc;
	implied = impliedc;
	values = valuesc;
	satisfied = satisfiedc;
	singletons = singletonsc;
	lits_left = lits_leftc;

	//check
	for(int i = 0; i < extra.size(); i++){
		auto c1 = added_clauses[i];

		if(extra[i]) continue; 
		//find the new implied literal
		int lit = 0;
		for(int lid = 0; lid < msSolver->clauses[c1].size() - 1; lid++){
			auto l = msSolver->clauses[c1][lid];
			if( l > 0 && !implied[l - 1] ){	lit = l; break;	}
			else if( l < 0 && !implied[(-1 * l) - 1] ){ lit = l; break; }
		}
		if(lit == 0) continue; // this singleton is already satisfied, i.e. there were two same singletons
		seed[c1] = true;

		int var = (lit > 0)? lit - 1 : (-1 * lit) - 1;
		int positive = lit > 0;
		implied[var] = true;
		values[var] = positive;
		violated = propagate_backbone(var, positive, lits_left, satisfied, top, singletons, seed, c1, influenced);
		if(!violated.empty()) break;
	}
	return violated;
}

void Master::implied_literal(int cl, Formula &implied, int &lit, int &var, bool &value){
	lit = 0;
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	for(int lid = 0; lid < msSolver->clauses[cl].size() - 1; lid++){
		auto l = msSolver->clauses[cl][lid];
		if( l > 0 && !implied[l - 1] ){	lit = l; break;	}
		else if( l < 0 && !implied[(-1 * l) - 1] ){ lit = l; break; }
	}	
	var = (lit > 0)? lit - 1 : (-1 * lit) - 1;
	value = lit > 0;
}

int Master::seek_conflict(vector<int> &unit, vector<bool> &unit_value, vector<int> &singletons, int from, int to, Formula &seed, Formula &top, Formula &implied, vector<int> &added){
	int lit, var;
	bool value;
	for(int i = from; i < to; i++){
		auto s = singletons[i];
		if(!top[s]) continue;
		implied_literal(s, implied, lit, var, value);
		if(lit == 0) continue; //this singleton is already satisfied, there is no unit literal in it
		if(unit[var] >= 0 && unit_value[var] != value ){
			if(!seed[s]){
				if(!seed[unit[var]]){ seed[unit[var]] = true; added.push_back(unit[var]); }
				return s;
			}else if(!seed[ unit[var] ]){
				if(!seed[s]){ seed[s] = true; added.push_back(s); }
				seed[s] = true;
				return unit[var];
			}else return -1;
		}
		unit[var] = s;
		unit_value[var] = value;
	}
	return -2;
}

vector<int> Master::backbone_find_seed_beta(Formula &seed, Formula implied, Formula values, Formula satisfied, vector<int> singletons, vector<int> lits_left, Formula &top){
	// unit[i] = s means that the value of variable i is implied by the singleton s if s >= 0
	// otherwise, the value of variable i is not implied yet.
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	vector<int> unit(msSolver->vars, -1);
	vector<bool> unit_value(msSolver->vars, false);
	vector<int> tmp_added;
	int conflicting = seek_conflict(unit, unit_value, singletons, 0, singletons.size(), seed, top, implied, tmp_added);
	if(conflicting >= -1)	
		return vector<int> (1,conflicting);


	Formula seed_ = seed;
	Formula implied_ = implied;
	Formula values_ = values;
	Formula satisfied_ = satisfied;
	vector<int> singletons_ = singletons;
	vector<int> lits_left_ = lits_left;
	vector<int> unit_ = unit;
	vector<bool> unit_value_ = unit_value;


	int lit, var;
	bool value;
	vector<int> violated;
	vector<vector<int>> influencing(dimension, vector<int>());
	vector<int> added;
	int iter = 0;
	while(!singletons.empty()){
		int singleton = singletons.back();
		singletons.pop_back();
		if(!top[singleton]) continue;
		if(!seed[singleton]) added.push_back(singleton);
		seed[singleton] = true;
		
		implied_literal(singleton, implied, lit, var, value);
		if(lit == 0) continue;
		implied[var] = true;
		values[var] = value;
		int singletons_size = singletons.size();
		violated = propagate_backbone(var, value, lits_left, satisfied, top, singletons, seed, singleton, influencing[singleton], true);
		if(!violated.empty()) break;
		
		int conflicting = seek_conflict(unit, unit_value, singletons, singletons_size, singletons.size(), seed, top, implied, added);
		if(conflicting >= -1){
			violated.push_back(conflicting);
			break;
		}	
	}
	if(violated.empty() || added.empty())
		return violated;
	//slicing
	vector<int> must_stay;
	for(auto c: influencing[ added[added.size() - 1] ])
		if(seed[c])
			must_stay.push_back(c);

	for(int i = added.size() - 2; i >= 0; i--){
		auto c = added[i];
		bool keep = false;
		for(auto influenced: influencing[c]){
			if(seed[influenced] || count(violated.begin(), violated.end(), influenced) > 0){
				keep = true;
				break;
			}
			if(keep) break;
		}
		if(!keep && count(must_stay.begin(), must_stay.end(), c) == 0 ) seed[c] = false;
	}
	return violated;
}

int Master::get_backbones_bones(Formula &seed, Formula &implied, Formula &values, vector<int> &variables_map, vector<int> &variables_map_inv){
	vector<vector<int>> cls;
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	for(int i = 0;  i < dimension; i++){
		if(seed[i]){
			bool sat = false;
			vector<int> lits;
			for(auto lit: msSolver->clauses[i]){
				int var = (lit > 0)? lit : (-1 * lit);
				int phase = (lit > 0)? 1: -1;
				bool value = phase == 1;
				if(var <= msSolver->vars){
					if( implied[var - 1] && values[var - 1] != value) continue;
					else if( implied[var - 1] && values[var - 1] == value){ sat = true; break; }
					else lits.push_back(variables_map[var] * phase);
				}
			}
			if(!sat){
				cls.push_back(lits);
			}
		}
	}

	Bones b = Bones();
	vector<int> backbone_lits = b.get_backbones(cls);
	for(auto lit: backbone_lits){
		int var = (lit > 0)? lit : (-1 * lit);
		int phase = (lit > 0)? 1: -1;
		int olit = variables_map_inv[var] * phase;
		int ovar = olit * phase - 1;
		implied[ovar] = true;
		values[ovar] = phase > 0; 
	}
	return 0;	
}

int Master::get_backbones(Formula &seed, Formula &implied, Formula &values, vector<int> &variables_map, vector<int> &variables_map_inv){

	//export the seed
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	vector<vector<int>> cls;
	for(int i = 0;  i < dimension; i++){
		if(seed[i]){
			bool sat = false;
			vector<int> lits;
			for(auto lit: msSolver->clauses[i]){
				int var = (lit > 0)? lit : (-1 * lit);
				int phase = (lit > 0)? 1: -1;
				bool value = phase == 1;
				if(var <= msSolver->vars){
					if( implied[var - 1] && values[var - 1] != value) continue;
					else if( implied[var - 1] && values[var - 1] == value){ sat = true; break; }
					else lits.push_back(variables_map[var] * phase);
				}
			}
			if(!sat){
				cls.push_back(lits);
			}
		}
	}

	stringstream ss;
	ss << "p cnf " << variables_map_inv.size() << " " << cls.size() << "\n";
	for(auto &cl: cls){
		for(auto l: cl)
			ss << l << " ";
		ss << "0" << "\n";
	}

	ofstream file;
	stringstream outputf;
	outputf << "minibones_export" << hash << "_" << count_ones(seed) << ".cnf"; 
	string filename = outputf.str();
	file.open(filename);
	file << ss.str();
	file.close();

	//compute the backbone (invoke minibones)
	stringstream cmd;
	cmd << "./minibones -e -i -c 100 " << filename << " 2> /dev/null";
	string output = exec(cmd.str().c_str());	
	stringstream oss(output);
	string line;
	while (getline(oss, line, '\n'))
	{
		if (line[0] != 'v')
			continue;
		istringstream is(line);	
		is.ignore(1,' '); //ignore the initial 'v'
		int lit;
		while(is >> lit){
			if(lit == 0) continue;
			int var = (lit > 0)? lit : (-1 * lit);			
			int phase = (lit > 0)? 1: -1;
			int olit = variables_map_inv[var] * phase;
			int ovar = olit * phase - 1;
			implied[ovar] = true;
			values[ovar] = phase > 0; 
		}
	}	
	remove(filename.c_str());
	return 0;
}

void Master::backbone_build_literal_map(Formula &seed, vector<int> &variables_map, vector<int> &variables_map_inv){
	//build the variable ids map
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);//
	vector<bool> variables (msSolver->vars + 1, false);
	for(int i = 0;  i < dimension; i++){
		if(seed[i]){
			for(auto lit: msSolver->clauses[i]){ //there are no trailing zeros in clause[i]		
				int var = (lit > 0)? lit : (-1 * lit);
				if(var <= msSolver->vars)
					variables[var] = true;
			}
		}
	}
	variables_map.resize(variables.size(), -1);
	variables_map_inv.resize(1,0);
	int current = 1;
	for(int i = 0; i < variables.size(); i++){
		if(variables[i]){
			variables_map[i] = current;
			variables_map_inv.push_back(i);
			current++;
		}
	}
}

void Master::backbone_simplify(Formula &seed, int c, Formula &implied, Formula &values){
	vector<int> lits_left (dimension,0);
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	for(int i = 0 ; i < dimension; i++){
		lits_left[i] = msSolver->clauses[i].size() - 1; // -1 because there is the activation literal
	}
	vector<bool> satisfied(dimension, false);
	vector<int> singletons;
	vector<int> influenced;
	for(int lid = 0; lid < msSolver->clauses[c].size() - 1; lid++){
		auto lit = msSolver->clauses[c][lid];
		int var = (lit > 0)? lit - 1 : (-1 * lit) - 1;					
		bool value = !(lit > 0);
		implied[var] = true;
		values[var] = value;
		propagate_backbone(var, value, lits_left, satisfied, seed, singletons, seed, c, influenced);
	}

	int iters = 0;
	while(!singletons.empty()){
		c = singletons.back();
		singletons.pop_back();		
		if(!seed[c]) continue;

		for(int lid = 0; lid < msSolver->clauses[c].size() - 1; lid++){
			auto lit = msSolver->clauses[c][lid];
			if(lit == 0) continue;		
			int var = (lit > 0)? lit - 1 : (-1 * lit) - 1;			
			if(var > msSolver->vars || implied[var]) continue;
			bool value = lit > 0;
			implied[var] = true;
			values[var] = value;
			vector<int> viol = propagate_backbone(var, value, lits_left, satisfied, seed, singletons, seed, c, influenced);
		}
		iters++;
	}
}

void Master::backbone_check_reminder(Formula &implied, Formula &values, Formula &top, Formula &m1){
	if(dimension > 50000) return;
	vector<int> assumptions;
/*	for(int i = 0; i < implied.size(); i++){
		if(implied[i]){
			int lit = ( values[i] )? i + 1 : (i + 1) * -1;
			assumptions.push_back(lit);
		}
	}
*/	//bool sat = satSolver->solve(top, assumptions);
	bool sat = is_valid(top, false, false);
	if(sat){
		block_down(top);
	}
	return;
}

int Master::backbone_mus_rotation(MUS &m1, Formula &top){
	MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
	vector<int> local_muses;
	int vars = msSolver->vars;	
	vector<int> violated;
	vector<int> variables_map;
	vector<int> variables_map_inv;
	backbone_build_literal_map(m1.bool_mus, variables_map, variables_map_inv);
	int iterations = 0;
	Formula original_top = top;
	bool extended = false;
	for(int d = 0; d < m1.without_crits.size(); d++){
		int c = m1.without_crits[ (d + muses.size()) % m1.without_crits.size() ];
		if(explorer->is_critical(c, top)) continue;
		if(iterations++ > 10) break;

		extended = false;
		if(explorer->mus_intersection[c]){
			original_top = top;
			top = Formula(dimension, true);
			top[c] = false;
			extended = true;
		}

		int iter_rot = 0;
		bool found_seed = false;		
		Formula seed = m1.bool_mus;
		seed[c] = false;
		top[c] = false;
		Formula implied(vars, false);
		Formula values(vars, false);
		vector<int> lits_left (dimension,0);
		vector<bool> satisfied(dimension, false);
		vector<int> singletons;
		violated  = backbone_init(seed, implied, values, lits_left, satisfied, singletons, top, variables_map, variables_map_inv, c);//) continue;
		for(auto c2: violated){
			Formula seed_c2 = seed; 
			seed_c2[c2] = true;
			MUS mus = shrink_formula(seed_c2);
			rotated_muses++;
			mark_MUS(mus);
			iter_rot++;
			top[c2] = false; //these muses have the form (m1 - {c}) \cup {c2}, thus removing c2 is the only option how to get "local" unexplored seed
		}
		Formula impliedc = implied;
		Formula valuesc = values;
		int inner_iter = 0;
		while(true){
			seed = m1.bool_mus;
			seed[c] = false;
			violated =  backbone_find_seed_beta( seed, implied, values, satisfied, singletons, lits_left, top );
			for(auto c2: violated){
				Formula seed_c2 = seed;
				if(c2 >= 0){
					seed_c2[c2] = true;
				}
				MUS mus = shrink_formula(seed_c2);
				rotated_muses++;
				mark_MUS(mus);
				iter_rot++;
				if(c2 >= 0 && mus.bool_mus[c2]) top[c2] = false;
				else{//find another c3 in overlap and block it
					for(auto c3: mus.without_crits){ if(!m1.bool_mus[c3]){ top[c3] = false; break; } }
					break; // since this MUS do not contain c2, we cannot be sure that the following seeds are unexplored
				}
			}
			if(violated.empty()) break;//todo: rotate the violated ones, find MHS of local MUSes, update top, and find new violated :]
			implied = impliedc;
			valuesc = values;
		}
		if(iter_rot == 0){ 
			backbone_check_reminder(implied, values, top, m1.bool_mus);
		}
		top[c] = true;
		if(extended){
			for(int i = 0; i < dimension; i++)
				if(!top[i]) original_top[i] = false;
			top = original_top;
		}
	
	}
	return 0;
}

