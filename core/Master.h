#ifndef REMUS_MASTER_H
#define REMUS_MASTER_H

#include "satSolvers/MSHandle.h"
#include "Explorer.h"
#ifndef NOSMT
	#include "satSolvers/Z3Handle.h"
#endif
#ifndef NOLTL 
	#include "satSolvers/SpotHandle.h"
	#include "satSolvers/NuxmvHandle.h"
#endif
#include "types_h.h"
#include <set>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <ctime>
#include <chrono>	
#include <unordered_map>

using namespace std;

class Master{
public:
	int dimension; //the number of constraints
	string algorithm; //MUS enumeration algorithm to be used
	int isValidExecutions;
	bool verbose; //TODO: make it int based 
	string output_file;
	bool validate_mus_c;
	float mus_approx;
	float crits_treshold;
	bool model_rotation;
	bool get_implies;
	bool criticals_rotation;
	string domain;
	bool useMatchmaker;
	bool useBackbone;
	int hash;
	int unex_unsat;
	int unex_sat;
	chrono::high_resolution_clock::time_point initial_time;
	vector<MUS> muses;
	Explorer* explorer;
	SatSolver* satSolver; 
	string sat_solver;

	Master(string filename, string alg);
	~Master();
	bool is_valid(Formula &f, bool core = false, bool grow = false);
	void block_down(Formula f);
	void block_up(MUS& f);
	void mark_MUS(MUS& m, bool block = true);
	MUS& shrink_formula(Formula& f, Formula crits = Formula());
	Formula grow_formula(Formula& f);
	void write_mus_to_file(MUS& f);
	void validate_mus(Formula &f);
	void enumerate();

	//reMUS algorithm functions
	int depthMUS;
	float dim_reduction;
	int current_depth;
	void find_all_muses_duality_based_remus(Formula subset, Formula crits, int depth);
	void extend_mus(Formula &top, Formula &mus, int dMUS = 1);

	//TOME algorithm functions
	void find_all_muses_tome();
	pair<Formula, Formula> local_mus(Formula bot, Formula top, int diff);

	//MARCO algorithm functions
	void marco_base();

	//dualMUS heuristic functions
	int scope_limit;

	//Matchmaker heuristic
	int rotated_muses;
	int recursive_rotation_delta(MUS m1, Formula &top, int depth);
	vector<int> minimal_hitting_set(vector<int> local_muses);

	//Backbone heuristic
	int backbone_mus_rotation(MUS &m1, Formula &top);	
	vector<int> propagate_backbone(int i, bool positive, vector<int> &lits_left, vector<bool> &satisfied, 
			Formula &top, vector<int> &singletons, Formula &seed, int c, vector<int> &influenced, 
			bool mark_influenced = false);
	vector<int> backbone_find_seed_beta(Formula &seed, Formula implied, Formula values, Formula satisfied, 
			vector<int> singletons, vector<int> lits_left, Formula &top);
	vector<int> backbone_init(Formula &seed, Formula &implied, Formula &values, vector<int> &lits_left, 
			Formula &satisfied, vector<int> &singletons, Formula &top,
			vector<int> &variables_map, vector<int> &variables_map_inv, int c);
	int get_backbones_bones(Formula &seed, Formula &implied, Formula &values, vector<int> &variables_map, 
			vector<int> &variables_map_inv);
	void backbone_build_literal_map(Formula &seed, vector<int> &variables_map, 
			vector<int> &variables_map_inv);
	void backbone_simplify(Formula &seed, int c, Formula &implied, Formula &values);
	void backbone_check_remainer(Formula &implied, Formula &values, Formula &top, Formula &m1);
	void implied_literal(int cl, Formula &implied, int &lit, int &var, bool &value);
	int seek_conflict(vector<int> &unit, vector<bool> &unit_value, vector<int> &singletons, int from, 
			int to, Formula &seed, Formula &top, Formula &implied, vector<int> &added);
};
#endif
