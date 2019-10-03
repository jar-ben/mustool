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
#include "counting/XorExplorer.h"

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
	bool useMixedHeuristic;
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
	void mark_MUS(MUS& m, bool block = true, bool verbose = true);
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

	//DAA algorithm functions
	void daa_base();

	//dualMUS heuristic functions
	int scope_limit;

	int approxMC2Core(float tresh, int &nCells);
	int counting(float e, float d);
	int bsat(float tresh);
	int bsat_xor(float tresh, XorExplorer &xe);
	vector<vector<CMSat::Lit>> blocksDown;
	vector<vector<CMSat::Lit>> blocksUp;
	int logSATSearch(vector<vector<int>> &As, int tresh, int mPrev);
	void initialBooster();
	Formula xor_shrink(Formula &t, XorExplorer &xe, int &limit);
};
#endif
