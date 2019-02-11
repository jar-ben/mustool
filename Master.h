#include "MSHandle.h"
#include "Explorer.h"
#include "Z3Handle.h"
#include "SpotHandle.h"
#include "NuxmvHandle.h"
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

struct OverlapHash{
	public:
		std::size_t operator()(std::vector<int> const& over) const{
			std::size_t ret = over.size();
			for(auto &c: over){
				ret ^= c + 0x9e3779b9 + (ret << 6) + (ret >> 2);
			}
			return ret;
		}
};

inline size_t key(int i, int j){ return (size_t) i << 32 | (unsigned int) j; }

using namespace std;

typedef vector<int> Clause;
typedef vector<Formula> Path;
typedef Path::iterator PathIter;

class Master{
public:
	int dimension; //dimension of the input instance
	int variant; //MUS enumeration algorithm to be used
	int isValidExecutions;
	bool verbose; //TODO: make it int based 
	int depthMUS;
	int extend_top_variant;
	float dim_reduction;
	string output_file;
	bool hsd;
	int hsd_succ;
	int block_down_mus;
	bool validate_mus_c;
	float mus_approx;
	bool verify_approx;
	int overapproximated_muses;
	int rightly_approximated_muses;	
	int total_difference_of_approximated_muses;
	int current_depth;
	int current_dimension;
	bool model_rotation;
	int explored_extensions;
	int spurious_muses;
	bool get_implies;
	int skipped_shrinks;
	float crits_treshold;
	bool criticals_rotation;
	string domain;

	int hash;

	int unex_unsat;
	int unex_sat;
	int known_unex_unsat;
	int extracted_unex;
	int scope_limit;

	chrono::high_resolution_clock::time_point initial_time;
	
	vector<MUS> muses;
	Explorer* explorer;
	SatSolver* satSolver; 
	string sat_solver;

	Master(string filename, int var = 1, bool vis = false, string sat_solver = "");
	~Master();
	bool is_valid(Formula &f, bool core = false, bool grow = false);
	void block_down(Formula f);
	void block_up(MUS& f);
	void mark_MUS(MUS& m, bool block = true);
	MUS& shrink_formula(Formula& f, Formula crits = Formula());
	Formula grow_formula(Formula& f);
	void write_mus_to_file(MUS& f);

	void validate_mus(Formula &f);

	//reMUS algorithm functions
	void find_all_muses_duality_based_remus(Formula subset, Formula crits, int depth);
	void find_all_muses_duality_based_mono_remus(Formula subset, Formula crits, int depth);
	void extend_mus(Formula &top, Formula &mus, int dMUS = 1);

	//TOME algorithm functions
	void find_all_muses_tome();
	pair<Formula, Formula> local_mus(Formula bot, Formula top, int diff);

	//MARCO algorithm functions
	void marco_base();
	void find_all_muses_marco_explicit();
	void find_all_muses_marco_complement();
	int explicit_seeds;

	void find_all_muses_marco_recursive_rotation();
	int recursive_rotation(MUS &m1, vector<int> &top_int, vector<int> &local_muses, vector<vector<int>> &scope, int &depth, vector<int> &scope_indicator);
	int recursive_rotation_init(MUS &m1, Formula top);
	bool use_edge_muses;


	void update_rot_muses(int c1, int c2, int rot_it, vector<int>& scope_indicator);
	int recursive_rotation_gama(MUS &m1, vector<int> &top_int, int &depth, vector<int>& scope_indicator, int rot_id);
	int recursive_rotation_delta(MUS m1, Formula &top, int depth);
	vector<unordered_map<int, vector<int>>> rot_muses;
	vector<int> minimal_hitting_set(vector<int> local_muses);


	int mus_rotation(Formula &top, MUS& mus);
	int rotated_muses;
	int duplicated;
	int max_round;

	bool use_backbone();
	int backbone_mus_rotation(MUS &m1, Formula &top);	
	int backbone_mus_rotation_all(MUS &m1, Formula &top);
	vector<int> propagate_backbone(int i, bool positive, vector<int> &lits_left, vector<bool> &satisfied, Formula &top, vector<int> &singletons, Formula &seed, int c, vector<int> &influenced, bool mark_influenced = false);
	vector<int> backbone_find_seed(Formula &seed, Formula &implied, Formula &values, Formula satisfied, vector<int> singletons, vector<int> lits_left, Formula &top);
	vector<int> backbone_find_seed_beta(Formula &seed, Formula implied, Formula values, Formula satisfied, vector<int> singletons, vector<int> lits_left, Formula &top);
	vector<int> backbone_init(Formula &seed, Formula &implied, Formula &values, vector<int> &lits_left, Formula &satisfied, vector<int> &singletons, Formula &top,
			vector<int> &variables_map, vector<int> &variables_map_inv, int c);
	int get_backbones(Formula &seed, Formula &implied, Formula &values, vector<int> &variables_map, vector<int> &variables_map_inv);
	
	int get_backbones_bones(Formula &seed, Formula &implied, Formula &values, vector<int> &variables_map, vector<int> &variables_map_inv);

	void backbone_build_literal_map(Formula &seed, vector<int> &variables_map, vector<int> &variables_map_inv);
	void backbone_simplify(Formula &seed, int c, Formula &implied, Formula &values);

	bool backbone_check_reminder(Formula &implied, Formula &values, Formula &top, Formula &m1);

	vector<bool> backbone_crit_value;
	vector<bool> backbone_crit_used;


	void implied_literal(int cl, Formula &implied, int &lit, int &var, bool &value);
	int seek_conflict(vector<int> &unit, vector<bool> &unit_value, vector<int> &singletons, int from, int to, Formula &seed, Formula &top, Formula &implied, vector<int> &added);

	bool check_mus_via_hsd(Formula &f);
	void enumerate();

	vector<Formula> used_backbones;
	int used_backbones_count;
	int original_backbones_count;

	void build_mapping(Formula &top, Formula &impl, Formula &val, Formula &core, vector<vector<int>> &reduced, vector<int> &map1, vector<int> &map2, int &max_var, vector<int> &implied_clauses);
};
