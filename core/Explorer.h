#ifndef EXPLORER_H
#define EXPLORER_H

#include "satSolvers/MSHandle.h"
#include "custom_minisat/Solver.h"
#include <string>
#include <map>
#include <ctime>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <utility>
#include <queue>
#include <tuple>
#include "types_h.h"
#include <unordered_map> 


struct CmpSeeds
{
    bool operator()(const triple& s1, const triple& s2)
    {
        return (std::get<4>(s1) - (std::get<1>(s1) + std::get<2>(s1))) > (std::get<4>(s2) - (std::get<1>(s2) + std::get<2>(s2)));
    }
};


class Explorer{
public:
	bool verbose;
	int dimension;
	CustomMinisat::Solver* solver;
	CustomMinisat::Solver* botSolver;
	CustomMinisat::Solver* topSolver;

	SatSolver* satSolver;
		
	int bot_calls, top_calls;
	int vars;
	int added_clauses;
	bool use_implications;
	bool test_rotation_unex;
	int rotation_pairs;
	int mus_size_limit;

	std::vector<std::vector<int>> core; //the set of all critical constraints. for each c in c, c[0] = id of constraints. 
	

	//std::vector<bool> get_unex_impl_sat(std::vector<bool> top, std::vector<bool> &searchspace);

	Explorer(int vars, bool verb = false);
	~Explorer();
	bool block_up(MUS& mus, bool block_unex = true);
	bool block_up_unex(std::vector<int> &formula);
        bool block_down(std::vector<bool> formula);
	std::vector<bool> get_unexplored(uint8_t polarity = 0, bool rnd_pol = false);
	std::vector<bool> get_bot_unexplored(std::vector<bool>);
	std::vector<bool> get_top_unexplored(std::vector<bool>);
	std::vector<bool> solve(std::vector<int> assumptions);
	bool checkValuation(std::vector<bool> valuation);
	bool checkValuationTopSolver(std::vector<bool> valuation);
	bool checkValuationBotSolver(std::vector<bool> valuation);
	bool validateTop(std::vector<bool> valuation);
	std::vector<bool> get_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_top_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_bot_unexplored(std::vector<int> assumptions);
	int get_implies(std::vector<bool>& implied, std::vector<bool>& formula);

	//data structures for model rotation
	int rotation_edge_uses;
	int rotated_unsat_seeds;
	float minimum_duration;
	float maximum_duration;
	float total_duration;
	int durations;
	int total_iters;

	std::vector<std::unordered_map<int, bool>> hitting_pairs;
	int digits;

	
	bool is_hitting_pair(int c1, int c2);
	bool is_backbone_pair(int c2, std::vector<bool> &implied, std::vector<bool>& values);
	std::vector<bool> get_unexplored(std::vector<bool> top, std::vector<bool>& mus);
	bool maximize_unex(std::vector<bool> &unex, int c1 = -1, int c2 = -1);
	int found_conflicting;
	int reduced_unsat;
	int mus_based_round;

	std::vector<std::pair<std::vector<bool>, int>> init_rotation(Formula &top, MUS &m2);
	std::vector<std::pair<std::vector<bool>, int>> init_rotation2(std::vector<bool> &top, std::vector<bool> &m2);
	std::vector<std::pair<std::vector<bool>, int>> init_rotation3(std::vector<bool> &top, std::vector<bool> &m2);

	std::vector<bool> init_rotation_single(Formula &top, MUS &m1);

	//data structures for maintaining criticals
	std::vector<bool> critical; //critical[c] = true iff c forms singleton MCS
	std::vector<bool> maybe_critical; //critical[c] = true iff c is a part of some MCS
	std::vector<bool> mus_intersection;
	std::vector<bool> mus_union;

	std::vector<std::vector<int>> mcses; //serves for storing mcses
	std::vector<std::vector<int>> parent_mcses; //parent_mcses[c] holds ids of MCSes stored in mcses
	bool is_critical(int c, std::vector<bool> &subset); // returns true if c is critical for subset
	int criticals;

	//data structures for maintaining availables - the opposites of criticals
	std::vector<bool> maybe_available; //available[c] = true iff c is a part of some MUS
	std::vector<MUS> muses;
	std::vector<std::vector<int>> parent_muses; //parent_muses[c] holds ids of MUSes stored in muses_int
	bool is_available(int c, std::vector<bool> &subset); // returns true if c is available for subset (can be added to it)


	//rotation solver beta
	std::vector<float> mus_shrink_duration;
	bool easy_to_shrink(MUS &mus);

	std::vector<bool> get_seed_explicit();
	bool is_unexplored(std::vector<bool> seed);
	bool is_unexplored(std::vector<bool> &m1, std::vector<bool> &m2, int c1, int c2);
	std::vector<int> rotation_start; //rotation_start[i] = j --> there is no rotation pair for mus i within the first 1....j muses

	//seed priority queue
	std::priority_queue<triple, std::vector<triple>, CmpSeeds> seeds_queue;
	bool seed_priority_queue;
	void add_seeds_to_queue(MUS& m1);


	bool rotate_co(std::vector<bool> &top);
	std::vector<bool> get_co_unexplored();

	int sat_rotation(std::vector<bool>& subset, std::vector<bool>& model, std::vector<std::vector<int>> &clauses, int &limit);
	bool check_maximality(std::vector<bool> &unex);
	int sat_rotated;

};

#endif
