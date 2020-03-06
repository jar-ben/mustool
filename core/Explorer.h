#ifndef EXPLORER_H
#define EXPLORER_H

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

class Explorer{
public:
	int dimension;
	int vars;
	CustomMinisat::Solver* solver;
	CustomMinisat::Solver* botSolver;
	CustomMinisat::Solver* topSolver;
	int calls;

	Explorer(int vars, bool verb = false);
	~Explorer();
	bool block_up(std::vector<bool> formula);
        bool block_down(std::vector<bool> formula);
	std::vector<bool> get_unexplored(uint8_t polarity = 0, bool rnd_pol = false);
	std::vector<bool> get_bot_unexplored(std::vector<bool>);
	std::vector<bool> get_top_unexplored(std::vector<bool>);
	std::vector<bool> get_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_top_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_bot_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_bot_unexplored_inside(std::vector<bool> subset);
	std::vector<bool> get_top_unexplored_inside(std::vector<bool> subset);
	std::vector<bool> get_bot_unexplored_containing(std::vector<bool> subset);
	std::vector<bool> get_unexplored(std::vector<bool> top, std::vector<bool>& mus);
	bool isUnexplored(std::vector<bool> valuation);
	bool isUnexploredSat(std::vector<bool> valuation);
	bool isUnexploredUnsat(std::vector<bool> valuation);
	int getImplied(std::vector<bool>& implied, std::vector<bool>& formula);
	int getConflicts(std::vector<bool>& conflicts, std::vector<bool>& formula);
	
	std::vector<bool> mus_intersection;
	std::vector<bool> mus_union;

	//data structures for maintaining criticals
	std::vector<bool> critical; //critical[c] = true iff c forms singleton MCS
	std::vector<bool> maybe_critical; //critical[c] = true iff c is a part of some MCS
	std::vector<bool> testedForCriticality;
	std::vector<std::vector<int>> mcses; //serves for storing mcses
	std::vector<std::vector<int>> parent_mcses; //parent_mcses[c] holds ids of MCSes stored in mcses
	bool is_critical(int c, std::vector<bool> &subset); // returns true if c is critical for subset
	int criticals;

	//data structures for maintaining availables - the opposites of criticals
	std::vector<bool> maybe_available; //available[c] = true iff c is a part of some MUS
	std::vector<std::vector<int>> muses;
	std::vector<std::vector<int>> parent_muses; //parent_muses[c] holds ids of MUSes stored in muses_int
	bool is_available(int c, std::vector<bool> &subset); // returns true if c is available for subset (can be added to it)
};

#endif
