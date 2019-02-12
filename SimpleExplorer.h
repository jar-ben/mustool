#ifndef SEXPLORER_H
#define SEXPLORER_H

#include "MSHandle.h"
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



class SimpleExplorer{
public:
	bool verbose;
	int dimension;
	CustomMinisat::Solver* solver;
	int vars;
	int added_clauses;

	SimpleExplorer(int vars, bool verb = false);
	~SimpleExplorer();
	bool block_up(std::vector<bool> mus);
        bool block_down(std::vector<bool> formula);
	std::vector<bool> get_unexplored(uint8_t polarity = 0, bool rnd_pol = false);
	std::vector<bool> get_bot_unexplored(std::vector<bool>);
	std::vector<bool> get_top_unexplored(std::vector<bool>);
	std::vector<bool> get_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_top_unexplored(std::vector<int> assumptions);
	std::vector<bool> get_bot_unexplored(std::vector<int> assumptions);
	bool check(std::vector<bool> f);
};

#endif
