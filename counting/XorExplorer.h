#ifndef XOREXPLORER_H
#define XOREXPLORER_H

#include <string>
#include <map>
#include <ctime>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <utility>
#include <queue>
#include <tuple>
#include "core/types_h.h"
#include <unordered_map> 
#include <cryptominisat5/cryptominisat.h>

class XorExplorer{
public:
	int dimension;
	int vars;
	int iter;

	CMSat::SATSolver* solver;
	XorExplorer(int vars, bool verb = false);
	~XorExplorer();
	bool block_up(MUS& mus);
        bool block_down(std::vector<bool> formula);
	Formula solve();
	void add_xor(int m, std::vector<std::vector<int>> As, std::vector<bool> alfa);
	bool temporary_block_up(Formula);
};

#endif
