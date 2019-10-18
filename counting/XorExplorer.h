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
	int xorSize;

	CMSat::SATSolver* solver;
	XorExplorer(int vars, std::vector<std::vector<CMSat::Lit>> &blocksDown, std::vector<std::vector<CMSat::Lit>> &blocksUp);
	~XorExplorer();
	std::vector<CMSat::Lit> block_up(Formula &f);
	std::vector<CMSat::Lit> block_down(Formula &f);
	Formula solve(std::vector<CMSat::Lit> assumptions = std::vector<CMSat::Lit>());
	Formula solve_subset(Formula &f);
	void add_xor(int m, std::vector<std::vector<int>> &As);
	bool block(Formula);
	bool check(Formula &mus);
	void checkUniformity();
	std::vector<std::vector<uint32_t>> xors;
	std::string export_xors();
};

#endif
