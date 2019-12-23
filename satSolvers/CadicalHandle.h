#ifndef CADICALHANDLE_H
#define CADICALHANDLE_H

#include "BooleanSolver.h"
#include "cadical.hpp"
#include <string>
#include <vector>
#include <unordered_map>

class CadicalHandle: public BooleanSolver{
public:

	CaDiCaL::Solver* solver;
	std::string toString(std::vector<bool> &mus);
	CadicalHandle(std::string filename);
	~CadicalHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);

	int vars;
	std::unordered_map<std::string, int> clauses_unique_map;
	std::vector<std::string> clauses_str;
	bool add_clause(std::vector<int> clause);
	bool parse(std::string filename);
	std::vector<int> convert_clause(std::string clause);
};

#endif
