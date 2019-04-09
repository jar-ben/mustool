#ifndef SIMPLESAT_H
#define SIMPLESAT_H

#include "custom_minisat/Solver.h"
#include <string>
#include <map>
#include <vector>
#include "SatSolver.h"
#include <unordered_map>

class SimpleSAT{
public:
	CustomMinisat::Solver* solver;

	SimpleSAT(std::vector<std::vector<int>> cls, int v);
	~SimpleSAT();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
	void get_unsat_core(std::vector<bool> &core);

	int vars;
	int checks;
	int dimension;
	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
};
#endif
