#include "custom_minisat/Solver.h"
#include <string>
#include <map>
#include <vector>


class MSSSolver{
public:
	Minisat::Solver* solver;
	std::vector<std::vector<int>> clauses;
	int vars;
	int unex_clauses;

	MSSSolver();
	~MSSSolver();
	bool solve(){	return solver->solve();	}
	std::vector<bool> solve(std::vector<bool> &controls);
	bool block_down(std::vector<bool> cl);
	bool block_up(std::vector<bool> cl);


	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
	int nVars(){	return solver->nVars(); }
	int nClauses(){	return solver->nClauses(); }
	bool parse_dimacs(std::string filename);

};
