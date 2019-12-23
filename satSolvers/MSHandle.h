#ifndef MSHANDLE_H
#define MSHANDLE_H

#include "custom_minisat/Solver.h"
#include "satSolvers/BooleanSolver.h"
#include <string>
#include <map>
#include <vector>
#include <unordered_map>


class MSHandle: public BooleanSolver{
public:
	CustomMinisat::Solver* solver;
	std::map<std::vector<int>,int> clauses_map;
	std::vector<std::string> clauses_str;
	std::unordered_map<std::string, int> clauses_unique_map;
	int vars;

	MSHandle(std::string filename);
	~MSHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
	bool solve(std::vector<bool> &f, std::vector<int> conflicts, bool shrink = false, bool grow = false);

	std::string toString(std::vector<bool> &mus);

	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
	bool parse(std::string filename);

	std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	std::vector<bool> shrink_muser(std::string input, int hash2);

	// helper function for muser2/dmuser manipulation
	void export_formula_crits(std::vector<bool> f, std::string filename, std::vector<bool> crits);
	std::vector<bool> import_formula_crits(std::string filename);

	// helper function for cmp manipulation
	std::vector<int> export_formula_wcnf(std::vector<bool> f, std::vector<bool> &conflicts, std::string filename);
	std::vector<bool> grow_cmp(std::vector<bool> &f, std::vector<bool> &conflicts);
	std::vector<bool> grow(std::vector<bool> &f, std::vector<bool> conflicts = std::vector<bool>());
	std::vector<bool> grow_uwrmaxsat(std::vector<bool> &f, std::vector<bool> &conflicts);
	std::vector<bool> satisfied(std::vector<int> &valuation);
	std::vector<std::vector<bool>> growMultiple(std::vector<bool> &f, std::vector<bool> conflicts = std::vector<bool>(), int limit = 1);

	//model rotation
	std::vector<bool> get_model();
	std::vector<bool> model_extension(std::vector<bool> subset, std::vector<bool> model);
};


#endif
