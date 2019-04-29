#ifndef MSHANDLE_H
#define MSHANDLE_H

#include "custom_minisat/Solver.h"
#include "satSolvers/SatSolver.h"
#include <string>
#include <map>
#include <vector>
#include <unordered_map>


class MSHandle: public SatSolver{
public:
	CustomMinisat::Solver* solver;
	std::vector<std::vector<int>> clauses;
	std::map<std::vector<int>,int> clauses_map;
	std::vector<std::string> clauses_str;
	std::unordered_map<std::string, int> clauses_unique_map;
	int vars;

	MSHandle(std::string filename);
	~MSHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
	void get_unsat_core(std::vector<bool> &core);

	std::string toString(std::vector<bool> &mus);

	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
	bool parse_dimacs(std::string filename);

	std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	std::vector<bool> shrink_muser(std::string input, int hash2);
	std::vector<bool> shrink_mcsmus(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());

	// helper function for muser2/dmuser manipulation
	void export_formula_crits(std::vector<bool> f, std::string filename, std::vector<bool> crits);
	std::vector<bool> import_formula_crits(std::string filename);

	//model rotation
	int rotated_crits;	
	std::vector<std::vector<int>> hitmap_pos;
	std::vector<std::vector<int>> hitmap_neg;
	std::vector<bool> flip_edges_computed;
	std::vector<std::vector<std::vector<int>>> flip_edges; // flip_edges[i][j][k] - i-th clause, j-th literal in the clause, k-th edge from i-th clause under j-th literal in the flip grap
	std::vector<std::vector<int>> flip_edges_flatten;
	void compute_flip_edges(int c);
	int model_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& subset,
		std::vector<bool>& model, std::vector<std::vector<bool>>& model_extensions);
	std::vector<bool> get_model();
	std::vector<bool> model_extension(std::vector<bool> subset, std::vector<bool> model);

	void criticals_rotation(std::vector<bool>& criticals, std::vector<bool> subset);
	int critical_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& model, std::vector<int>& new_criticals);

	int get_implied(std::vector<bool>& controls, std::vector<bool>& implied, std::vector<bool>& values);
};


#endif
