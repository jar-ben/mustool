#ifndef GLHANDLE_H
#define GLHANDLE_H

#include "glucose/Solver.h"
#include "glucose/SimpSolver.h"
#include "satSolvers/SatSolver.h"
#include <string>
#include <map>
#include <vector>
#include <unordered_map>


class GlucoseHandle: public SatSolver{
public:
	//CustomGlucose::Solver* solver;
	CustomGlucose::SimpSolver* solver;
	std::vector<std::vector<int>> clauses;
	std::map<std::vector<int>,int> clauses_map;
	std::vector<std::string> clauses_str;
	std::unordered_map<std::string, int> clauses_unique_map;
	int vars;

	GlucoseHandle(std::string filename);
	~GlucoseHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);

	std::string toString(std::vector<bool> &mus);

	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
	bool parse(std::string filename);

	std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());

	// helper function for muser2/dmuser manipulation
	void export_formula_crits(std::vector<bool> f, std::string filename, std::vector<bool> crits);

	//model rotation
	int rotated_crits;	
	std::vector<std::vector<int>> hitmap_pos;
	std::vector<std::vector<int>> hitmap_neg;
	std::vector<bool> flip_edges_computed;
	std::vector<std::vector<std::vector<int>>> flip_edges; // flip_edges[i][j][k] - i-th clause, j-th literal in the clause, k-th edge from i-th clause under j-th literal in the flip grap
	std::vector<std::vector<int>> flip_edges_flatten;
	void compute_flip_edges(int c);
	std::vector<bool> get_model();
	std::vector<bool> model_extension(std::vector<bool> subset, std::vector<bool> model);
	void criticals_rotation(std::vector<bool>& criticals, std::vector<bool> subset);
};


#endif
