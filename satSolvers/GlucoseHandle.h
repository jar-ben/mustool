#ifndef GLHANDLE_H
#define GLHANDLE_H

#include "glucose/Solver.h"
#include "glucose/SimpSolver.h"
#include "satSolvers/BooleanSolver.h"
#include <string>
#include <map>
#include <vector>
#include <unordered_map>


class GlucoseHandle: public BooleanSolver{
public:
	//CustomGlucose::Solver* solver;
	CustomGlucose::SimpSolver* solver;
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
	std::vector<bool> get_model();
	std::vector<bool> model_extension(std::vector<bool> subset, std::vector<bool> model);
};


#endif
