#ifndef MSHANDLE_H
#define MSHANDLE_H

#include "custom_minisat/Solver.h"
#include <string>
#include <map>
#include <vector>
#include "SatSolver.h"
#include <unordered_map>

//#include "muser2_api.hh"

class MSHandle: public SatSolver{
public:
	CustomMinisat::Solver* solver;
//	std::vector<std::vector<int>> clauses;
	std::map<std::vector<int>,int> clauses_map;
	std::vector<std::string> clauses_str;
	std::unordered_map<std::string, int> clauses_unique_map;
		

	MSHandle(std::string filename);
	~MSHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
	bool solve_glucose(std::vector<bool> &f, bool shrink = false, bool grow = false);
	bool solve_lingeling(std::vector<bool> &f);
	bool solve(std::vector<bool> &seed, std::vector<int> &assumptions);
	void get_unsat_core(std::vector<bool> &core);


	bool add_clause(std::vector<int> clause);
	bool add_unit(int lit);
	int nVars(){	return solver->nVars(); }
	int nClauses(){	return solver->nClauses(); }
	bool parse_dimacs(std::string filename);
	void add_flip_edge(std::vector<int> cl);
	void compute_flip_edges(int c);

	std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	std::vector<bool> custom_shrink(std::vector<bool> f);
	std::vector<bool> shrink_model_extension(std::vector<bool> &f, std::vector<bool> &crits, std::vector<std::vector<bool>>& extensions);

	std::vector<bool> shrink_muser(std::string input, int hash2);
	std::vector<bool> shrink_dmuser(std::string input, int hash2);
	std::vector<bool> shrink_mcsmus(std::string input, int hash2);


	std::vector<bool> shrink_api_muser(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());

	// helper function for muser2/dmuser manipulation
	void export_formula(std::vector<bool> f, std::string filename);
	void export_formula_crits(std::vector<bool> f, std::string filename, std::vector<bool> crits);
	std::vector<bool> import_formula(std::string filename);
	std::vector<bool> import_formula_crits(std::string filename);
	int muser_output(std::string filename);

	//model rotation
	void build_flip_graph();
	int model_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& subset,
		std::vector<bool>& model, std::vector<std::vector<bool>>& model_extensions);
	int model_rotation_beta(std::vector<bool>& criticals, int critical, std::vector<bool>& subset,
		std::vector<bool>& model, std::vector<std::vector<bool>>& model_extensions);
	std::vector<bool> get_model();
	std::vector<bool> model_extension(std::vector<bool> subset, std::vector<bool> model);

	void criticals_rotation(std::vector<bool>& criticals, std::vector<bool> subset);
	int critical_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& model, std::vector<int>& new_criticals);

	int get_implied(std::vector<bool>& controls, std::vector<bool>& implied, std::vector<bool>& values);
};


#endif
