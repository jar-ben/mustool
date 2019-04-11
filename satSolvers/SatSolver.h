#ifndef SATSOLVER_H
#define SATSOLVER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

class SatSolver{
public:
	int dimension;
	std::string source;
	std::vector<std::vector<int>> clauses;
	std::vector<std::string> clauses_string;
	std::map<std::string, int> clauses_string_map;
	std::string shrink_alg;
	std::string sat_solver;
	int checks;
	int rotated_crits;
	int vars;
	int hash;

	//data structures for model rotation
	std::vector<std::vector<int>> hitmap_pos;
	std::vector<std::vector<int>> hitmap_neg;

	std::vector<bool> flip_edges_computed;
	std::vector<std::vector<std::vector<int>>> flip_edges; // flip_edges[i][j][k] - i-th clause, j-th literal in the clause, k-th edge from i-th clause under j-th literal in the flip graph
	std::vector<std::vector<int>> flip_edges_flatten;
	virtual void compute_flip_edges(int c){}

	SatSolver(std::string filename): source(filename), shrink_alg("default"), checks(0){ rotated_crits = 0; }
	~SatSolver(){}
	virtual bool solve(std::vector<bool> &f, bool core = false, bool grow = false) = 0;
	virtual int get_implied(std::vector<bool>& controls, std::vector<bool>& implied, std::vector<bool>& values){}
	virtual std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	virtual std::vector<bool> grow(std::vector<bool> &f);
	virtual int model_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& subset,
		std::vector<bool>& model, std::vector<std::vector<bool>>& model_extensions){}
	virtual std::vector<bool> get_model(){}
	virtual bool validate_model(std::vector<bool> model, std::vector<bool> subset, std::vector<bool> original){}
	virtual void export_formula(std::vector<bool> mus, std::string file){}
	virtual	std::vector<bool> shrink_model_extension(std::vector<bool> &f, std::vector<bool> &crits, std::vector<std::vector<bool>>& extensions){}
	virtual void criticals_rotation(std::vector<bool>& criticals, std::vector<bool> subset){}
	virtual int critical_rotation(std::vector<bool>& criticals, int critical, std::vector<bool>& model, std::vector<int>& new_criticals){};
	virtual bool solve(std::vector<bool> &seed, std::vector<int> &assumptions){};
	virtual void get_unsat_core(std::vector<bool> &core){};
};


#endif

