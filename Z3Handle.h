#include <string>
#include <map>
#include <vector>
#include "z3++.h"
#include "SatSolver.h"



class Z3Handle: public SatSolver{
public:
	z3::context ctx;
	z3::solver* s;	
	std::vector<z3::expr> controls;
	std::vector<z3::expr> implications;
	std::vector<z3::expr> clauses;
	std::map<z3::expr, int> controls_map;

	Z3Handle(std::string);
	~Z3Handle();
	std::vector<bool> shrink(std::vector<bool> &formula, std::vector<bool> crits = std::vector<bool>());	
	std::vector<bool> shrink_hsmtmuc(std::vector<bool> &formula, std::vector<bool> crits = std::vector<bool>());
	std::vector<bool> custom_shrink(std::vector<bool> &formula, std::vector<bool> crits = std::vector<bool>());
	bool solve(std::vector<bool> &formula, bool core = false, bool grow = false);

	// helper function for hsmtmuc manipulation
	void export_formula(std::vector<bool> f, std::string filename, std::vector<bool> crits);
	std::vector<bool> import_formula(std::string filename);


};
