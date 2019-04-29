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
	void addExpr(z3::expr &e);
	void parseMultiple(Z3_ast_vector &a_vec, int size);
	void parseAndFormula(Z3_ast_vector &a_vec);
	std::vector<bool> shrink(std::vector<bool> &formula, std::vector<bool> crits = std::vector<bool>());	
	bool solve(std::vector<bool> &formula, bool core = false, bool grow = false);
	std::string toString(std::vector<bool> &mus);
};
