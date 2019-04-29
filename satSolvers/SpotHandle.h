#ifndef NOLTL
#include "SatSolver.h"
#include <iostream>
#include <fstream>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>
#include <spot/tl/simplify.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/emptiness.hh>
#include <spot/twaalgos/gtec/gtec.hh>

class SpotHandle: public SatSolver{
public:
	std::vector<spot::formula> formulae;

	std::string toString(std::vector<bool> &mus);
	SpotHandle(std::string filename);
	~SpotHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
};
#endif
