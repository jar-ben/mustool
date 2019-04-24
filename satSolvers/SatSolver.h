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
	std::string shrink_alg;
	std::string sat_solver;
	int checks;
	int hash;
	int shrinks;
	int grows;

	SatSolver(std::string filename): shrink_alg("default"), checks(0), shrinks(0), grows(0){}
	~SatSolver(){}
	virtual void exportMUS(std::vector<bool> mus, std::string outputFile) = 0;
	virtual bool solve(std::vector<bool> &f, bool core = false, bool grow = false) = 0;
	virtual std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	virtual std::vector<bool> grow(std::vector<bool> &f);
};


#endif

