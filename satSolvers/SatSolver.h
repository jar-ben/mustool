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
	std::string grow_alg;
	std::string sat_solver;
	int checks;
	int hash;
	int shrinks;
	int grows;
	int exported_muses;
	int exported_msses;
	std::string mcslsArgs;
	int growStrategy;

	SatSolver(std::string filename): shrink_alg("default"), checks(0), shrinks(0), grows(0), exported_muses(0){}
	~SatSolver(){}
	virtual std::string toString(std::vector<bool> &mus) = 0;
	virtual bool solve(std::vector<bool> &f, bool core = false, bool grow = false) = 0;
	virtual bool solve(std::vector<bool> &f, std::vector<int> conflicts, bool core = false, bool grow = false);
	virtual std::vector<bool> shrink(std::vector<bool> &f, std::vector<bool> crits = std::vector<bool>());
	virtual std::vector<bool> grow2(std::vector<bool> &f, std::vector<std::vector<bool>> &unsat, std::vector<bool> conflicts = std::vector<bool>());
	virtual std::vector<bool> grow(std::vector<bool> &f, std::vector<bool> conflicts = std::vector<bool>());
	virtual std::vector<std::vector<bool>> growMultiple(std::vector<bool> &f, std::vector<bool> conflicts = std::vector<bool>(), int limit = 1);
	void exportMUS(std::vector<bool> mus, std::string outputFile);
	void exportMSS(std::vector<bool> mss, std::string outputFile);
};


#endif

