#ifdef UMCSMUS

#include <errno.h>
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "satSolvers/mcsmus_handle.h"
#include "mcsmus/minisat-wrapper.hh"
#include "mcsmus/glucose-wrapper.hh"
#include "mcsmus/minisolver.hh"
#include "mcsmus/parseutils.hh"
#include "mcsmus/options.hh"
#include "mcsmus/mcsmus.hh"
#include "mcsmus/dimacs.hh"
#include "mcsmus/system.hh"

using namespace mcsmus;

std::vector<Lit> intToLit(std::vector<int> cls){
	std::vector<Lit> lits;
	//we are taking the clasue from the MSHandle class, and each clause contain additional control literal
	//here we have to exclude the control literal
	//TODO: we should not store the control literals in MSHandle.clauses, get rid of them
	for(int i = 0; i < cls.size() - 1; i++){
		auto l = cls[i];
		int var = abs(l)-1;
		lits.push_back( mkLit(var, l<0) );
	}
	return lits;
}

std::vector<bool> shrink_mcsmus(std::vector<bool> &f, std::vector<std::vector<int>> &clauses, std::vector<bool> crits){
	setX86FPUPrecision();
	Wcnf wcnf;
	std::unique_ptr<BaseSolver> s;
	s = make_glucose_simp_solver();
	MUSSolver mussolver(wcnf, s.get());

	
	Control* control = getGlobalControl();;
	control->verbosity = 0;
	control->status_interval = 10;
	mussolver.single_mus = true;
	mussolver.initializeControl();
	control->force_exit_on_interrupt();

	//add the clauses
	std::vector<int> constraintGroupMap;
	int cnt = 0;
	for(int i = 0; i < f.size(); i++){
		if(f[i]){
			if(crits[i]){
				wcnf.addClause(intToLit(clauses[i]), 0);
			}else{
				constraintGroupMap.push_back(i);
				cnt++;
				wcnf.addClause(intToLit(clauses[i]), cnt);
			}
		}
	}
	
	std::vector<Lit> mus_lits;
	wcnf.relax();
	mussolver.find_mus(mus_lits, false);

	std::vector<bool> mus(f.size(), false);
	for(auto b : mus_lits){
		int gid = wcnf.bvarIdx(b);
		if(gid > 0){
			mus[constraintGroupMap[gid - 1]] = true;
		}	
	}
	for(int i = 0; i < f.size(); i++){
		if(crits[i]){
			mus[i] = true;
		}
	}
	return mus;
}
#endif
