#include <errno.h>
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "satSolvers/MSHandle.h"
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

vector<Lit> do_solve(mcsmus::Wcnf& w, mcsmus::MUSSolver& S){
	Control* control = getGlobalControl();;
	control->verbosity = 0;
	control->status_interval = 10;
	S.single_mus = true;
	S.initializeControl();
	control->force_exit_on_interrupt();
	w.relax();
	vector<Lit> mus;
	S.find_mus(mus, false);
	return mus;
}

std::vector<bool> MSHandle::shrink_mcsmus(std::vector<bool> &f, std::vector<bool> crits){
	setX86FPUPrecision();
	Wcnf wcnf;
	std::unique_ptr<BaseSolver> s;
	s = make_glucose_simp_solver();
	MUSSolver mussolver(wcnf, s.get());
	int cnt = 0;
	for(int i = 0; i < dimension; i++){
		if(f[i]){
			cnt++;
			wcnf.addClause(intToLit(clauses[i]), cnt);
		}
	}

	std::vector<Lit> mus_lits = do_solve(wcnf, mussolver);
	std::vector<int> f_int;
	for(int i = 0 ; i < dimension; i++){
		if(f[i]){
			f_int.push_back(i);
		}
	}
	std::vector<bool> mus(dimension, false);
	for (auto b : mus_lits){
		mus[f_int[wcnf.bvarIdx(b) - 1]] = true;
	}
	return mus;
}

