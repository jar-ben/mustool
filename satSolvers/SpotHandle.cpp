#ifndef NOLTL
#include "SpotHandle.h"

using namespace std;


SpotHandle::SpotHandle(string filename):SatSolver(filename){
	ifstream file;
	file.open(filename, ifstream::in);
	string tf;
	while (getline(file, tf)){
		formulae.push_back(spot::parse_formula(tf));
	}
	dimension = formulae.size();
}

SpotHandle::~SpotHandle(){
}

string SpotHandle::toString(vector<bool> &mus){
	stringstream result;
	for(int i = 0; i < dimension; i++){
		if(mus[i]){
			result << formulae[i] << "\n";
		}
	}
	return result.str();
}

bool SpotHandle::solve(std::vector<bool> &f, bool shrink, bool grow){
	checks++;
	vector<spot::formula> form;
	for(int i = 0; i < f.size(); i++)
		if(f[i])
			form.push_back(formulae[i]);
	spot::formula con = spot::formula::And(form);
	spot::tl_simplifier simp;
	spot::translator trans(&simp);
	spot::twa_graph_ptr a = trans.run(con);
	spot::couvreur99_check checker(a);
	auto res = checker.check();
	bool result = (bool) res;
	return result;
}
#endif
