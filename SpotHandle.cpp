#include "SpotHandle.h"

using namespace std;


SpotHandle::SpotHandle(string filename):SatSolver(filename){
	ifstream file;
	file.open(filename, ifstream::in);
	std::string tf;
	while (std::getline(file, tf)){
		formulae.push_back(spot::parse_formula(tf));
		clauses_string.push_back(tf);
	}
	dimension = formulae.size();
}

SpotHandle::~SpotHandle(){
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

