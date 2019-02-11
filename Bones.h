#include "minibones/src/ToolConfig.hh"
#include "minibones/src/CoreBased.hh"
#include "minibones/src/LitSet.hh"
#include <iostream>
#include "types_h.h"

class Bones{
public:
	ToolConfig config;
	Bones(){ set_options(); }
	~Bones(){}

	void set_options();	
	std::vector<int> get_backbones(std::vector<std::vector<int>> &cls);
	
};

Lit itoLitB(int i){
	bool sign = i < 0;
	int var = (sign)? -i-1 : i-1;
	return (sign) ? ~mkLit(var) : mkLit(var);
}

void Bones::set_options(){
      config.set_backbone_insertion(1);
      config.set_use_core_based(true);
      config.set_chunk_size(100);
      config.set_input_file_name("test_backbone.cnf");
}

std::vector<int> Bones::get_backbones(std::vector<std::vector<int>> &cls){
	CNF clauses;
	for(auto c: cls){
		vector<Lit> literals;
		for(auto l: c){
			literals.push_back( itoLitB(l) );
		}
		clauses.push_back(LitSet(literals));
	}
	
	/*	CoreBased *pcorebased = new CoreBased(config, output, reader.get_max_id(), reader.get_clause_vector());
	CoreBased& corebased=*pcorebased;
	if (!corebased.initialize()) {//unsatisfiable
		cout << "unsat backbone input" << endl;
		delete pcorebased;
		exit(1);
	}
	corebased.run();
	std::cout << "backbone computed" << endl;
	delete pcorebased;
*/	return vector<int>();
}

