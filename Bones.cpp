#ifdef MINIBONES_SRC
#include "Bones.h"

using namespace minibones;

Var itoVar(int i){
	bool sign = i < 0;
	return (sign)? -i-1 : i-1;
}

Lit itoLitB(int i){
	return ( i<0 ) ? ~mkLit(itoVar(i)) : mkLit(itoVar(i));
}

Bones::~Bones(){ }

Bones::Bones(){ set_options(); }

void Bones::set_options(){
      config.set_backbone_insertion(1);
      config.set_use_core_based(true);
      config.set_chunk_size(100);
}

vector<int> backbone_literals(CoreBased& c, int max_var){
	vector<int> backbones;
	for(int i = 1; i < max_var; i++){
		Lit v = itoLitB(i);
		Var v_var = itoVar(i);
		if(c.is_backbone(v)){
			int sign = (c.backbone_sign(v_var))? 1 : -1;
			backbones.push_back( i * sign);
		}
	}
	return backbones;
}

std::vector<int> Bones::get_backbones(std::vector<std::vector<int>> &cls){
	CNF clauses;
	Minisat::Var max_id = 0;
	for(auto c: cls){
		vector<Lit> literals;
		for(auto l: c){
			literals.push_back( itoLitB(l) );
			if(itoVar(l) + 1 > max_id)
				max_id = itoVar(l) + 1;
		}
		clauses.push_back(LitSet(literals));
	}
	ostream& output=std::cout;
	CoreBased *pcorebased = new CoreBased(config, output, max_id, clauses);
	CoreBased& corebased=*pcorebased;
	if (!corebased.initialize()) {//unsatisfiable
		cout << "unsat backbone input" << endl;
		delete pcorebased;
		exit(1);
	}
	corebased.run();
	vector<int> backbones = backbone_literals(corebased, max_id);
	delete pcorebased;
	return backbones; 
}
#endif
