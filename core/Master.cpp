#include "Master.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


Master::Master(string filename, string alg){
        isValidExecutions = 0;
	algorithm = alg;
	sat_solver = "default";
	if(ends_with(filename, "smt2")){
		#ifdef NOSMT
			std::cout << "working with smt is currently not enabled, plese build the tool with the flag USESMT=YES, e.g. 'make USESMT=YES'" << std::endl;
			exit(1);
		#else
		satSolver = new Z3Handle(filename);
	       	#endif	
		domain = "smt";
	}
	else if(ends_with(filename, "cnf")){
		satSolver = new groupMSHandle(filename);
		domain = "sat";
	}
	else if(ends_with(filename, "ltl")){
		#ifdef NOLTL
			std::cout << "working with ltl is currently not enabled, plese build the tool with the flag USELTL=YES, e.g. 'make USELTL=YES'" << std::endl;
			exit(1);
		#else
		if(sat_solver == "nuxmv")
			satSolver = new NuxmvHandle(filename); 
		else
			satSolver = new SpotHandle(filename);
		#endif
		domain = "ltl";
	}
	else
		print_err("wrong input file");
	dimension = satSolver->dimension;	
	cout << "Dimension:" << dimension << endl;
        explorer = new Explorer(dimension);	
        verbose = false;
	depthMUS = 0;
	dim_reduction = 0.5;
	output_file = "";
	validate_mus_c = false;
	current_depth = 0;
	unex_sat = unex_unsat = 0;
	scope_limit = 100;	
        hash = random_number();
	satSolver->hash = hash;
}

Master::~Master(){
	delete satSolver;
	delete explorer;
}

void Master::write_mus_to_file(MUS& f){
	satSolver->exportMUS(f.bool_mus, output_file);
}

// mark formula and all of its supersets as explored
void Master::block_up(MUS& formula){
	explorer->block_up(formula);
}

// mark formula and all of its subsets as explored
void Master::block_down(Formula formula){
        explorer->block_down(formula);
}

// check formula for satisfiability
// core and grow controls optional extraction of unsat core and model extension (replaces formula)
bool Master::is_valid(Formula &formula, bool core, bool grow){
	bool sat = satSolver->solve(formula, core, grow); 
	if(sat && model_rotation){
		groupMSHandle *msSolver = static_cast<groupMSHandle*>(satSolver);
		Formula model = msSolver->get_model();
		vector<bool> crits (dimension, false);
		for(int i = 0; i < dimension; i++){
			if(!formula[i] && !explorer->is_available(i, formula)){
				crits[i] = true;
			}
		}
		for(int i = 0; i < dimension; i++){
			if(crits[i]){
				vector<Formula> model_extensions;
				int rotated = msSolver->model_rotation(crits, i, formula, model, model_extensions);
				int fresh = 0;
				for(auto &ext: model_extensions){
					if(explorer->isUnexploredSat(ext)){
						block_down(ext);
						fresh++;
					}
				}
				cout << "rotated: " << rotated << ", extensions: " << model_extensions.size() << ", fresh: " << fresh << endl;
			}
		}
	}
	return sat;
}

//verify if f is a MUS
void Master::validate_mus(Formula &f){
	if(is_valid(f))
		print_err("the mus is SAT");
	if(!explorer->isUnexplored(f))
		print_err("this mus has been already explored");
	for(int l = 0; l < f.size(); l++)
		if(f[l]){
			f[l] = false;
			if(!is_valid(f))
				print_err("the mus has an unsat subset");	
			f[l] = true;
		}	
}

MUS& Master::shrink_formula(Formula &f, Formula crits){
	int f_size = count_ones(f);
	chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();
//	cout << "shrinking dimension: " << f_size << endl;
	if(algorithm == "tome")
		is_valid(f,true,true); //get core before shrink in the case of TOME
	f_size = count_ones(f);
	if(crits.empty()) crits = explorer->critical;
	if(get_implies){ //get the list of known critical constraints	
		explorer->getImplied(crits, f);	
		cout << "criticals before rot: " << count_ones(crits) << endl;	
		if(criticals_rotation){
			int before = count_ones(crits);
			groupMSHandle *msSolver = static_cast<groupMSHandle*>(satSolver);
			msSolver->criticals_rotation(crits, f);
			cout << "rotated: " << (count_ones(crits) - before) << endl;
		}
		float ones_crits = count_ones(crits);		 
		if(f_size == ones_crits){ // each constraint in f is critical for f, i.e. it is a MUS 
			muses.push_back(MUS(f, -1, muses.size(), f_size)); //-1 duration means skipped shrink
			return muses.back();
		}		
		if((ones_crits/f_size) > crits_treshold){
			if(!is_valid(crits, false, false)){ 
				muses.push_back(MUS(crits, -1, muses.size(), f_size));//-1 duration means skipped shrink
				return muses.back();
			}
		}
	}

	Formula mus = satSolver->shrink(f, crits);
	chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>( end_time - start_time ).count() / float(1000000);
	muses.push_back(MUS(mus, duration, muses.size(), f_size));
	return muses.back();
}

//grow formula into a MSS
Formula Master::grow_formula(Formula &f){
	return satSolver->grow(f);
}


void Master::mark_MUS(MUS& f, bool block_unex, bool verbose){	
	if(validate_mus_c) validate_mus(f.bool_mus);		
	explorer->block_up(f);

	if(verbose){
		chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
		cout << "Found MUS #" << muses.size() <<  ", mus dimension: " << f.dimension;
		cout << ", checks: " << satSolver->checks << ", time: " << duration;
		cout << ", unex sat: " << unex_sat << ", unex unsat: " << unex_unsat << ", criticals: " << explorer->criticals;
		cout << ", intersections: " << std::count(explorer->mus_intersection.begin(), explorer->mus_intersection.end(), true);
		cout << ", union: " << std::count(explorer->mus_union.begin(), explorer->mus_union.end(), true) << ", dimension: " << dimension;
		cout << ", seed dimension: " << f.seed_dimension << ", duration: " << f.duration;
		cout << ", shrinks: " << satSolver->shrinks;
		cout << endl;
	}

	if(output_file != "")
		write_mus_to_file(f);
}

void Master::enumerate(){
	initial_time = chrono::high_resolution_clock::now();
	cout << "running algorithm: " << algorithm << endl;
	Formula whole(dimension, true);
	if(is_valid(whole))
		print_err("the input instance is satisfiable");

	if(algorithm == "remus"){
		find_all_muses_duality_based_remus(Formula (dimension, true), Formula (dimension, false), 0);
	}
	else if(algorithm == "tome"){
		find_all_muses_tome();
	}
	else if(algorithm == "marco"){
		marco_base();
	}
	else if(algorithm == "daa"){
		daa_base();
	}
	else if(algorithm == "counter"){
		counting(0.8, 0.2);
	}
	else if(algorithm == "union"){
		union_base();
	}
	else{
		print_err("invalid algorithm chosen");
	}
	return;
}

void Master::union_base(){
	Formula seed(dimension, true);
	while(true){
		MUS mus = shrink_formula(seed);
		mark_MUS(mus);
		ofstream file;
		file.open("exp.cnf", std::ofstream::out);
		Formula whole(dimension, true);				
		file << satSolver->toString(whole);
		for(int i = 0; i < dimension; i++)
			if(explorer->mus_union[i])
				file << "c " << i << endl;
		file.close();
		stringstream cmd;
		cmd << "python cadet.py exp.cnf";
		string result = exec(cmd.str().c_str());
		if (result.find("UNSAT") != string::npos) {
			vector<string> res = split(result);
			assert(res[0] == "UNSAT:");
			for(int i = 0; i < dimension; i++){
				seed[i] = !explorer->mus_union[i];
			}
			for(int i = 1; i < res.size(); i++){
				seed[stoi(res[i])] = true;
			}
		}else if(result.find("SAT") != string::npos){
			break;
		}else{
			print_err("unexpected return value from cadet.py");
		}
	}
}
