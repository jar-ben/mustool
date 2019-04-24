#include "Master.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


Master::Master(string filename, int var, bool vis, string s_solver){
        isValidExecutions = 0;
        variant = var;
	sat_solver = s_solver;
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
		#ifdef NOSAT
			std::cout << "working with sat is currently not enabled, plese build the tool with the flag USESAT=YES, e.g. 'make USESAT=YES'" << std::endl;
			exit(1);
		#else
		satSolver = new MSHandle(filename);
		#endif
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
	explorer->satSolver = satSolver;
        verbose = false;
	depthMUS = 0;
	dim_reduction = 0.5;
	output_file = "";
	validate_mus_c = false;
	current_depth = 0;
	unex_sat = unex_unsat = 0;
	rotated_muses = 0;
	scope_limit = 100;	
        hash = random_number();
	satSolver->hash = hash;
	useBackbone = useMatchmaker = false;
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
	return satSolver->solve(formula, core, grow);
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
	cout << "shrinking dimension: " << f_size << endl;
	chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();
	if(crits.empty()) crits = explorer->critical;
	if(get_implies){ //get the list of known critical constraints	
		explorer->getImplied(crits, f);	
		cout << "criticals before rot: " << count_ones(crits) << endl;	
		if(criticals_rotation){
			MSHandle *msSolver = static_cast<MSHandle*>(satSolver);
			msSolver->criticals_rotation(crits, f);
		}
		cout << "criticals after rot: " << count_ones(crits) << endl;
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


void Master::mark_MUS(MUS& f, bool block_unex){	
	if(validate_mus_c) validate_mus(f.bool_mus);		
	explorer->block_up(f);

	chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
        cout << "Found MUS #" << muses.size() <<  ", mus dimension: " << f.dimension;
	cout << ", checks: " << satSolver->checks << ", time: " << duration;
	cout << ", unex sat: " << unex_sat << ", unex unsat: " << unex_unsat << ", criticals: " << explorer->criticals;
	cout << ", intersections: " << std::count(explorer->mus_intersection.begin(), explorer->mus_intersection.end(), true);
	cout << ", rotated MUSes: " << rotated_muses;
	cout << ", union: " << std::count(explorer->mus_union.begin(), explorer->mus_union.end(), true) << ", dimension: " << dimension;
	cout << ", seed dimension: " << f.seed_dimension << ", duration: " << f.duration;
	cout << ", shrinks: " << satSolver->shrinks;
	cout << endl;

	if(output_file != "")
		write_mus_to_file(f);
}

void Master::enumerate(){
	initial_time = chrono::high_resolution_clock::now();
	cout << "running algorithm variant " << variant << endl;
	Formula whole(dimension, true);
	if(is_valid(whole))
		print_err("the input instance is satisfiable");

	switch (variant) {
		case 1:
			find_all_muses_duality_based_remus(Formula (dimension, true), Formula (dimension, false), 0);
			break;						
		case 2:
			find_all_muses_tome();
			break;
		case 3:
			marco_base();
			break;
		default:
			print_err("invalid algorithm chosen");
			break;
	}
	return;
}

