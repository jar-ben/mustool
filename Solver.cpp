#include "Solver.h"
#include "misc.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <random>


Solver::Solver(string filename, int var, bool vis, string s_solver){
        isValidExecutions = 0;
        variant = var;
	sat_solver = s_solver;
	if(ends_with(filename, "smt2")){
		satSolver = new Z3Handle(filename); 
		domain = "smt";
	}
	else if(ends_with(filename, "cnf")){
		satSolver = new MSHandle(filename);
		domain = "sat";
	}
	else if(ends_with(filename, "ltl")){
		if(sat_solver == "nuxmv")
			satSolver = new NuxmvHandle(filename); 
		else
			satSolver = new SpotHandle(filename);
		domain = "ltl";
	}
	else{
		print_err("wrong input file");
		exit(1);
	}
	dimension = satSolver->dimension;	
	cout << "Dimension:" << dimension << endl;
        explorer = new Explorer(dimension);	
	explorer->satSolver = satSolver;
        verbose = false;
	depthMUS = 0;
	extend_top_variant = 0;
	dim_reduction = 0.5;
	output_file = "";
	validate_mus_c = false;
	hsd = false;
	hsd_succ = 0;
	verify_approx = false;
	total_difference_of_approximated_muses = 0;
	overapproximated_muses = rightly_approximated_muses = 0;
	current_depth = 0;
	explored_extensions = 0;
	spurious_muses = 0;
	skipped_shrinks = 0;
	unex_sat = unex_unsat = known_unex_unsat = extracted_unex = 0;
	rotated_muses = 0;
	duplicated = 0;
	max_round = 0;
	explicit_seeds = 0;

	scope_limit = 100;	

	use_edge_muses = true;
	rot_muses.resize(dimension, unordered_map<int, vector<int>>());

	backbone_crit_value.resize(dimension, false);
	backbone_crit_used.resize(dimension, false);

        std::random_device rd;
        std::mt19937 eng(rd());
        std::uniform_int_distribution<> distr(1, 1000000000);
        hash = distr(eng);
	satSolver->hash = hash;
	used_backbones_count = 0;
	original_backbones_count = 0;
}

Solver::~Solver(){
	delete satSolver;
	delete explorer;
}

void Solver::write_mus_to_file(MUS& f){
	if(satSolver->clauses_string.size() != dimension){
		print_err("clauses_string error");
		return;
	}
	ofstream outfile;
	outfile.open(output_file, std::ios_base::app);
	for(auto c: f.int_mus)
		outfile << satSolver->clauses_string[c] << endl;
	outfile << endl;
	outfile.close(); 
}

// mark formula and all of its supersets as explored
void Solver::block_up(MUS& formula){
	explorer->block_up(formula);
}

// mark formula and all of its subsets as explored
void Solver::block_down(Formula formula){
        explorer->block_down(formula);
}

// experimental function
// allows output overapproximated MUSes 
bool Solver::check_mus_via_hsd(Formula &f){
	Formula bot = explorer->get_bot_unexplored(f);
	int size_bot = count_ones(bot);
	int size_f = count_ones(f);
	if(verbose) std::cout << "HSD: " << (size_f - size_bot) << ", hsd_succ: " << hsd_succ << endl << endl; 
	cout << "hsdd: " << (size_f - size_bot) << endl;
	bool result =  false;

	if(bot == f){
		result = true;
		cout << "guessed MUS" << endl;
	}

	else if((size_f - size_bot) < (dimension * mus_approx)){
		result = true;
		cout << "approximated MUS" << endl;
	}

	if(result){
		hsd_succ++;
		if(verify_approx){		
			Formula mus = satSolver->shrink(f);
			int f_size = count_ones(f);
			int mus_size = count_ones(mus);
			int distance = (f_size - mus_size);
			if(distance == 0){
				rightly_approximated_muses++;
				cout << "MUS approximation: success" << endl;				
			}
			else{
				overapproximated_muses++;
				total_difference_of_approximated_muses += distance;	
				cout << "MUS approximation: overapproximated, difference in size (over real MUS): " << distance << endl;
			}

			if(verbose){
				cout << "Total number of approximations: " << (rightly_approximated_muses + overapproximated_muses) << endl;
				cout << "Right approximations: " << rightly_approximated_muses << endl;
				cout << "Overapproximations: " << overapproximated_muses << endl;
				if(overapproximated_muses > 0){
					cout << "ratio right approximations / overapproximations: " << (rightly_approximated_muses / (float) overapproximated_muses) << endl;
					cout << "Average diff of overapproximation: " << (total_difference_of_approximated_muses / (float) overapproximated_muses) << endl;
				}	
			}
		}
	}
	return result;	
}

// check formula for satisfiability
// core and grow controls optional extraction of unsat core and model extension (replaces formula)
bool Solver::is_valid(Formula &formula, bool core, bool grow){
	return satSolver->solve(formula, core, grow);
}

//verify if f is a MUS
bool Solver::validate_mus(Formula &f){
	if(is_valid(f)){
		print_err("the mus is SAT");
		exit(1);
	}
	if(!explorer->checkValuation(f)){
		print_err("this mus has been already explored");
		exit(1);
		duplicated++;
	}	
	for(int l = 0; l < f.size(); l++)
		if(f[l]){
			f[l] = false;
			if(!is_valid(f)){
				print_err("the mus has unsat subset");	
				satSolver->export_formula(f, "spurious_mus");
				exit(1);
				spurious_muses++;				
				return false;
			}
			f[l] = true;
		}	
	return true;
}

MUS& Solver::shrink_formula(Formula &f, Formula crits){
	int f_size = count_ones(f);
	chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();
	if(crits.empty()) crits = explorer->critical;
	if(get_implies){ //get the list of known critical constraints	
		explorer->get_implies(crits, f);	
		cout << "criticals before rot: " << count_ones(crits) << endl;	
		if(criticals_rotation) satSolver->criticals_rotation(crits, f);
		cout << "criticals after rot: " << count_ones(crits) << endl;
		float ones_crits = count_ones(crits);		 
		if(f_size == ones_crits){ // each constraint in f is critical for f, i.e. it is a MUS 
			skipped_shrinks++; 
			muses.push_back(MUS(f, -1, muses.size(), f_size)); //-1 duration means skipped shrink
			return muses.back();
		}		
		if((ones_crits/f_size) > crits_treshold){
			if(!is_valid(crits, false, false)){ 
				skipped_shrinks++; 
				muses.push_back(MUS(crits, -1, muses.size(), f_size));//-1 duration means skipped shrink
				return muses.back();
			}
		}
	}

	if(use_core_backbones && count_ones(explorer->critical) > explorer->core.size())
		core_backbones();
	Formula mus = satSolver->shrink(f, crits);
	chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>( end_time - start_time ).count() / float(1000000);
	muses.push_back(MUS(mus, duration, muses.size(), f_size));
	return muses.back();
}

//grow formula into a MSS
Formula Solver::grow_formula(Formula &f){
	return satSolver->grow(f);
}


void Solver::mark_MUS(MUS& f, bool block_unex){	
	if(validate_mus_c) validate_mus(f.bool_mus);		
	explorer->block_up(f, block_unex);

	chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>( now - initial_time ).count() / float(1000000);
        cout << "Found MUS #" << muses.size() <<  ", mus dimension: " << f.dimension;
	cout << ", checks: " << satSolver->checks << ", time: " << duration;
//	cout << ", sat rotated: " << explorer->sat_rotated;
	cout << ", unex sat: " << unex_sat << ", unex unsat: " << unex_unsat << ", criticals: " << explorer->criticals;
	cout << ", intersections: " << std::count(explorer->mus_intersection.begin(), explorer->mus_intersection.end(), true);
	cout << ", rotated MUSes: " << rotated_muses << ", explicit seeds: " << explicit_seeds;
	cout << ", union: " << std::count(explorer->mus_union.begin(), explorer->mus_union.end(), true) << ", dimension: " << dimension;
	cout << ", seed dimension: " << f.seed_dimension << ", duration: " << f.duration;
       	cout << ", used b: " << used_backbones_count << ", original b: " << original_backbones_count << endl;
	cout << ((original_backbones_count > 0)? (used_backbones_count/float(original_backbones_count)) : 0 ) << endl;

	if(output_file != "")
		write_mus_to_file(f);
}


int Solver::backbone_dual_mus(MUS &m1, Formula &top){
	if(!explorer->easy_to_shrink(m1)) return 0;
	int rotated = 0;
	int mid_limit = muses.size() - scope_limit;
	if(mid_limit < 0) mid_limit = 0;
	int vars = satSolver->vars;	

	vector<int> top_int;
	for(int i = 0; i < dimension; i++)
		if(!top[i]) top_int.push_back(i);

	vector<int> indicator(muses.size(), 0);
	for(int mid = muses.size() - 1; mid > mid_limit; mid--){
		auto &m = muses[mid];
		for(auto c2: top_int)
			if(m.bool_mus[c2])
				indicator[mid]++;
	}

	vector<vector<int>> local_muses_overlaps;
	vector<int> local_muses_ids;
	unordered_set<vector<int>, OverlapHash> overlaps;	
	vector<bool> already_used_mus(muses.size(), false);

	vector<int> variables_map;
	vector<int> variables_map_inv;
	backbone_build_literal_map(m1.bool_mus, variables_map, variables_map_inv);

	//just for experimental statistics
	int round = 0;

	for(auto c1: m1.without_crits){
		if(explorer->mus_intersection[c1]) continue;

		//compute backbones
		Formula implied(vars, false);
		Formula values(vars, false);
		vector<int> lits_left (dimension,0);
		vector<bool> satisfied(dimension, false);
		vector<int> singletons;
		Formula seed = m1.bool_mus; seed[c1] = false;
		backbone_simplify(seed, c1, implied, values);	
		get_backbones(seed, implied, values, variables_map, variables_map_inv);

		vector<int> pairs;	
		for(auto c: top_int){
			if(explorer->is_backbone_pair(c, implied, values)){ pairs.push_back(c); }
		}		
		for(auto c2: pairs){			
			
			Formula tmp = top;
			tmp[c1] = false;

			for(auto mid: explorer->parent_muses[c2]){
				if(indicator[mid] != 1) continue;// || already_used_mus[mid]) continue;
				auto &m2 = muses[mid];
				if(m2.bool_mus[c1] || !explorer->easy_to_shrink(m2)) continue;
				round++;
//				cout << "round: " << round << endl;
 				vector<int> m1_overlap;
				for(auto o: m2.without_crits)
					if(o != c2 && !m1.bool_mus[o])
						m1_overlap.push_back(o);
	//			if(overlaps.find(m1_overlap) != overlaps.end()) continue;
				overlaps.insert(m1_overlap);
				bool unexplored = true;					
				for(int i = 0; i < local_muses_overlaps.size(); i++){
					if(muses[local_muses_ids[i]].bool_mus[c1]) continue; //hitten by c1
					auto &overlap = local_muses_overlaps[i];
					bool hit = false;
					for(auto l: overlap){
						if(!m2.bool_mus[l]){
							hit = true; break; 
						}
					}
					if(!hit){ unexplored = false; break; }
				}
				if(!unexplored) continue;	
				seed = m1.bool_mus; seed[c1] = false;			
				for(auto o: m1_overlap) seed[o] = true;

				MUS m = shrink_formula(seed);
				rotated_muses++;
				mark_MUS(m);
				indicator.push_back(0);
				rotated++;
				already_used_mus[mid] = true;
				//local mus overlap
				vector<int> overlap;
				for(auto l: m.without_crits)
					if(!m1.bool_mus[l])
						overlap.push_back(l);
				local_muses_overlaps.push_back(overlap);
				local_muses_ids.push_back(muses.size() - 1);
			}
		}		
	}
//	if(rotated > 0)
//		cout << "overlaps: " << overlaps.size() << ", rotated: " << rotated << ", " << (float(rotated) / overlaps.size()) <<  endl;
	return rotated;
}


void Solver::enumerate(){
	initial_time = chrono::high_resolution_clock::now();
	cout << "running algorithm variant " << variant << endl;
	Formula whole(dimension, true);
	if(is_valid(whole)){
		cout << "the input instance is satisfiable" << endl;
		exit(1);
	}

	switch (variant) {
		case 1:
		case 10:
		case 11:
		case 12:
			find_all_muses_duality_based_remus(Formula (dimension, true), Formula (dimension, false), 0);
			break;						
		case 2:
		case 21:
		case 22:
			find_all_muses_tome();
			break;
		case 30:
		case 31:
		case 32:
		case 33:
		case 34:
			marco_base();
			break;
		default:
			print_err("invalid algorithm chosen");
			break;
	}
	return;
}

