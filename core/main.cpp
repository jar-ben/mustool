#include <ctime>
#include <time.h>
#include "Master.h"
#include <getopt.h>
#include <algorithm>
#include <csignal>
#include <experimental/filesystem>
namespace fs = std::filesystem;

string add_line_breaks(string s, int len = 50){
	int pos = len;
	while(1){
		if(pos > s.size() - 1)
			break;
		while(pos < s.size() && s[pos] != ' ')
			pos++;
		if(pos < s.size())
			s.replace(pos, 1, "\n\t\t");		
		pos += len;
	}
	return s;
}

void print_help(){
	int args = 14;
		string help[args][4] = 
		{
			{"variant", "v", "VARIANT", "Algorithm to be used: 1 - ReMUS, 2 - TOME, 3 - MARCO. (default: 1)"},
			{"output-file", "o", "OUTPUT_FILE", "Output file. (default: None)"},
			{"verbose", "i", "", "Verbose output. (default: None)"},
			{"sat-solver", "n", "", "Choose satisfiability solver. Depends on the constraint domain. In the SAT domain: minisat (default). In the SMT domain: z3 (default). In the LTL domain: spot (default). [Currently, we support only one satisfiability solver in each domain.]"},
			{"shrink-alg", "m", "ALG", "Choose algorithm for performing the shrink procedure. Depends on the constraint domain. In the SAT domain: muser (default), dmuser, custom. In the SMT domain: custom (default). In the LTL domain: 1 - custom (default)"},
			{"max-recursion-depth", "d", "DEPTH", "ReMUS only. Set the limit on the depth of recursion calls(maximal level of nesting). (default: 6)"},
			{"dim-reduction-factor", "r", "FACTOR", "ReMUS only. Factor of dimensionality reduction. For example, 0.8 means reduce the dimension by 20 percent in every recursion call. (default: 0.9)"},

			{"verify-muses", "c", "","Experimental option. Verify each outputted MUS for being a MUS."},
			{"block-down-muses", "a", "", "Experimental option."},
			{"mus-approximation", "s", "", "Experimental option. Enable MUS approximation."},
			{"mus-approximation-treshold", "e", "TRE", "Experimental option. Threshold used to determine whether to use MUS approximation or not. (default: 0)"},
			{"verify-approximated-muses", "f", "", "Experimental option. Check whether are approximated MUSes actually MUSes."},
			{"heuristic", "w", "", "Only for SAT domain. Enable either Backbone (-w 1) or Matchmaker (-w 2) heuristic."},
			{"help", "h", "", "Prints this help message."}
        	};
	cout << "usage: ./mvc <arguments> <input_file_name>" << endl;
	cout << "Available arguments:" << endl;
	int line_len;
	for(auto &opt: help){
		if(line_len > 50){
			cout << endl << "             ";
			line_len = 0;
		}
		cout << "[-" << opt[1] << ((opt[2].size() > 0)? " " : "") << opt[2] << "] ";
		line_len += 5 + opt[2].size();
	}
	cout << endl << endl;

	cout << "mandatory arguments:" << endl;
	cout << "\t input_file_name \t\t input file, allowed formats: .cnf, .smt2, .ltl (see example input files in ./examples/)" << endl << endl;

	cout << "optional arguments:" << endl;


	for(auto &opt: help){
		cout << "\t-" << opt[1] << ", --" << opt[0] << " " << opt[2] << endl << "\t\t" << add_line_breaks(opt[3]) << endl << endl;
	}

	cout << endl << "example usage: ./mvc -v 1 ./examples/test1.cnf" << endl << endl; 

}

void signal_handler(int signal){
	if(signal == 15){
		cout << "the computation was interupted due to a timeout" << endl;
		for (const auto & entry : fs::directory_iterator("./"))
			std::cout << entry.path() << std::endl;
	}
	exit(1);
}

int main(int argc, char *argv[]){

	std::signal(SIGTERM, signal_handler);

	int opt;
        int variant = 1;
        bool verbose = false;
        string input, output = "";
	int depthMUS = 6;
	float dim_reduction = 0.9;
	int extend_top_variant = 1;
	bool hsd = false;
	bool block_down_mus = false;
	bool validate_mus = false;
	float mus_approx = 0;
	bool verify_approx = false;	
	string sat_solver = "default";
	string shrink_alg = "default";
	bool visualize = false;
	bool model_rotation = false;
	bool get_implies = false;
	float crits_treshold = 0.9;
	bool criticals_rotation = false;
	bool unex_impl = false;
	bool test_rotation_unex = false;
	int mus_size_limit = 10000;
	int scope_limit = 1000000;
	bool useBackbone = false;
	bool useMatchmaker = false;

	int c;
	while(1){
		static struct option long_options[] =
		{
			{"variant",     required_argument,       0, 'v'},
			{"output-file",  required_argument,       0, 'o'},
			{"verbose",  no_argument,       0, 'i'},
			{"sat-solver", required_argument,	0, 'n'},
			{"shrink-alg",  required_argument,       0, 'm'},
			{"max-recursion-depth",  required_argument,       0, 'd'},
			{"extend-top-variant",  required_argument,       0, 'j'},
			{"dim-reduction-coef",  required_argument,       0, 'r'},
			{"mus-approximation",  no_argument,       0, 's'},
			{"block-down-muses",  no_argument,       0, 'a'},
			{"verify-muses",  no_argument,       0, 'c'},
			{"mus-approximation-treshold",  required_argument,       0, 'e'},
			{"verify-approximated-muses",  no_argument,       0, 'f'},
			{"help",  no_argument,       0, 'h'},
			{"visualize",  no_argument,       0, 't'},
			{"model-rotation", no_argument, 	0, 'x'},
			{"heuristic", required_argument, 0, 'w'},
			{0, 0, 0, 0}
        	};

		int option_index = 0;
		c = getopt_long (argc, argv, "l:v:o:im:d:j:b:r:sace:fhn:xygt:zkqp:uw:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;
		switch (c){
			case 'l':
				mus_size_limit = atoi(optarg);
				break;
			case 'a':
				block_down_mus = true;
				break;
			case 'k':
				test_rotation_unex = true;
				break;
			case 'z':
				unex_impl = true;
				break;
			case 'c':
				validate_mus = true;
				break;
			case 'e':
				mus_approx = stof(optarg);
				break;
			case 'f':
				verify_approx = true;
				break;
			case 's':
				hsd = true;
				break;
			case 'j':
				extend_top_variant = atoi(optarg);
				break;
			case 'r':
				dim_reduction = stof(optarg);
				break; 
                        case 'v':
                                variant = atoi(optarg);
                                break;
                        case 'o':
                                output = optarg;
                                break;
                        case 'i':
                                verbose = true;
                                break;
			case 'm':
				shrink_alg = optarg;
				break;
			case 'n':
				sat_solver = optarg;
				break;

			case 'd':
				depthMUS = atoi(optarg);
				break;
			case 'p':
				scope_limit = atoi(optarg);
				break;
			case 'w':
				useBackbone = atoi(optarg) == 1;
				useMatchmaker = atoi(optarg) == 2;
				break;
			case 'x':
				model_rotation = true;
				break;
			case 'y':
				criticals_rotation = true;
				break;
			case 'g':
				get_implies = true;
				break;
			case 't':
				crits_treshold = stof(optarg);
				break;
                        case 'h':
				print_help();
				return 0;
			}
	}


        if(optind >= argc) {
		return 0; //no input file specified
        }

        input = argv[optind];

	if(output != ""){
		std::ofstream ofs;
		ofs.open(output, std::ofstream::out | std::ofstream::trunc);
		ofs.close();
	}

	Master solver(input, variant, visualize, sat_solver);
	solver.output_file = output;
	solver.verbose = verbose;
	solver.depthMUS = (depthMUS >= 0)? depthMUS : solver.dimension;
	solver.dim_reduction = dim_reduction;
	solver.extend_top_variant = extend_top_variant;
	solver.hsd = hsd;
	solver.block_down_mus = block_down_mus;
	solver.validate_mus_c = validate_mus;
	solver.mus_approx = mus_approx;
	solver.verify_approx = verify_approx;
	solver.satSolver->sat_solver = sat_solver;
	solver.satSolver->shrink_alg = shrink_alg;
	solver.model_rotation = model_rotation;
	solver.get_implies = get_implies;
	solver.crits_treshold = crits_treshold;
	solver.criticals_rotation = criticals_rotation;
	solver.explorer->use_implications = unex_impl;
	solver.explorer->test_rotation_unex = test_rotation_unex;
	solver.explorer->mus_size_limit = mus_size_limit;
	solver.scope_limit = scope_limit;
	solver.useBackbone = useBackbone;
	solver.useMatchmaker = useMatchmaker;
	solver.enumerate();

	
	cout << "Enumeration completed" << endl;
       	cout << "Number of MUSes: " << solver.muses.size() << endl;

	return 0;
}

