#include <ctime>
#include <time.h>
#include "Master.h"
#include <getopt.h>
#include <algorithm>
#include <csignal>
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraint.h"

void signal_handler(int signal){
	if(signal == 15){
		cout << "the computation was interupted due to a timeout" << endl;
	}
	exit(1);
}

int main(int argc, char *argv[]){

	std::signal(SIGTERM, signal_handler);

	int scope_limit = 1000000;
	string input = "examples/bf1355-228.cnf";


	try{
		TCLAP::CmdLine cmd("Domain Agnostic MUS Enumeration Tool (DAMUSET), Jaroslav Bendik, 2019.", ' ', "");
		vector<string> allowedAlgs {"remus", "tome", "marco"};
		TCLAP::ValuesConstraint<string> allowedVals(allowedAlgs);
		TCLAP::ValueArg<string> algorithm("a","algorithm","MUS enumeration algorithm to be used.",false,"remus",&allowedVals);
		cmd.add(algorithm);

		TCLAP::ValueArg<std::string> output("o","output-file","A file where the identified MUSes will be exported to.",false,"","string");
		cmd.add(output);
		TCLAP::SwitchArg verbose("v","verbose","Verbose output", cmd, false);
		TCLAP::ValueArg<std::string> shrink("s","shrink","Specifies the shrinking algorithm (single MUS extraction subroutine).",false,"default","string");
		cmd.add(shrink);
		TCLAP::ValueArg<int> recursionDepthLimit("","max-recursion-depth","Affects only the algorithm ReMUS. Sets the depth limit on recursion calls in the algorithm.",false,6,"N+ or -1 for unlimited.");
		cmd.add(recursionDepthLimit);
		TCLAP::ValueArg<float> reductionCoeff("","dimension-reduction-coefficient","Affects only the algorithm ReMUS. Sets the dimension reduction coefficient used for the recursion calls of the algorithm.",false,0.9,"float 0-1");
		cmd.add(reductionCoeff);
		TCLAP::SwitchArg verify("c","verify-muses","Used for testing purposes. Verify that the outputted MUSes are indeed MUSes.", cmd, false);
		TCLAP::SwitchArg getImplied("g","get-implied","Based on already found correction sets (satisfiable subsets), determines some critical constraints of seeds before shrinking and thus may speed up (or even completely avoid) the shrinking.", cmd, false);

		TCLAP::SwitchArg backbone("","backbone","Available only in the SAT domain. Use the heuristic Backbone for finding seeds without explicitely checking subsets for satisfiability.", cmd, false);
		TCLAP::SwitchArg matchmaker("","matchmaker","Available only in the SAT domain. Use the heuristic Matchmaker for finding seeds without explicitely checking subsets for satisfiability.", cmd, false);
		TCLAP::SwitchArg modelRotation("","model-rotation","Available only in the SAT domain. Every time a subset is checked for satisfiability and find to be satisfiable, the corresponding model is rotated to identify additional satisfiable subsets. This technique is similar to the famous (recursive) model rotation for finding aditional critical constraints, i.e. singleton correction sets. In our case, we identify arbitrary correction sets.", cmd, false);

		cmd.parse(argc, argv);

		if(output.getValue() != ""){
			std::ofstream ofs;
			ofs.open(output.getValue(), std::ofstream::out | std::ofstream::trunc);
			ofs.close();
		}

		Master solver(input, algorithm.getValue());
		solver.output_file = output.getValue();
		solver.verbose = verbose.getValue();
		solver.depthMUS = (recursionDepthLimit.getValue() >= 0)? recursionDepthLimit.getValue() : solver.dimension;
		solver.dim_reduction = reductionCoeff.getValue();
		solver.validate_mus_c = verify.getValue();
		solver.satSolver->shrink_alg = shrink.getValue();
		solver.model_rotation = modelRotation.getValue();
		solver.get_implies = getImplied.getValue();
		solver.useBackbone = backbone.getValue();
		solver.useMatchmaker = matchmaker.getValue();
		
		solver.scope_limit = scope_limit;
		solver.criticals_rotation = true; //criticals_rotation;
		
		solver.enumerate();

		
		cout << "Enumeration completed" << endl;
		cout << "Number of MUSes: " << solver.muses.size() << endl;

	}catch (TCLAP::ArgException &e){
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}

