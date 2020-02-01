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

	try{
		TCLAP::CmdLine cmd("domain agnostic MUS enumeration Tool (MUST), Jaroslav Bendik, 2019.", ' ', "");
		vector<string> allowedAlgs {"remus", "tome", "marco", "comarco", "duremus", "unibase", "unibase2", "unimus", "counimus"};
		TCLAP::ValuesConstraint<string> allowedVals(allowedAlgs);
		TCLAP::ValueArg<string> algorithm("a","algorithm","MUS enumeration algorithm to be used.",false,"remus",&allowedVals);
		cmd.add(algorithm);

		vector<string> allowedSolvers {"minisat", "glucose", "cadical", "default"};
		TCLAP::ValuesConstraint<string> allowedSolversVals(allowedSolvers);
		TCLAP::ValueArg<string> satsolver("","sat-solver","SAT-solver to be used.",false,"default",&allowedSolversVals);
		cmd.add(satsolver);

		TCLAP::ValueArg<std::string> output("o","output-file","A file where the identified MUSes will be exported to.",false,"","string");
		cmd.add(output);

		TCLAP::ValueArg<std::string> mcslsArgs("","mcsls-args","",false,"","string");
		cmd.add(mcslsArgs);

		TCLAP::SwitchArg conflictsNegation("","conflicts-negation","Negate known conflicting clauses during satsolver.solve() calls.", cmd, false);
		TCLAP::SwitchArg mssRotation("","mss-rotation","Use mss-rotation technique", cmd, false);
		TCLAP::SwitchArg verbose("v","verbose","Verbose output", cmd, false);
		vector<string> allowedShrinks {"default", "muser"};
		TCLAP::ValuesConstraint<string> allowedValsShrink(allowedShrinks);
		TCLAP::ValueArg<std::string> shrink("s","shrink","Specifies the shrinking algorithm (single MUS extraction subroutine). In the SMT and LTL domain, only the default one is supported. In SAT domain, you can opt between default (implemented as mcsmus) and muser.",false,"default",&allowedValsShrink);
		cmd.add(shrink);
		
		vector<string> allowedGrows {"default", "cmp", "uwr", "combined", "mcsls", "fixpoint"};
		TCLAP::ValuesConstraint<string> allowedValsGrow(allowedGrows);
		TCLAP::ValueArg<std::string> grow("","grow","Specifies the growing algorithm (single MSS/MCS extraction subroutine). In the SMT and LTL domain, only the default one is supported. In SAT domain, you can opt between default (naive) and cmp.",false,"default",&allowedValsGrow);
		cmd.add(grow);

		TCLAP::ValueArg<int> cmpStrategy("","cmp-strategy","CMP grow strategy.",false,1,"1,2,3 or 4");
		cmd.add(cmpStrategy);
		TCLAP::ValueArg<int> recursionDepthLimit("","max-recursion-depth","Affects only the algorithm ReMUS. Sets the depth limit on recursion calls in the algorithm.",false,100,"N+ or -1 for unlimited.");
		cmd.add(recursionDepthLimit);
		TCLAP::ValueArg<float> reductionCoeff("","dimension-reduction-coefficient","Affects only the algorithm ReMUS. Sets the dimension reduction coefficient used for the recursion calls of the algorithm.",false,0.9,"float 0-1");
		cmd.add(reductionCoeff);
		TCLAP::SwitchArg verify("c","verify-muses","Used for testing purposes. Verify that the outputted MUSes are indeed MUSes.", cmd, false);
		TCLAP::SwitchArg getImplied("g","get-implied","Based on already found correction sets (satisfiable subsets), determines some critical constraints of seeds before shrinking and thus may speed up (or even completely avoid) the shrinking. Use this flag to disable it (not recommended).", cmd, true);

		TCLAP::SwitchArg criticalsRotation("","criticals-rotation","Available only in the SAT domain and used only if the flag -g is not set. Allows to find additional critical constraints based on the already found ones. Use this flag to turn the feature on.", cmd, false);

		TCLAP::UnlabeledValueArg<string>  input( "input_file", "Input file, either .cnf, .smt2, or .ltl. See the ./examples/.", true, "", "input_file"  );
		cmd.add(input);
		cmd.parse(argc, argv);

		//clear the output file
		if(output.getValue() != ""){
			std::ofstream ofs;
			ofs.open(output.getValue(), std::ofstream::out | std::ofstream::trunc);
			ofs.close();
		}

		Master solver(input.getValue(), algorithm.getValue(), satsolver.getValue());
		solver.output_file = output.getValue();
		solver.verbose = verbose.getValue();
		solver.depthMUS = (recursionDepthLimit.getValue() >= 0)? recursionDepthLimit.getValue() : solver.dimension;
		solver.dim_reduction = reductionCoeff.getValue();
		solver.validate_mus_c = verify.getValue();		
		solver.conflicts_negation = conflictsNegation.getValue();
		std::string shr = shrink.getValue();
		if(solver.domain != "sat") shr = "default";
		solver.satSolver->shrink_alg = shr;
		solver.satSolver->grow_alg = grow.getValue();
		solver.get_implies = getImplied.getValue();
		solver.criticals_rotation = criticalsRotation.getValue(); //criticals_rotation;
		solver.satSolver->growStrategy = cmpStrategy.getValue();
		solver.satSolver->mcslsArgs = mcslsArgs.getValue();
		solver.mss_rotation = mssRotation.getValue();
		solver.enumerate();
		
		cout << "Enumeration completed" << endl;
		cout << "Number of MUSes: " << solver.muses.size() << endl;

	}catch (TCLAP::ArgException &e){
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	return 0;
}

