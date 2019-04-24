#include <errno.h>
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "satSolvers/MSHandle.h"

#include "mcsmus/minisat-wrapper.hh"
#include "mcsmus/glucose-wrapper.hh"
#include "mcsmus/minisolver.hh"
#ifdef USE_LINGELING
#include "mcsmus/lingeling-solver.hh"
#endif

#include "mcsmus/parseutils.hh"
#include "mcsmus/options.hh"
#include "mcsmus/mcsmus.hh"
#include "mcsmus/dimacs.hh"
#include "mcsmus/system.hh"

//=================================================================================================
// Main:

using namespace mcsmus;

IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", 0, IntRange(0, INT32_MAX));
IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", 0, IntRange(0, INT32_MAX));
BoolOption   strictp("MAIN", "strict", "Validate DIMACS header during parsing.", false);

BoolOption   print_core("MAIN", "print-core", "Print core as CNF to <input>.min.cnf", false);
BoolOption   print_gmus("MAIN", "print-gmus", "Print gMUS as GCNF to <input>.min.gcnf", false);
BoolOption   print_core_idx("MAIN", "print-core-indices", "Print indices of core clauses to stdout", false);
BoolOption   all_cores("MAIN", "all-cores", "Find all cores", false);
BoolOption   some_cores("MAIN", "some-cores", "Find more than one core (no guarantee of completeness)", false);
BoolOption   all_mcs_times("MAIN", "all-mcs-times", "Print times for each MCS to stdout", false);
BoolOption   all_mus_times("MAIN", "all-mus-times", "Print times for each MUS to stdout", false);
BoolOption   all_times("MAIN", "all-times", "Implies -all-mcs-times and -all-mus-times", false);

BoolOption   opt_print_mcses("MAIN", "print-mcses", "Print each MCS as a line of clause indices to <input>.mcs", false);
BoolOption   opt_print_muses("MAIN", "print-muses", "Print each MUS as a line of clause indices to <input>.mus", false);
BoolOption   opt_print_msses_cnf("MAIN", "print-msses-cnf", "Print each MSS as a CNF to <input>.mss.<num>.cnf", false);
BoolOption   opt_print_muses_cnf("MAIN", "print-muses-cnf", "Print each MUS as a CNF to <input>.mus.<num>.cnf", false);

#ifdef USE_LINGELING
StringOption opt_backend("MAIN", "backend", "SAT engine to use [minisat[-simp]/glucose[-simp]/lingeling", "glucose-simp");
#else
StringOption opt_backend("MAIN", "backend", "SAT engine to use [minisat[-simp]/glucose[-simp]", "glucose-simp");
#endif

IntOption    opt_status_interval   ("MAIN", "status-interval", "How often we print the current status [0 = never]",  10, IntRange(0, INT32_MAX));

/* basename("/path/to/file.ext", ".ext") returns "file" */
std::string basename2(std::string fn, std::string extension)
{
    using std::string;
    size_t dotpos = fn.rfind(extension);
    if( dotpos < fn.size() - extension.size() )
        dotpos = string::npos;
    size_t dirpos = fn.rfind('/');
    if( dirpos == string::npos)
        dirpos = 0;
    else
        ++dirpos;
    return fn.substr(dirpos, dotpos-dirpos);
}

void do_solve2(
		std::vector<bool> &formula,
		std::vector<bool> &crits,
     		std::string inf,
       		mcsmus::Wcnf& w,
              mcsmus::MUSSolver& S)
{
    double initial_time = cpuTime();

    Control* control = getGlobalControl();;

    control->verbosity = verb;
    control->status_interval = opt_status_interval;

    S.single_mus = true;

    bool find_all = bool(all_cores) || bool(some_cores);
    if( find_all ) {
        S.single_mus = false;
        S.all_muses = true;
    }


    using cs_cb = mcsmus::MUSSolver::cs_callback;
    using mus_cb = mcsmus::MUSSolver::mus_callback;

    using mcsmus::Lit;
    using mcsmus::lbool;

    int num_muses{0}, num_mcses{0};

    bool print_mcses = opt_print_mcses,
        print_muses = opt_print_muses,
        print_muses_cnf = opt_print_muses_cnf,
        print_msses_cnf = opt_print_msses_cnf;

    std::ofstream mcs_out, mus_out;

    std::vector<int> mcs_indices;
    cs_cb printmcses = [&](std::vector<Lit>const & mcs, model_t const&,
                           bool) {
        mcs_indices.clear();
        for(Lit l : mcs)
            mcs_indices.push_back(w.bvarIdx(l));
        sort(begin(mcs_indices), end(mcs_indices));
        for(int cl : mcs_indices)
            mcs_out << cl << ' ';
        mcs_out << "0" << std::endl;
    };
    cs_cb printmssescnf = [&](std::vector<Lit> const& mcs, model_t const&, bool) {
        std::ofstream ofs( basename2(inf, ".cnf") + ".mss."
                           + std::to_string(num_mcses) + ".cnf" );
        S.print_mss(mcs, ofs);
    };
    std::vector<cs_cb> mcscb;

    if( print_mcses ) {
        mcs_out.open( basename2(inf, ".cnf")+".mcs" );
        mcscb.push_back(printmcses);
    }

    if( print_msses_cnf )
        mcscb.push_back(printmssescnf);

    std::vector<int> mus_indices;
    mus_cb printmuses = [&](std::vector<Lit> const& mus, bool) {
        mus_indices.clear();
        for(Lit l : mus)
            mus_indices.push_back(w.bvarIdx(l));
        sort(begin(mus_indices), end(mus_indices));
        for(int cl : mus_indices)
            mus_out << cl << ' ';
        mus_out << "0" << std::endl;
    };
    mus_cb printmusescnf = [&](std::vector<Lit> const& mus, bool) {
        std::ofstream ofs( basename2(inf, ".cnf") + ".mus."
                           + std::to_string(num_muses) + ".cnf" );
        S.print_mus(mus, ofs);
    };
    std::vector<mus_cb> muscb;

    if( print_muses ) {
        mus_out.open( basename2(inf, ".cnf")+".mus" );
        muscb.push_back(printmuses);
    }

    if( print_muses_cnf )
        muscb.push_back(printmusescnf);

    if( all_mcs_times || all_times ) {
        mcscb.push_back([&](std::vector<Lit> const&, model_t const&,
                            bool) {
                printf("S %d : %f\n", num_mcses, cpuTime());
                fflush(stdout);
            });
    }

    if( all_mus_times || all_times ) {
        muscb.push_back([&](std::vector<Lit> const&, bool) {
                printf("U %d : %f\n", num_muses, cpuTime());
                fflush(stdout);
            });
    }

    if( !mcscb.empty() )
        S.register_cs_callback([&](std::vector<Lit> const& mcs, model_t const& model,
                                   bool minimal) {
            if(!minimal)
                return;
            ++num_mcses;
            for(auto& cb : mcscb)
                cb(mcs, model, minimal);
        });

    if( !muscb.empty() )
        S.register_mus_callback([&](std::vector<Lit> const& mus, bool minimal) {
                if( !minimal )
                    return;
                ++num_muses;
                for(auto& cb : muscb)
                    cb(mus, minimal);
        });

    S.initializeControl();
    control->force_exit_on_interrupt();

    // Try to set resource limits:
    if (cpu_lim != 0) limitTime(cpu_lim);
    if (mem_lim != 0) limitMemory(mem_lim);


    gzFile in = gzopen(inf.c_str(), "rb");

    auto filetype = parse_auto_detect(in, w, strictp);
    if( filetype == input_type::CNF_FILE )
        printf("c CNF on input -- MUS problem\n");
    else
        printf("c GCNF on input -- GMUS problem\n");
    gzclose(in);
    double parsed_time = cpuTime();

    if (control->verbosity > 0)
        S.printHeader(parsed_time-initial_time);

    control->notify_on_interrupt();

    vector<Lit> mus;
    w.relax();

    bool ret;
    if( !find_all )
        ret = S.find_mus(mus, false);
    else {
        if( all_cores )
            ret = S.find_all_mus();
        else
            ret = S.find_some_mus();
    }

    if (control->verbosity >= 1)
        S.printFooter();

    if (control->verbosity > 0){
        S.printStats();
        printf("\n"); }

    if( find_all ) {
        std::cout << (ret ? "ALL MUSES GENERATED\n" : "INCOMPLETE GENERATION\n");
        std::cout << S.num_muses << " MUSes found\n";
        return;
    }

    std::cout << (ret ? "MINIMAL\n" : "INCOMPLETE MINIMIZATION\n");
    std::cout << "MUS size " << mus.size() << "/" << (w.groups().size()-1) << "\n";
}

std::vector<bool> MSHandle::shrink_mcsmus(std::vector<bool> &f, std::vector<bool> crits){
    try {
        setUsageHelp("USAGE: %s [options] <input-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        setX86FPUPrecision();

    std::stringstream exp;
	int hash = 123;
    exp << "./f_" << hash << ".cnf";
    std::string inf = exp.str();
	export_formula_crits(f, exp.str(), crits);

        Wcnf wcnf;
        std::unique_ptr<BaseSolver> s;
        s = make_glucose_simp_solver();

        MUSSolver mussolver(wcnf, s.get());
        do_solve2(f, crits, inf, wcnf, mussolver);

        return std::vector<bool>();
    } catch (std::bad_alloc&){
        printf("INDETERMINATE\n");
        exit(0);
	return std::vector<bool>();
    }
}
