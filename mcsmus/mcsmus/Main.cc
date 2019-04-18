/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <errno.h>
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>

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
std::string basename(std::string fn, std::string extension)
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

void do_solve(int argc, char **argv,
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

    std::string inf{argv[1]};

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
        std::ofstream ofs( basename(inf, ".cnf") + ".mss."
                           + std::to_string(num_mcses) + ".cnf" );
        S.print_mss(mcs, ofs);
    };
    std::vector<cs_cb> mcscb;

    if( print_mcses ) {
        mcs_out.open( basename(inf, ".cnf")+".mcs" );
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
        std::ofstream ofs( basename(inf, ".cnf") + ".mus."
                           + std::to_string(num_muses) + ".cnf" );
        S.print_mus(mus, ofs);
    };
    std::vector<mus_cb> muscb;

    if( print_muses ) {
        mus_out.open( basename(inf, ".cnf")+".mus" );
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

    if (argc == 1)
        printf("Reading from standard input... Use '--help' for help.\n");

    gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
    if (in == NULL)
        printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

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

    if (print_core) {
        using std::string;
        string out{((argc > 1) ? argv[1] : "")};
        string outf = basename(out, ".cnf") + ".min.cnf";
        std::cout << "printing core to \'" << outf << "\'\n";
        std::ofstream ofs(outf);
        S.print_mus(mus, ofs);
    }
    if( print_gmus ) {
        using std::string;
        string out{((argc > 1) ? argv[1] : "")};
        string outf = basename(out, ".gcnf") + ".min.gcnf";
        std::cout << "printing gMUS to \'" << outf << "\'\n";
        std::ofstream ofs(outf);
        S.print_gmus(mus, ofs);
    }
    if (print_core_idx) {
        std::sort(begin(mus), end(mus));
        std::cout << "v ";
        for (auto b : mus)
            std::cout << w.bvarIdx(b) << " ";
        std::cout << "0\n";
    }
}

int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        setX86FPUPrecision();

        // Extra options:
        //

        parseOptions(argc, argv, true);

        if( argc < 2 ) {
            printUsageAndExit(argc, argv);
        }

        std::string backend((const char*)opt_backend);
        Wcnf wcnf;
        std::unique_ptr<BaseSolver> s;
        if (backend == "minisat")
            s = make_minisat_solver();
        else if (backend == "minisat-simp")
            s = make_minisat_simp_solver();
        else if (backend == "glucose")
            s = make_glucose_solver();
        else if (backend == "glucose-simp")
            s = make_glucose_simp_solver();
#ifdef USE_LINGELING
        else if (backend == "lingeling")
            s = make_lingeling_solver();
#endif
        else {
            std::cout << "Invalid backend \"" << backend << "\"\n";
            exit(1);
        }

        MUSSolver mussolver(wcnf, s.get());
        do_solve(argc, argv, wcnf, mussolver);

#ifdef NDEBUG
        exit(0);
#else
        return 0;
#endif
    } catch (std::bad_alloc&){
        printf("INDETERMINATE\n");
        exit(0);
    }
}
