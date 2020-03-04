/***********[mcsmus.hh]
Copyright (c) 2012-2014 Jessica Davies, Fahiem Bacchus
Copyright (c) 2014-2015 George Katsirelos

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********/

#ifndef MCSMUS_MCSMUS_HH
#define MCSMUS_MCSMUS_HH

#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <algorithm>

#include "mcsmus/basesolver.hh"
#include "minisat/core/mcsmus_SolverTypes.h"
#include "minisat/utils/mcsmus_System.h"
#include "mcsmus/algorithm.hh"

namespace mcsmus {

using std::unique_ptr;
using std::vector;
using std::map;

using Minisat::cpuTime;

using ClsIdx = int;
const ClsIdx ClsIdx_Undef = -1;
class MUSSolver;
struct hash_lit;

class Wcnf
{
public:
    using CRef = Minisat::CRef;
    using Clause = Minisat::Clause;

private:
    Minisat::ClauseAllocator ca;

    // indexed by toInt(external literal), holds CRefs to the hard
    // clauses
    vector<vector<CRef>> hard_occurrences;
    // same, but holds the bvars of the groups in which it occurs
    vector<vector<Lit>> bvar_occurrences;
    // index by bvarIdx, holds CRefs to the soft clauses
    vector<vector<CRef>> groups_;

    // no cheap way to compute this (would have to accumulate sizes in
    // groups_), so store it explicitly
    std::size_t num_clauses{0};

    MUSSolver* solver;

    int firstBVar{-1}, lastBVar{-2};

    Var maxvar{var_Undef};
    std::set<int> distinct_vars;

    // used in getGroupLiterals, here to avoid reallocations
    vector<bool> gather_seen;
public:
    vector<vector<CRef>> const& groups() const { return groups_; }

    // note this is just a convenient container. It is not related to
    // the clauses used during solving
    Clause& getClause(CRef cr) { return ca[cr]; }
    std::vector<CRef> const& getGroupClauses(Lit exbv)
    {
        return groups_[bvarIdx(exbv)];
    }
    std::vector<CRef> const& getGroupClauses(Var exbv)
    {
        return getGroupClauses(mkLit(exbv));
    }

    vector<Lit> const& softBvarOccurrences(Lit p)
    {
        if( bvar_occurrences.size() <= static_cast<unsigned>(toInt(p)) )
            bvar_occurrences.resize(toInt(p)+1);
        return bvar_occurrences[toInt(p)];
    }
    vector<CRef> const& hardOccurrences(Lit p)
    {
        if( hard_occurrences.size() <= static_cast<unsigned>(toInt(p)) )
            hard_occurrences.resize(toInt(p)+1);
        return hard_occurrences[toInt(p)];
    }

    vector<Lit> getGroupLiterals(Lit exbv);
    vector<Var> getGroupVariables(Lit exbv);

    void setSolver(MUSSolver* S);

    // if grpidx == 0, clause is hard (assumed to be implicitly part
    // of every MUS we report). Otherwise it is an optional clause,
    // assigned to that group. If any group has cardinality > 1, the
    // problem is a group-MUS problem. This matches the semantics of
    // the GCNF format.
    void addClause(std::vector<Lit> const& lits, int grpidx);

    int nVars() const { return (lastBVar > 0) ? (lastBVar+1) : (maxvar+1); }
    int nXVars() const { return distinct_vars.size(); }
    int nClauses() const { return num_clauses; }

    void relax() {
        firstBVar = maxvar+1;
        for(auto group : groups_) {
            Var bvar = ++maxvar;
            lastBVar = bvar;
            if( bvarIdx(bvar) == 0 ) // hard clauses here
                continue;
            for( auto cr : group ) {
                for(Lit l : getClause( cr ) ) {
                    if( bvar_occurrences.size() <= unsigned(toInt(l)) )
                        bvar_occurrences.resize(toInt(l)+1);
                    bvar_occurrences[toInt(l)].push_back(mkLit(bvar));
                }
            }
        }
    }

    bool isBvar(Lit l) const { return isBvar(var(l)); }
    bool isBvar(Var v) const { return firstBVar <= v && v <= lastBVar; }

    int bvarIdx(Var bv) const;
    int bvarIdx(Lit bl) const;
    Var varFromBvarIdx(int bvidx) const;

    // return a range containing all bvars corresponding to soft
    // groups (i.e., omits bvar 0, which is the group of hard clauses)
    auto bvars() const { return int_range<int>{firstBVar+1, lastBVar + 1}; }
};

inline bool modelSatisfies(model_t& model, Wcnf& theWcnf, Minisat::CRef cr)
{
    auto& cl = theWcnf.getClause(cr);
    return std::any_of(
        begin(cl), end(cl), [&](Lit l) { return model.value(l) == l_True; });
}

/* The wrapper provides a bunch of typedefs that allow us to use
   different solvers. */
class MUSSolver {
public:
    using Solver = BaseSolver;
    using CRef = Wcnf::CRef;
protected:
    Wcnf& theWcnf;
    BaseSolver* solver;
    unique_ptr<Solver> hssolver; // for all muses only, remains empty for single mus

    Control* control;
public:
  MUSSolver(Wcnf& f, BaseSolver* solver);
  virtual ~MUSSolver() {}

  void initializeControl();

  Solver& satsolver() { return *solver; }

  void relax(vector<Lit>& torelax) { ensureSoftCls(torelax); }

  bool musBudget(vector<Lit>& conflict, bool require_global_mcses, int64_t propBudget) {
    /*******************************************************************
     Reduce a conflict (a set of blits) to a MUS. Try to remove
     elements of conflict in the order given (try removing conflict[0]
     before conflict[1] etc.). So if there is a priority for
     reductions pass a sorted conflict vector.

     Set total propagation budget for computing the MUS.  (-1 = no budget).

     Return reduced conflict vector. Return true if the returned
     vector is a MUS, false if only a partial reduction was done.
    *******************************************************************/
    stime = cpuTime();
    prevTotalTime = totalTime;
    auto val = mus_(conflict, require_global_mcses, propBudget);
    totalTime += cpuTime() - stime;
    stime = -1;
    solves++;
    return val;
  }

    bool musBudget(vector<Lit>& conflict, bool require_global_mcses,
                   double timeLimit)  {
    int64_t propBudget, props;
    props = getProps();
    if(props > 0 && totalTime > 0)
      propBudget = props/totalTime * timeLimit;
    else
      propBudget = 1024*128;

    stime = cpuTime();
    prevTotalTime = totalTime;
    auto val = mus_(conflict, require_global_mcses, propBudget);
    solves++;
    if(props == 0 && !val) {
      auto moreProps = int64_t(getProps()/totalTime*timeLimit - propBudget);
      if (moreProps > propBudget*.5)
        //try again
        val =  mus_(conflict, require_global_mcses, moreProps);
      solves++;
    }
    totalTime += cpuTime() - stime;
    stime = -1;
    return val;
  }

  // minimize the core described in conflict. Returns that core in
  // the same vector. Return true if minimal, false if interrupted
  bool mus(vector<Lit>& conflict, bool require_global_mcses)
    { return  musBudget(conflict, require_global_mcses,
                        static_cast<int64_t>(-1)); }

  // find all cores that are a subset of candidates. Return true if
  // it generated all cores, false if interrupted
  bool allmus(vector<Lit>& candidates);

  // assuming the formula is unsat if all bvars are false, find a
  // minimal core. The result is in mus (ignored as input). Returns
  // true if core is minimal (may be interrupted by resource limits)
  bool find_mus(vector<Lit>& mus, bool require_global_mcses);

  // check if the input is unsat if only the bvars in opt are
  // false. If unsat, minimize the core. Returns <l_False,true> if
  // unsat and opt contains the mus, <l_False, false> if unsat but mus
  // is not minimal (because interrupted or exceeded budget), <l_True,
  // true> if sat, <l_Undef, false> if interrupted or exceeded budget
  // before determining SAT/UNSATness. require_global_mcses is as
  // elsewhere
  std::pair<lbool, bool> solve_mus(vector<Lit>& opt, bool require_global_mcses);

  // as above, but with assumptions. Note the mus generated is not a
  // mus of the entire formula, only of F|_{assumptions}
  std::pair<lbool, bool> solve_mus(vector<Lit>& opt,
      vector<Lit> const& assumptions, bool require_global_mcses);

  // find all cores. Returns true if it
  // generated all cores (may be interrupted)
  bool find_all_mus();

  // find some cores, by making the choice of critical clause in
  // mishmash into a choice point. This algorithm does not guarantee
  // that it will find all muses. Additionally we perform clause set
  // refinement which may remove even more muses. Returns true if not
  // interrupted.
  bool find_some_mus();

  // find an MCS under the current assumptions. Returns l_True if it
  // found an MCS, l_Undef if it was interrupted, l_False if no MCS
  // exists (UNSAT under the assumptions). The second version can do
  // it in a different solver. This is used in allmus for doing
  // destructive updates of the formula. Note that this makes no
  // attempt to map Lits and Vars between mcssolver and this->solver.
  lbool findMCS(vector<Lit> const& conflict, vector<Lit> const& all,
      bool require_global_mcs, bool block_mcs, vector<Lit>& cs);
  lbool findMCS(Solver& mcssolver, vector<Lit> const& conflict,
      vector<Lit> const& all, bool require_global_mcs, bool block_mcs,
      vector<Lit>& cs);

  // same as findMCS, but cs is assumed to contain a CS already which
  // we try to minimize
  lbool findMCS_hint(vector<Lit> const& conflict, vector<Lit> const& all,
      bool require_global_mcs, bool block_mcs, vector<Lit>& cs,
      model_t const& model);
  lbool findMCS_hint(Solver& mcssolver, vector<Lit> const& conflict,
      vector<Lit> const& all, bool require_global_mcs, bool block_mcs,
      vector<Lit>& cs, model_t const& model);

  // as above, but CS may not be minimal and it cannot refine the conflict
  lbool findCS(Solver& solver, vector<Lit> const& conflict, vector<Lit>& cs);

  // print mus as dimacs cnf
  void print_mus(vector<Lit> const& mus, std::ostream& outf);

  // print group MUS as gcnf
  void print_gmus(vector<Lit> const& mus, std::ostream& outf);

  // print MSS (complement of the given MCS) as dimacs cnf
  void print_mss(vector<Lit> const& mcs, std::ostream& outf);

  /* register function to be called when a new CS or US is found.

     A CS is a conjunction of positive literals such that if they are
     all true, the formula is satisfiable.

     A MUS is a disjunction of positive literals such that if they are
     all false, the formula is unsatisfiable.

     Both are reported as vector<Lit>s. A CS also gets the model that
     witnesses satisfiability.

     Note that the CSes given here may not be MCSes if finding a
     single MUS. This is because we only minimize with respect to the
     current set of maybe clauses, which excludes all the redundant
     clauses.

     Either may be non-minimal if the solver is interrupted, either by
     user input or by exceeding the budget.

     The arguments of the functions are (signatures in the typedefs
     below):

     CS -> (cs, model, isMinimal)
     MUS -> (mus, isMinimal)
  */
  using cs_callback =  std::function<void(std::vector<Lit> const&,
                                          model_t const&, bool)>;
  using mus_callback =  std::function<void(std::vector<Lit> const&, bool)>;
  void register_cs_callback(cs_callback f) { report_cs = f; }
  void register_mus_callback(mus_callback f) { report_mus = f; }

  // accessing the current model and forced literals at the root
  vector<Lit> getForced(int index);
  void updateForced(vector<Lit>& frc);
  void addClause(const vector<Lit>& lts);
  lbool curVal(Var x) const;
  lbool curVal(Lit x) const;

  // internal to external variable mapping. This is only for directly
  // accessing the satsolver, all other public functions of this class
  // perform the internal <-> external translation automatically.  In
  // most applications every internal literal of the SatSolver is
  // associated with an external literal on creation.  So this array
  // function is safe...i.e., won't add var_Undef to output vector. An
  // array version of ex2in is typically not safe in this way so is
  // not provided.
  void in2ex(const Minisat::vec<Lit> &from, vector<Lit> &to) const;
  void in2ex(const vector<Lit> &from, vector<Lit> &to) const;
  void in2ex(const Minisat::vec<Var> &from, vector<Var> &to) const;
  void in2ex(const model_t& from, model_t& to) const;

  Var ex2in(Var v) const;
  Lit ex2in(Lit lt) const;
  Lit in2ex(Lit lt) const;
  Var in2ex(Var v) const;

  void ensure_mapping(Var ex);
  void ensure_mapping(Lit lt);

  //stats
  //Use startTimer and elapTime to accumulate time over a set of muser calls.
  //nSolves() counts number of calls since startTimer
  //solveTime for time of most recent call.
  void startTimer() { timer = totalTime; prevSolves = solves; prevSatSolves = satSolves; }
  double elapTime() { return totalTime-timer; }
  double solveTime() { return totalTime-prevTotalTime; }
  double total_time() {
    if (stime >=0) {
        totalTime += cpuTime() - stime;
      stime = -1;
    }
    return totalTime;
  }
  int nSolvesSinceTimer() { return (int) solves-prevSolves; }
  int nSolves() { return (int) solves; }
  int nSatSolvesSinceTimer() { return (int) satSolves - prevSatSolves; }
  int nSatSolves() { return (int) satSolves; }

    virtual void printStats() const;
    virtual void printHeader(double parsetime) const;
    virtual void printFooter() const;

    std::tuple<int, int, int, int> progress() const
    {
        return std::make_tuple(
            crits_size, maybe_size, redundant_size, current_cs_size);
    }

    // parameters
    /* if true, we'll find a single mus and exit. enables some
       otherwise unsafe optimizations */
    bool single_mus{false};
    /* if true, we willll find all muses. Difference from single_mus:
       single_mus promises that the external caller will make a single
       call to mus(). !single_mus says the external caller may make
       many calls to mus(). all_muses means we'll find all muses from
       a single call to allmus() and implies !single_mus:
    */
    bool all_muses{false};
    /* if true, choose next clause to test at random */
    bool shuffle_conflict{false};
    /* if true, post exactly-one constraint instead of at-most-one in
       addMoreCrits */
    bool exactly_one{false};
    /* if true, do model rotation */
    bool mr{true};
    /* as it says: if true, do recursive model rotation */
    bool recursive_mr{true};
    /* if true do haifamuc's eager rmr */
    bool eager_rmr{true};
    /* if true, try flipping variables to get smaller correction sets
       before calling the sat solver. */
    bool mcs_rotate{true};
    /* if true, detect and falsify correction set backbone literals
       during CS minimization */
    bool mcs_backbone{false};
    /* if true (default), use addMoreCrits in mus_ to find transition clauses */
    bool amc{true};
    /* if true, use f_b^eq, not f_b */
    bool fbeq{true};
    /* if true, find a (potentially non-minimal) correction set before each set */
    bool use_cs{true};
    /* if true, use solver.conflict from the last unsat instance of
       MCS to refine the set of maybes. (if the MCS has size > 1,
       otherwise there is no UNSAT instance to use) */
    bool mcs_ref{false};
    /* which algorithm we will use to find a MUS */
    enum algorithm_t { ALG_DELETE, ALG_MISHMASH };
    algorithm_t algorithm{ALG_MISHMASH};
    /* which algorithm we use to find all MUSes. */
    enum allmus_algorithm_t { ALL_ALG_GUIDED, ALL_ALG_UNGUIDED };
    allmus_algorithm_t allmus_algorithm;
    /* heuristic to choose a clause from an MCS. Could be extended to CS and test. */
    enum clause_selection_t { FRONT, BACK, SMALLEST, MOST_ACTIVE };
    clause_selection_t clause_selection{MOST_ACTIVE};

    /* Should we demand MCSes when finding all MUSes? */
    bool allmus_find_mcses{true};

    /* Whether we use the backtracking algorithm for incomplete MUS
       enumeration during all-MUS enumeration */
    bool allmus_use_bt{true};

    /* amc bounding stuff */
    /* do we bound amc at all? */
    bool bound_amc{true};
    /* initial upper bound on number of conflicts */
    uint64_t bound_amc_cf{100};
    /* decay factor on upper bound (multiply by that each time we exceed bounds) */
    double bound_amc_decay{0.9};
    /* additive bump factor (add that each time we successfully solve the instance) */
    int bound_amc_bump{20};
    /* lower bound on the upper bound: do not get it to zero */
    uint64_t bound_amc_lb{10};
    /* upper bound on the upper bound: cannot increase beyond this */
    uint64_t bound_amc_ub{10000};

    /* force unsat checks after some time with no checks */
    bool force_unsat{true};
    /* timeout for force unsat checks */
    double force_unsat_timeout{100.0};
    /* Remove learnts, reset activities after some time with no clauses moved to redundants */
    bool stall_reset{false};
    /* timeout for resets */
    double stall_reset_timeout{500.0};
    /* if true and the number of groups added by MR over the last few
       iterations falls far below the average, reset the solver's
       heuristics (phase, learned clauses, activities) */
    bool mr_trigger_reset{true};
    double mr_trigger_ratio{0.1};
    unsigned mr_trigger_window{100};

    /* control whether we print acticity information (quite verbose) */
    bool print_activity_line{false};

    // Stats
    double timer{0}, totalTime{0}, prevTotalTime{0}, stime{-1};
    int solves{0}, prevSolves{0}, satSolves{0}, prevSatSolves{0};
    uint64_t prevConflicts{0};
    int test_clauses_add{0}, model_rotate_clauses{0};
    int amo_clauses{0}, amo_clauses_removed{0}, refine_clauses{0};
    int mcsadded{0}, mcsremoved{0}, mcs_refine_removed{0};
    int csadded{0}, csunsatremoved{0};
    double prevTime{0.0};
    mutable double last_activity_time{0.0};

    int mcs_rotate_removals{0};
    int mcs_backbone_literals{0};

    uint64_t num_muses{0}, num_mcses{0}, num_singleton_mcses{0};

    uint64_t easymcses{0}, trivialmcses{0};

    uint64_t mr_trigger_resets{0}, mr_trigger_blocks{0};

    bool isMinableCritical(vector<Lit> const& unknown, vector<Lit> const& crits, Lit c); //JB: used for mining critical clauses
    void addMinableBlockDown(vector<int> &block); //JB: used for mining critical clauses
    int minedCriticals{0}; //JB: used for mining critical clauses
    bool conflictMining{false}; //Control variable (switch) that enables the mining (in default false, i.s. the mining is disabled)

 protected:

  vector<vector<int>> blockDowns; //JB: added to store map^+ (information from Explorer)

  //Stats reporting
  bool status_timeout{false};

  enum solve_type { TEST_SOLVE, CS_SOLVE, MCS_SOLVE, AMC_SOLVE, RS_SOLVE,
                    HS_SOLVE };
  vector<std::string> solve_type_strings{ "test", "cs", "mcs", "amc", "rs", "hs" };
  vector<std::string> solve_type_desc{ "Literal test", "Correction Set",
          "Minimal Correction Set", "Singleton Correction Set",
          "Relaxation Search", "Hitting Set" };
  struct solve_stats {
      vector<uint64_t> conflicts_sat, conflicts_unsat, conflicts_undef;
      vector<double> times_sat, times_unsat, times_undef;
  };
  struct iteration_stats {
      int64_t satsolves{0}, clauses_determined{0};
      uint64_t conflicts{0}, propagations{0};
      double time{0};
  };
  vector<iteration_stats> delete_iteration_stats, mishmash_iteration_stats;
  vector<double> mr_times, mcs_rotate_times;
  solve_stats solvestats_total;
  map<solve_type, solve_stats> solvestats;
  void recordSolve(lbool result, solve_type st);
  void recordSolve(Solver& rssolver, lbool result, solve_type st);
  void reportSolveStats(const char *name, solve_stats const& st) const;

  // keep track of model rotation success: number of clauses added
  // for each iteration, average over all iterations and average
  // over the latest window size
  vector<uint64_t> mr_added;
  double mr_added_average{0.0}, mr_added_moving_average{0.0};
  // if we use the moving_avg/avg ratio to reset heuristics, we can
  // also block resets for a number of iterations, so that we do not
  // reset all the time after we hit the ratio
  int mr_block_reset{0};

  // set of assumptions under which to MUS. see solve_mus for more
  vector<Lit> assumptions;

  // current status of all clauses
  vector<uint8_t> crits, maybe, redundant, assumed;
  int crits_size, maybe_size, redundant_size, assumed_size;
  int current_cs_size{0};

  //External to internal variable mapping.
  vector<Var>in2ex_map;
  vector<Var>ex2in_map;
  int firstbvar{-1};

  // last model discovered by MCS (needed for rotation and reporting,
  // if the MCS is not singleton)
  model_t mcs_witness;

  // reporting MCSes and MUSes to the user with callbacks
  cs_callback report_cs;
  mus_callback report_mus;

  // the following are member variables (instead of local variables in
  // some member function) only to avoid reallocations
  std::vector<Lit> reportsingleton_cs;
  std::vector<Lit> external_csmus; /* CSes and MUSes have the same
                                      representation, so it can be
                                      either, depending on if we are
                                      in reportCS() or reportMUS() */
  model_t external_model;
  std::vector<bool> reportcs_seen; // a bitset
  std::vector<uint8_t> minimizecs_incs, minimizecs_ass;
  std::vector<uint8_t> eager_rmr_K;
  std::vector<Lit> eager_rmr_K_toclear;

  // When singleton MCSes are discovered during allMUSing, the
  // hssolver is not at dlvl 0, so we cannot add them as clauses. So
  // we store them here and add them when it goes back to dlvl 0
  // again. We also use them in the backtracking algorithm, replaying
  // them everytime we backtrack
  std::vector<Lit> allmus_singleton_mcses;

  // When allMUSing with the backtracking-based algorithm, adding an
  // MCS-blocking clause to hssolver may cause it to backtrack
  // arbitrarily (similar to what happens with singleton MCSes in the
  // general case). So we detect such clauses and instead of adding
  // them to hssolver, we store them here, along with the level (lvl)
  // to which they would backtrack hssolver. Until we get back to that
  // level, we add their effect (cl[0]) manually and when we do get
  // back to it, we add the clause and remove it from this list.
  struct deferred_clause {
      int lvl;
      std::vector<Lit> cl;
  };
  std::vector< deferred_clause > allmus_deferred_clauses;

  // report a CS to the user. Note that CSes have to be extended to
  // include any violated redundant clauses before being passed on to
  // report_cs. Ugh: mixing vector and vec
  void reportCS(std::vector<Lit> const& cs,
                model_t const& model, bool minimal);
  // save the (internal) caller the trouble of constructing
  // singleton CS (in MR, aMC etc)
  void reportSingletonCS(Lit b, model_t const& model, bool minimal);
  // report a MUS to the user.
  void reportMUS(std::vector<Lit> const& mus, bool minimal);

  //computation routines
  bool inSolver(Lit lt) const;
  bool inSolver(Var v) const;

  // external interface, meaning conflict contains external variables
  bool mus_(vector<Lit>& conflict, bool require_global_mcses, int64_t propBudget);
  // internal interface: conflict, knowncrits, all, contain internal
  // variables. knowncrits is a starting point of critical clauses. If
  // require_global_mcses is true, it will find CSes (and minimize
  // them) with respect to all the bvars, not just those in conflict.
  bool mus_internal(vector<Lit>& conflict, vector<Lit> const& knowncrits,
                    vector<Lit> all, bool require_global_mcses,
                    int64_t propBudget);
  // the standard deletion-based algorithm
  bool mus_delete(vector<Lit>& conflict, vector<Lit> const& all,
                  bool require_global_mcses, int64_t propBudget);
  // the mcsmus algorithm (which can be read as ``mishmash'' with some
  // imagination)
  bool mus_mishmash(vector<Lit>& conflict, vector<Lit> const& all,
                    bool require_global_mcses, int64_t propBudget);
  // perform one iteration of the "delete" algorithm (test a literal,
  // do core refinement, model rotation etc). 'conflict' is the
  // current set of undecided clauses (not known to be either
  // redundant or critical) and cs is a cache. Returns false iff
  // interrupted.
  bool mus_delete_iteration(vector<Lit>& conflict, vector<Lit>& cs);
  // perform one iteration of mishmash. Arguments 'conflict' and 'cs'
  // as well as return value are as above. 'all' is the set of all
  // optional clauses. If 'require_global_mcses' is true, this
  // iteration will find an MCS of the entire formula, not of hard U
  // crits U unkown.
  bool mus_mishmash_iteration(vector<Lit>& conflict, vector<Lit>& cs,
      vector<Lit> const& all, bool require_global_mcses);

  // prove unsatisfiability of the current formula and refine set of
  // unknowns using the core returned by the solver. rrsolver is the
  // solver to use (its current set of assumptions must match the set
  // of critical clauses) and conflict is the set of clauses of
  // unknown status. Returns false if interrupted
  bool refute_and_refine(Solver& rrsolver, vector<Lit>& conflict);

  bool allmus_unguided(vector<Lit>& candidates);
  bool allmus_guided(vector<Lit>& candidates);

  /* Construct a MUS from a given maximal hitting set and previously
     discovered MCSes (and new MCSes too) */
  bool allmus_find_mus(vector<Lit> const& all_candidates,
                       vector<Lit>& cand_crits, vector<Lit>& cand_maybes);


  // find some muses that are a subset of candidates. See
  // find_some_mus() documentation for what ``some'' means. It
  // constructs (somemus)/uses (somemus_) hssolver to store MCSes it
  // discovers and reuses them for finding critical clauses. The
  // candidates are *internal* variables. somemus_() assumes hssolver
  // is already setup and that the appropriate cs_callback has been
  // installed. somemus() does all that and then calls somemus_()
  bool somemus(vector<Lit>& candidates);
  bool somemus_(vector<Lit> const& all_candidates, vector<Lit>& candidates,
                Solver const& fork_from);

  // block a MUS in hssolver. The MUS is given as a set of
  // *negative* literals (because it comes from solver.assumptions,
  // as in reportMUS)
  void blockMUS(vector<Lit> const& MUS);

  // for allmus methods: make sure whatever is implied by
  // restricting ourselves to cand_maybes is mirrored in
  // crits/maybes. This requires that hssolver is in a fresh
  // decision level, i.e., not at level 0.
  //
  // in terms of implementation, if hssolver has a model, move
  // everything assigned to false to redundant and make that literal
  // false in hssolver. Then unit propagate. Anything that is true
  // is critical.
  void syncWithHS(vector<Lit> const& all_candidates,
                  vector<Lit>& cand_maybes);

  void ensureSoftCls(vector<Lit> const& conflict);
  void ensureSoftCls(Lit lt);
  void ensureXVars(Lit blit); // map the xvars only

  uint64_t getProps() const { return solver->getStats().propagations; }

  vector<Lit> addAmoUnk(vector<Lit>& unknowns);
  vector<Lit> addExactlyOne(vector<Lit>& unknowns);
  bool addMoreCrits(vector<Lit>& conflict);

  int modelRotate(Solver& solver, vector<Lit>& conflict, Lit violated);
  template<typename RotatePred, typename MarkCritFunc>
  int modelRotate_(Solver& solver, vector<Lit>& conflict,
                   Lit violated, model_t& model,
                   RotatePred can_rotate, MarkCritFunc mark_crit);

  /* If we flip the literals ls in the given model, do we still satisfy all
     clauses in cls? */
  template<class ...Lits>
  bool canFlip(vector<CRef> const& cls, model_t const& model,
               Lits...  ls);
  /* If we flip the literals ls in the given model, do we still satisfy all
     clauses in which p occurs? */
  template<class ...Lits>
  bool canFlip(Lit p, model_t const& model,
               Lits...  ls);

  /* Helper for canFlip: if we flip these literals, is literal 'inclause' true? */
  bool flipSatisfies(model_t const& model, Lit inclause,
                     Lit l);
  template<class ...Lits>
  bool flipSatisfies(model_t const& model, Lit inclause,
                     Lit l, Lits... ls);

  // Use rcsolver.conflict to remove clauses from the set of
  // maybes. Attribute any reduction to stat. Call f for each new
  // redundant clause. Returns the number of clauses removed. It only
  // removes clauses that satisfy the predicate p. The second version
  // uses true/nop
  template<typename Pred, typename Func>
  int refineConflict(Solver& rcsolver, vector<Lit> &conflict, int &stat, Pred p,
                     Func f);
  int refineConflict(Solver& rcsolver, vector<Lit> &conflict, int &stat);

  /* minimizeCS: minimizes an existing CS

     conflict [in] : the set over which we minimize.
     cs [in/out]: a not necessarily minimal cs on input, minimized on out
  */
  bool minimizeCS(Solver& mcssolver,
                  vector<Lit> const& conflict, vector<Lit>& cs,
                  bool block_mcs);
  bool minimizeCS_rs(Solver& mcssolver,
                     vector<Lit> const& conflict, vector<Lit>& cs,
                     bool block_mcs);
  bool minimizeCS_cl(Solver& mcssolver,
                     vector<Lit> const& conflict, vector<Lit>& cs,
                     bool block_mcs);

  // helper for minimizeCS_cl: for each literal contained in one of
  // the clauses in 'cs', try to flip it. If all clauses in
  // 'assumptions' are still satisfied after the flip, we can remove
  // all the clauses in 'cs' that this flip satisfies.
  void mcsRotate(
      Solver& mcssolver, vector<Lit> const& assumptions, vector<Lit>& cs);

  // helper for minimizeCS_cl: find backbone literals of a
  // correction set, i.e., literals $l$ that appear in all clauses
  // of the CS. Since F is SAT, F+CS is UNSAT and setting $l$ to
  // true would satisfy all clauses in CS, it follows that all
  // models of F falsify l. So add it to the assumptions during this
  // minimization phase.
  void mcsFindBackboneLits(Solver& mcssolver, vector<Lit> const& cs);

  void preProcessConflict(vector<Lit>& conflict);

  void moveToCrits(Lit crit); // crit must be a literal to be pushed into assumptions. SIDE EFFECT: enqueues crit
  void markRedundant(Lit red); // red must be a literal that appears in conflict. SIDE EFFECT: red no longer a dvar

  /* Wrapper for SAT solving */
  lbool satsolve(Solver& solver, solve_type st);

  // reporting
  std::ostream& report_activity(std::string const& activity) const;
  void report_activityln(std::string const& activity) const;
  void printStatusLine();

    // debugging: verify that assumptions \cup conflict is unsat
    void verify_unsat(vector<Lit> const& assumptions, vector<Lit> const& conflict);
    void verify_unsat(Solver& origsolver, vector<Lit> const& assumptions, vector<Lit> const& conflict);
    // verify that mcs is an MCS of assumptions \cup
    // conflict. Implemented as testing that assumptions \cup
    // (conflicts \setminus mcs) \cup (\lor_i c_i \in mcs) is unsat
    void verify_mcs(vector<Lit> const& assumptions, vector<Lit> const& conflict,
                    vector<Lit> const& mcs);
    // verify that assumptions \cup conflict is unsat, (assumptions
    // \cup conflict) \setminus {crit} is sat. crit is the positive
    // literal of the bvar.
    void verify_critical(vector<Lit> const& assumptions, vector<Lit> const& conflict,
                         Lit crit);
}; //Class Muser

/***************************************************
 * Inline functions
 **************************************************/
inline int Wcnf::bvarIdx(Var v) const
{
    return v - firstBVar;
}

inline int Wcnf::bvarIdx(Lit bl) const
{
    return bvarIdx(var(bl));
}

inline auto MUSSolver::curVal(Var x) const -> lbool {
    Var in = ex2in(x);
    if(in == var_Undef)
        return l_Undef;
    else
        return solver->value(in);
}

inline auto MUSSolver::curVal(Lit x) const -> lbool {
    Lit in = ex2in(x);
    if(in == lit_Undef)
        return l_Undef;
    else
        return solver->value(in);
}

inline auto MUSSolver::in2ex(Lit lt) const -> Lit
{
  if(var(lt) >= (int) in2ex_map.size() ||
     in2ex_map[var(lt)] == var_Undef)
    return lit_Undef;
  else
    return mkLit(in2ex_map[var(lt)], sign(lt));
}

inline auto MUSSolver::in2ex(Var v) const -> Var
{
  if(v >= (int) in2ex_map.size())
    return var_Undef;
  else
    return in2ex_map[v];
}

inline void MUSSolver::in2ex(const Minisat::vec<Lit> &from, vector<Lit> &to) const
{
  to.clear();
  for(int i = 0; i < from.size(); i++)
    to.push_back(in2ex(from[i]));
}

inline void MUSSolver::in2ex(const vector<Lit> &from, vector<Lit> &to) const
{
    to.clear();
    for(Lit l : from)
        to.push_back(in2ex(l));
}

inline void MUSSolver::in2ex(const Minisat::vec<Var> &from, vector<Var> &to) const
{
  to.clear();
  for(int i = 0; i < from.size(); i++)
      to.push_back(in2ex(from[i]));
}

inline void MUSSolver::in2ex(model_t const &from, model_t& to) const
{
    to.model.clear();
    to.model.resize(theWcnf.nVars(), l_Undef);
    for(unsigned i = 0; i != from.size(); ++i)
        if( in2ex(i) >= 0 )
            to[in2ex(i)] = from[i];
}

inline auto MUSSolver::ex2in(Var v) const -> Var
{
    if(v >= (int) ex2in_map.size())
      return var_Undef;
    else
      return ex2in_map[v];
}

inline auto MUSSolver::ex2in(Lit lt) const -> Lit
{
  if(var(lt) >= (int) ex2in_map.size() ||
     ex2in_map[var(lt)] == var_Undef)
    return lit_Undef;
  else
      return mkLit(ex2in_map[var(lt)], sign(lt));
}

inline bool MUSSolver::inSolver(Lit lt) const
{
  return inSolver(var(lt));
}

inline bool MUSSolver::inSolver(Var v) const
{
  return ex2in(v) != var_Undef;
}

} //namesspace
#endif
