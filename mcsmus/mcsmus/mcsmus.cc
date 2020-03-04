/***********[mcsmus.cc]
Copyright (c) 2012-2013 Jessica Davies, Fahiem Bacchus
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
#include <stdexcept> //added by Jaroslav Bendik
#include "mcsmus/algorithm.hh"
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <cassert>

#include "mcsmus/mcsmus.hh"
#include "mcsmus/options.hh"
#include "mcsmus/system.hh"

#include "mcsmus/debug.hh"

using std::cout;
using std::unordered_set;
using std::unordered_map;

/***********************************************************************/

#define ACT(x) if( !print_activity_line ); else report_activity(x)

/***********************************************************************/

namespace mcsmus {

static const char* _cat = "MUS";
static BoolOption opt_shuffle_conflict(_cat, "shuffle-conflict", "Shuffle conflict at each test", false);
static BoolOption opt_exactly_one(_cat, "exactly-one", "Post exactly-one instead of at-most-one in AMC", false);
static BoolOption opt_mr(_cat, "mr", "Do model rotation", true);
static BoolOption opt_recursive_mr(_cat, "recursivemr", "Do recursive model rotation", true);
static BoolOption opt_eager_rmr(_cat, "eagerrmr", "Do eager model rotation (like HaifaMUC)", true);
static BoolOption opt_mcs_rotate(_cat, "mcs-rotate", "Do model rotation to find smaller CSes", true);
static BoolOption opt_mcs_backbone(_cat, "mcs-backbone", "Falsify CS backbone literals during minimization", false);
static BoolOption opt_amc(_cat, "amc", "Call addMoreCrits", false);
static BoolOption opt_fbeq(_cat, "fbeq", "Use f_b^eq, not f_b", true);
static BoolOption opt_cs(_cat, "cs", "Use correction sets as heuristics", true);
static BoolOption opt_mcs_ref(_cat, "mcs-ref", "Use MCS UNSAT instances for clause refinement", false);
static StringOption opt_alg(_cat, "alg", "Set MUS algorithm to use [delete/mishmash]", "mishmash");
static StringOption opt_allmus_alg(_cat, "all-alg", "Set ALLMUS algorithm to use [unguided/guided]", "guided");
static IntOption opt_cs_heuristic(_cat, "cs-heuristic", "Clause selection heuristic from MCS [0=first,1=last,2=smallest, 3=most active]", 3, IntRange(0,3));
static BoolOption opt_act_line(_cat, "act-line", "Print activity information (verbose)", false);

static BoolOption opt_allmus_find_mcses(_cat, "allmus-find-mcses", "Require MCSes only when finding all MUSes", true);
static BoolOption opt_allmus_use_bt(_cat, "allmus-use-bt", "Use backtracking algorithm for incomplete MUS enumeration", true);

static BoolOption opt_bound_amc(_cat, "bound-amc", "Place an upper bound on conflicts in each addMoreCrits solve", true);
static IntOption opt_bound_amc_initcf(_cat, "bound-amc-init-cf", "Initial upper bound for bound-amc", 100, IntRange(0, INT32_MAX));
static DoubleOption opt_bound_amc_decay(_cat, "bound-amc-decay", "Decay factor for bound-amc", .9, DoubleRange(0.0, true, 1.0, true));
static IntOption opt_bound_amc_bump(_cat, "bound-amc-bump", "Bump factor for bound-amc", 1, IntRange(0, INT32_MAX));
static IntOption opt_bound_amc_ub(_cat, "bound-amc-ub", "Upper bound for bound-amc", 10000, IntRange(0, INT32_MAX));

static BoolOption opt_force_unsat(_cat, "force-unsat", "Force unsat checks every so often", true);
static DoubleOption opt_force_unsat_timeout(_cat, "force-unsat-timeout", "... every this often", 100.0);
static BoolOption opt_stall_reset(_cat, "stall-reset", "Reset learnts and activities when stall detected", false);
static DoubleOption opt_stall_reset_timeout(_cat, "stall-check-timeout", "Stall detected after this much time with no removals", 500.0);

static BoolOption opt_mr_trigger_reset(_cat, "mr-trigger-reset",
    "Reset learnts and activities when MR stalls", true);
static DoubleOption opt_mr_trigger_ratio(_cat, "mr-trigger-ratio",
    "Reset when recent MR average falls below ratio*overall average", 0.1);
static IntOption opt_mr_trigger_window(_cat, "mr-trigger-window",
    "Size of window (# of iterations) to compute recent MR average", 100);

/***********************************************************************/

struct hash_lit {
    size_t operator()(Lit p) const { return toInt(p); }
};

/***********************************************************************/

std::vector<Lit> Wcnf::getGroupLiterals(Lit bv)
{
    // we can only be called after relax()
    assert(lastBVar > 0);
    gather_seen.resize(2 * (lastBVar + 1));
    vector<Lit> rv;
    for (auto cr : getGroupClauses(bv))
        for (Lit l : getClause(cr))
            if (!gather_seen[toInt(l)]) {
                rv.push_back(l);
                gather_seen[toInt(l)] = true;
            }
    for (Lit l : rv)
        gather_seen[toInt(l)] = false;
    return rv;
}

std::vector<Var> Wcnf::getGroupVariables(Lit bv)
{
    // we can only be called after relax()
    assert(lastBVar > 0);
    gather_seen.resize(2 * (lastBVar + 1));
    vector<Var> rv;
    for (auto cr : getGroupClauses(bv))
        for (Lit l : getClause(cr)) {
            Lit posl = mkLit(var(l));
            if (!gather_seen[toInt(posl)]) {
                rv.push_back(var(l));
                gather_seen[toInt(posl)] = true;
            }
        }
    for (Var v : rv)
        gather_seen[toInt(mkLit(v))] = false;
    return rv;
}

void Wcnf::addClause(std::vector<Lit> const& lits, int grpid)
{
    bool optional = (grpid > 0);
    maxvar = std::accumulate(begin(lits), end(lits), maxvar,
        [](int m, Lit l) { return std::max(m, var(l)); });
    for_each(begin(lits), end(lits), [&](Lit l) {
            distinct_vars.insert(var(l));
        });
    // we store soft clauses locally anyway, and additionally hard
    // clauses for model rotation. Note we could be cleverer here and
    // try to keep up with the solver's state.
    Minisat::vec<Lit> mlits;
    for(Lit l : lits)
        mlits.push(l);

    CRef cr = ca.alloc(mlits, false);
    if (groups_.size() <= unsigned(grpid))
        groups_.resize(grpid + 1);
    groups_[grpid].push_back(cr);
    ++num_clauses;
    if (!optional) {
        solver->addClause(lits);
    }

    // update occurrences
    if (optional) {
        if (firstBVar < 0 )
            return;
        Var bvar = firstBVar + grpid;
        lastBVar = std::max(bvar, lastBVar);
        for (Lit l : lits) {
            if (bvar_occurrences.size() <= unsigned(toInt(l)))
                bvar_occurrences.resize(toInt(l) + 1);
            bvar_occurrences[toInt(l)].push_back(mkLit(bvar));
        }
        return;
    }
    vector<vector<CRef>>& occ = hard_occurrences;
    for (Lit l : lits) {
        if (occ.size() <= unsigned(toInt(l)))
            occ.resize(toInt(l) + 1);
        occ[toInt(l)].push_back(cr);
    }
}

void Wcnf::setSolver(MUSSolver* S) { solver = S; }

auto Wcnf::varFromBvarIdx(int bvidx) const -> Var
{
    return bvidx + firstBVar;
}

/***********************************************************************/

MUSSolver::MUSSolver(Wcnf& f, BaseSolver *psolver) :
    theWcnf(f),
    solver(psolver),
    control(getGlobalControl()),
    // options
    shuffle_conflict(opt_shuffle_conflict),
    exactly_one(opt_exactly_one),
    mr(opt_mr),
    recursive_mr(opt_recursive_mr),
    eager_rmr(opt_eager_rmr),
    mcs_rotate(opt_mcs_rotate),
    mcs_backbone(opt_mcs_backbone),
    amc(opt_amc),
    fbeq(opt_fbeq),
    use_cs(opt_cs),
    mcs_ref(opt_mcs_ref),
    clause_selection( clause_selection_t(int(opt_cs_heuristic)) ), // awkwaaaard
    allmus_find_mcses( opt_allmus_find_mcses ),
    allmus_use_bt( opt_allmus_use_bt ),
    bound_amc(opt_bound_amc),
    bound_amc_cf(opt_bound_amc_initcf),
    bound_amc_decay(opt_bound_amc_decay),
    bound_amc_bump(opt_bound_amc_bump),
    bound_amc_ub(opt_bound_amc_ub),
    force_unsat(opt_force_unsat),
    force_unsat_timeout(opt_force_unsat_timeout),
    stall_reset(opt_stall_reset),
    stall_reset_timeout(opt_stall_reset_timeout),
    mr_trigger_reset(opt_mr_trigger_reset),
    mr_trigger_ratio(opt_mr_trigger_ratio),
    mr_trigger_window(opt_mr_trigger_window),
    print_activity_line(opt_act_line)
{
    theWcnf.setSolver(this);
    std::string alg { (const char*)opt_alg };
    if( alg == "delete" )
        algorithm = ALG_DELETE;
    else if( alg == "mishmash" )
        algorithm = ALG_MISHMASH;
    else {
        cout << "Unknown algorithm \"" << opt_alg << "\"\n";
        exit(1);
    }

    std::string allalg{ (const char*)opt_allmus_alg };
    if( allalg == "guided" )
        allmus_algorithm = ALL_ALG_GUIDED;
    else if( allalg == "unguided" )
        allmus_algorithm = ALL_ALG_UNGUIDED;
}

void MUSSolver::initializeControl()
{
    control->fn_printStats = [this]() { printStats(); };
    control->fn_printHeader = [this](double pt) { printHeader(pt); };
    control->fn_printFooter = [this]() { printFooter(); };
    control->fn_printStatusLine = [this]() { printStatusLine(); };
}

void MUSSolver::ensure_mapping(const Var ex)
{
  if (ex >= (int) ex2in_map.size())
    ex2in_map.resize(ex+1,var_Undef);

  if(ex2in_map[ex] == var_Undef) {
    Var in {solver->newVar()};
    ex2in_map[ex] = in;
    if (in >= (int) in2ex_map.size())
      in2ex_map.resize(in+1, var_Undef);
    in2ex_map[in] = ex;
  }
}

inline void MUSSolver::ensure_mapping(const Lit l)
{
  ensure_mapping(var(l));
}

auto MUSSolver::getForced(int index) -> vector<Lit>
{
  static vector<Lit> forced {};
  updateForced(forced);

  vector<Lit> tmp;
  for(size_t i = index; i < forced.size(); i++)
    tmp.push_back(in2ex(forced[i]));
  return tmp;
}

void MUSSolver::updateForced(vector<Lit>& frc) {
  vector<Lit> trail = solver->entailed();
  auto limit = trail.size();
  unsigned i {0};

  if(frc.size() > 0) {
    i = frc.size() - 1;
    while(i < limit && trail[i++] != frc.back());
  }

  for( ; i < limit; i++)
    if(in2ex(trail[i]) != lit_Undef)
      frc.push_back(trail[i]);
}

void MUSSolver::addClause(const vector<Lit>& lts)
{
  //if a b-var is referred to make sure we add the equivalent softcls.
  vector<Lit> ps;
  for(auto lt: lts) {
    if(theWcnf.isBvar(lt))
      ensureSoftCls(lt);
    ensure_mapping(lt);
    ps.push_back(ex2in(lt));
  }
  solver->addClause(ps);
}

void MUSSolver::ensureSoftCls(vector<Lit> const& conflict) {
    //Ensure that all soft clauses referenced by conflict are in the solver.
    //This means that we incrementally add only the soft clauses of the theory
    //that are actually involved in minimizations.

    // we do x vars first to ensure that bvars get mapped above
    for(auto blit : conflict) {
        if (inSolver(blit))
            continue;
        ensureXVars(blit);
    }
    if (firstbvar < 0)
        firstbvar = solver->nUsedVars();
    for(auto blit : conflict)
        ensureSoftCls(blit);
}

void MUSSolver::ensureXVars(Lit blit) {
    for(auto cref: theWcnf.getGroupClauses(blit))
        for(Lit l : theWcnf.getClause(cref) )
            if (!theWcnf.isBvar(l))
                ensure_mapping(l);
}

void MUSSolver::ensureSoftCls(Lit blit)
{
  //Ensure that all soft clauses referenced by conflict are in the solver.
  //This means that we incrementally add only the soft clauses of the theory
  //that are actually involved in minimizations.
  if(inSolver(blit))
    return;

  if(!theWcnf.isBvar(blit))
    cout << " ERROR conflict does not contain only bvars\n";

  Lit b = mkLit(var(blit)); // note that blit might be ~b.
  // however: negative literals in the conflict do not make sense if
  // we do not use fbeq
  assert(!sign(blit) || fbeq);

  ensure_mapping(b);
  Lit inbPos = ex2in(b);
  auto& group = theWcnf.getGroupClauses(b);
  vector<Lit> ps;
  for (auto cr : group) {
      ps.clear();
      ps.push_back(inbPos);
      for (Lit lt : theWcnf.getClause(cr)) {
          ensure_mapping(lt);

          Lit inLt = ex2in(lt);
          ps.push_back(inLt);

          // disable fbeq for groups of size greater than 1
          if (fbeq && group.size() == 1)
              solver->addClause({~inbPos, ~inLt}); // full equivalence
      }
      solver->addClause(ps);
  }

  solver->setFrozen(var(inbPos), true);
}

void MUSSolver::preProcessConflict(vector<Lit>& conflict) {
  size_t j {0};
  for(auto lt : conflict) {
    auto val = solver->value(ex2in(lt));
    if(val == l_Undef)
      conflict[j++] = ex2in(lt);
    else if(val == l_True) {
      conflict = vector<Lit> { ex2in(lt) };
      return;
    }
  }
  conflict.resize(j);
}

void MUSSolver::moveToCrits(Lit crit)
{
    DOUT << "moving " << crit << " to crits\n";

    DBG{
        vector<Lit> conflict;
        for(int i = 0; i != solver->nUsedVars(); ++i)
            if( unsigned(toInt(mkLit(i))) < maybe.size() &&
                maybe[toInt(mkLit(i))] )
                conflict.push_back(mkLit(i));
        verify_critical(solver->assumptions(), conflict, ~crit);
    }

    assert(maybe[toInt(~crit)]);
    assert(!crits[toInt(crit)]);
    assert(!redundant[toInt(~crit)]);
    assert(maybe_size > 0);

    solver->assumptions().push_back(crit);
    crits[toInt(crit)] = 1;
    ++crits_size;
    maybe[toInt(~crit)] = 0;
    --maybe_size;
    if( single_mus ) {
        solver->addClause({crit});
    }
}

void MUSSolver::markRedundant(Lit red)
{
    DOUT << "moving " << red << " to redundants\n";

    DBG{
        vector<Lit> conflict;
        for(int i = 0; i != solver->nUsedVars(); ++i)
            if( unsigned(toInt(mkLit(i))) < maybe.size() &&
                maybe[toInt(mkLit(i))] && red != mkLit(i))
                conflict.push_back(mkLit(i));
        verify_unsat(solver->assumptions(), conflict);
    }

    assert(maybe[toInt(red)]);
    assert(!crits[toInt(~red)]);
    assert(!redundant[toInt(red)]);
    assert(maybe_size > 0);

    if( single_mus ) {
        solver->setDecisionVar(var(red), false);
        solver->setRedundantBvar(var(red));
    }
    maybe[toInt(red)] = 0;
    redundant[toInt(red)] = 1;
    --maybe_size;
    ++redundant_size;
}

auto MUSSolver::satsolve(Solver& thesolver, solve_type st) -> lbool
{
    prevConflicts = thesolver.getStats().conflicts;
    lbool val = thesolver.solve();
    recordSolve(thesolver, val, st);
    if( val != l_True )
        thesolver.cancelUntil(0);
    else
        mcs_witness = thesolver.model();
    return val;
}

bool MUSSolver::mus_(vector<Lit>& conflict,
               bool require_global_mcses,
               int64_t propBudget)
{
    assert(solver->assumptions().size() == 0);
    vector<Lit> dummy;

    ensureSoftCls(conflict);
    preProcessConflict(conflict);

    vector<Lit> all;
    for(Var exv : theWcnf.bvars())
        all.push_back( mkLit(ex2in(exv)) );
    return mus_internal(conflict, dummy, all, require_global_mcses, propBudget);
}

bool MUSSolver::mus_internal(vector<Lit>& conflict,
                       vector<Lit> const& knowncrits,
                       vector<Lit> all,
                       bool require_global_mcses,
                       int64_t propBudget)
{
    solver->assumptions() = assumptions;
    solver->setBvars(firstbvar);

    crits.clear();
    maybe.clear();
    redundant.clear();
    assumed.clear();
    crits.resize(2*solver->nUsedVars());
    maybe.resize(2*solver->nUsedVars());
    redundant.resize(2*solver->nUsedVars());
    assumed.resize(solver->nUsedVars());
    fill(begin(crits), end(crits), 0);
    fill(begin(maybe), end(maybe), 0);
    fill(begin(redundant), end(redundant), 0);
    fill(begin(assumed), end(assumed), 0);

    for (auto l : assumptions)
        assumed[var(l)] = 1;

    erase_if(all, [&](Lit l) { return assumed[var(l)]; });
    erase_if(conflict, [&](Lit l) { return assumed[var(l)]; });

    for( auto l: conflict )
        maybe[toInt(l)] = 1;
    for( auto l: knowncrits )
        maybe[toInt(l)] = 1;

    for (auto l : all) {
        if (!maybe[toInt(l)] && !assumed[var(l)]) {
            redundant[toInt(l)] = 1;
            ++redundant_size;
        }
    }

    crits_size = 0;
    maybe_size = conflict.size() + knowncrits.size();
    redundant_size = all.size() - conflict.size() - knowncrits.size();

    if(conflict.size() <= 1) {
        if(conflict.size() == 0)
            DOUT << "c mcsmus found false conflict\n";
        if(conflict.size() == 1) {
            DOUT << "c mcsmus found unit conflict on input\n";
            conflict[0] = in2ex(conflict[0]);
        }
        return true;
    }

    for(auto l: knowncrits) {
        assert(maybe[toInt(l)] == 1);
        moveToCrits(~l);
    }

    switch(algorithm) {
    case ALG_DELETE:
        return mus_delete(conflict, all, require_global_mcses, propBudget);
    case ALG_MISHMASH:
        return mus_mishmash(conflict, all, require_global_mcses, propBudget);
    }
    return false; // unreachable, but gcc cannot tell
}

bool MUSSolver::mus_delete_iteration(vector<Lit>& conflict, vector<Lit>& cs)
{
    cs.clear();
    Lit test = lit_Undef;
    if( use_cs ) {
        cs.clear();
        auto val = satsolve(*solver, CS_SOLVE);
        if( val == l_Undef ) {
            for(size_t i = 0; i < conflict.size(); i++)
                solver->assumptions().push_back(~conflict[i]); //act like all remaining conflict lits are critical
            conflict.clear();
            maybe_size = 0;
            crits_size = solver->assumptions().size(); // we do not have to update maybe and crits, we are exiting
            reportMUS(solver->assumptions(), false);
            return false;
        } else if( val == l_False ) {
            // we are done
            ACT("cs") << " val = " << val << std::endl;
            conflict.clear();
            csunsatremoved = maybe_size;
            maybe_size = 0;
            reportMUS(solver->assumptions(), true);
            return true;
        } else { // l_True
            auto model = solver->model();
            for(Lit l : conflict) {
                if(model[var(l)] == l_True)
                    cs.push_back(l);
            }
            assert( !cs.empty() );
            test = cs.back();
            ACT("cs") << " val = " << val << " CS size = " << cs.size()
                           << std::endl;
            auto it = find(begin(conflict), end(conflict), test);
            conflict.erase(it);
            if( cs.size() == 1 )
                ++csadded;
        }
    } else {
        test = conflict.back();
        conflict.pop_back();
    }
    int ncrits = solver->assumptions().size();
    for(Lit l : conflict)
      solver->assumptions().push_back(~l);

    lbool val;
    if( cs.size() == 1 )
        val = l_True;
    else {
        val = satsolve(*solver, TEST_SOLVE);
    }

    ACT("test") << " test = " << test << " val = " << val << std::endl;

    solver->shrinkAssumptions(solver->assumptions().size() - ncrits);

    if(val == l_Undef) { //timed out
      solver->assumptions().push_back(~test);
      for(size_t i = 0; i < conflict.size(); i++)
        solver->assumptions().push_back(~conflict[i]); //act like all remaining conflict lits are critical
      conflict.clear();
      maybe_size = 0;
      crits_size = solver->assumptions().size(); // we do not have to update maybe and crits, we are exiting
      return false;
    }
    else if(val == l_False) { //redundant
        markRedundant(test);
        refineConflict(*solver, conflict, refine_clauses);
        if( cs.size() == 2 ) {
            // the other clause just became critical, so we can add it
            Lit crit = ~cs.front();
            moveToCrits(crit);
            // XXX: is the model we found in the findCS() call still
            // valid after the UNSAT later solve? Let's say yes.
            reportSingletonCS(~crit, solver->model(), false); //XXX: unclear
            auto it = find(begin(conflict), end(conflict), ~crit);
            conflict.erase(it);
            ++csadded;
            ACT("cs") << " from size 2 CS" << std::endl;
        }
    }
    else { // l_True:
      ++test_clauses_add;
      moveToCrits(~test);
      reportSingletonCS(test, solver->model(), false); // XXX: not minimal unless no redundants
      if( mr )
          modelRotate(*solver, conflict, test);
      if( amc )
          return addMoreCrits(conflict);
    }
    return true;
}

bool MUSSolver::mus_delete(vector<Lit>& conflict, vector<Lit> const& /*all*/,
                           bool require_global_mcses,
                           int64_t propBudget)
{
  if( require_global_mcses ) {
      throw std::runtime_error("amc cannot find global MCSes");
  }
  //reduce conflict to a MUS, return true if we succeeded.
  //don't use more than probBudget props to do so.
  bool isMus {true};

  solver->setPropagationBudget( (propBudget < 0) ? -1
                                : solver->getStats().propagations + propBudget );

  if( shuffle_conflict) {
      std::random_device rd;
      std::mt19937 rng(rd());
      shuffle(begin(conflict), end(conflict), rng);
  }

  vector<Lit> cs; // a (potentially non-minimal) correction set

  while(conflict.size() > 0) {
      double start_time = cpuTime();
      auto start_confl = solver->getStats().conflicts;
      auto start_prop = solver->getStats().propagations;
      auto start_solves = satSolves;
      auto start_maybe = maybe_size;
      auto dummy = on_scope_exit([&]() {
          delete_iteration_stats.push_back(
              {satSolves - start_solves, start_maybe - maybe_size,
                  solver->getStats().conflicts - start_confl,
                  solver->getStats().propagations - start_prop,
                  cpuTime() - start_time});
      });

      isMus = mus_delete_iteration(conflict, cs);

      assert(int(conflict.size()) == maybe_size);
      assert(int(solver->assumptions().size()) == crits_size);
  }

  assert(conflict.size() == 0);
  for(Lit l : solver->assumptions() ) {
      conflict.push_back(in2ex(~l));
      assert(theWcnf.isBvar(in2ex(l)));
  }

  return isMus;
}

bool MUSSolver::refute_and_refine(Solver& rrsolver, vector<Lit>& conflict)
{
    for (Lit l : conflict)
        rrsolver.assumptions().push_back(~l);
    auto val = satsolve(rrsolver, TEST_SOLVE);
    
    if(val == l_True){ // Added by Jaroslav Bendik, when I use mcsmus inside my MUS enumeration tool, var(lt) is sometimes <0 at this point which causes a segfault 
    	    throw std::runtime_error("Unexpected Failure in mcsmus");
    }
    assert(val != l_True);
    rrsolver.shrinkAssumptions(
        rrsolver.assumptions().size() - crits_size - assumptions.size());
    if (val == l_Undef)
        return false;
    int removed = refineConflict(*solver, conflict, refine_clauses);
    ACT("test") << " removed " << removed << "\n";
    return true;
}

bool MUSSolver::mus_mishmash_iteration(vector<Lit>& conflict, vector<Lit>& cs,
    vector<Lit> const& all, bool require_global_mcses)
{
    if (mr_trigger_reset) {
        if (mr_added_moving_average < mr_added_average * mr_trigger_ratio
            && mr_block_reset == 0) {
            ACT("reset") << "@" << mr_added.size() << " iterations, "
                         << "average " << mr_added_average << ", last "
                         << mr_trigger_window << " average "
                         << mr_added_moving_average << std::endl;
            ++mr_trigger_resets;
            solver->resetHeuristics();
            mr_block_reset = mr_trigger_window;
        }
        if (mr_added_moving_average < mr_added_average * mr_trigger_ratio
            && (mr_block_reset > 0)) {
            ++mr_trigger_blocks;
            ACT("block reset") << " will keep blocking for " << mr_block_reset
                               << " more iterations\n";
        }
    }

    if (mr_block_reset > 0)
        --mr_block_reset;

    cs.clear();
    auto val = findMCS(conflict, all, require_global_mcses,
        require_global_mcses || single_mus, cs);
    if (val == l_Undef) {
        for (size_t i = 0; i < conflict.size(); i++)
            solver->assumptions().push_back(~conflict[i]); // act like all
                                                           // remaining conflict
                                                           // lits are critical
        conflict.clear();
        maybe_size = 0;
        crits_size = solver->assumptions().size(); // we do not have to update
                                                   // maybe and crits, we are
                                                   // exiting
        reportMUS(solver->assumptions(), false);
        return false;
    } else if (val == l_False) {
        // we are done
        ACT("cs") << " val = " << val << std::endl;
        conflict.clear();
        csunsatremoved = maybe_size;
        maybe_size = 0;
        return true;
    }

    DBG {
        verify_mcs(solver->assumptions(), conflict, cs);
    }

    // l_True
    bool added = false;
    Lit addedlit = lit_Undef;

    switch (clause_selection) {
    case FRONT:
        break;
    case BACK:
        reverse(begin(cs), end(cs));
        break;
    case SMALLEST:
        nth_element(begin(cs), begin(cs), end(cs), [&](Lit x, Lit y) {
            auto& gx = theWcnf.getGroupClauses(in2ex(x));
            auto& gy = theWcnf.getGroupClauses(in2ex(y));
            if (gx.size() == 1 && gy.size() == 1) {
                auto& cx = theWcnf.getClause(gx[0]);
                auto& cy = theWcnf.getClause(gy[0]);
                return cx.size() < cy.size();
            } else
                return gx.size() < gy.size();
        });
        break;
    case MOST_ACTIVE:
        nth_element(begin(cs), begin(cs), end(cs), [&](Lit x, Lit y) {
            return solver->varActivity(var(x)) > solver->varActivity(var(y));
        });
    }

    for (Lit l : cs) {
        assert(theWcnf.isBvar(in2ex(l)));
        assert(require_global_mcses || maybe[toInt(l)]);
        assert(maybe[toInt(l)] || redundant[toInt(l)]);
        if (!maybe[toInt(l)])
            continue;
        if (!added) {
            moveToCrits(~l);
            // no call to reportCS here, findMCS() did it already
            addedlit = ~l;
            ++mcsadded;
            added = true;
        } else {
            ++mcsremoved;
            markRedundant(l);
        }
    }
    assert(added);
    erase_if(conflict, [&](Lit q) { return !maybe[toInt(q)]; });
    ACT("MCS") << " MCS size " << cs.size() << std::endl;

    if (mr) {
        auto added = modelRotate(*solver, conflict, ~addedlit);
        auto& avg{mr_added_average};
        auto& mvavg{mr_added_moving_average};
        auto wnd{mr_trigger_window};
        mr_added.push_back(added);
        avg = (avg * (mr_added.size() - 1) + added) / mr_added.size();
        if (mr_added.size() > wnd)
            mvavg = (mvavg * wnd - mr_added[mr_added.size() - wnd - 1] + added)
                / wnd;
        else
            mvavg = avg;
    }

    // reprove unsatisfiability to refine core
    if (cs.size() > 1) {
        bool success = refute_and_refine(*solver, conflict);
        if (!success) { //JB: this probably means that refute_and_refine was interrupted (timeouted), i.e. the reportMUS below is an overapproximated MUS
            for (Lit l : conflict)
                solver->assumptions().push_back(~l);
            conflict.clear();
            maybe_size = 0;
            crits_size = solver->assumptions().size();
            reportMUS(solver->assumptions(), false);
            return false;
        }

        DBG {
            verify_unsat(solver->assumptions(), conflict);
        }
    }

    if (amc)
        return addMoreCrits(conflict);
    return true;
}


// Added by JB
// Add a clause from Explorer (from the overal MUS enumeration tool) that represents a mcs
void MUSSolver::addMinableBlockDown(vector<int> &block){
	blockDowns.push_back(block);
}

// Added by JB
// checks if the clause unknown[c] is minable critical for F = (unknown \cup crits)
bool MUSSolver::isMinableCritical(vector<Lit> const& unknown, vector<Lit> const& crits, Lit c){
	vector<bool> f(theWcnf.nClauses(), false);
	for(auto l: unknown){
		f[theWcnf.bvarIdx(in2ex(l)) - 1] = true;
	}
	for(auto l: crits){
		f[theWcnf.bvarIdx(in2ex(l)) - 1] = true;
	}
	f[theWcnf.bvarIdx(in2ex(c)) - 1] = false;

	for(auto &block: blockDowns){
		bool found = false;
		for(auto cl: block){
			if(f[cl]){
				found = true;
				break;
			}			
		}
		if(!found) return true;// the block (mcs) is not hit by F - {unknown[c]}
	}
	return false;
}

bool MUSSolver::mus_mishmash(vector<Lit>& conflict, vector<Lit> const& all,
    bool require_global_mcses, int64_t propBudget)
{
    bool isMus {true};

    solver->setPropagationBudget( (propBudget < 0) ? -1
                                  : solver->getStats().propagations + propBudget );

    vector<Lit> cs; // a (potentially non-minimal) correction set

    double beginTime = cpuTime();
    double lastRemoved{beginTime}, lastTest{beginTime};

    while(conflict.size() > 0) {
        DBG {
            verify_unsat(solver->assumptions(), conflict);
        }

        double start_time = cpuTime();
        auto start_confl = solver->getStats().conflicts;
        auto start_prop = solver->getStats().propagations;
        auto start_solves = satSolves;
        auto start_maybe = maybe_size;
        auto dummy = on_scope_exit([&]() {
            mishmash_iteration_stats.push_back(
                {satSolves - start_solves, start_maybe - maybe_size,
                    solver->getStats().conflicts - start_confl,
                    solver->getStats().propagations - start_prop,
                    cpuTime() - start_time});
        });

        auto prev_redundant{redundant_size};
        isMus = mus_mishmash_iteration(conflict, cs, all, require_global_mcses);

	//JB: Critical Clauses Extension based on information from Explorer (the overall MUST tool)
	if(conflictMining){
		vector<Lit> newCrits;
		for(auto c: conflict){
			if(isMinableCritical(conflict, solver->assumptions(), c)){			
				minedCriticals++;
				newCrits.push_back(c);
			}else{
			}
		}
		for(auto &c: newCrits){
			moveToCrits(~c);
			conflict.erase(find(begin(conflict), end(conflict), c));
		}
		if(!conflict.empty()){
			isMus = true;
		}
	}
	//End of added code by JB

        auto num_removed{redundant_size - prev_redundant};
        double current_time = cpuTime();
        if (num_removed > 0)
            lastRemoved = current_time;
        else {
            if (force_unsat && current_time - lastTest > force_unsat_timeout) {
                bool success = refute_and_refine(*solver, conflict);
                if (!success) {
                    for (Lit l : conflict)
                        solver->assumptions().push_back(~l);
                    conflict.clear();
                    maybe_size = 0;
                    crits_size = solver->assumptions().size();
                    reportMUS(solver->assumptions(), false);
                    isMus = false;
                }
                // see if refute_and_refine did anything
                num_removed = redundant_size - prev_redundant;
                if (stall_reset && num_removed == 0
                    && current_time - lastRemoved > stall_reset_timeout) {
                    std::cout << "Stalled, resetting learnts and activities\n";
                    solver->resetHeuristics();
                }
                DBG {
                    verify_unsat(solver->assumptions(), conflict);
                }
            }
        }

        assert(int(conflict.size()) == maybe_size);
        assert(solver->assumptions().size() == crits_size + assumptions.size());
    }

    DBG {
        verify_unsat(solver->assumptions(), conflict);
    }

    reportMUS(solver->assumptions(), true);


    assert(conflict.size() == 0);
    for (Lit l : solver->assumptions()) {
        if (!assumed[var(l)]) {
            conflict.push_back(in2ex(~l));
            assert(theWcnf.isBvar(in2ex(l)));
        }
    }

    return isMus;
}

template<typename Pred, typename Func>
int MUSSolver::refineConflict(Solver& rcsolver,
                        vector<Lit>& conflict, int& stat, Pred p,
                        Func f)
{
    int numremoved{0};
    auto& rccore = rcsolver.core();
    for(Lit l : conflict)
        if(!rccore.has(l) && p(l)) {
            markRedundant(l);
            f(l);
            ++numremoved;
        }
    erase_if( conflict, [&](Lit l) {
            return !rccore.has(l) && p(l); });
    ACT("refine") << " removed " << numremoved << " clauses\n";
    stat += numremoved;
    return numremoved;
}

int MUSSolver::refineConflict(Solver& rcsolver,
                                      vector<Lit>& conflict, int& stat)
{
    return refineConflict(rcsolver,conflict, stat,
                          [](Lit) { return true; }, [](Lit) {});
}

bool MUSSolver::allmus(vector<Lit>& candidates)
{
    all_muses = true;

    hssolver = solver->clone_blank();

    ensureSoftCls(candidates);
    preProcessConflict(candidates);

    crits.clear();
    maybe.clear();
    redundant.clear();
    assumed.clear();
    crits.resize(2*solver->nUsedVars());
    maybe.resize(2*solver->nUsedVars());
    redundant.resize(2*solver->nUsedVars());
    assumed.resize(solver->nUsedVars());
    fill(begin(crits), end(crits), 0);
    fill(begin(maybe), end(maybe), 0);
    fill(begin(redundant), end(redundant), 0);
    fill(begin(assumed), end(assumed), 0);

    for (auto l : assumptions)
        assumed[var(l)] = 1;

    for(auto l: candidates)
        maybe[toInt(l)] = 1;

    crits_size = solver->assumptions().size();
    maybe_size = candidates.size();
    redundant_size = 0;

    solver->assumptions().clear();
    solver->setBvars(firstbvar);

    for(Lit l : candidates) {
        while( hssolver->nUsedVars() <= var(l) )
            hssolver->newVar();
        hssolver->setVarPolarity(var(l), l_True);
        hssolver->setFrozen(var(l), true);
    }

    switch( allmus_algorithm ) {
    case ALL_ALG_UNGUIDED:
        return allmus_unguided(candidates);
    case ALL_ALG_GUIDED:
        return allmus_guided(candidates);
    }
    return false; // unreachable
}

bool MUSSolver::allmus_unguided(vector<Lit>& candidates)
{
    lbool val = l_True;
    vector<Lit> mcs, mus;
    vector<Lit> ps;

    // set polarity to always prefer false
    for(Lit l : candidates)
        hssolver->setVarPolarity(var(l), l_True);

    do {
        mcs.clear();
        solver->assumptions().clear();
        for(Lit l : mus)
            solver->assumptions().push_back(~l);
        lbool mcsval;
        if( allmus_find_mcses ) {
            // we call with require_global_mcs = false, because
            // candidates itself is also global
            mcsval = findMCS(candidates, candidates, false, true, mcs);
        } else
            mcsval = findCS(*solver, candidates, mcs);
        if( mcsval == l_Undef ) { // interrupted
            cout << "Interrupted after "
                 << num_muses << " MUSes, "
                 << num_mcses << " MCSes, "
                 << num_singleton_mcses << " singleton MCSes\n";
            return false;
        }
        if( mcsval == l_True ) {
            ps.clear();
            for(Lit l : mcs)
                ps.push_back(l);
            hssolver->addClause(ps);
            if( mcs.size() == 1 ) {
                ++num_singleton_mcses;
                solver->enqueue(~mcs[0]);
            }
            ++num_mcses;
        } else {
            // found a MUS, report it
            ++num_muses;
            reportMUS(solver->assumptions(), true);
            // block it
            ps.clear();
            for(Lit l : mus)
                ps.push_back(~l);
            hssolver->addClause(ps);
        }
        val = satsolve(*hssolver, HS_SOLVE);
        if( val == l_True ) {
            mus.clear();
            for(Lit l : candidates)
                if( hssolver->model().value(var(l)) == l_True )
                    mus.push_back( l );
            ACT("hs") << " val = " << val << " HS size = " << mus.size() << "\n";
        }
    } while( val == l_True );
    cout << num_muses << " MUSes found\n";
    return true;
}

namespace {
    // compute the level to which we have to backtrack to in order to
    // add ps to solver. As a side-effect, ps is sorted so that ps[0]
    // is the effect of the clause at btlvl *if* ps[1] is false at
    // that level.
    int computeClauseBtlvl(BaseSolver& solver, vector<Lit>& ps)
    {
        // same computation as MiniSolver::addClause_internal
        auto pivot = std::partition(begin(ps), end(ps),
            [&](Lit p) { return solver.value(p) != l_False; });
        auto ppos = std::distance(begin(ps), pivot);
        if (ppos < 2) {
            std::partial_sort(
                pivot, begin(ps) + 2, end(ps), [&](Lit l1, Lit l2) {
                    return solver.levelOf(var(l1)) > solver.levelOf(var(l2));
                });
            int btlvl = solver.levelOf(var(ps[1]));
            if (ppos == 0)
                btlvl = std::max(
                    0, std::min(btlvl, solver.levelOf(var(ps[0])) - 1));
            return btlvl;
        }
        return solver.decisionLevel();
    }
}

bool MUSSolver::allmus_guided(vector<Lit>& candidates)
{
    lbool val = l_True;
    vector<Lit> ps;

    // set polarity to always prefer true
    for(Lit l : candidates)
        hssolver->setVarPolarity(var(l), l_False);

    // singleton MCSes cannot be added as clauses as we discover them,
    // because hssolver is not at dlvl 0. We store them here and add
    // them when it reaches level 0 again.
    vector<Lit>& singletons = allmus_singleton_mcses;

    auto prevcb = report_cs;
    auto cb = [&](std::vector<Lit> const& mcs,
                  model_t const& model,
                  bool minimal) {
        if( prevcb )
            prevcb(mcs, model, minimal);

        if( allmus_use_bt && !minimal && mcs.size() > 100 ) {
            // do not record very large non-minimal CSes, although we
            // have no choice for allmus_find_mus (it is a limitation
            // of that implementation that we have to have all CSes as
            // clauses)
            return;
        }

        ps.clear();
        for(Lit l : mcs)
            ps.push_back(ex2in(l));
        if( ps.size() == 1 ) {
            // it may be false already if it is a singleton MCS and we
            // have discovered all MUSes and UP all discovered MCSes
            // and MUSes makes it false.
            assert(hssolver->value(ps[0]) != l_False
                || hssolver->levelOf(var(ps[0])) == 0);
            hssolver->enqueue(ps[0]);
            singletons.push_back(ps[0]);
        } else {
            for(Lit l : ps) {
                assert(!sign(l));
                assert(hssolver->value(l) != l_True);
            }
            // either add the clause to hssolver or defer it for when
            // we backtrack to the appropriate level.
            int btlvl{computeClauseBtlvl(*hssolver, ps)};
            if (btlvl < hssolver->decisionLevel()) {
                // we cannot be adding a falsified clause. That would
                // imply marking as redundant all clauses of an MCS.
                assert(hssolver->value(ps[0]) != l_False && hssolver->value(ps[1]) == l_False);
                allmus_deferred_clauses.push_back({btlvl, ps});
                hssolver->enqueue(ps[0]);
            } else {
                hssolver->addClause(ps);
            }
        }
        if( mcs.size() == 1 ) {
            ++num_singleton_mcses;
            solver->addClause({ex2in(~mcs[0])});
        }
        ++num_mcses;
    };
    register_cs_callback(cb);

    vector<Lit> cand_crits, cand_maybes; // candidates
    vector<Lit> dummy_mcs; // we never use the result of findMCS()
    do {
        // all singleton MCSes were unit clauses that we could not add
        // while we were at dlvl 1 -- add them now
        assert(hssolver->decisionLevel() == 0);
        for(Lit l : singletons)
            hssolver->addClause({l});
        singletons.clear();

        control->checkPrintStatusLine();

        val = satsolve(*hssolver, HS_SOLVE);
        // a weakness in Minisolver, because solving does not
        // explicitly return to dlvl 0. It is normally transparent,
        // but we mix in newDecisionLevel/cancelUntil and value()
        hssolver->cancelUntil(0);
        if (val == l_Undef) {
            cout << "HS interrupted\n";
            return false;
        }
        if( val != l_True ) {
            cout << "HS is unsat\n";
            break;
        }
        cand_maybes.clear();
        for(Lit l : candidates)
            if( hssolver->model().value(var(l)) == l_True )
                cand_maybes.push_back( l );
        ACT("hs") << " val = " << val << " HS size = " << cand_maybes.size() << "\n";

        crits.clear();
        maybe.clear();
        redundant.clear();
        crits.resize(2*solver->nUsedVars());
        maybe.resize(2*solver->nUsedVars());
        redundant.resize(2*solver->nUsedVars());
        fill(begin(crits), end(crits), 0);
        fill(begin(maybe), end(maybe), 0);
        fill(begin(redundant), end(redundant), 0);

        for( auto l : candidates )
            maybe[toInt(l)] = 1;
        crits_size = 0;
        maybe_size = candidates.size();
        redundant_size = 0;

        /* If the problem is SAT, find an MCS and try again */
        assert( solver->assumptions().empty() );
        for(Lit l : cand_maybes)
            solver->assumptions().push_back(~l);
        dummy_mcs.clear();
        val = findCS(*solver, candidates, dummy_mcs);
        solver->assumptions().clear();
        if( val == l_Undef )
            return false;
        else if( val == l_True ) {
            // because we find a maximal solution and it is sat, it is
            // an MCS
            reportCS(dummy_mcs, solver->model(), true);
            continue;
        }

        ACT("cs") << "\n";

        /* val == l_False, construct a MUS */
        bool minimized{false};
        if( allmus_use_bt ) {
            assert(cand_crits.size() == 0);
            minimized = somemus_(candidates, cand_maybes, *solver);
            cand_maybes.clear();
        } else {
            minimized = allmus_find_mus(candidates, cand_crits, cand_maybes);
            assert( solver->assumptions().size() > 0 );
        }
        if( !minimized )
            return false;
        assert( allmus_deferred_clauses.empty() );
        assert( cand_maybes.empty() );
        /* no call to reportMUS here, it is done by allmus_find_mus/somemus_ */

        // block it
        if( !allmus_use_bt ) {
            ++num_muses;
            blockMUS(solver->assumptions());
            solver->assumptions().clear();
        }
        // XXX: bumping redundant variables was previously done here
        // but is unsupported under the new interface
    } while( true ); // we never actually exit this way
    cout << num_muses << " MUSes found\n";
    return true;
}

bool MUSSolver::allmus_find_mus(vector<Lit> const& all_candidates,
                          vector<Lit>& cand_crits,
                          vector<Lit>& cand_maybes)
{
    /* dlvl 0 - singleton MCSes, singleton redundants (could happen
       after some time)

       dlvl 1 - candidates \setminus cand_maybes are set to
       false. Propagation may set some variables to true, so they are
       local singleton MCSes. After each iteration of the do-while
       loop, (a) at least one bvar clause is added to cand_crits and
       made true here and (b) some bvars are removed from cand_maybes
       and made false here. We never get a conflict. At the end, all
       the variables that are true form a MUS (but we use the
       moveToCrits machinery, so the MUS is also stored in
       solver.assumption, which is what we use)
    */
    vector<Lit> localmcs;

    cand_crits.clear();
    hssolver->newDecisionLevel();
    syncWithHS(all_candidates, cand_maybes);
    for(Lit l : solver->assumptions())
        cand_crits.push_back( ~l );

    auto forked_solver = solver->fork();
    // ensure all bvars exist in the fork. Although we freeze bvars,
    // if a clause contains pure literals, we may eliminate all
    // clauses that contain some bvar, so fork_solver will not copy it
    for(Lit l : all_candidates)
        while( forked_solver->nUsedVars() <= var(l) )
            forked_solver->newVar();
    forked_solver->setBvars(firstbvar);
    ACT("fork") << " #assumptions " << solver->assumptions().size()
                << std::endl;

    vector<Lit> cs;
    auto prevcs_size = 0u;
    int iteration = 0;
    while( !cand_maybes.empty() ) {
        ++iteration;

        for(Lit l : cand_maybes )
            assert( hssolver->value(l) == l_Undef );
        if( iteration == 1 || prevcs_size > 1 ) {
            // the simplest way to do clause refinement is to treat it
            // as another source of non-critical clauses. Note that we
            // do this on the first iteration, to exploit the fact
            // that we just proved the exact same thing to be unsat in
            // allmus_guided(). Therefore the first test() of each MUS
            // is guaranteed to be almost immediate.
            mcs_witness.clear(); // disables model rotation
            forked_solver->assumptions().clear();
            forked_solver->assumptions() = solver->assumptions();
            for(Lit l : cand_maybes)
                forked_solver->assumptions().push_back(~l);
            if( iteration > 1 ) {
                std::sort(begin(forked_solver->assumptions()), end(forked_solver->assumptions()),
                          [&](Lit l1, Lit l2) {
                              return forked_solver->varActivity(var(l1)) > forked_solver->varActivity(var(l2));
                          });
            } else {
                std::sort(begin(forked_solver->assumptions()), end(forked_solver->assumptions()),
                          [&](Lit l1, Lit l2) {
                              return solver->varActivity(var(l1)) > solver->varActivity(var(l2));
                          });
            }
            auto val = satsolve(*forked_solver, TEST_SOLVE);
            assert(val != l_True);
            forked_solver->assumptions().clear();
            assert( val != l_True );
            if( val == l_Undef )
                return false;
            cs.clear();
            auto &core = forked_solver->core();
            for(Lit l : cand_maybes)
                if(!core.has(l))
                    cs.push_back(l);
            ACT("test") << " val = " << val << " may remove " << cs.size()
                        << ", " << forked_solver->getStats().conflicts << " conflicts"
                        << ", " << forked_solver->getStats().decisions << " decisions"
                        << std::endl;
            prevcs_size = 0;
        } else {
            forked_solver->assumptions().clear();
            for(Lit l : cand_crits)
                forked_solver->assumptions().push_back(~l);
            mcs_witness.clear();
            lbool val = findMCS(*forked_solver, cand_maybes, cand_maybes, false, true, cs);
            if( val == l_Undef ) {
                return false;
            } else if( val == l_False ) {
                cs = cand_maybes;
            } else {
                ACT("MCS") << " MCS size " << cs.size() << std::endl;
                prevcs_size = cs.size();
            }
        }

        for(Lit l : cs) {
            if( hssolver->value(l) != l_Undef )
                continue;
            hssolver->enqueue(~l);
            auto res = hssolver->propagate();
            assert( res == prop_result::ok );
        }

        for(Lit l : cs) {
            assert(hssolver->value(l) != l_Undef);
            if( hssolver->value(l) == l_True ) {
                moveToCrits(~l);
                cand_crits.push_back(l);
            } else {
                if( !redundant[toInt(l)] )
                    markRedundant(l);
            }
        }

        ACT("up") << std::endl;

        // model rotation only applies if the last model we found
        // violated exactly one of the criticals.
        if( mr && forked_solver->model().size() > 0 ) {
            int numviol = 0;
            Lit lviol = lit_Undef;
            for(Lit l : cs ) {
                if( hssolver->value(l) == l_True ) {
                    ++numviol;
                    lviol = l;
                }
            }
            if( numviol == 1 ) {
                // make sure crits/maybes are in sync with HS solver
                // before calling modelRotate()
                for(Lit l : cand_maybes ) {
                    if( hssolver->value(l) == l_True && !crits[toInt(~l)]) {
                        moveToCrits(~l);
                        cand_crits.push_back(l);
                    }
                }

                erase_if( cand_maybes, [&](Lit l) { return hssolver->value(var(l)) != l_Undef; });
                for(Lit l : cand_maybes )
                    assert( maybe[toInt(l)] );

                modelRotate(*forked_solver, cand_maybes, lviol);
            }
        }

        // gather the results of unit propagation to mark clauses as
        // critical (true in the HS solver) ...
        for(Lit l : cand_maybes ) {
            if( hssolver->value(l) == l_True && !crits[toInt(~l)]) {
                moveToCrits(~l);
                cand_crits.push_back(l);
            }
        }

        // ... or redundant (false, not caused by unit propagation, only
        // explicitly enqueued by us)
        erase_if( cand_maybes, [&](Lit l) { return hssolver->value(var(l)) != l_Undef; });
        for(Lit l : cand_maybes )
            assert( maybe[toInt(l)] );

        // now fix all crits and redundants in forked_solver
        for(Lit l : all_candidates) {
            if( crits[toInt(~l)] && forked_solver->value(l) != l_False ) {
                // this variable may actually already be true, which
                // means we have reached the end now. In the next
                // iteration we will do a "test" solve and it will be unsat
                forked_solver->enqueue(~l);
            }
            if( !maybe[toInt(l)] )
                forked_solver->setRedundantBvar(var(l));
        }

    }

    hssolver->cancelUntil(0);

    reportMUS(solver->assumptions(), true);

    return true;
}

bool MUSSolver::somemus(vector<Lit>& candidates)
{
    assert(solver->assumptions().size()==0);

    all_muses = true;

    hssolver = solver->clone_blank();

    ensureSoftCls(candidates);
    preProcessConflict(candidates);

    crits.resize(2*solver->nUsedVars());
    maybe.resize(2*solver->nUsedVars());
    redundant.resize(2*solver->nUsedVars());
    fill(begin(crits), end(crits), 0);
    fill(begin(maybe), end(maybe), 0);
    fill(begin(redundant), end(redundant), 0);

    for(auto l: candidates)
        maybe[toInt(l)] = 1;

    crits_size = 0;
    maybe_size = candidates.size();
    redundant_size = 0;

    solver->assumptions().clear();
    solver->setBvars(firstbvar);

    for(Lit l : candidates) {
        while( hssolver->nUsedVars() <= var(l) )
            hssolver->newVar();
        hssolver->setVarPolarity(var(l), l_True);
        hssolver->setFrozen(var(l), true);
    }

    for(Lit l : candidates)
        solver->assumptions().push_back(~l);
    auto val = satsolve(*solver, TEST_SOLVE);
    solver->assumptions().clear();
    if( val == l_True )
        throw "SAT instance\n";
    else if( val == l_Undef )
        return false;

    vector<Lit> all(candidates);
    return somemus_(all, candidates, *solver);
}

bool MUSSolver::somemus_(vector<Lit> const& all_candidates,
                   vector<Lit>& candidates,
                   Solver const& fork_from)
{
    // As we discover each non-singleton MCS, we store it here. On
    // backtrack, we pop the last element of the choice point, mark it
    // redundant *at the previous decision level* and then choose the
    // last element as the critical. In this way each choice point is
    // always an MCS of the working formula.
    using vsize = vector<Lit>::size_type;
    vector<vector<Lit>> cp;
    vector<vector<Lit>> cands;
    vector<Lit> crit_trail;
    vector<vsize> crit_trail_lim;
    vector<Lit> red_trail;
    vector<vsize> red_trail_lim;
    vector<int> assumptions_trail_lim;

    vector<vector<Lit>> muses; // keep everything we've generated here
                               // so we can block it afterwards

    bool rv{true};

    assert(hssolver->decisionLevel() == 0);
    hssolver->newDecisionLevel();
    syncWithHS(all_candidates, candidates);
    assert(unsigned(crits_size) == solver->assumptions().size());

    std::unique_ptr<Solver> forked_solver{ fork_from.clone() };
    forked_solver->assumptions() = solver->assumptions();
    for(Lit l : solver->assumptions())
        forked_solver->addClause({l});
    ACT("fork") << " #assumptions " << solver->assumptions().size()
                << std::endl;

    vector<Lit> cs;

    if( forked_solver->core().size() > 0 ) {
        // we have already proved UNSAT, so use it to refine the
        // conflict.
        int nr = refineConflict(*forked_solver, candidates, refine_clauses,
                                [](Lit) { return true; },
                                [&](Lit l) {
                                    forked_solver->setFrozen(var(l), false);
                                    forked_solver->setDecisionVar(var(l), false);
                                    assert(hssolver->value(l) != l_True);
                                    if( hssolver->value(l) != l_Undef )
                                        hssolver->enqueue(~l);
                                });
        ACT("refine") << " removed " << nr << "\n";

        hssolver->propagate();
        auto hstrail = hssolver->forced();
        auto added = 0u;
        for(Lit l : hstrail) {
            if( sign(l) )
                continue;
            if( crits[toInt(~l)] )
                continue;
            ++added;
            moveToCrits(~l);
            forked_solver->assumptions().push_back(~l);
            if( forked_solver->value(~l) != l_False )
                forked_solver->enqueue(~l);
            else
                forked_solver->addClause({~l});
        }
        ACT("up") << " after initial refine, added "
                  << added << std::endl;
    }


    lbool val{l_Undef};
    do {
        assert(forked_solver->assumptions().size() == solver->assumptions().size());
        assert(std::equal(begin(solver->assumptions()), end(solver->assumptions()),
                          begin(forked_solver->assumptions())));
        val = findMCS(*forked_solver,
                      candidates, candidates, false, true, cs);
        if( val == l_Undef ) {
            // interrupted, do not emit any more MUSes
            rv = false;
            break;
        } else if( val == l_False ) {
            // no more MCSes. emit a MUS
            if (control->verbosity >= 2) {
                cout << "reporting MUS at dlvl " << cp.size() << "\n";
                cout << "\tsize " << solver->assumptions().size() << "\n";
            }
            ++num_muses;
            reportMUS(solver->assumptions(), true);
            muses.push_back(solver->assumptions());

            if (control->verbosity >= 2) {
                cout << forked_solver->getStats().conflicts << " conflicts "
                     << forked_solver->getStats().decisions << " decisions "
                     << forked_solver->getStats().propagations
                     << " propagations\n";
            }

            // ... backtrack
            if( cp.empty() ) {
                // we are done here
                break;
            }

            if (control->verbosity >= 2)
                cout << "Backtracking to level " << (cp.size() - 1) << ", "
                     << (cp.back().size() - 1) << " choices remaining"
                     << " from " << cp.back() << ", "
                     << allmus_deferred_clauses.size() << " deferred clauses"
                     << "\n";

            hssolver->cancelUntil(cp.size()); // hssolver is always one level deeper than cp.size()

            while( crit_trail.size() > crit_trail_lim.back() ) {
                Lit l = crit_trail.back();
                crit_trail.pop_back();
                maybe[toInt(~l)] = 1;
                crits[toInt(l)] = 0;
                ++maybe_size;
                --crits_size;
            }
            crit_trail_lim.pop_back();

            while(red_trail.size() > red_trail_lim.back() ) {
                Lit l = red_trail.back();
                red_trail.pop_back();
                maybe[toInt(l)] = 1;
                redundant[toInt(l)] = 0;
                ++maybe_size;
                --redundant_size;
            }
            red_trail_lim.pop_back();

            forked_solver->shrinkAssumptions(solver->assumptions().size() - assumptions_trail_lim.back());
            solver->shrinkAssumptions(solver->assumptions().size() - assumptions_trail_lim.back());
            assumptions_trail_lim.pop_back();
            assert((unsigned)crits_size == solver->assumptions().size());

            candidates = move(cands.back());
            cands.pop_back();

            // reapply the effects of deferred clauses or potentially
            // add them to hssolver
            for(auto& dcl : allmus_deferred_clauses) {
                Lit l = dcl.cl[0];
                if( dcl.lvl == hssolver->decisionLevel() ) {
                    hssolver->addClause(dcl.cl);
                } else {
                    assert(hssolver->value(dcl.cl[0]) != l_False
                           && hssolver->value(dcl.cl[1]) == l_False);
                    hssolver->enqueue(dcl.cl[0]);
                }
                if( !crits[toInt(~l)] ) {
                    moveToCrits(~l);
                    forked_solver->assumptions().push_back(~l);
                    forked_solver->addClause({~l});
                    crit_trail.push_back(~l);
                }
            }
            erase_if(allmus_deferred_clauses, [&](deferred_clause& cl) {
                    return cl.lvl == hssolver->decisionLevel(); });

            // ... re-fork from the original solver
            forked_solver = std::unique_ptr<Solver>{fork_from.clone()};
            forked_solver->assumptions() = solver->assumptions();
            mcs_witness.clear();
            for(Lit l : solver->assumptions() ) {
                forked_solver->addClause( {l} );
            }
            ACT("fork") << " #assumptions " << solver->assumptions().size()
                        << std::endl;

            // ... use the last choice point as our MCS
            cs = cp.back();
            cp.pop_back();

            Lit prevcrit = cs.back();
            markRedundant(prevcrit);
            forked_solver->setFrozen(var(prevcrit), false);
            forked_solver->setDecisionVar(var(prevcrit), false);
            red_trail.push_back(prevcrit);
            assert(hssolver->value(prevcrit) != l_True);
            if( hssolver->value(prevcrit) != l_Undef )
                hssolver->enqueue(~prevcrit);
            cs.pop_back();

            // replay all singleton MCSes, some of which may be new
            for(Lit l : allmus_singleton_mcses) {
                if( hssolver->value(l) != l_True ) {
                    assert(hssolver->value(l) != l_False);
                    hssolver->enqueue(l);
                }
                if( !crits[ toInt(~l) ] ) {
                    moveToCrits(~l);
                    forked_solver->assumptions().push_back(~l);
                    forked_solver->addClause({~l});
                    crit_trail.push_back(~l);
                }
            }

            if (control->verbosity >= 2) {
                cout << "after backtracking"
                     << " ncrits " << crits_size << " maybe " << maybe_size
                     << " redundant " << redundant_size << "\n";
                cout << "using cs " << cs << "\n";
            }
            val = l_True; // it is as if we found an MCS
        } else {
            //cout << "using new MCS " << cs << "\n";
        }

        // MCS found
        if( cs.size() > 1 ) {
            if (control->verbosity >= 2) {
                cout << "branching on an mcs of size " << cs.size()
                     << " at dlvl " << (cp.size() + 1) << " ncrits "
                     << crits_size << " maybe " << maybe_size << " redundant "
                     << redundant_size << "\n";
                cout << "\t" << cs << "\n";
            }
            // new decision level
            hssolver->newDecisionLevel();
            cp.push_back(cs);
            cands.push_back(candidates);
            crit_trail_lim.push_back(crit_trail.size());
            red_trail_lim.push_back(red_trail.size());
            assumptions_trail_lim.push_back(solver->assumptions().size());
            assert(unsigned(hssolver->decisionLevel()) == cp.size()+1);
        }

        assert(forked_solver->decisionLevel() == 0);

        Lit newcrit = cs.back();
        for( Lit l : cs ) {
            if( l == newcrit )
                continue;
            markRedundant(l);
            forked_solver->setFrozen(var(l), false);
            forked_solver->setDecisionVar(var(l), false);
            red_trail.push_back(l);
            assert(hssolver->value(l) != l_True);
            if( hssolver->value(l) != l_Undef )
                hssolver->enqueue(~l);
        }
        moveToCrits(~newcrit);
        forked_solver->assumptions().push_back(~newcrit);
        if(forked_solver->value(~newcrit) != l_False) {
            forked_solver->enqueue(~newcrit);
        } else {
            // if false already, it is in the assumptions, so it will
            // cause a conflict in the next solve()
            if (control->verbosity >= 2)
                cout << "cs: " << ~newcrit
                     << " is false already in forked_solver\n";
            forked_solver->addClause({~newcrit});
        }
        //forked_solver->addClause(~newcrit);
        crit_trail.push_back(~newcrit);

        erase_if(candidates, [&](Lit q) { return !maybe[toInt(q)]; });
        ACT("MCS") << " MCS size " << cs.size() << std::endl;

        // note we may come here after backtracking and making a
        // different choice in the choice point. We then have no
        // model, so we have to check that model.size() > 0
        if( mr && forked_solver->model().size() > 0 ) {
            modelRotate(*forked_solver, candidates, newcrit);
            // modelRotate calls moveToCrits which puts new critical
            // clauses into solver->assumptions(). Copy those over to
            // forked_solver->assumptions (it is an invariant that
            // these two are always the same) and record them for
            // backtracking, as we do when we call moveToCrits above.
            for(unsigned q = forked_solver->assumptions().size();
                q != solver->assumptions().size(); ++q) {
                forked_solver->assumptions().push_back(solver->assumptions()[q]);
                if(forked_solver->value(solver->assumptions()[q]) != l_False)
                    forked_solver->enqueue(solver->assumptions()[q]);
                else if (control->verbosity >= 2)
                    cout << "mr: " << solver->assumptions()[q] << " is false already in forked_solver\n";
                crit_trail.push_back(solver->assumptions()[q]);
            }
            forked_solver->assumptions().clear();
            forked_solver->assumptions() = solver->assumptions();
        }

        if( cs.size() > 1 ) {
            // clause set refinement
            for(Lit l : candidates)
                forked_solver->assumptions().push_back(~l);
            val = satsolve(*forked_solver, TEST_SOLVE);
            assert(val != l_True);
            forked_solver->shrinkAssumptions(forked_solver->assumptions().size() - crits_size);
            if( val == l_Undef ) {
                rv = false;
                break;
            }
            int nr = refineConflict(*forked_solver, candidates, refine_clauses,
                                    [](Lit) { return true; },
                                    [&](Lit l) {
                                        forked_solver->setFrozen(var(l), false);
                                        forked_solver->setDecisionVar(var(l), false);
                                        red_trail.push_back(l);
                                        assert(hssolver->value(l) != l_True);
                                        if( hssolver->value(l) != l_Undef )
                                            hssolver->enqueue(~l);
                                    });
            ACT("test") << " removed " << nr << "\n";

            // unit propagation in hssolver
            auto forced = hssolver->forced();
            auto thead = forced.size();
            hssolver->propagate();
            forced = hssolver->forced();
            if( forced.size() > thead ) {
                for(auto q = thead; q != forced.size(); ++q) {
                    Lit l = forced[q];
                    assert( !sign(l) );
                    if( crits[toInt(~l)] )
                        continue;
                    moveToCrits(~l);
                    forked_solver->assumptions().push_back(~l);
                    if( forked_solver->value(~l) != l_False )
                        forked_solver->enqueue(~l);
                    else
                        forked_solver->addClause({~l});
                }
            }
            ACT("up") << " added "
                      << (thead - forced.size()) << std::endl;
        }
    } while( !cp.empty() || val == l_True );

    solver->assumptions().clear();
    hssolver->cancelUntil(0);

    for( vector<Lit>& mus : muses )
        blockMUS(mus);

    if (control->verbosity >= 2)
        cout << "partial MUS generation finished\n";

    return rv;
}


void MUSSolver::syncWithHS(vector<Lit> const& all_candidates,
                     vector<Lit>& cand_maybes)
{
    // limit ourselves to cand_maybe
    if( hssolver->model().size() > 0 )
        for(Lit l : all_candidates)
            if( hssolver->model().value(var(l)) == l_False ) {
                markRedundant(l);
                hssolver->enqueue( ~l );
            }
    auto lvl1result = hssolver->propagate();
    assert(lvl1result == prop_result::ok);
    // Find crits that are immediately implied at lvl 1
    for(Lit l : all_candidates)
        if( hssolver->value(var(l)) == l_True ) {
            moveToCrits(~l);
        }
    erase_if( cand_maybes, [&](Lit l) { return hssolver->value(var(l)) != l_Undef; });
}

void MUSSolver::blockMUS(vector<Lit> const& mus)
{
    for( Lit l : mus )
        assert(sign(l));
    hssolver->addClause(mus);
}

bool MUSSolver::addMoreCrits(vector<Lit>& conflict)
{
  if(conflict.size() <= 1)
    return true;

  int critsFnd {0};

  uint64_t conflictsBefore{solver->getStats().conflicts};

  int64_t orig_conflict_budget{solver->conflictBudget()};

  if( bound_amc )
      solver->setConflictBudget(solver->getStats().conflicts + bound_amc_cf);

  //Insert removable most one constraint over conflicts to more criticals
  vector<Lit> clits = exactly_one ? addExactlyOne(conflict) : addAmoUnk(conflict);
  lbool critVal;

  solver->setFrozen(var(clits[0]), true);
  solver->assumptions().push_back(~clits[0]); //activate at-most-one
  while((critVal = satsolve(*solver, AMC_SOLVE)) == l_True) {
    auto model = solver->model();
    for(size_t j = 0; j < conflict.size(); j++) {
      if(model.value(conflict[j]) == l_True) {
        ++amo_clauses;
        ++critsFnd;
        moveToCrits(~conflict[j]);
        reportSingletonCS(conflict[j], model, false); // XXX: unless no redundants
        DOUT << "amc added in:" << ~conflict[j]
             << " ex: " << in2ex(~conflict[j]) << std::endl;
        ACT("amc") << " conflicts " << (solver->getStats().conflicts - conflictsBefore)
            << std::endl;
        conflictsBefore = solver->getStats().conflicts;
        assert(theWcnf.isBvar(in2ex(~conflict[j])));
        DBG {
            DOUT << "verifying that its clauses are violated\n";
            Lit exb = in2ex(~conflict[j]);
            auto& clauses = theWcnf.getGroupClauses(exb);
            for(auto cr : clauses) {
                auto& vc = theWcnf.getClause(cr);
                for(auto q : vc )
                    DOUT << q << (solver->model().value(ex2in(q))) << ' ';
                DOUT << "\n";
            }
            DOUT << "----\n";
        }
        if( mr )
            modelRotate(*solver, conflict, conflict[j]);
        break;
      }
    }

    if( bound_amc ) {
        bound_amc_cf += bound_amc_bump;
        bound_amc_cf = std::min(bound_amc_cf, bound_amc_ub);
        solver->setConflictBudget( solver->getStats().conflicts + bound_amc_cf );
    }
  }
  solver->cancelUntil(0);
  for(auto lt : clits) {
      if( !solver->isEliminated(var(lt)) )
          solver->releaseVar(lt);
  }

  solver->assumptions().clear();
  for(Var i = 0; i != solver->nUsedVars(); ++i) {
      int idx = toInt(mkLit(i, true));
      if( crits.size() > unsigned(idx) && crits[idx] ) {
          assert( theWcnf.isBvar(in2ex(mkLit(i))));
          solver->assumptions().push_back(mkLit(i,true));
      }
  }

  solver->setConflictBudget(orig_conflict_budget);

  if(critVal == l_Undef) {  //timed out
      // was it an amc-specific timeout or global?
      if( solver->withinBudget() ) {
          // amc-specific timeout
          bound_amc_cf = std::max(uint64_t(bound_amc_cf * bound_amc_decay),
                                  bound_amc_lb);
      } else {
          for( auto q : conflict )
              solver->assumptions().push_back(~q); //act like all remaining conflict lits are critical
          conflict.clear();
          maybe_size = 0;
          crits_size = solver->assumptions().size();
          return false;
      }
  }

  erase_if(conflict, [&](Lit l) { return maybe[toInt(l)] == 0; });
  if(conflict.size() > 0 && critVal != l_Undef ) {
      ACT("amc") << " val = " << critVal
                 << " conflicts " << (solver->getStats().conflicts - conflictsBefore)
                 << std::endl;
      markRedundant(conflict.back());
      conflict.pop_back();
      ++amo_clauses_removed;
  }

  ACT("amc") << " added " << critsFnd << " val = " << critVal
             << " cnfl bound = " << bound_amc_cf
             << std::endl;
  return true;
}

auto MUSSolver::addAmoUnk(vector<Lit>& unknowns) -> vector<Lit>
{
  //Add an at-most-one constraint to the solver.
  //  at most one unknown lit can be TRUE.
  //In:  unknowns.size() > 1
  //Out: clits[0] is the control lit that deletes all clauses
  //     clits[1]...are the other auxiliary variables
  //     added to encode the constraint.

    //    cout << "Adding AMO " << unknowns << "\n";

  vector<Lit> clits;

  for(size_t i = 0; i < unknowns.size(); i++) {
      Var x = solver->newVar();
      solver->setDecisionVar(x, false);
      clits.push_back(mkLit(x, false));
  }

  for(size_t i = 0; i < unknowns.size()-1 ; i++) {
      //unk_i --> clit_i+1
      //~clit_i+1 >> ~unk_i
      solver->addClause({~unknowns[i], clits[i+1], clits[0]});
  }

  for(size_t i = 1; i < clits.size()-1 ; i++) {
      //clit_i --> clit_i+1
      //~clit_i+1 --> ~clit_i
      solver->addClause( {~clits[i], clits[i+1], clits[0]} );
  }

  for(size_t i = 1; i < clits.size(); i++) {
      //clit_i --> ~unk_i
      //unk_i --> ~clit_i
      solver->addClause( {~clits[i], ~unknowns[i], clits[0]} );
  }

  return clits;
}

auto MUSSolver::addExactlyOne(vector<Lit>& unknowns) -> vector<Lit>
{
  //Add an exactly-one constraint to the solver.
  //  at most one unknown lit can be TRUE.
  //In:  unknowns.size() > 1
  //Out: clits[0] is the control lit that deletes all clauses
  //     clits[1]...are the other auxiliary variables
  //     added to encode the constraint.

  vector<Lit> clits;

  for(size_t i = 0; i <= unknowns.size(); i++) {
      Var x = solver->newVar();
      solver->setDecisionVar(x, false);
      clits.push_back(mkLit(x, false));
  }

  for(size_t i = 0; i < unknowns.size() ; i++) {
      //unk_i --> clit_i+1
      //~clit_i+1 >> ~unk_i
      solver->addClause({~unknowns[i], clits[i+1], clits[0]});
  }

  for(size_t i = 1; i < clits.size()-1 ; i++) {
      //clit_i --> clit_i+1
      //~clit_i+1 --> ~clit_i
      solver->addClause({ ~clits[i], clits[i+1], clits[0] });
  }

  for(size_t i = 1; i < clits.size(); i++) {
      //clit_i --> ~unk_i
      //unk_i --> ~clit_i
      solver->addClause({ ~clits[i], ~unknowns[i], clits[0] });
  }

  for(size_t i = 1; i != unknowns.size(); ++i) {
      // c_i+1 --> u_i+1 \/ c_i
      solver->addClause({ ~clits[i+1], unknowns[i], clits[i], clits[0] });
  }

  solver->addClause( {clits.back(), clits[0]} );

  return clits;
}

int MUSSolver::modelRotate(
    Solver& rotsolver, vector<Lit>& conflict, Lit inviolated)
{
    int numadded{0};
    auto start_time = cpuTime();
    auto mrmodel = mcs_witness;
    vector<Lit> newcrits;
    if (eager_rmr) {
        eager_rmr_K.resize(2 * solver->nUsedVars());
        numadded = modelRotate_(rotsolver, conflict, inviolated, mrmodel,
            [&](Lit b) { return eager_rmr_K[toInt(b)] == 0; },
            [&](Lit b) {
                if (!crits[toInt(~b)])
                    newcrits.push_back(~b);
                assert(eager_rmr_K[toInt(b)] == 0);
                eager_rmr_K[toInt(b)] = 1;
                eager_rmr_K_toclear.push_back(b);
            });
        for (Lit l : eager_rmr_K_toclear)
            eager_rmr_K[toInt(l)] = 0;
    } else
        numadded = modelRotate_(rotsolver, conflict, inviolated, mrmodel,
            [&](Lit b) -> bool { return maybe[toInt(b)]; },
            [&](Lit b) { newcrits.push_back(~b); });

    for (Lit l : newcrits)
        if (!assumed[var(l)])
            moveToCrits(l);

    erase_if(conflict, [&](Lit l) { return maybe[toInt(l)] == 0; });
    mr_times.push_back(cpuTime() - start_time);
    ACT("MR") << " added " << numadded << std::endl;
    return numadded;
}

bool MUSSolver::flipSatisfies(model_t const& model, Lit inclause, Lit l)
{
    return inclause == ~l
        || (inclause != l && (model.value(inclause) == l_True));
}

template <class... Lits>
bool MUSSolver::flipSatisfies(
    model_t const& model, Lit inclause, Lit l, Lits... ls)
{
    return inclause == ~l
        || (inclause != l && flipSatisfies(model, inclause, ls...));
}

template <class... Lits>
bool MUSSolver::canFlip(
    vector<CRef> const& cls, model_t const& model, Lits... ls)
{
    DOUT << "Checking " << cls.size() << " hard clauses\n";
    for (CRef cr : cls) {
        Wcnf::Clause const& c = theWcnf.getClause(cr);
        bool sat{false};
        for (Lit exl : c)
            if (flipSatisfies(model, ex2in(exl), ls...)) {
                sat = true;
                break;
            }
        if (!sat) {
            DOUT << "Violates hard clause " << c << "\n";
            return false;
        }
    }
    return true;
}

template <class... Lits>
bool MUSSolver::canFlip(Lit p, model_t const& model, Lits... ls)
{
    if (p == lit_Undef || ~p == lit_Undef)
        return true;
    return canFlip(theWcnf.hardOccurrences(in2ex(p)), model, ls...);
}


template <typename RotatePred, typename MarkCritFunc>
int MUSSolver::modelRotate_(Solver& rotsolver, vector<Lit>& conflict,
    Lit inviolated, model_t& model, RotatePred can_rotate_to_p,
    MarkCritFunc make_crit)
{
    assert(!sign(inviolated));
    assert(rotsolver.decisionLevel() == 0);
    if(var(inviolated) < 0) // Added by Jaroslav Bendik, when I use mcsmus inside my MUS enumeration tool, var(lt) is sometimes <0 at this point which causes a segfault 
	throw std::runtime_error("Unexpected Failure in mcsmus");
	    
    Lit violated = in2ex(inviolated);
    auto vc = theWcnf.getGroupVariables(violated);

    control->checkPrintStatusLine();

    DBG {
        DOUT << "Model rotate around b" << violated << ", in: " << inviolated
             << "\n";
        DOUT << "external clauses:\n";
        for (auto cr : theWcnf.getGroupClauses(violated))
            DOUT << "\t" << theWcnf.getClause(cr) << "\n";
        DOUT << "-----\n";
        DOUT << "again:\n";
        for (auto cr : theWcnf.getGroupClauses(violated)) {
            cout << "[";
            for(Lit l : theWcnf.getClause(cr))
                cout << l << model.value(ex2in(l)) << " ";
            cout << "]\n";
        }
        DOUT << "-----\n";
    }

    if (rotsolver.isEliminated(var(inviolated)))
        DOUT << "eliminated bvar\n";

    DBG
        for (auto q : vc) {
            DOUT << mkLit(q) << model.value(ex2in(q)) << ' ';
        }
    DOUT << "\n";
    DBG
        for (auto q : vc) {
            DOUT << mkLit(q) << rotsolver.value(ex2in(q)) << "("
                 << rotsolver.levelOf(ex2in(q)) << ") ";
        }
    DOUT << "\n";
    if (theWcnf.getGroupClauses(violated).size() == 1) {
        auto clref = theWcnf.getGroupClauses(violated)[0];
        auto& cl = theWcnf.getClause(clref);
        for (auto q : cl) {
            assert(model.value(ex2in(q)) == l_False);
            // we are at dlvl 0, so if it is true it means it got
            // propagated when we called moveToCrits on its clause. It may
            // have had other consequences and we cannot tell them apart
            // from things that happened earlier and we can safely ignore
            // (such as the if(value == l_False) continue; test below. So
            // bail.
            if (rotsolver.value(ex2in(q)) == l_True)
                return 0;
        }
    }

    /* Store rotations so we can recurse */
    unordered_map<Var, Lit> rotations;

    int numadded = 0;

    for (Var exflipvar : vc) {
        Var inflipvar{ex2in(exflipvar)};
        Lit inl{mkLit(inflipvar, model.value(inflipvar) == l_True)};
        Lit exl{in2ex(inl)};
        Lit infalselit = ~inl; // lit that WILL be false if flipped
        Lit exfalselit = ~exl;

        auto isSatisfied = [&](CRef cr) {
            auto& cl = theWcnf.getClause(cr);
            DBG {
                DOUT << "Checking external clause " << cl << "\n";
                // DOUT << 'b' << bvarlit << model.value(bvarlit) << ' ';
                for (Lit exol : cl) {
                    Lit ol = ex2in(exol);
                    DOUT << ol << (model[var(ol)] ^ sign(ol)) << ' ';
                }
                DOUT << "\n";
            }
            bool sat = std::any_of(begin(cl), end(cl), [&](Lit exol) {
                Lit ol = ex2in(exol);
                if (rotsolver.isEliminated(var(ol)))
                    DOUT << "ol " << ol << " is eliminated\n";
                bool lsat = (var(ol) != inflipvar && (model.value(ol) == l_True)
                    || ol == inl);
                if (lsat)
                    DOUT << "satisfied by " << ol << "\n";
                return lsat;
            });
            if (sat)
                DOUT << "satisfied\n";
            return sat;
        };

        if (rotsolver.value(inl) != l_Undef)
            continue; // cannot be flipped
        DOUT << "Checking flip of " << exl
             << " (currently: " << model.value(inl) << ")\n"
             << "\tinternal: " << inl << "\n";
        if (theWcnf.isBvar(exl)) // cannot flip bvars in MR
            continue;
        if (rotsolver.isEliminated(inflipvar))
            DOUT << "\tEliminated\n";
        DOUT << "\t" << theWcnf.hardOccurrences(infalselit).size()
             << " hard occurrences\n";
        DOUT << "soft occurrences: " << theWcnf.softBvarOccurrences(exfalselit) << "\n";
        DBG {
            vector<int> p;
            for(Lit b : theWcnf.softBvarOccurrences(exfalselit))
                p.push_back(theWcnf.bvarIdx(b));
            DOUT << "(groups: " << p << ")\n";
        }
        int numviolations = 0;
        Lit inrotateb = lit_Undef;
        bool has_ignored{false};

        // check first whether this group is satisfied by flipping this var
        auto& thiscls = theWcnf.getGroupClauses(violated);
        bool canrotatetoself{can_rotate_to_p(inviolated)};
        DOUT << PV(canrotatetoself) << "\n";
        if( !std::all_of(begin(thiscls), end(thiscls), isSatisfied) ) {
            DOUT << "rotating does not satisfly self group\n";
            // regardless, we allow to rotate to this new model, as it
            // may yield other rotations (with eager MR). Even in EMR,
            // however, it will not try to rotate on the same group
            // twice (the group is placed in K). This means that
            // alternate MR strategies may enter an infinite loop if
            // they are even more permissive that EMR. Also note that
            // this may only occur in groupMUS, not normal MUS instances.
            inrotateb = inviolated;
            ++numviolations;
            if( !canrotatetoself )
                continue;
        }

        for (Lit exbvarlit : theWcnf.softBvarOccurrences(exfalselit)) {
            auto& cls = theWcnf.getGroupClauses(exbvarlit);
            Lit bvarlit = ex2in(exbvarlit);

            DOUT << "Checking group " << theWcnf.bvarIdx(exbvarlit) << "\n";

            if (redundant[toInt(bvarlit)]) {
                has_ignored = true;
                continue;
            }

            bool canrotateto{can_rotate_to_p(bvarlit)};
            bool satisfied{std::all_of(begin(cls), end(cls), isSatisfied)};

            if (!satisfied) {
                inrotateb = bvarlit;
                ++numviolations;
                if (!canrotateto)
                    ++numviolations;
            }
            DOUT << PV(satisfied) << "," << PV(canrotateto) << "\n";
            if (numviolations > 1)
                break;
        }

        if (single_mus && has_ignored) {
            static bool p = false;
            if (!p && control->verbosity >= 2)
                cout << "missing cleanup of redundant/satisfied clauses\n";
            p = true;
            // erase_if( exoccurrences[toInt(exfalselit)], [&] (Lit l) {
            //         return redundant[toInt(ex2in(l))]; } );
        }

        if (!canFlip(infalselit, model, infalselit, inviolated, ~inrotateb)
            || !canFlip(inviolated, model, infalselit, inviolated, ~inrotateb)
            || !canFlip(~inrotateb, model, infalselit, inviolated, ~inrotateb))
            continue;

        if (numviolations == 0) {
            /* this is slightly complicated. In the single-mus case we
               cannot arrive here. In the many-mus case, e.g. in
               maxsat, we can block an MCS and hence make some subset
               S of the clauses unsat that would be sat if the mcs was
               not blocked. Model rotation could generate that MCS (or
               a superset) as a model of S and end up with 0
               violations. But it is expensive to detect this (we'd
               have to check the occurrences of all redundant clauses
               that are satisfied by flipping the current literal), so
               we just accept that this may happen and reject this
               flip. The only problem is we cannot just assert that
               numviolations > 0 anymore. */
            if (single_mus) {
                cout << "numviolations == 0. Balls\n";
                verify_unsat(solver->assumptions(), conflict);
                assert(0);
            }
            assert(!single_mus);
            break;
        }

        assert(numviolations > 0);
        if (numviolations > 1) {
            DOUT << "At least " << numviolations
                 << " violations or violates hard, skipping\n";
            continue;
        }
        assert(numviolations == 1);

        DBG
            if (!crits[toInt(~inrotateb)])
                verify_critical(solver->assumptions(), conflict, inrotateb);

        // numviolations == 1, we can rotate
        Lit b = inrotateb;
        DOUT << "Need in b" << b << " ex: b" << in2ex(b) << "\n";
        DOUT << "\t" << theWcnf.getGroupClauses(in2ex(b)) << "\n";
        if (!crits[toInt(~b)]) {
            ++numadded;
            ++model_rotate_clauses;
        }
        make_crit(b);
        // XXX: we could call reportSingletonCS here, but it takes too
        // much time in large instances (reportCS needs to find all
        // redundant clauses violated by the rotated model)

        if (recursive_mr)
            rotations[inflipvar] = b;
    }

    assert(recursive_mr || rotations.empty());
    for (auto& r : rotations) {
        Var inflipvar = r.first;
        Lit b = r.second;
        model[inflipvar] = model[inflipvar] ^ 1;
        model[var(b)] = model[var(b)] ^ 1;
        model[var(inviolated)] = model[var(inviolated)] ^ 1;
        numadded += modelRotate_(
            rotsolver, conflict, b, model, can_rotate_to_p, make_crit);
        model[inflipvar] = model[inflipvar] ^ 1;
        model[var(b)] = model[var(b)] ^ 1;
        model[var(inviolated)] = model[var(inviolated)] ^ 1;
    }

    return numadded;
}

auto MUSSolver::findCS(
    Solver& mcssolver, vector<Lit> const& conflict, vector<Lit>& cs) -> lbool
{
    cs.clear();
    auto cnfl_before = mcssolver.getStats().conflicts;
    auto val = satsolve(mcssolver, CS_SOLVE);
    if( val != l_True ) {
        ACT("cs") << " val = " << val << std::endl;
        return val;
    }

    // l_True
    for(Lit l : conflict) {
        if(mcs_witness.value(var(l)) == l_True)
            cs.push_back(l);
    }
    assert( !cs.empty() );
    ACT("cs") << " val = " << val << " CS size = " << cs.size()
              << " conflicts = " << (mcssolver.getStats().conflicts - cnfl_before)
              << std::endl;
    if( cs.size() == 1 )
        ++csadded;

    return val;
}

auto MUSSolver::findMCS(vector<Lit> const& conflict, vector<Lit> const& all,
                  bool require_global_mcs, bool block_mcs,
                  vector<Lit>& cs) -> lbool
{
    return findMCS(*solver, conflict, all, require_global_mcs, block_mcs, cs);
}

auto MUSSolver::findMCS(Solver &mcssolver,
                  vector<Lit> const& conflict, vector<Lit> const& all,
                  bool require_global_mcs, bool block_mcs,
                  vector<Lit>& cs) -> lbool
{
    current_cs_size = conflict.size();

    lbool val{findCS(mcssolver, conflict, cs)};
    if( val != l_True )
        return val;

    return findMCS_hint(mcssolver, conflict, all, require_global_mcs, block_mcs,
        cs, mcs_witness);
}

auto MUSSolver::findMCS_hint(vector<Lit> const& conflict,
    vector<Lit> const& all, bool require_global_mcs, bool block_mcs,
    vector<Lit>& cs, model_t const& model) -> lbool
{
    return findMCS_hint(
        *solver, conflict, all, require_global_mcs, block_mcs, cs, model);
}

auto MUSSolver::findMCS_hint(Solver& mcssolver, vector<Lit> const& conflict,
    vector<Lit> const& all, bool require_global_mcs, bool block_mcs,
    vector<Lit>& cs, model_t const& model) -> lbool
{
    current_cs_size = cs.size();
    mcs_witness = model;

    // note if require_global_mcs == true, we do two-phase
    // minimization: first minimal wrt conflict, then globally minimal

    bool minimized = true;
    if( cs.size() > 1  ) {
        minimized = minimizeCS(mcssolver, conflict, cs,
                               block_mcs && !require_global_mcs);
    }

    if( require_global_mcs ) {
        auto before = cs.size();
        DBG {
            cout << "local cs " << cs << "\n";
            bool foundr{false};
            for(Lit l : cs) {
                assert(!sign(l));
                assert(!crits[toInt(~l)]);
                assert(maybe[toInt(l)] || redundant[toInt(l)]);
                if( maybe[toInt(l)] )
                    cout << "M ";
                if( redundant[toInt(l)] ) {
                    foundr = true;
                    cout << "R ";
                }
            }
            cout << "\n";
            assert(!foundr);
        }
        cs.clear();
        for(Lit l : all)
            if( mcs_witness.value( l ) == l_True )
                cs.push_back(l);
        DBG if( cs.size() == before )
            cout << "global cs " << cs << ", no other clauses violated\n";
        if( !minimized )
            return l_Undef;
        if( cs.size() == before ) {
            if( block_mcs ) {
                // it was not blocked in minimizeCS, in case it was
                // locally but not globally minimal
                mcssolver.cancelUntil(0);
                if( &mcssolver == solver ) {
                    vector<Lit> ps;
                    for(Lit l : cs)
                        ps.push_back(~l);
                    theWcnf.addClause(ps, false); // this is so we also update
                                                  // the occlists, see comment
                                                  // in minimizeCS_cl
                } else {
                    vector<Lit> ps;
                    for(Lit l : cs)
                        ps.push_back(~l);
                    mcssolver.addClause(ps);
                }
            }
            reportCS(cs, mcs_witness, require_global_mcs);
            return l_True;
        }
        ACT("mcs") << " local MCS size " << before << " global " << cs.size()
                   << " minimizing global" << std::endl;
        minimized =  minimizeCS(mcssolver, all, cs, block_mcs);
        DBG {
            cout << "global cs " << cs << "\n";
            for(Lit l : cs) {
                assert(!sign(l));
                assert(!crits[toInt(~l)]);
                assert(maybe[toInt(l)] || redundant[toInt(l)]);
                if( maybe[toInt(l)] )
                    cout << "M ";
                if( redundant[toInt(l)] )
                    cout << "R ";
            }
            cout << "\n";
        }
    }

    reportCS(cs, mcs_witness, require_global_mcs);

    mcssolver.cancelUntil(0);

    current_cs_size = 0;

    if( !minimized )
        return l_Undef;

    return l_True;
}

bool MUSSolver::minimizeCS(Solver &mcssolver,
                     vector<Lit> const& conflict,
                     vector<Lit>& cs,
                     bool block_mcs)
{
    static const bool relaxation_search = false;
    if( relaxation_search )
        return minimizeCS_rs(mcssolver, conflict, cs, block_mcs);
    else
        return minimizeCS_cl(mcssolver, conflict, cs, block_mcs);
}

bool MUSSolver::minimizeCS_rs(Solver &mcssolver,
                        vector<Lit> const& conflict,
                        vector<Lit>& cs,
                        bool /*block_mcs*/)
{
    for(Lit l : cs) {
        mcssolver.setVarPriority(var(l), 1);
        mcssolver.setVarPolarity(var(l), l_True);
    }

    auto asize = mcssolver.assumptions().size();
    for(Lit l :conflict) {
        if(mcssolver.varPriority(var(l)) == 0)
            mcssolver.assumptions().push_back(~l);
    }

    auto val = satsolve(mcssolver, RS_SOLVE);
    mcssolver.shrinkAssumptions(mcssolver.assumptions().size() - asize);
    assert(mcssolver.assumptions().size()
        == static_cast<unsigned>(crits_size) + assumptions.size());

    for(Lit l : cs) {
        mcssolver.setVarPriority(var(l), 0);
        mcssolver.setVarPolarity(var(l), l_Undef);
    }

    assert( val != l_False);
    if( val == l_True ) {
        erase_if(cs, [&](Lit l) { return mcs_witness[var(l)] != l_True; });
        assert(!cs.empty());
    }
    return val == l_True;
}

bool MUSSolver::minimizeCS_cl(Solver &mcssolver,
                        vector<Lit> const& conflict,
                        vector<Lit>& cs,
                        bool block_mcs)
{
    auto& incs = minimizecs_incs;
    auto& ass = minimizecs_ass;
    incs.resize(2*solver->nUsedVars());
    ass.resize(2*solver->nUsedVars());

    for(Lit l : cs)
        incs[toInt(l)] = 1;
    for(Lit l : mcssolver.assumptions())
        ass[toInt(l)] = 1;

    Var blockingvar{var_Undef};

    auto asize = mcssolver.assumptions().size();

    for(Lit l : conflict)
        if( !incs[toInt(l)] && !ass[toInt(~l)] )
            mcssolver.assumptions().push_back(~l);

    for(Lit l : cs)
        incs[toInt(l)] = 0;
    for(Lit l : mcssolver.assumptions())
        ass[toInt(l)] = 0;

    if( !block_mcs ) {
        blockingvar = mcssolver.newVar();
        mcssolver.assumptions().push_back(~mkLit(blockingvar));
    }

    // ps: added to the solver
    // final_blocking_clause: add to the Wcnf for model rotation (in
    // non-single-mus mode only)
    vector<Lit> ps, final_blocking_clause;

    // for keeping stats on changes between successive models
    model_t last_model{mcs_witness};

    auto mcsprevconf = mcssolver.getStats().conflicts;

    lbool val;
    do {
        if (mcs_rotate)
            mcsRotate(mcssolver, mcssolver.assumptions(), cs);
        if(mcs_backbone)
            mcsFindBackboneLits(mcssolver, cs);
        if( cs.size() == 1 ) {
            //cout << "minimizing CS with size 1. hm\n";
            continue;
        }
        ps.clear();
        for(Lit l : cs)
            ps.push_back(~l);
        if( !block_mcs )
            ps.push_back(mkLit(blockingvar));
        final_blocking_clause = ps;
        mcssolver.addClause(ps);

        val = satsolve(mcssolver, MCS_SOLVE);

        if( val == l_True ) {
            int diff = 0;
            for (unsigned v = 0; v != mcs_witness.size(); ++v)
                if (!theWcnf.isBvar(in2ex(v)) && last_model.size() > v
                    && mcs_witness[v] != last_model[v])
                    ++diff;
            last_model = mcs_witness;
            current_cs_size = cs.size();

            for(Lit l : cs)
                if( mcs_witness[var(l)] == l_False )
                    mcssolver.assumptions().push_back(~l);
            erase_if(cs, [&](Lit l) { return mcs_witness[var(l)] == l_False; });
            ACT("cs") << " minimizing CS size " << cs.size()
                      << " conflicts = " << (mcssolver.getStats().conflicts - mcsprevconf)
                      << " diff = " << diff
                      << std::endl;
            if( diff == 1 )
                ++easymcses;
            // In the group MUS problem or when !fbeq, diff == 0 is a
            // possibility: we can set the x vars so that they satisfy
            // a group but that group's bvar is true, because it is
            // not an equivalence.
            if( diff == 0 )
                ++trivialmcses;
            assert(!cs.empty());
            mcsprevconf = mcssolver.getStats().conflicts;
        }
    } while( val == l_True && cs.size() > 1 );

    ACT("cs") << " val = " << val << ", finishing"
              << ", cs size = " << cs.size()
              << std::endl;

    mcssolver.cancelUntil(0);
    mcssolver.shrinkAssumptions(mcssolver.assumptions().size() - asize);
    if( !block_mcs )
        mcssolver.releaseVar(mkLit(blockingvar));
    else if( &mcssolver == solver && !single_mus ) {
        /* If we do not do this, model rotation may stop working
           because it may think it has found a model which in fact
           violates a blocking clause. This can only happen if
           !single_mus and block_mcs. Additionally, just to avoid
           trouble with the all-mus case, we check that &mcssolver ==
           &solver */
        if( final_blocking_clause.size() > 1 ) {
            vector<Lit> exfbc;
            in2ex(final_blocking_clause, exfbc);
            theWcnf.addClause(exfbc, false);
        }
    }

    return (val == l_False || cs.size() == 1);
}

void MUSSolver::mcsRotate(
    Solver& mcssolver, vector<Lit> const& assumptions, vector<Lit>& cs)
{
    auto start_time = cpuTime();
    auto removals_before = mcs_rotate_removals;

    control->checkPrintStatusLine();

    auto const& model = mcs_witness;

    // gather all literals in the CS
    vector<bool> exlitsbm(theWcnf.nVars());
    vector<Lit> exlits;
    for (Lit inbv : cs) {
        Lit exbv{ in2ex(inbv) };
        if (theWcnf.getGroupClauses(exbv).empty()) {
            cs.erase(find(begin(cs), end(cs), inbv));
            return;
        }
        if (theWcnf.getGroupClauses(exbv).size() > 1)
            // we cannot handle groups of size > 1 here
            return;
        auto cl = theWcnf.getGroupLiterals(exbv);
        for (Lit exl : cl)
            if (!exlitsbm[var(exl)]) {
                exlits.push_back(exl);
                exlitsbm[var(exl)] = true;
            }
    }

    // mark all the clauses in 'conflicts', they cannot be violated by flips
    vector<bool> exassumptionbm(theWcnf.nVars());
    for (Lit inl : assumptions) {
        Lit exl = in2ex(inl);
        if( exl == lit_Undef )
            continue;
        exassumptionbm[var(in2ex(inl))] = true;
    }

    // for each of the literals, try to flip it and keep this rotation
    // if it does not violate any hard clause or any of the clauses in
    // conflict
    for (Lit exfliplit : exlits) {
        Lit infliplit{ ex2in(exfliplit) };
        if (!canFlip(~infliplit, model, ~infliplit))
            continue;
        bool satallothers{true};
        for (Lit exbv : theWcnf.softBvarOccurrences(~exfliplit)) {
            if (!exassumptionbm[var(exbv)])
                continue;
            auto& grclauses = theWcnf.getGroupClauses(exbv);
            auto& excl = theWcnf.getClause(grclauses[0]);
            bool sat{ false };
            for (Lit exl : excl) {
                assert(exl != exfliplit);
                if (exl != ~exfliplit && model.value(ex2in(exl)) == l_True) {
                    sat = true;
                    break;
                }
            }
            if (!sat) {
                satallothers = false;
                break;
            }
        }
        if (satallothers) {
            // keep the flip, remove clauses that it satisfies from cs
            erase_if(cs, [&](Lit inb) {
                Lit exb{ in2ex(inb) };
                auto& grclauses = theWcnf.getGroupClauses(exb);
                auto& cl = theWcnf.getClause(grclauses[0]);
                if (std::find(begin(cl), end(cl), exfliplit) != end(cl)) {
                    ++mcs_rotate_removals;
                    mcssolver.assumptions().push_back(~inb);
                    exassumptionbm[var(exb)] = true;
                    mcs_witness.model[var(inb)] = l_False;
                    return true;
                }
                return false;
            });
            mcs_witness.model[var(infliplit)]
                = mcs_witness.model[var(infliplit)] ^ 1;
            mcssolver.initializePhase(
                var(infliplit), sign(infliplit) ? l_True : l_False);
        }
    }

    mcs_rotate_times.push_back(cpuTime() - start_time);
    ACT("MCSrot") << ", " << (mcs_rotate_removals - removals_before)
                  << " clauses removed\n";
}

void MUSSolver::mcsFindBackboneLits(Solver& mcssolver, vector<Lit> const& cs)
{
    control->checkPrintStatusLine();

    // count occurrences for all literals in the CS
    vector<size_t> exlitscnt(theWcnf.nVars());
    vector<Lit> exlits;
    for (Lit inbv : cs) {
        Lit exbv{ in2ex(inbv) };
        auto& grclauses = theWcnf.getGroupClauses(exbv);
        if (grclauses.size() > 1)
            return;
        auto& cl = theWcnf.getClause(grclauses[0]);
        for (Lit exl : cl) {
            if (!exlitscnt[var(exl)])
                exlits.push_back(exl);
            ++exlitscnt[var(exl)];
        }
    }

    // we do not want to add it again if it is in assumptions already
    for (Lit l : mcssolver.assumptions())
        exlitscnt[var(in2ex(l))] = 0;

    for (Lit exl : exlits) {
        if (exlitscnt[var(exl)] == cs.size()) {
            ++mcs_backbone_literals;
            mcssolver.assumptions().push_back(ex2in(~exl));
        }
    }
}

void MUSSolver::reportCS(vector<Lit> const& cs,
                   model_t const& model,
                   bool minimal)
{
    using namespace std::placeholders;
    using std::ref;

    if( !report_cs )
        return;
    vector<Lit>& extended = external_csmus;
    in2ex(cs, extended);
    for(Lit l : extended)
        assert(!sign(l));

    std::vector<bool>& seen = reportcs_seen;
    seen.resize(solver->nUsedVars());
    for(Lit l : cs)
        seen[var(l)] = 1;
    DOUT << "ReportCS " << cs << "\n";

    in2ex(model, external_model);

    // extend with violated redundant clauses.
    // if we only had models generated by calls to the solver, we
    // would just need to check the value of the bvar, but model
    // rotation generates models with incorrect values for the bvars
    for (Var exv : theWcnf.bvars()) {
        Var v = ex2in(exv);
        assert(theWcnf.bvarIdx(exv) >= 0);
        // XXX: is this line below incorrect?
        // if( !redundant[toInt(mkLit(v))] )
        //  continue;
        if (seen[v])
            continue;

        auto& cls = theWcnf.getGroupClauses(exv);
        bool sat = std::all_of(begin(cls), end(cls),
            std::bind(modelSatisfies, ref(external_model), ref(theWcnf), _1));
/*        if(! ((!sat || model[v] == l_False) && (sat || model[v] == l_True))) {
            cout << "bvar = " << mkLit(exv) << " group = " << theWcnf.bvarIdx(exv) << "\n";
            cout << "interval bvar = " << mkLit(v) << "\n";
            cout << "value(bvar) = " << model.value(v) << "\n";
            cout << "soft clause sat = " << sat << "\n";
            cout << "group size = " << cls.size() << "\n";
            for(auto cr : cls) {
                auto& cl = theWcnf.getClause(cr);
                cout << cl << "\n";
                cout << "lit values: ";
                for(Lit l : cl)
                    cout << model.value(ex2in(l)) << " ";
                cout << "\n";
                cout << "var values: ";
                for(Lit l : cl)
                    cout << model.value(ex2in(var(l))) << " ";
                cout << "\n";
                vector<Lit> incl(begin(cl), end(cl));
                for(Lit& l : incl)
                    l = ex2in(l);
                cout << "in: " << incl << "\n";
            }
        }
*/        //assert(((!sat || model[v] == l_False) && (sat || model[v] == l_True)));
        if (sat)
            continue;

        if (minimal || crits[toInt(~mkLit(v))]) {
            cout << "Oops: Also violates " << mkLit(v) << std::endl;
            assert(0);
        } else
            DOUT << "Also violates in" << mkLit(v) << " ex" << in2ex(mkLit(v))
                 << std::endl;
        extended.push_back(in2ex(mkLit(v)));
        assert(!minimal);
    }
    if (extended.size() == 1)
        minimal = true;

    for(Lit l : cs)
        seen[var(l)] = 0;

    report_cs(extended, external_model, minimal);
}

void MUSSolver::reportSingletonCS(Lit b, model_t const& model,
                            bool minimal)
{
    assert(!sign(b));
    if( !report_cs )
        return;
    reportsingleton_cs.clear();
    reportsingleton_cs.push_back(b);
    reportCS(reportsingleton_cs, model, minimal);
}

void MUSSolver::reportMUS(vector<Lit> const& assumptions, bool minimal)
{
    if( !report_mus )
        return;
    vector<Lit> mus;
    mus.reserve(assumptions.size());
    for(Lit l : assumptions) {
        assert(sign(l));
        mus.push_back(~l);
    }
    in2ex(mus, external_csmus);
    report_mus(external_csmus, minimal);
}

bool MUSSolver::find_mus(vector<Lit>& themus, bool require_global_mcses)
{
    themus.clear();
    auto bv = theWcnf.bvars();
    themus.resize(bv.size());
    std::transform(begin(bv), end(bv), begin(themus),
                   [](Var x) { return mkLit(x); });
    return mus(themus, require_global_mcses);
}

auto MUSSolver::solve_mus(vector<Lit>& themus, vector<Lit> const& exassumptions,
    bool require_global_mcses) -> std::pair<lbool, bool>
{
    solver->cancelUntil(0);
    ensureSoftCls(themus);
    ensureSoftCls(exassumptions);

    assumptions.clear();
    for(Lit l: exassumptions)
        assumptions.push_back(ex2in(l));
    preProcessConflict(themus);

    solver->assumptions().clear();
    for (Lit l : assumptions)
        solver->assumptions().push_back(l);
    for (Lit l : themus)
        solver->assumptions().push_back(~l);
    lbool r1 = satsolve(*solver, TEST_SOLVE);
    ACT("test") << " val = " << r1 << std::endl;
    solver->assumptions().clear();
    if( r1 == l_True )
        return std::make_pair(r1, true);
    else if( r1 == l_Undef )
        return std::make_pair(r1, false);

    // unsat
    themus.clear();
    for(Lit l : solver->core().toVec())
        themus.push_back(in2ex(l));
    if( themus.empty() )
        return std::make_pair(l_False, true);
    return std::make_pair(l_False, mus(themus, require_global_mcses));
}

auto MUSSolver::solve_mus(vector<Lit>& themus, bool require_global_mcses)
    -> std::pair<lbool, bool>
{
    return solve_mus(themus, {}, require_global_mcses);
}

bool MUSSolver::find_all_mus()
{
    vector<Lit> themus;
    auto bv = theWcnf.bvars();
    themus.resize(bv.size());
    std::transform(begin(bv), end(bv), begin(themus),
                   [](Var x) { return mkLit(x); });
    return allmus(themus);
}

bool MUSSolver::find_some_mus()
{
    vector<Lit> themus;
    auto bv = theWcnf.bvars();
    themus.resize(bv.size());
    std::transform(begin(bv), end(bv), begin(themus),
                   [](Var x) { return mkLit(x); });
    return somemus(themus);
}

void MUSSolver::recordSolve(lbool res, solve_type st)
{
    recordSolve(*solver, res, st);
}

void MUSSolver::recordSolve(Solver& rssolver, lbool res, solve_type st)
{
    ++satSolves;
    double time = cpuTime();
    if( res == l_True ) {
        solvestats[st].conflicts_sat.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats[st].times_sat.push_back(time - prevTime);
        solvestats_total.conflicts_sat.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats_total.times_sat.push_back(time - prevTime);
    } else if( res == l_False ) {
        solvestats[st].conflicts_unsat.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats[st].times_unsat.push_back(time - prevTime);
        solvestats_total.conflicts_unsat.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats_total.times_unsat.push_back(time - prevTime);
    } else {
        solvestats[st].conflicts_undef.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats[st].times_undef.push_back(time - prevTime);
        solvestats_total.conflicts_undef.push_back(rssolver.getStats().conflicts - prevConflicts);
        solvestats_total.times_undef.push_back(time - prevTime);
    }
    prevConflicts = rssolver.getStats().conflicts;
    prevTime = time;
}

void MUSSolver::print_mus(vector<Lit> const& mus, std::ostream& ofs)
{
    int nv{0};
    int nc{0};
    auto compute_max = [&](Lit b) {
        for (auto cr : theWcnf.getGroupClauses(b)) {
            ++nc;
            for (auto l : theWcnf.getClause(cr))
                nv = std::max(nv, var(l));
        }
    };

    auto b0{mkLit(theWcnf.varFromBvarIdx(0))};

    compute_max(b0);
    for_each(begin(mus), end(mus), compute_max);

    auto print_clauses = [&](Lit b) {
        for (auto cr : theWcnf.getGroupClauses(b)) {
            for (auto l : theWcnf.getClause(cr))
                ofs << (1 - 2 * sign(l)) * (var(l) + 1) << ' ';
            ofs << "0\n";
        }
    };

    ofs << "p cnf " << (nv + 1) << ' ' << nc << "\n";
    print_clauses(b0);
    for_each(begin(mus), end(mus), print_clauses);
}

void MUSSolver::print_gmus(vector<Lit> const& mus, std::ostream& ofs)
{
    int nv{0};
    int nc{0};
    int ng{0};
    auto compute_max = [&](Lit b) {
        ng = std::max(ng, theWcnf.bvarIdx(b));
        for (auto cr : theWcnf.getGroupClauses(b)) {
            ++nc;
            for (auto l : theWcnf.getClause(cr))
                nv = std::max(nv, var(l));
        }
    };

    auto b0{mkLit(theWcnf.varFromBvarIdx(0))};

    compute_max(b0);
    for_each(begin(mus), end(mus), compute_max);

    auto print_clauses = [&](Lit b) {
        for (auto cr : theWcnf.getGroupClauses(b)) {
            ofs << "{" << theWcnf.bvarIdx(b) << "} ";
            for (auto l : theWcnf.getClause(cr))
                ofs << (1 - 2 * sign(l)) * (var(l) + 1) << ' ';
            ofs << "0\n";
        }
    };

    ofs << "p gcnf " << (nv + 1) << ' ' << nc << ' ' << ' ' << ng << "\n";
    print_clauses(b0);
    for_each(begin(mus), end(mus), print_clauses);
}

void MUSSolver::print_mss(vector<Lit> const& mcs, std::ostream& ofs)
{
    std::set<Lit> mcset(begin(mcs), end(mcs));
    int nv{0}, nc{0};
    for(auto b : theWcnf.bvars()) {
        if( mcset.find(mkLit(b)) != mcset.end() )
            continue;
        for(auto cr : theWcnf.getGroupClauses(mkLit(b))) {
            for(auto l : theWcnf.getClause(cr))
                nv = std::max(nv, var(l));
            ++nc;
        }
    }
    ofs << "p cnf " << (nv+1) << ' ' << nc << "\n";
    for(auto b : theWcnf.bvars()) {
        if( mcset.find(mkLit(b)) != mcset.end() )
            continue;
        for(auto cr : theWcnf.getGroupClauses(mkLit(b))) {
            for(auto l : theWcnf.getClause(cr))
                ofs << (1-2*sign(l))*(var(l)+1) << ' ';
            ofs << "0\n";
        }
    }
}

template<typename Num>
struct series_stats {
    int n{};
    Num min{}, max{}, sum{};
    Num median{}, perc90{}; // 50th and 90th percentile
};

namespace{

    struct identity {
        template<typename T>
        auto operator()(T&& t) { return t; }
    };

    /* helpers for printing stats */
    template <typename T, typename Proj = identity,
        typename Num = typename std::remove_reference<
            typename std::result_of<Proj(T)>::type>::type>
    series_stats<Num> computeSeriesStats(
        std::vector<T> const& vp, Proj proj = Proj{})
    {
        series_stats<Num> s{};

        vector<Num> v;
        v.reserve(vp.size());
        for (auto&& e : vp)
            v.push_back(std::invoke(proj, e));

        s.n = vp.size();
        s.sum = accumulate(begin(v), end(v), Num{});
        if (!vp.empty()) {
            nth_element(begin(v), begin(v), end(v));
            nth_element(begin(v), begin(v) + v.size() / 2, end(v));
            nth_element(begin(v), begin(v)+9*v.size()/10, end(v));
            nth_element(begin(v), prev(end(v)), end(v));
            s.min=v[0];
            s.max=v.back();
            s.median=v[v.size()/2];
            s.perc90=v[9*v.size()/10];
        }
        return s;
    }

    template<typename Num>
        void reportSeriesStats(const char *cat,
                               const char *type,
                               const char *format,
                               const char *avgformat,
                               bool printn,
                               series_stats<Num> const& ss)
    {
        printf("%21s :", cat);
        if( printn )
            printf("%5d, ", ss.n);
        else
            printf("       ");
        printf(format, ss.sum);
        printf(" %s, ", type);
        printf(avgformat, ss.sum/double(ss.n));
        printf(" avg, ");
        printf(format, ss.median);
        printf(" median, ");
        printf(format, ss.perc90);
        printf(" 90th %%, ");
        printf(format, ss.min);
        printf(" min, ");
        printf(format, ss.max);
        printf(" max\n");
    }
}

void MUSSolver::reportSolveStats(const char *type, solve_stats const& st) const
{
    auto cfsat = computeSeriesStats(st.conflicts_sat),
        cfunsat = computeSeriesStats(st.conflicts_unsat),
        cfundef = computeSeriesStats(st.conflicts_undef);
    auto tmsat = computeSeriesStats(st.times_sat),
        tmunsat = computeSeriesStats(st.times_unsat),
        tmundef = computeSeriesStats(st.times_undef);

    printf("%s\n", type);

    vector<uint64_t> cf(st.conflicts_sat);
    cf.insert(cf.end(), begin(st.conflicts_unsat), end(st.conflicts_unsat));
    cf.insert(cf.end(), begin(st.conflicts_undef), end(st.conflicts_undef));
    vector<double> tm(st.times_sat);
    tm.insert(tm.end(), begin(st.times_unsat), end(st.times_unsat));
    tm.insert(tm.end(), begin(st.times_undef), end(st.times_undef));
    auto cfall = computeSeriesStats(cf);
    auto tmall = computeSeriesStats(tm);

    reportSeriesStats("solver invocations", "conflicts", "%" PRIu64, "%.0f", true, cfall);
    reportSeriesStats("", "sec", "%.3f", "%.3f", false, tmall);
    printf("                      :       %.0f leaves/sec (conflicts+solutions)\n",
           (cfall.sum + cfsat.n)/double(tmall.sum));

    if( cfsat.n ) {
        reportSeriesStats("   SAT", "conflicts", "%" PRIu64, "%.0f", true, cfsat);
        reportSeriesStats("", "sec", "%.3f", "%.3f", false, tmsat);
        printf("                      :       %.0f leaves/sec (conflicts+solutions)\n",
               (cfsat.sum + cfsat.n)/double(tmsat.sum));
    }
    if( cfunsat.n ) {
        reportSeriesStats("   UNSAT", "conflicts", "%" PRIu64, "%.0f", true, cfunsat);
        reportSeriesStats("", "sec", "%.3f", "%.3f", false, tmunsat);
        printf("                      :       %.0f conflicts/sec\n",
               (cfunsat.sum)/double(tmunsat.sum));
    }
    if( cfundef.n ) {
        reportSeriesStats("   interrupted", "conflicts", "%" PRIu64, "%.0f", true, cfundef);
        reportSeriesStats("", "sec", "%.3f", "%.3f", false, tmundef);
        printf("                      :       %.0f conflicts/sec\n",
               (cfundef.sum)/double(tmundef.sum));
    }
}

void MUSSolver::printStats() const
{
    reportSolveStats("total", solvestats_total);
    for(auto& st : solvestats)
        reportSolveStats(solve_type_desc[st.first].c_str(), st.second);

    reportSeriesStats("Model Rotations      ", "sec", "%.3f", "%.3f", true,
        computeSeriesStats(mr_times));
    reportSeriesStats("MCS-MUS Iterations   ", "instances", "%" PRIu64, "%.3f",
        true, computeSeriesStats(mishmash_iteration_stats,
                          &iteration_stats::satsolves));
    reportSeriesStats("", "sec", "%.3f", "%.3f", false,
        computeSeriesStats(mishmash_iteration_stats, &iteration_stats::time));
    reportSeriesStats("", "conflicts", "%" PRIu64, "%.3f", false,
        computeSeriesStats(mishmash_iteration_stats,
                          &iteration_stats::conflicts));
    reportSeriesStats("", "clauses", "%" PRIu64, "%.3f", false,
        computeSeriesStats(mishmash_iteration_stats,
                          &iteration_stats::clauses_determined));
    reportSeriesStats("Delete Iterations    ", "instances", "%" PRIu64, "%.3f",
        true, computeSeriesStats(
                          delete_iteration_stats, &iteration_stats::satsolves));
    reportSeriesStats("", "sec", "%.3f", "%.3f", false,
        computeSeriesStats(delete_iteration_stats, &iteration_stats::time));
    reportSeriesStats("", "conflicts", "%" PRIu64, "%.3f", false,
        computeSeriesStats(
                          delete_iteration_stats, &iteration_stats::conflicts));
    reportSeriesStats("", "clauses", "%" PRIu64, "%.3f", false,
        computeSeriesStats(delete_iteration_stats,
                          &iteration_stats::clauses_determined));

    auto mcsstatsi = solvestats.find(MCS_SOLVE);
    if( mcsstatsi != solvestats.end() ) {
        auto& mcsstats = mcsstatsi->second;
        auto mcssolved = mcsstats.times_sat.size()
            + mcsstats.times_unsat.size()
            + mcsstats.times_undef.size();
        printf("Easy MCS instances    : %" PRIu64 " (%.2f %%)\n", easymcses,
               easymcses*100.0/mcssolved);
        printf("Trivial MCS instances : %" PRIu64 " (%.2f %%)\n", trivialmcses,
               trivialmcses*100.0/mcssolved);
        printf("MCS rotations         : %d (%.2f %%)\n", mcs_rotate_removals,
            mcs_rotate_removals * 100.0 / mcssolved);
        printf("MCS backbone literals : %d (%.1f/iteration)\n",
            mcs_backbone_literals,
            mcs_backbone_literals / double(mishmash_iteration_stats.size()));
    }
    printf("MR triggered resets   : %" PRIu64" (%" PRIu64" blocked)\n",
           mr_trigger_resets, mr_trigger_blocks);
    if( !all_muses ) {
        printf("Clauses\n");
        printf("   added by MR        : %d\n", model_rotate_clauses);
        printf("   added by aMC       : %d\n", amo_clauses);
        printf("   added by test      : %d\n", test_clauses_add);
        printf("   added by CS        : %d\n", csadded);
        printf("   added by MCS       : %d\n", mcsadded);
        printf("   removed by aMC     : %d\n", amo_clauses_removed);
        printf("   removed by test    : %d\n", refine_clauses);
        printf("   removed by MCS     : %d\n", mcsremoved);
        printf("   removed by MCS ref : %d\n", mcs_refine_removed);
        printf("   removed by CS unsat: %d\n", csunsatremoved);
    }


    Solver const& s(*solver); // short alias

    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    auto st = s.getStats();
    printf("restarts              : %" PRIu64"\n", st.starts);
    printf("leaves                : %-12" PRIu64"   (%.0f /sec)\n", st.conflicts+st.solutions   , (st.conflicts+st.solutions)   /cpu_time);
    printf("  conflicts           : %-12" PRIu64"   (%.0f /sec)\n", st.conflicts   , st.conflicts   /cpu_time);
    printf("  solutions           : %-12" PRIu64"   (%.0f /sec)\n", st.solutions   , st.solutions   /cpu_time);
    printf("decisions             : %-12" PRIu64"   (%4.2f %% random) (%.0f /sec)\n", st.decisions, (float)st.rnd_decisions*100 / (float)st.decisions, st.decisions   /cpu_time);
    printf("propagations          : %-12" PRIu64"   (%.0f /sec)\n", st.propagations, st.propagations/cpu_time);
    printf("conflict literals     : %-12" PRIu64"   (%4.2f %% deleted)\n", st.tot_literals, (st.max_literals - st.tot_literals)*100 / (double)st.max_literals);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("Preprocessing time    : %g s (%.2f %%)\n",
           st.pre_time,
           st.pre_time/cpu_time*100.);
    printf("CPU time              : %g s\n", cpu_time);
}

std::ostream& MUSSolver::report_activity(std::string const& act) const
{
    double this_activity_time = cpuTime();
    cout << std::setw(9) << this_activity_time
         << std::setw(9) << (this_activity_time - last_activity_time)
         << std::setw(8) << act
         << " ncrits " << crits_size
         << " maybe " << maybe_size
         << " redundant " << redundant_size;
    last_activity_time = this_activity_time;
    return cout;
}

void MUSSolver::report_activityln(std::string const& act) const
{
    report_activity(act) << std::endl;
}

void MUSSolver::printStatusLine()
{
    Solver const& s(*solver); // short alias
    double time = cpuTime();
    if(all_muses) {
        assert(hssolver);
        Solver const& h(*hssolver);
        printf("| %8.2f |%9" PRIu64" %9" PRIu64" %15" PRIu64" |"
               "%9d %9d %9" PRIu64" |"
               "%9d %9d %9" PRIu64" |"
               "%9d %9d %9d |\n",
               time, num_muses, num_mcses, num_singleton_mcses,
               s.nVars(),
               s.nClauses(), s.getStats().conflicts+s.getStats().solutions,
               h.nVars() - theWcnf.nXVars(),
               h.nClauses(), h.getStats().conflicts,
               crits_size, maybe_size, redundant_size);
        fflush(stdout);
        return;
    }
    printf("| %8.2f | %10lu | %9d %9d %9d | %9d %9d %9d | %7d %8d %6.0f | "
           "%6.3f %% |\n",
        time, mishmash_iteration_stats.size() + delete_iteration_stats.size(),
        crits_size, maybe_size, redundant_size,
        (int)(s.getStats().conflicts + s.getStats().solutions),
        (int)s.getStats().conflicts, (int)s.getStats().solutions,
        (int)s.nVars(), s.nClauses(),
        (double)s.getStats().learnts_literals / s.getStats().num_learnts,
        100.0 * (crits_size + redundant_size)
            / double(crits_size + redundant_size + maybe_size));
    fflush(stdout);
}

void MUSSolver::printHeader(double parse_time) const
{
    if( all_muses ) {
        printf("=======================================[ Problem Statistics ]==================================================\n");
        printf("|                                                                                                             |\n");
        printf("|  Number of variables:  %12d                                                                         |\n", theWcnf.nXVars());
        printf("|  Number of clauses:    %12d                                                                         |\n",
               theWcnf.nClauses());
        printf("|  Parse time:           %12.2f s                                                                       |\n", parse_time);
        printf("|                                                                                                             |\n");
        printf("=======================================[ Search Statistics ]==================================================================================\n");
        printf("|          |            Progress                |          MCS Solver          |           HS Solver          |         Current MUS          |\n");
        printf("|     Time |    MUSes     MCSes Singleton MCSes |     Vars   Clauses    Leaves |     Vars   Clauses Conflicts | Critical     Maybe Redundant |\n");
        printf("==============================================================================================================================================\n");
        return;
    }
    printf("=======================================[ Problem Statistics ]=================================================================\n");
    printf("|                                                                                                                            |\n");
    printf("|  Number of variables:  %12d                                                                                        |\n", theWcnf.nXVars());
    printf("|  Number of clauses:    %12d                                                                                        |\n",
           theWcnf.nClauses());
    printf("|  Parse time:           %12.2f s                                                                                      |\n", parse_time);
    printf("|                                                                                                                            |\n");
    printf("====================================================[ Search Statistics ]=====================================================\n");
    printf("|          |            |            Core               |            Leaves             |        SAT instance     | Progress |\n");
    printf("|     Time | Iterations |  Critical     Maybe Redundant |     Total Conflicts Solutions |    Vars  Clauses Lit/Cl |          |\n");
    printf("==============================================================================================================================\n");
}

void MUSSolver::printFooter() const
{
    if( all_muses )
        printf("===============================================================================================================\n");
    else
        printf("==============================================================================================================================\n");
}

void MUSSolver::verify_unsat(vector<Lit> const& assumptions,
                       vector<Lit> const& conflict)
{
    verify_unsat(*solver, assumptions, conflict);
}

void MUSSolver::verify_unsat(Solver& origsolver,
                       vector<Lit> const& assumptions,
                       vector<Lit> const& conflict)
{
    // do it in a different solver so as not to perturb the current
    // solver's state and hide the bug
    cout << "verifying unsatness\n";
    auto s = origsolver.clone_input();
    for(Lit l : assumptions) {
        s->addClause( {l} );
    }
    for(Lit l : conflict) {
        Lit unit = ~l;
        if( s->value(unit) == l_False )
            return; // unsat already
        s->addClause( {unit} );
    }
    lbool itsval = s->solve();
    ACT("verify") << " val = " << itsval
                  << " conflicts = " << s->getStats().conflicts
                  << std::endl;
    if(itsval != l_False) std::cout << "varify_unsat failed (it is sat)" << std::endl; //added by JB 
    assert(itsval == l_False);
}

void MUSSolver::verify_mcs(vector<Lit> const& assumptions,
                     vector<Lit> const& conflict,
                     vector<Lit> const& mcs)
{
    // do it in a different solver so as not to perturb the current
    // solver's state and hide the bug
    cout << "verifying MCSness of " << mcs << "\n";
    auto s = solver->clone_input();
    for(Lit l : assumptions) {
        s->addClause( { mkLit(var(l), sign(l)) } );
    }
    std::set<Lit> csset(begin(mcs), end(mcs));
    for(Lit l : conflict) {
        if( csset.find(l) != csset.end() )
            continue;
        s->addClause( { mkLit( var(l), !sign(l)) } );
    }
    vector<Lit> mcscl;
    for(Lit l : mcs)
        mcscl.push_back(mkLit(var(l), !sign(l)));
    s->addClause(mcscl);
    lbool itsval = s->solve();
    ACT("verify") << " val = " << itsval
                  << " conflicts = " << s->getStats().conflicts
                  << std::endl;
    if(itsval != l_False) std::cout << "varify_mcs failed (it is not a mcs)" << std::endl; //added by JB 
    assert(itsval == l_False);
}

void MUSSolver::verify_critical(vector<Lit> const& assumptions,
                          vector<Lit> const& conflict,
                          Lit crit)
{
    cout << "verifying criticality of " << crit << "\n";
    cout << "[but first ";
    verify_unsat(assumptions, conflict);
    cout << "]" << std::endl;
    cout << "[and then ";
    vector<Lit> mcs;
    mcs.push_back(crit);
    verify_mcs(assumptions, conflict, mcs);
    cout << "]" << std::endl;
}

} // namespace mcsmus
