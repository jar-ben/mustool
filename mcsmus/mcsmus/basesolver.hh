/***********[basesolver.hh]
Copyright (c) 2015 George Katsirelos

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

#ifndef MCSMUS_BASESOLVER_HH
#define MCSMUS_BASESOLVER_HH

#include "minisat/core/mcsmus_SolverTypes.h"
#include "mcsmus/control.hh"

#include <memory>

namespace mcsmus {

using Var = Minisat::Var;
const Var var_Undef = Minisat::var_Undef;

using Lit = Minisat::Lit;
const Lit lit_Undef = Minisat::lit_Undef;
const Lit lit_Error = Minisat::lit_Error;
using Minisat::mkLit;

using lbool = Minisat::lbool;
const lbool l_True = Minisat::l_True;
const lbool l_False = Minisat::l_False;
const lbool l_Undef = Minisat::l_Undef;

using LSet = Minisat::LSet;

struct Statistics {
    uint64_t solves, starts;

    // current problem size
    uint64_t dec_vars, num_clauses;

    // tree stats
    uint64_t decisions, rnd_decisions,
        propagations, conflicts;

    // learnt clause database size
    uint64_t num_learnts, clauses_literals,
        learnts_literals, max_literals,
        tot_literals;

    // solutions found
    uint64_t solutions;

    // time: total and for preprocessing
    double time, pre_time;
};

struct model_t
{
    std::vector<lbool> model;
    lbool value(Var x) const { return model[x]; }
    lbool value(Lit x) const { return value(var(x)) ^ sign(x); }

    lbool& operator[](Var x) { return model[x]; }
    lbool operator[](Var x) const { return model[x]; }

    model_t() {}
    explicit model_t(std::vector<lbool> const& m) : model(m) {}

    std::vector<lbool>::size_type size() const { return model.size(); }
    bool empty() const { return model.empty(); }
    void clear() { model.clear(); }
};

enum class prop_result { ok, failed };

class BaseSolver {
public:
    virtual ~BaseSolver() {}

    // --------------------------------------------------
    // Cloning

    // Create a clone of the solver
    virtual std::unique_ptr<BaseSolver> clone() const = 0;

    // create a fork: a clone where all assignments are made into unit
    // clauses.
    virtual std::unique_ptr<BaseSolver> fork() const = 0;

    // create an empty solver
    virtual std::unique_ptr<BaseSolver> clone_blank() const = 0;

    // a copy of the solver with no learnt clauses
    virtual std::unique_ptr<BaseSolver> clone_input() const = 0;

    // --------------------------------------------------
    // problem setup

    // add a clause. No need to explicitly add variables
    virtual void addClause(std::vector<Lit> const& ps) = 0;
    virtual void addClause(std::vector<Lit> const&& ps) = 0;

    // Ask the solver for a free variable id.
    virtual Var newVar() = 0;

    // Release a variable id by making the literal l true at the root,
    // hence removing it from all clauses. This id may in the future
    // be given back to the user by newVar()
    virtual void releaseVar(Lit l) = 0;

    // current and original size of the problem. nVars() is a
    // statistic: it is allowed to report the effective number of
    // variables, ignoring eliminated variables, either by BVE or by
    // finding unit clauses or any other method. nUsedVars() warrants
    // that there is no variable with greater index in the
    // solver. *** This is different from minisat ***
    virtual int nVars() const = 0;
    virtual int nUsedVars() const = 0;
    virtual int nClauses() const = 0;

    // --------------------------------------------------
    // unit propagation/backtracking interface

    // unit propagate. Returns 'ok' or 'failed'
    virtual prop_result propagate() = 0;

    // value of a variable/literal in the current state of the solver
    virtual lbool value(Var x) const = 0;
    virtual lbool value(Lit l) const = 0;
    // Get the level a variable was assigned in. Easily deducible
    // outside the solver, so given as a convenience here
    virtual int levelOf(Var x) const = 0;

    // The decision made at a specific level (assuming
    // 1<=lvl<=decisionLevel()). Can be lit_Undef, as there are cases
    // (e.g., assumptions) where a decision level is empty
    virtual Lit decisionAtLevel(int lvl) const = 0;

    // current set of literals forced at level 0
    virtual std::vector<Lit> entailed() const = 0;

    // current set of forced literals, *given the currrent state of
    // the solver*. Therefore it is only useful in conjunction with
    // the 3 functions below (enqueue, newDecisionLevel, cancelUntil)
    virtual std::vector<Lit> forced() const = 0;
    virtual std::vector<Lit> forcedAt(int lvl) const = 0;

    // assert a literal
    virtual void enqueue(Lit l) = 0;

    // decision level manipulation
    virtual void newDecisionLevel() = 0;
    virtual void cancelUntil(int lvl) = 0;
    virtual int decisionLevel() const = 0;

    // --------------------------------------------------
    // assumptions. can be freely mainpulated. shrink is a convenience
    // to match minisat's interface and means shrink *by* num elements
    // (not *to* a number of elements)

    virtual std::vector<Lit>& assumptions() = 0;
    virtual void shrinkAssumptions(int num) = 0;

    // --------------------------------------------------
    // solve. Return l_True for SAT, l_False for UNSAT, l_Undef if
    // interrupted (e.g. for exceeding limits). This will undo
    // everything done with newDecisionLevel etc.
    virtual lbool solve() = 0;

    // accessing more info about the solve: model. Assumes the last
    // call was SAT, otherwise returns an empty model
    virtual model_t const& model() const = 0;

    // accessing more info about the solve: core. Assumes the call was
    // UNSAT, otherwise throws
    virtual LSet const& core() const = 0;

    // --------------------------------------------------
    // branching heuristic

    // decision var: if false, the solver does not have to branch on
    // this variable (but it is not forbidden to, either)
    virtual void setDecisionVar(int var, bool dvar) = 0;
    virtual bool isDecisionVar(int var) = 0;

    // priority: vars with higher priority are guaranteed to be
    // branched on before those with lower priority. Default is 0.
    virtual void setVarPriority(int var, int priority) = 0;
    virtual int varPriority(int var) = 0;

    // polarity: l_Undef lets the solver decide. With l_True
    // (l_False), it is always instantiated to false (true). Default
    // is l_Undef.
    virtual void setVarPolarity(int var, lbool pol) = 0;
    virtual lbool varPolarity(int var) = 0;

    // initial phase: as with polarity, but it can be overwritten by
    // the solver. 'phase' is interpreted same as polarity. Useful for
    // example to hand a (near-)solution to the solver.
    virtual void initializePhase(int var, lbool phase) = 0;

    // activity: get the activity of a variable. This is as heuristic
    // as the activity itself, so the solver can decide to give back
    // arbitrary values
    virtual double varActivity(int var) const = 0;

    // reset heuristics. For the ill-defined case where the solver has
    // ``stalled'', or as a crude way to get solution (or core)
    // diversity. May clear all learnt clauses, VSIDS scores and
    // phases, or nothing at all.
    virtual void resetHeuristics() = 0;

    // --------------------------------------------------
    // incremental interface

    // blocking variables must have the largest indices. Here, set the
    // id of the first one. This is an optimization only, so we can
    // use a non-bvar as an assumption (as long as it is frozen)
    virtual void setBvars(int firstbvar) = 0;

    // tell the solver that a bvar will never be used as an assumption
    // again. This is different from releaseVar(), as it does not
    // release the id 'bvar'
    virtual void setRedundantBvar(int bvar) = 0;

    // --------------------------------------------------
    // blocking variable elimination. A frozen variable cannot be
    // eliminated. Once a variable is eliminated, it cannot be frozen
    virtual void setFrozen(Var var, bool f) = 0;
    virtual bool isFrozen(Var var) const = 0;
    virtual bool isEliminated(Var var) const = 0;

    // --------------------------------------------------
    // Budget. These are cumulative, across all solves, not for each
    // individual solving episode. If an episode exceeds its budget,
    // we can increase the budget and try again. If any of the budgets
    // is <0, it means no limit.

    virtual int64_t conflictBudget() const = 0;
    virtual int64_t propagationBudget() const = 0;
    virtual void setConflictBudget(int64_t cb) = 0;
    virtual void setPropagationBudget(int64_t pb) = 0;
    virtual bool withinBudget() const = 0;

    // --------------------------------------------------
    // The derived class should use the singleton mcsmus::Control
    // instance for printing, timeouts etc. This sets the control
    // object to use this solver's printing functions
    virtual void initializeControl() const = 0;

    // --------------------------------------------------
    // Statistics
    virtual Statistics getStats() const = 0;
};

}

#endif
