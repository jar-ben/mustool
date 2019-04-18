/***********[minisat-wrapper-common.hh]
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

#ifndef MCSMUS_MINISAT_WRAPPER_COMMON_HH
#define MCSMUS_MINISAT_WRAPPER_COMMON_HH

#include "mcsmus/basesolver.hh"
#include "mcsmus/minisat-wrapper.hh"
#include "mcsmus/glucose-wrapper.hh"
#include "mcsmus/system.hh"
#include "minisat/mtl/mcsmus_Vec.h"

#include <algorithm>

namespace mcsmus
{

/* Implement the BaseSolver interface for all minisat variants:
   [Minisat|Glucose]::[Simp]Solver */
template <typename Solver> class MiniSolver : public BaseSolver
{
public:
    MiniSolver()
    {
        solver.ovr.fn_analyzeFinal = [this](Solver& s, Lit p,
            Minisat::LSet& out_cnfl) { analyzeFinal(s, p, out_cnfl); };
    }
    MiniSolver(MiniSolver const& other)
        : solver(other.solver),
          last_model(other.last_model),
          last_result(other.last_result),
          next_assumptions(other.next_assumptions)
    {
        solver.ovr.fn_analyzeFinal = [this](Solver& s, Lit p,
            Minisat::LSet& out_cnfl) { analyzeFinal(s, p, out_cnfl); };
    }
    virtual ~MiniSolver() {}

    /* This *does not* override anything in BaseSolver, but instead
       gives direct access to the underlying minisat|glucose
       solver. */
    Solver& getSolver() { return solver; }

    virtual std::unique_ptr<BaseSolver> clone() const override
    {
        return std::unique_ptr<BaseSolver>(new MiniSolver<Solver>(*this));
    }

    virtual std::unique_ptr<BaseSolver> fork() const override
    {
        auto* other = new MiniSolver<Solver>();
        auto to_p = std::unique_ptr<BaseSolver>(other);
        auto& to = other->solver;
        other->next_assumptions = next_assumptions;

        for (Lit l : next_assumptions) {
            auto v = var(l);
            while (to.nVars() <= v)
                to.newVar();
            to.addClause(l);
        }

        int end = solver.trail_lim.size() > 0 ? solver.trail_lim[0]
                                              : solver.trail.size();
        for (int q = 0; q != end; ++q) {
            auto v = var(solver.trail[q]);
            while (to.nVars() <= v)
                to.newVar();
            to.addClause({ solver.trail[q] });
        }

        for (auto cr : solver.clauses) {
            Minisat::vec<Lit> ps;
            auto& c = solver.getClause(cr);
            if (solver.isRemoved(cr))
                continue;
            bool sat = false;
            int maxcvar{ var_Undef };
            for (Lit l : c) {
                if (solver.value(l) == l_False)
                    continue;
                if (solver.value(l) == l_True)
                    sat = true;
                assert(!solver.isEliminated(var(l)));
                if (var(l) < to.nVars()) {
                    if (to.value(l) == l_False)
                        continue;
                    if (to.value(l) == l_True)
                        sat = true;
                }
                auto v = var(l);
                maxcvar = std::max(maxcvar, v);
                ps.push(l);
            }
            if (sat)
                continue;
            while (to.nVars() <= maxcvar)
                to.newVar();
            to.addClause(ps);
        }

        for (auto cr : solver.learnts) {
            Minisat::vec<Lit> ps;
            auto& c = solver.getClause(cr);
            if (solver.isRemoved(cr))
                continue;
            bool sat = false;
            int maxcvar{ var_Undef };
            for (Lit l : c) {
                if (solver.value(l) == l_False)
                    continue;
                if (solver.value(l) == l_True)
                    sat = true;
                assert(!solver.isEliminated(var(l)));
                if (var(l) < to.nVars()) {
                    if (to.value(l) == l_False)
                        continue;
                    if (to.value(l) == l_True)
                        sat = true;
                }
                auto v = var(l);
                maxcvar = std::max(maxcvar, v);
                ps.push(l);
            }
            if (sat)
                continue;
            while (to.nVars() <= maxcvar)
                to.newVar();
            if (ps.size() == 0) {
                to.addClause(ps); // just to mark it as failed
                return to_p;
            } else if (ps.size() == 1) {
                to.enqueue(ps[0]);
                to.propagate();
            } else {
                auto tocr = to.ca.alloc(ps, true);
                to.learnts.push(tocr);
                to.attachClause(tocr);
            }
        }

        // copy the frozen/eliminated status of each variable and get a
        // copy of the eliminated clauses
        for (Var v = 0; v != to.nVars(); ++v) {
            to.setFrozen(v, solver.isFrozen(v));
            if (solver.isEliminated(v))
                to.markEliminatedUnsafe(v);
        }
        solver.elimclauses.copyTo(to.elimclauses);

        return to_p;
    }

    virtual std::unique_ptr<BaseSolver> clone_blank() const override
    {
        return std::unique_ptr<BaseSolver>(new MiniSolver<Solver>());
    }

    virtual std::unique_ptr<BaseSolver> clone_input() const override
    {
        auto* other = new MiniSolver<Solver>();
        auto to_p = std::unique_ptr<BaseSolver>(other);
        auto& to = other->solver;
        for (Lit l : solver.assumptions) {
            auto v = var(l);
            while (to.nVars() <= v)
                to.newVar();
            to.addClause(l);
        }

        int end = solver.trail_lim.size() > 0 ? solver.trail_lim[0]
                                              : solver.trail.size();
        for (int q = 0; q != end; ++q) {
            auto v = var(solver.trail[q]);
            while (to.nVars() <= v)
                to.newVar();
            to.addClause({ solver.trail[q] });
        }

        for (auto cr : solver.clauses) {
            Minisat::vec<Lit> ps;
            auto& c = solver.getClause(cr);
            if (solver.isRemoved(cr))
                continue;
            bool sat = false;
            int maxcvar{ var_Undef };
            for (Lit l : c) {
                if (solver.value(l) == l_False)
                    continue;
                if (solver.value(l) == l_True)
                    sat = true;
                // XXX: it can actually happen some times that the
                // following assert fails. I believe it is a bug but I
                // have not tried to find it. Instead, we treat the
                // clause as satisfied, which should be safe
                if( solver.isEliminated(var(l)) ) {
                    sat = true;
                    break;
                }
                assert(!solver.isEliminated(var(l)));
                if (var(l) < to.nVars()) {
                    if (to.value(l) == l_False)
                        continue;
                    if (to.value(l) == l_True)
                        sat = true;
                }
                auto v = var(l);
                maxcvar = std::max(maxcvar, v);
                ps.push(l);
            }
            if (sat)
                continue;
            while (to.nVars() <= maxcvar)
                to.newVar();
            to.addClause(ps);
        }

        // copy the frozen/eliminated status of each variable and get a
        // copy of the eliminated clauses
        for (Var v = 0; v != to.nVars(); ++v) {
            to.setFrozen(v, solver.isFrozen(v));
            if (solver.isEliminated(v))
                to.markEliminatedUnsafe(v);
        }
        solver.elimclauses.copyTo(to.elimclauses);

        solver.assumptions.copyTo(to.assumptions);

        return to_p;
    }

    virtual void addClause(std::vector<Lit> const& ps) override
    {
        addClause_internal(begin(ps), end(ps));
    }
    virtual void addClause(std::vector<Lit> const&& ps) override
    {
        addClause_internal(begin(ps), end(ps));
    }

    virtual Var newVar() override { return solver.newVar(); }
    virtual void releaseVar(Lit l) override { solver.releaseVar(l); }

    virtual int nVars() const override { return solver.nVars(); }
    virtual int nUsedVars() const override { return solver.nVars(); }
    virtual int nClauses() const override { return solver.nClauses(); }

    virtual prop_result propagate() override
    {
        auto cr = solver.propagate();
        if (cr == Minisat::CRef_Undef)
            return prop_result::ok;
        else
            return prop_result::failed;
    }

    virtual lbool value(Var x) const override { return solver.value(x); }
    virtual lbool value(Lit l) const override { return solver.value(l); }
    virtual int levelOf(Var x) const override { return solver.level(x); }
    virtual Lit decisionAtLevel(int lvl) const override
    {
#ifdef _GLIBCXX_DEBUG
        assert(lvl>=1 && lvl <= solver.decisionLevel());
#endif
        auto pos = solver.trail_lim[lvl - 1];
        if (lvl < solver.decisionLevel() && pos == solver.trail_lim[lvl])
            return lit_Undef;
        if (lvl == solver.decisionLevel() && pos == solver.trail.size())
            return lit_Undef;
        return solver.trail[pos];
    }
    virtual std::vector<Lit> entailed() const override
    {
        return std::vector<Lit>(begin(solver.trail),
            (solver.trail_lim.size() == 0
             ? end(solver.trail)
             : begin(solver.trail) + solver.trail_lim[0]));
    }

    virtual std::vector<Lit> forced() const override
    {
        return std::vector<Lit>(begin(solver.trail), end(solver.trail));
    }

    virtual std::vector<Lit> forcedAt(int lvl) const override
    {
#ifdef _GLIBCXX_DEBUG
        assert(lvl>=0 && lvl <= solver.decisionLevel());
#endif
        int b{0}, e{solver.trail.size()};
        if( lvl != 0 )
            b = solver.trail_lim[lvl-1];
        if( lvl < solver.decisionLevel() )
            e = solver.trail_lim[lvl];
        return std::vector<Lit>(begin(solver.trail)+b, begin(solver.trail)+e);
    }

    virtual void enqueue(Lit l) override { solver.enqueue(l); }

    virtual void newDecisionLevel() override { solver.newDecisionLevel(); }
    virtual void cancelUntil(int lvl) override { solver.cancelUntil(lvl); }
    virtual int decisionLevel() const override
    {
        return solver.decisionLevel();
    }

    virtual std::vector<Lit>& assumptions() override
    {
        return next_assumptions;
    }
    virtual void shrinkAssumptions(int num) override
    {
        next_assumptions.resize(next_assumptions.size() - num);
    }

    virtual lbool solve() override
    {
        auto startTime = cpuTime();
        // compute maximum safe bt level
        unsigned i;
        for (i = 0; i != static_cast<unsigned>(solver.decisionLevel())
                 && i != next_assumptions.size()
                 && (next_assumptions[i] == decisionAtLevel(i + 1)
                        || (decisionAtLevel(i + 1) == lit_Undef
                               && solver.value(next_assumptions[i]) == l_True
                               && static_cast<unsigned>(solver.level(
                                      var(next_assumptions[i]))) < i + 1));
             ++i)
            ;
        solver.cancelUntil(i);
        // sync assumptions
        solver.assumptions.clear();
        for (Lit l : next_assumptions)
            solver.assumptions.push(l);
        // actual solve
        lbool rv = solver.solve_nocancel_();
        if (rv == l_True)
            last_model.model
                = std::vector<lbool>(begin(solver.model), end(solver.model));
        else
            last_model.clear();
        time += cpuTime() - startTime;
        return rv;
    }

    virtual model_t const& model() const override { return last_model; }

    virtual LSet const& core() const override { return solver.conflict; }

    virtual void setDecisionVar(int var, bool dvar) override
    {
        solver.setDecisionVar(var, dvar);
    }
    virtual bool isDecisionVar(int var) override
    {
        return solver.decision[var];
    }

    virtual void setVarPriority(int var, int priority) override
    {
        solver.priority[var] = priority;
        solver.order_heap.update(var);
    }
    virtual int varPriority(int var) override { return solver.priority[var]; }

    virtual void setVarPolarity(int var, lbool pol) override
    {
        solver.user_pol[var] = pol;
    }
    virtual lbool varPolarity(int var) override { return solver.user_pol[var]; }

    virtual void initializePhase(int var, lbool phase) override
    {
        if (phase == l_Undef)
            return;
        bool b = phase == l_True;
        solver.polarity[var] = b;
    }

    virtual double varActivity(int var) const override
    {
        return solver.activity[var];
    }

    virtual void resetHeuristics() override
    {
        using namespace Glucose;
        solver.cancelUntil(0);
        for (CRef cr : solver.learnts)
            solver.removeClause(cr);
        solver.learnts.clear();
        for (Var v = 0; v != solver.nVars(); ++v) {
            solver.activity[v] = 0;
            solver.polarity[v] = true;
        }
        solver.simplify(true);
    }

    virtual void setBvars(int firstbvar) override
    {
        solver.setBvars(firstbvar);
    }

    virtual void setRedundantBvar(int bvar) override
    {
        solver.setRedundantBvar(bvar);
    }

    virtual void setFrozen(Var var, bool f) override
    {
        solver.setFrozen(var, f);
    }
    virtual bool isFrozen(Var var) const override
    {
        return solver.isFrozen(var);
    }
    virtual bool isEliminated(Var var) const override
    {
        return solver.isEliminated(var);
    }

    virtual int64_t conflictBudget() const override
    {
        return solver.conflict_budget;
    }
    virtual int64_t propagationBudget() const override
    {
        return solver.propagation_budget;
    }
    virtual void setConflictBudget(int64_t cb) override
    {
        solver.conflict_budget = cb;
    }
    virtual void setPropagationBudget(int64_t pb) override
    {
        solver.propagation_budget = pb;
    }
    virtual bool withinBudget() const override { return solver.withinBudget(); }

    virtual void initializeControl() const override
    {
        Control* control = getGlobalControl();
        control->fn_printStats = [&]() { solver.printStats(); };
        control->fn_printHeader = [&](double t) { solver.printHeader(t); };
        control->fn_printFooter = [&]() { solver.printFooter(); };
        control->fn_printStatusLine = [&]() { solver.printStatusLine(); };
    }

    virtual Statistics getStats() const override
    {
        return { solver.solves, solver.starts, solver.dec_vars,
            (uint64_t)solver.nClauses(), solver.decisions, solver.rnd_decisions,
            solver.propagations, solver.conflicts, (uint64_t)solver.nLearnts(),
            solver.clauses_literals, solver.learnts_literals,
            solver.max_literals, solver.tot_literals, solver.solutions, time,
            solver.preprocessingTime() };
    }

protected:
    void analyzeFinal(Solver& afsolver, Lit p, Minisat::LSet& out_conflict)
    {
        // Changes from original: stop resolving backwards when we hit an
        // assumption literal, bump activity of literals we meet

        out_conflict.clear();
        out_conflict.insert(p);

        if (afsolver.decisionLevel() == 0)
            return;

        std::vector<uint8_t> seen(afsolver.nVars() + 1);

        seen[var(p)] = 1;

        Minisat::LSet assumps;
        for (int i = 0; i < afsolver.assumptions.size(); i++) {
            assumps.insert(afsolver.assumptions[i]);
        }

        for (int i = afsolver.trail.size() - 1; i >= afsolver.trail_lim[0];
             i--) {
            Var x = var(afsolver.trail[i]);
            if (seen[x]) {
                if (afsolver.reason(x) == Minisat::CRef_Undef
                    || assumps.has(afsolver.trail[i])) {
                    // note assumptions has ~p in it and ~p was
                    // forced. So x == p won't trigger.
                    assert(afsolver.level(x) > 0);
                    assert(x != var(p));
                    out_conflict.insert(~afsolver.trail[i]);
                    afsolver.varBumpActivity(var(afsolver.trail[i]));
                } else {
                    auto& c = afsolver.getClause(afsolver.reason(x));
                    for (int j = 0; j < c.size(); j++)
                        if (var(c[j]) != x && afsolver.level(var(c[j])) > 0) {
                            seen[var(c[j])] = 1;
                            afsolver.varBumpActivity(var(c[j]));
                        }
                }
                seen[x] = 0;
            }
        }

        seen[var(p)] = 0;
    }

    template <typename Iterator, typename EndIter>
    void addClause_internal(Iterator b, EndIter e)
    {
        auto& ps = addClause_tmp;
        ps.clear();
        for (; b != e; ++b) {
            ps.push(*b);
            while (solver.nVars() <= var(*b))
                solver.newVar();
        }

        if (ps.size() < 2) {
            solver.cancelUntil(0);
            solver.addClause_(ps);
            return;
        }

        // not at dlvl 0, we may have to backjump
        auto pivot = std::partition(begin(ps), end(ps),
            [&](Lit p) { return solver.value(p) != l_False; });
        auto ppos = std::distance(begin(ps), pivot);
        if (ppos < 2) {
            std::partial_sort(
                pivot, begin(ps) + 2, end(ps), [&](Lit l1, Lit l2) {
                    return solver.level(var(l1)) > solver.level(var(l2));
                });
            int btlvl = solver.level(var(ps[1]));
            if (ppos == 0)
                btlvl = std::max(
                    0, std::min(btlvl, solver.level(var(ps[0])) - 1));
            solver.cancelUntil(btlvl);
        }
        if (solver.decisionLevel() == 0) {
            solver.addClause_(ps);
            return;
        } else {
            auto cl = solver.addClauseR_unsafe_(ps);
            if (solver.value(ps[1]) == l_False)
                solver.enqueue(ps[0], cl.second);
        }
    }

protected:
    Solver solver;

    // information on last solve
    model_t last_model;
    lbool last_result;

    // we let the user manipulate this and sync with
    // solver.assumptions before each solve, performing minimal
    // backtracking
    std::vector<Lit> next_assumptions;

    // amazingly minisat does not explicitly keep this
    double time{ 0.0 };

    // class member to avoid allocations in addClause_internal
    Minisat::vec<Lit> addClause_tmp;
};

inline std::unique_ptr<MiniSolver<MinisatWrapper>> make_minisat_solver()
{
    // return make_unique<MiniSolver<MinisatWrapper>>();
    return std::unique_ptr<MiniSolver<MinisatWrapper> >(
        new MiniSolver<MinisatWrapper>());
}

inline
std::unique_ptr<MiniSolver<MinisatSimpWrapper>> make_minisat_simp_solver()
{
    // return make_unique<MiniSolver<MinisatSimpWrapper>>();
    return std::unique_ptr<MiniSolver<MinisatSimpWrapper> >(
        new MiniSolver<MinisatSimpWrapper>());
}

inline std::unique_ptr<MiniSolver<GlucoseWrapper>> make_glucose_solver()
{
    // return make_unique<MiniSolver<GlucoseWrapper>>();
    return std::unique_ptr<MiniSolver<GlucoseWrapper> >(
        new MiniSolver<GlucoseWrapper>());
}

inline
std::unique_ptr<MiniSolver<GlucoseSimpWrapper>> make_glucose_simp_solver()
{
    // return make_unique<MiniSolver<GlucoseSimpWrapper>>();
    return std::unique_ptr<MiniSolver<GlucoseSimpWrapper> >(
        new MiniSolver<GlucoseSimpWrapper>());
}
}

#endif
