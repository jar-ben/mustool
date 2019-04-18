/***********[lingeling-solver.hh]
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

#ifndef MCSMUS_LINGELING_SOLVER_H
#define MCSMUS_LINGELING_SOLVER_H

#include "mcsmus/basesolver.hh"
#include "mcsmus/system.hh"

extern "C" {
#include "lglib.h"
}

namespace mcsmus
{

class LingelingSolver : public BaseSolver
{
public:
    LingelingSolver();
    LingelingSolver(LingelingSolver const& rhs);
    virtual ~LingelingSolver();

    virtual std::unique_ptr<BaseSolver> clone() const override;

    // note our notion of forking is different from lingeling's. See
    // comments in basesolver.hh vs liglib.h
    virtual std::unique_ptr<BaseSolver> fork() const override;

    virtual std::unique_ptr<BaseSolver> clone_blank() const override;

    // hard to do in lingeling. so it's just clone
    virtual std::unique_ptr<BaseSolver> clone_input() const override;

    virtual void addClause(std::vector<Lit> const& ps) override
    {
        for (Lit l : ps)
            lgladd(lgl, l2int(l));
        lgladd(lgl, 0);
    }
    virtual void addClause(std::vector<Lit> const&& ps) override
    {
        for (Lit l : ps)
            lgladd(lgl, l2int(l));
        lgladd(lgl, 0);
    }

    virtual Var newVar() override;
    virtual void releaseVar(Lit l) override;

    virtual int nVars() const override { return lglnvars(lgl); }
    virtual int nUsedVars() const override { return lglmaxvar(lgl); }
    virtual int nClauses() const override { return lglnclauses(lgl); }

    virtual prop_result propagate() override
    {
        return lglinconsistent(lgl) ? prop_result::failed : prop_result::ok;
    }

    virtual lbool value(Var x) const override;
    virtual lbool value(Lit l) const override;
    virtual int levelOf(Var /*x*/) const override { return nVars(); }
    virtual Lit decisionAtLevel(int /*lvl*/) const override
    {
        return lit_Undef;
    }

    virtual std::vector<Lit> entailed() const override;
    virtual std::vector<Lit> forced() const override;
    virtual std::vector<Lit> forcedAt(int) const override;

    // unimplemented in lingeling
    virtual void enqueue(Lit /*l*/) override { assert(0); }
    virtual void newDecisionLevel() override { assert(0); }
    virtual void cancelUntil(int /*lvl*/) override {}
    virtual int decisionLevel() const override { return 0; }

    virtual std::vector<Lit>& assumptions() override
    {
        return next_assumptions;
    }
    virtual void shrinkAssumptions(int num) override
    {
        next_assumptions.resize(next_assumptions.size() - num);
    }

    virtual lbool solve() override;
    virtual model_t const& model() const override;
    virtual LSet const& core() const override;

    virtual void setDecisionVar(int /*var*/, bool /*dvar*/) override {}
    virtual bool isDecisionVar(int /*var*/) override { return true; }

    virtual void setVarPriority(int /*var*/, int /*priority*/) override
    {
        assert(0);
    }
    virtual int varPriority(int /*var*/) override { return 0; }

    virtual void setVarPolarity(int var, lbool pol) override
    {
        if (pol == l_Undef)
            lglresetphase(lgl, l2int(var));
        else
            lglsetphase(lgl, pol == l_True ? -l2int(var) : l2int(var));
    }
    virtual lbool varPolarity(int /*var*/) override
    {
        return l_Undef; // unknown unless we save it
    }

    virtual void initializePhase(int /*var*/, lbool /*phase*/) override {}

    virtual double varActivity(int /*var*/) const override { return 0.0; }

    virtual void resetHeuristics() override { lglflushcache(lgl); }

    virtual void setBvars(int /*firstbvar*/) override {}

    virtual void setRedundantBvar(int /*bvar*/) override {}

    virtual void setFrozen(Var var, bool f) override
    {
        if (f)
            lglfreeze(lgl, l2int(var));
        else
            lglmelt(lgl, l2int(var));
    }
    virtual bool isFrozen(Var var) const override
    {
        return lglfrozen(lgl, l2int(var));
    }
    virtual bool isEliminated(Var var) const override
    {
        return !lglusable(lgl, l2int(var));
    }

    virtual int64_t conflictBudget() const override { return cflbudget; }
    virtual int64_t propagationBudget() const override { return propbudget; }
    virtual void setConflictBudget(int64_t cb) override { cflbudget = cb; }
    virtual void setPropagationBudget(int64_t pb) override { propbudget = pb; }
    virtual bool withinBudget() const override;

    virtual void initializeControl() const override {}

    virtual Statistics getStats() const override;

private:
    LGL* lgl;
    std::vector<Lit> next_assumptions, last_assumptions;

    model_t model_cache;
    LSet core_cache;

    uint64_t solves{ 0 }, solutions{ 0 };

    int64_t cflbudget{-1}, propbudget{-1};

    int maxvar() const { return lglmaxvar(lgl); }
    int l2int(Lit l) const { return (var(l) + 1) * (1 - 2 * sign(l)); }
    int l2int(Var x) const { return x + 1; }
    Lit int2l(int x) const { return x > 0 ? mkLit(x - 1) : ~mkLit(x - 1); }
};

inline std::unique_ptr<BaseSolver> make_lingeling_solver()
{
    return std::unique_ptr<BaseSolver>(new LingelingSolver());
}
}

#endif
