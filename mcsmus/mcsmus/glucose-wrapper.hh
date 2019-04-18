#ifndef MCSMUS_GLUCOSE_WRAPPER_H
#define MCSMUS_GLUCOSE_WRAPPER_H

#include "glucose/core/SolverTypes.h"
#include "glucose/core/Solver.h"
#include "glucose/simp/SimpSolver.h"
#include "minisat/utils/mcsmus_System.h"
#include "mcsmus/basesolver.hh"
#include <functional>
#include <vector>

namespace mcsmus {

    template<class SolverType>
    class glucose_overrides
    {
    public:
        std::function<void(SolverType&, Glucose::ClauseAllocator&,
                           Glucose::ClauseAllocator&)> fn_relocAll;
        std::function<void(SolverType&, Lit, LSet&)> fn_analyzeFinal;
    };

    class GlucoseWrapper : public Glucose::Solver
    {
    public:
        GlucoseWrapper() = default;
        GlucoseWrapper(GlucoseWrapper const& other) : Solver(other),
                                                      ovr(other.ovr)
        {
            update_members_with_references();
        }

        Glucose::Clause& getClause(Glucose::CRef cr) { return ca[cr]; }
        Glucose::Clause const& getClause(Glucose::CRef cr) const { return ca[cr]; }

        // promote to public
        using Solver::propagate;
        using Solver::removeClause;
        using Solver::enqueue;
        using Solver::newDecisionLevel;
        using Solver::cancelUntil;
        using Solver::decisionLevel;
        using Solver::assumptions;

        using Solver::level;
        using Solver::reason;

        using Solver::propagation_budget;
        using Solver::conflict_budget;
        using Solver::withinBudget;

        using Solver::solve_;
        using Solver::solve_nocancel_;

        using Solver::decision;
        using Solver::activity;
        using Solver::priority;
        using Solver::user_pol;
        using Solver::polarity;
        using Solver::order_heap;
        using Solver::varBumpActivity;
        using Solver::rebuildOrderHeap;

        using Solver::trail;
        using Solver::trail_lim;

        /* overrides */
        virtual void printStats() const override;
        virtual void printHeader(double parsetime) const override;
        virtual void printFooter() const override;
        virtual void printStatusLine() const override;

        virtual void relocAll(Glucose::ClauseAllocator& to) override;
        virtual void analyzeFinal(Glucose::Lit p, Glucose::LSet& out_conflict) override;

        glucose_overrides<GlucoseWrapper> ovr;

        // handling clauses moved to redundants
        virtual void removeUseless(Glucose::vec<Glucose::CRef>& cs) override;
        virtual void removeUselessVars() override;

        void setBvars(int firstbvar)
        {
            setIncrementalMode();
            initNbInitialVars(firstbvar-1);
        }
        void setRedundantBvar(int bvar)
        {
            redundant.resize(std::max(size_t(2*(bvar+1)), redundant.size()));
            redundant[toInt(mkLit(bvar))] = 1;
        }
        bool isEliminated(Glucose::Var /*v*/) const { return false; }

        // support for VE-aware forking, so nops here. elimclauses is
        // not used anywhere
        void setFrozen(Glucose::Var, bool) {}
        bool isFrozen(Glucose::Var) const { return false; }
        void markEliminatedUnsafe(Glucose::Var) { }
        Glucose::vec<uint32_t> elimclauses;

        // preprocessing time
        double preprocessingTime() const { return 0.0; }

        // allow us to only do extendModel on demand
        Glucose::lbool solve_nocancel_noextend_() { return solve_nocancel_(); }
        void extendModel() {}

        // for forking (including debugging)
        using Solver::clauses;
        using Solver::learnts;
        using Solver::ca;
        using Solver::attachClause;
        bool isRemoved(Glucose::CRef cr) const { return ca[cr].mark() == 1; }
    private:
        std::vector<uint8_t> redundant;
    };

    class GlucoseSimpWrapper : public Glucose::SimpSolver {
    public:
        GlucoseSimpWrapper() = default;
        GlucoseSimpWrapper(GlucoseSimpWrapper const& other) :
            SimpSolver(other), ovr(other.ovr)
        {
            update_members_with_references();
        }

        Glucose::Clause& getClause(Glucose::CRef cr) { return ca[cr]; }
        Glucose::Clause const& getClause(Glucose::CRef cr) const { return ca[cr]; }

        // promote to public
        using SimpSolver::propagate;
        using SimpSolver::removeClause;
        using SimpSolver::enqueue;
        using SimpSolver::newDecisionLevel;
        using SimpSolver::cancelUntil;
        using SimpSolver::decisionLevel;
        using SimpSolver::assumptions;

        using SimpSolver::level;
        using SimpSolver::reason;

        using SimpSolver::propagation_budget;
        using SimpSolver::conflict_budget;
        using SimpSolver::withinBudget;

        using SimpSolver::solve_;
        using SimpSolver::solve_nocancel_;

        using SimpSolver::decision;
        using SimpSolver::activity;
        using SimpSolver::priority;
        using SimpSolver::user_pol;
        using SimpSolver::polarity;
        using SimpSolver::order_heap;
        using SimpSolver::varBumpActivity;
        using SimpSolver::rebuildOrderHeap;

        using SimpSolver::trail;
        using SimpSolver::trail_lim;

        /* overrides */
        virtual void printStats() const override;
        virtual void printHeader(double parsetime) const override;
        virtual void printFooter() const override;
        virtual void printStatusLine() const override;

        virtual void relocAll(Glucose::ClauseAllocator& to) override;

        virtual void analyzeFinal(Glucose::Lit p, Glucose::LSet& out_conflict) override;

        glucose_overrides<GlucoseSimpWrapper> ovr;

        // handling clauses moved to redundants
        void setRedundantArray(std::vector<uint8_t> const* red);
        virtual void removeUseless(Glucose::vec<Glucose::CRef>& cs) override;
        virtual void removeUselessVars() override;

        // functionality missing from minisat
        void setBvars(int firstbvar)
        {
            setIncrementalMode();
            initNbInitialVars(firstbvar-1);
            for(int i = firstbvar; i != nVars(); ++i)
                setFrozen(i, true);
        }
        void setRedundantBvar(int bvar)
        {
            setFrozen(bvar, false);
            updateElimHeap(bvar);
            setDecisionVar(bvar, false);
        }

        // support for VE-aware forking
        bool isFrozen(Glucose::Var v) const { return frozen[v]; }
        void markEliminatedUnsafe(Glucose::Var v) { eliminated[v] = true; }
        using SimpSolver::elimclauses;

        // preprocessing time
        double preprocessingTime() const { return elim_time; }

        // allow us to only do extendModel on demand
        using SimpSolver::extendModel;

        // for forking (including debugging)
        using SimpSolver::clauses;
        using SimpSolver::learnts;
        using SimpSolver::ca;
        using SimpSolver::attachClause;
        bool isRemoved(Glucose::CRef cr) const { return ca[cr].mark() == 1; }
    };

    class GlucoseSimpSolver : public BaseSolver {
        GlucoseSimpWrapper solver;
    public:
        void resetHeuristics() override;
    };
}


#endif // __GLUCOSE_WRAPPER_H
