#ifndef MCSMUS_MINISAT_WRAPPER_H
#define MCSMUS_MINISAT_WRAPPER_H

#include "minisat/core/mcsmus_SolverTypes.h"
#include "minisat/core/mcsmus_Solver.h"
#include "minisat/simp/mcsmus_SimpSolver.h"
#include "minisat/utils/mcsmus_System.h"
#include <functional>

namespace mcsmus {

    template<typename SolverType>
    struct minisat_overrides {
        std::function<void(SolverType&,
                           Minisat::ClauseAllocator&, Minisat::ClauseAllocator&)>
             fn_relocAll;
        std::function<void(SolverType&, Minisat::Lit, Minisat::LSet&)> fn_analyzeFinal;
    };

    class MinisatWrapper : public Minisat::Solver {
    public:
        MinisatWrapper() = default;
        MinisatWrapper(MinisatWrapper const& other) : Solver(other),
                                                      ovr(other.ovr)
        {
            update_members_with_references();
        }

        Minisat::Clause& getClause(Minisat::CRef cr) { return ca[cr]; }
        Minisat::Clause const& getClause(Minisat::CRef cr) const { return ca[cr]; }

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

        using Solver::trail;
        using Solver::trail_lim;

        /* overrides */
        virtual void printStats() const override;
        virtual void printHeader(double parsetime) const override;
        virtual void printFooter() const override;
        virtual void printStatusLine() const override;

        virtual void relocAll(Minisat::ClauseAllocator& to) override;
        virtual void analyzeFinal(Minisat::Lit p, Minisat::LSet& out_conflict) override;

        /* the above overrides simply call these functions */
        minisat_overrides<MinisatWrapper> ovr;

        // handling clauses moved to redundants
        virtual void removeUseless(Minisat::vec<Minisat::CRef>& cs) override;
        virtual void removeUselessVars() override;

        // functionality missing from minisat
        void setBvars(int /*firstbvar*/) {}
        void setRedundantBvar(int bvar)
        {
            redundant.resize(std::max(size_t(2*(bvar+1)), redundant.size()));
            redundant[toInt(Minisat::mkLit(bvar))] = 1;
        }
        bool isEliminated(Minisat::Var /*v*/) const { return false; }

        // support for VE-aware forking, so nops here. elimclauses is
        // not used anywhere
        void setFrozen(Minisat::Var, bool) {}
        bool isFrozen(Minisat::Var) const { return false; }
        void markEliminatedUnsafe(Minisat::Var) { }
        Minisat::vec<uint32_t> elimclauses;

        // preprocessing time
        double preprocessingTime() const { return 0.0; }

        // allow us to only do extendModel on demand
        Minisat::lbool solve_nocancel_noextend_() { return solve_nocancel_(); }
        void extendModel() {}

        // for forking (including debugging)
        using Solver::clauses;
        using Solver::learnts;
        using Solver::ca;
        using Solver::attachClause;
        using Solver::isRemoved;
    private:
        std::vector<uint8_t> redundant;
    };

    class MinisatSimpWrapper : public Minisat::SimpSolver {
    public:
        MinisatSimpWrapper() = default;
        MinisatSimpWrapper(MinisatSimpWrapper const& other) :
            SimpSolver(other),
            ovr(other.ovr)
        {
            update_members_with_references();
        }

        Minisat::Clause& getClause(Minisat::CRef cr) { return ca[cr]; }
        Minisat::Clause const& getClause(Minisat::CRef cr) const { return ca[cr]; }

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

        using SimpSolver::trail;
        using SimpSolver::trail_lim;

        /* overrides */
        virtual void printStats() const override;
        virtual void printHeader(double parsetime) const override;
        virtual void printFooter() const override;
        virtual void printStatusLine() const override;

        virtual void relocAll(Minisat::ClauseAllocator& to) override;

        virtual void analyzeFinal(Minisat::Lit p, Minisat::LSet& out_conflict) override;

        /* the above overrides simply call these functions */
        minisat_overrides<MinisatSimpWrapper> ovr;

        // handling clauses moved to redundants
        virtual void removeUseless(Minisat::vec<Minisat::CRef>& cs) override;
        virtual void removeUselessVars() override;

        // functionality missing from minisat
        void setBvars(int firstbvar)
        {
            for(int i = firstbvar; i != nVars(); ++i)
                setFrozen(i, true);
        }
        void setRedundantBvar(int bvar)
        {
            setFrozen(bvar, false);
            updateElimHeap(bvar);
        }

        // support for VE-aware forking
        bool isFrozen(Minisat::Var v) const { return frozen[v]; }
        void markEliminatedUnsafe(Minisat::Var v) { eliminated[v] = true; }
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
        using SimpSolver::isRemoved;
    };

}

#endif // __MINISAT_WRAPPER_H
