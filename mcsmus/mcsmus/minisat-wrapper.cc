#include <vector>

#include "mcsmus/minisat-wrapper.hh"
#include "mcsmus/options.hh"

namespace mcsmus {

void MinisatWrapper::printStats() const
{
    if( control->fn_printStats )
        control->fn_printStats();
}

void MinisatWrapper::printHeader(double parsetime) const {
    if( control->fn_printHeader )
        control->fn_printHeader(parsetime);
}

void MinisatWrapper::printFooter() const {
    if( control->fn_printFooter )
        control->fn_printFooter();
}

void MinisatWrapper::printStatusLine() const {
    if( control->fn_printStatusLine )
        control->fn_printStatusLine();
    else
        Solver::printStatusLine();
}

void MinisatWrapper::relocAll(Minisat::ClauseAllocator& to) {
    if( ovr.fn_relocAll )
        ovr.fn_relocAll(*this, ca, to);
    relocAll_(to);
}

void MinisatWrapper::removeUseless(Minisat::vec<Minisat::CRef>& cs) {
    erase_if(cs, [&](Minisat::CRef cr) {
            Minisat::Clause &cl = getClause(cr);
            if( cl.mark() == 1 ) return true;
            for(Minisat::Lit l : cl) {
                Minisat::Var v = var(l);
                Minisat::Lit lp = Minisat::mkLit(v);
                if( unsigned(toInt(lp)) < redundant.size()
                    && redundant[toInt(lp)] ) {
                    assert(cl.mark() == 0);
                    removeClause(cr);
                    return true;
                }
            }
            return false;
        });
}

void MinisatWrapper::removeUselessVars() { }

void MinisatWrapper::analyzeFinal(Minisat::Lit p, Minisat::LSet& out_conflict)
{
    if( ovr.fn_analyzeFinal )
        ovr.fn_analyzeFinal(*this, p, out_conflict);
    else
        Solver::analyzeFinal(p, out_conflict);
}


void MinisatSimpWrapper::printStats() const
{
    if( control->fn_printStats )
        control->fn_printStats();
}

void MinisatSimpWrapper::printHeader(double parsetime) const {
    if( control->fn_printHeader )
        control->fn_printHeader(parsetime);
}

void MinisatSimpWrapper::printFooter() const {
    if( control->fn_printFooter )
        control->fn_printFooter();
}

void MinisatSimpWrapper::printStatusLine() const {
    if( control->fn_printStatusLine )
        control->fn_printStatusLine();
    else
        SimpSolver::printStatusLine();
}

void MinisatSimpWrapper::relocAll(Minisat::ClauseAllocator& to) {
    if( ovr.fn_relocAll )
        ovr.fn_relocAll(*this, ca, to);
    SimpSolver::relocAll(to); // TODO: GAH!
}

void MinisatSimpWrapper::removeUseless(Minisat::vec<Minisat::CRef>&) {
    // nothing here, varelim should take care of it
}

void MinisatSimpWrapper::removeUselessVars() {
    // nothing
}

void MinisatSimpWrapper::analyzeFinal(Minisat::Lit p, Minisat::LSet& out_conflict)
{
    if( ovr.fn_analyzeFinal )
        ovr.fn_analyzeFinal(*this, p, out_conflict);
    else
        SimpSolver::analyzeFinal(p, out_conflict);
}


}
