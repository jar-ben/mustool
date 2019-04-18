#include "mcsmus/glucose-wrapper.hh"
#include "mcsmus/options.hh"

namespace mcsmus {

void GlucoseWrapper::printStats() const
{
    if( control->fn_printStats )
        control->fn_printStats();
}

void GlucoseWrapper::printHeader(double parsetime) const {
    if( control->fn_printHeader )
        control->fn_printHeader(parsetime);
}

void GlucoseWrapper::printFooter() const {
    if( control->fn_printFooter )
        control->fn_printFooter();
}

void GlucoseWrapper::printStatusLine() const {
    if( control->fn_printStatusLine )
        control->fn_printStatusLine();
    else
        Solver::printStatusLine();
}

void GlucoseWrapper::relocAll(Glucose::ClauseAllocator& to) {
    if( ovr.fn_relocAll )
        ovr.fn_relocAll(*this, ca, to);
    relocAll_(to);
}

void GlucoseWrapper::removeUseless(Glucose::vec<Glucose::CRef>& cs) {
    erase_if(cs, [&](Glucose::CRef cr) {
            Glucose::Clause &cl = ca[cr];
            if( cl.mark() == 1 ) return true;
            for(Glucose::Lit l : cl) {
                Glucose::Var v = var(l);
                Glucose::Lit lp = Glucose::mkLit(v);
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

void GlucoseWrapper::removeUselessVars() { }

void GlucoseWrapper::analyzeFinal(Glucose::Lit p, Glucose::LSet& out_conflict)
{
    if( ovr.fn_analyzeFinal )
        ovr.fn_analyzeFinal(*this, p, out_conflict);
    else
        Solver::analyzeFinal(p, out_conflict);
}

void GlucoseSimpWrapper::printStats() const
{
    if( control->fn_printStats )
        control->fn_printStats();
}

void GlucoseSimpWrapper::printHeader(double parsetime) const {
    if( control->fn_printHeader )
        control->fn_printHeader(parsetime);
}

void GlucoseSimpWrapper::printFooter() const {
    if( control->fn_printFooter )
        control->fn_printFooter();
}

void GlucoseSimpWrapper::printStatusLine() const {
    if( control->fn_printStatusLine )
        control->fn_printStatusLine();
    else
        Solver::printStatusLine();
}

void GlucoseSimpWrapper::relocAll(Glucose::ClauseAllocator& to) {
    if( ovr.fn_relocAll )
        ovr.fn_relocAll(*this, ca, to);
    SimpSolver::relocAll(to);
}

void GlucoseSimpWrapper::removeUseless(Glucose::vec<Glucose::CRef>&) {
    //empty
}

void GlucoseSimpWrapper::removeUselessVars() {
    //empty
}

void GlucoseSimpWrapper::analyzeFinal(Glucose::Lit p, Glucose::LSet& out_conflict)
{
    if( ovr.fn_analyzeFinal )
        ovr.fn_analyzeFinal(*this, p, out_conflict);
    else
        Solver::analyzeFinal(p, out_conflict);
}

}
