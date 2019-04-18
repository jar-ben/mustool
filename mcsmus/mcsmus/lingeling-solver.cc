/***********[lingeling-solver.cc]
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

#ifdef USE_LINGELING
#include "mcsmus/lingeling-solver.hh"

using namespace mcsmus;
using namespace std;

int lingeling_control_helper(void* v)
{
    Control* c = getGlobalControl();
    c->checkPrintStatusLine();
    if (c->asynch_interrupt)
        return true;
    LingelingSolver* l = reinterpret_cast<LingelingSolver*>(v);
    return !l->withinBudget();
}

LingelingSolver::LingelingSolver()
{
    lgl = lglinit();
    lglseterm(lgl, lingeling_control_helper, reinterpret_cast<void*>(this));
}

LingelingSolver::LingelingSolver(LingelingSolver const& rhs)
    : next_assumptions(rhs.next_assumptions)
{
    lgl = lglclone(rhs.lgl);
}

LingelingSolver::~LingelingSolver() { lglrelease(lgl); }

unique_ptr<BaseSolver> LingelingSolver::clone() const
{
    return unique_ptr<BaseSolver>(new LingelingSolver(*this));
}

unique_ptr<BaseSolver> LingelingSolver::fork() const
{
    auto rv = clone();
    for (Lit l : next_assumptions)
        rv->addClause({ l });
    return rv;
}

unique_ptr<BaseSolver> LingelingSolver::clone_blank() const
{
    return unique_ptr<BaseSolver>(new LingelingSolver());
}

unique_ptr<BaseSolver> LingelingSolver::clone_input() const { return clone(); }

Var LingelingSolver::newVar() { return lglincvar(lgl) - 1; }

void LingelingSolver::releaseVar(Lit l)
{
    lgladd(lgl, l2int(l));
    lgladd(lgl, 0);
}

lbool LingelingSolver::solve()
{
    ++solves;
    for (Lit l : next_assumptions)
        lglassume(lgl, l2int(l));
    last_assumptions = next_assumptions;
    int ret = lglsat(lgl);
    model_cache.clear();
    core_cache.clear();
    switch (ret) {
    case LGL_UNKNOWN:
        return l_Undef;
    case LGL_SATISFIABLE:
        ++solutions;
        for (int i = 0, end = maxvar(); i != end; ++i)
            model_cache.model.push_back(
                lglderef(lgl, l2int(i)) > 0 ? l_True : l_False);
        return l_True;
    case LGL_UNSATISFIABLE:
        for (Lit l : last_assumptions) {
            int ic = lglfailed(lgl, l2int(l));
            if (!ic)
                continue;
            if (ic > 0)
                core_cache.insert(~l);
            else
                core_cache.insert(l);
        }
        return l_False;
    default:
        assert(0);
    }
}

lbool LingelingSolver::value(Var x) const
{
    int f = lglfixed(lgl, l2int(x));
    if (f == 0)
        return l_Undef;
    if (f > 0)
        return l_True;
    return l_False;
}

lbool LingelingSolver::value(Lit l) const { return value(var(l)) ^ sign(l); }

void forced_helper(void* r, int unit)
{
    vector<Lit>* f = reinterpret_cast<vector<Lit>*>(r);
    f->push_back(unit > 0 ? mkLit(unit - 1) : ~mkLit(unit - 1));
}

vector<Lit> LingelingSolver::forced() const
{
    vector<Lit> f;
    lglutrav(lgl, reinterpret_cast<void*>(&f), forced_helper);
    return f;
}

vector<Lit> LingelingSolver::entailed() const
{
    // same as forced, because that's all we have access to
    return forced();
}

vector<Lit> LingelingSolver::forcedAt(int lvl) const
{
    if( lvl == 0 )
        return forced();
    return {};
}

model_t const& LingelingSolver::model() const { return model_cache; }

LSet const& LingelingSolver::core() const { return core_cache; }

bool LingelingSolver::withinBudget() const
{
    if (cflbudget >= 0 && lglgetconfs(lgl) >= cflbudget)
        return false;
    if (propbudget >= 0 && lglgetprops(lgl) >= propbudget)
        return false;
    return true;
}

Statistics LingelingSolver::getStats() const
{
    return { solves, solves /* lingeling does not report # restarts*/,
        static_cast<uint64_t>(nVars()), static_cast<uint64_t>(nClauses()),
        static_cast<uint64_t>(lglgetdecs(lgl)), 0,
        static_cast<uint64_t>(lglgetprops(lgl)),
        static_cast<uint64_t>(lglgetconfs(lgl)), 0, 0, 0, 0,
        0, /* no info on learnt clauses */
        solutions, lglsec(lgl), 0 };
}

#endif
