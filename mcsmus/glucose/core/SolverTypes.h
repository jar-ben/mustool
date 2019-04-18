/***********************************************************************************[SolverTypes.h]
 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France

Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose are exactly the same as Minisat on which it is based on. (see below).

---------------
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#ifndef Glucose_SolverTypes_h
#define Glucose_SolverTypes_h

#include <assert.h>

#include "minisat/mtl/mcsmus_IntTypes.h"
#include "minisat/mtl/mcsmus_Alg.h"
#include "minisat/mtl/mcsmus_Vec.h"
#include "minisat/mtl/mcsmus_Map.h"
#include "minisat/mtl/mcsmus_Alloc.h"
#include "minisat/mtl/mcsmus_IntMap.h"
#include "minisat/mtl/mcsmus_Queue.h"
#include "minisat/mtl/mcsmus_Heap.h"

#include "minisat/core/mcsmus_SolverTypes.h"

namespace Glucose {

//=================================================================================================
// Variables, literals, lifted booleans, clauses:

template<typename T>
using vec = Minisat::vec<T>;
template<typename K, typename T>
using IntMap = Minisat::IntMap<K, T>;

template<typename T>
using Queue = Minisat::Queue<T>;

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

using MkIndexLit = Minisat::MkIndexLit;

template<typename T>
using VMap = Minisat::VMap<T>;
template<typename T>
using LMap = Minisat::LMap<T>;

using LSet = Minisat::LSet;

template<typename T>
using RegionAllocator = Minisat::RegionAllocator<T>;

template<class K, class Comp, class MkIndex = Minisat::MkIndexDefault<K> >
using Heap = Minisat::Heap<K, Comp, MkIndex>;

//=================================================================================================
// Clause -- a simple class for representing a clause:

class Clause;
typedef RegionAllocator<uint32_t>::Ref CRef;

class Clause {
    struct {
      unsigned mark      : 2;
      unsigned learnt    : 1;
      unsigned has_extra : 1;
      unsigned reloced   : 1;
      unsigned canbedel  : 1;
      unsigned lbd       : 26;
      unsigned size      : 32;
      unsigned szWithoutSelectors : 32;

    }                            header;
    union { Lit lit; float act; uint32_t abs; CRef rel; } data[0];

    friend class ClauseAllocator;

    // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
    template<class V>
    Clause(const V& ps, bool use_extra, bool learnt) {
        header.mark      = 0;
        header.learnt    = learnt;
        header.has_extra = use_extra;
        header.reloced   = 0;
        header.size      = ps.size();
        header.lbd = 0;
        header.canbedel = 1;
        for (int i = 0; i < ps.size(); i++)
            data[i].lit = ps[i];

        if (header.has_extra){
          if (header.learnt)
                data[header.size].act = 0;
            else
                calcAbstraction(); }
    }

public:
    void calcAbstraction() {
        assert(header.has_extra);
        uint32_t abstraction = 0;
        for (int i = 0; i < size(); i++)
            abstraction |= 1 << (var(data[i].lit) & 31);
        data[header.size].abs = abstraction;  }


    int          size        ()      const   { return header.size; }
    void         shrink      (int i)         { assert(i <= size()); if (header.has_extra) data[header.size-i] = data[header.size]; header.size -= i; }
    void         pop         ()              { shrink(1); }
    bool         learnt      ()      const   { return header.learnt; }
    bool         has_extra   ()      const   { return header.has_extra; }
    uint32_t     mark        ()      const   { return header.mark; }
    void         mark        (uint32_t m)    { header.mark = m; }
    const Lit&   last        ()      const   { return data[header.size-1].lit; }

    bool         reloced     ()      const   { return header.reloced; }
    CRef         relocation  ()      const   { return data[0].rel; }
    void         relocate    (CRef c)        { header.reloced = 1; data[0].rel = c; }

    // NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
    //       subsumption operations to behave correctly.
    Lit&         operator [] (int i)         { return data[i].lit; }
    Lit          operator [] (int i) const   { return data[i].lit; }
    operator const Lit* (void) const         { return (Lit*)data; }

    float&       activity    ()              { assert(header.has_extra); return data[header.size].act; }
    uint32_t     abstraction () const        { assert(header.has_extra); return data[header.size].abs; }

    Lit          subsumes    (const Clause& other) const;
    void         strengthen  (Lit p);
    void         setLBD(int i)  {header.lbd = i;}
    // unsigned int&       lbd    ()              { return header.lbd; }
    unsigned int        lbd    () const        { return header.lbd; }
    void setCanBeDel(bool b) {header.canbedel = b;}
    bool canBeDel() {return header.canbedel;}
    void setSizeWithoutSelectors   (unsigned int n)              {header.szWithoutSelectors = n; }
    unsigned int        sizeWithoutSelectors   () const        { return header.szWithoutSelectors; }

};

inline Lit const *begin(Clause const&c) { return (Lit const*)c; }
inline Lit const *end(Clause const&c) { return begin(c)+c.size(); }

inline
std::ostream& operator<<(std::ostream& os, Clause& c)
{
    os << "cl[";
    bool first = true;
    for( auto&& t: c ) {
        if( first )
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}


//=================================================================================================
// ClauseAllocator -- a simple class for allocating memory for clauses:


const CRef CRef_Undef = RegionAllocator<uint32_t>::Ref_Undef;
class ClauseAllocator : public RegionAllocator<uint32_t>
{
    static int clauseWord32Size(int size, bool has_extra){
        return (sizeof(Clause) + (sizeof(Lit) * (size + (int)has_extra))) / sizeof(uint32_t); }
 public:
    bool extra_clause_field;

    ClauseAllocator(uint32_t start_cap) : RegionAllocator<uint32_t>(start_cap), extra_clause_field(false){}
    ClauseAllocator() : extra_clause_field(false){}

    void moveTo(ClauseAllocator& to){
        to.extra_clause_field = extra_clause_field;
        RegionAllocator<uint32_t>::moveTo(to); }

    template<class Lits>
    CRef alloc(const Lits& ps, bool learnt = false)
    {
        assert(sizeof(Lit)      == sizeof(uint32_t));
        assert(sizeof(float)    == sizeof(uint32_t));
        bool use_extra = learnt | extra_clause_field;

        CRef cid = RegionAllocator<uint32_t>::alloc(clauseWord32Size(ps.size(), use_extra));
        new (lea(cid)) Clause(ps, use_extra, learnt);

        return cid;
    }

    // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
    Clause&       operator[](Ref r)       { return (Clause&)RegionAllocator<uint32_t>::operator[](r); }
    const Clause& operator[](Ref r) const { return (Clause&)RegionAllocator<uint32_t>::operator[](r); }
    Clause*       lea       (Ref r)       { return (Clause*)RegionAllocator<uint32_t>::lea(r); }
    const Clause* lea       (Ref r) const { return (Clause*)RegionAllocator<uint32_t>::lea(r); }
    Ref           ael       (const Clause* t){ return RegionAllocator<uint32_t>::ael((uint32_t*)t); }

    void free(CRef cid)
    {
        Clause& c = operator[](cid);
        RegionAllocator<uint32_t>::free(clauseWord32Size(c.size(), c.has_extra()));
    }

    void reloc(CRef& cr, ClauseAllocator& to)
    {
        Clause& c = operator[](cr);

        if (c.reloced()) { cr = c.relocation(); return; }

        cr = to.alloc(c, c.learnt());
        c.relocate(cr);

        // Copy extra data-fields:
        // (This could be cleaned-up. Generalize Clause-constructor to be applicable here instead?)
        to[cr].mark(c.mark());
        if (to[cr].learnt())        {
          to[cr].activity() = c.activity();
          to[cr].setLBD(c.lbd());
          to[cr].setSizeWithoutSelectors(c.sizeWithoutSelectors());
          to[cr].setCanBeDel(c.canBeDel());
        }
        else if (to[cr].has_extra()) to[cr].calcAbstraction();
    }
};


//=================================================================================================
// OccLists -- a class for maintaining occurence lists with lazy deletion:

template<class K, class Vec, class Deleted, class MkIndex = Minisat::MkIndexDefault<K> >
using OccLists = Minisat::OccLists<K, Vec, Deleted, MkIndex>;

//=================================================================================================
// CMap -- a class for mapping clauses to values:

template<typename T>
using CMap = Minisat::CMap<T>;

/*_________________________________________________________________________________________________
|
|  subsumes : (other : const Clause&)  ->  Lit
|
|  Description:
|       Checks if clause subsumes 'other', and at the same time, if it can be used to simplify 'other'
|       by subsumption resolution.
|
|    Result:
|       lit_Error  - No subsumption or simplification
|       lit_Undef  - Clause subsumes 'other'
|       p          - The literal p can be deleted from 'other'
|________________________________________________________________________________________________@*/
inline Lit Clause::subsumes(const Clause& other) const
{
    //if (other.size() < size() || (extra.abst & ~other.extra.abst) != 0)
    //if (other.size() < size() || (!learnt() && !other.learnt() && (extra.abst & ~other.extra.abst) != 0))
    assert(!header.learnt);   assert(!other.header.learnt);
    assert(header.has_extra); assert(other.header.has_extra);
    if (other.header.size < header.size || (data[header.size].abs & ~other.data[other.header.size].abs) != 0)
        return lit_Error;

    Lit        ret = lit_Undef;
    const Lit* c   = (const Lit*)(*this);
    const Lit* d   = (const Lit*)other;

    for (unsigned i = 0; i < header.size; i++) {
        // search for c[i] or ~c[i]
        for (unsigned j = 0; j < other.header.size; j++)
            if (c[i] == d[j])
                goto ok;
            else if (ret == lit_Undef && c[i] == ~d[j]){
                ret = c[i];
                goto ok;
            }

        // did not find it
        return lit_Error;
    ok:;
    }

    return ret;
}

inline void Clause::strengthen(Lit p)
{
    remove(*this, p);
    calcAbstraction();
}

//=================================================================================================
}


#endif
