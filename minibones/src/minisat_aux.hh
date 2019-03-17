/* 
 * File:   minisat_aux.hh
 * Author: mikolas
 *
 * Created on October 21, 2010, 2:10 PM
 */

#ifndef MINISAT_AUX_HH
#define	MINISAT_AUX_HH
#include "minisat/core/Solver.h"
#include "auxiliary.hh"
#include <ostream>
using std::ostream;
using std::vector;
using Minisat::Solver;
using Minisat::lbool;
using Minisat::l_False;
using Minisat::l_True;
using Minisat::l_Undef;
using Minisat::mkLit;
using Minisat::sign;
using Minisat::var;
using Minisat::vec;
using Minisat::Lit;
using Minisat::Var;
typedef std::vector<Lit> LiteralVector;

inline bool is_true(Lit l, const vec<lbool>& model) { 
  const Var v = var(l);
  assert (v<model.size());
  return (model[v]==l_False)==sign(l);
}

ostream& print_model(ostream& out, const vec<lbool>& lv);
void     print(const vec<Lit>& lv);
ostream& print(const Lit& l);
ostream& operator << (ostream& outs, Lit lit);

ostream& print(ostream& out, const vec<Lit>& lv);
inline ostream& operator << (ostream& out, const vec<Lit>& lv) { return print(out,lv); }
inline ostream& operator << (ostream& out, const Minisat::LSet& lv) { return print(out,lv.toVec()); }

inline size_t literal_index(Lit l) { 
  assert(var(l) > 0);
  const size_t v = (size_t) var(l);
  return sign(l) ? v<<1 : ((v<<1)+1);
}

inline Lit index2literal(size_t l) { 
  const bool positive = (l&1);
  const Var  variable = l>>1;
  return positive ? mkLit(variable) : ~mkLit(variable);
}

class Lit_equal {
public:
  inline bool operator () (const Lit& l1,const Lit& l2) const { return l1==l2; }
};

class Lit_hash {
public:
  inline size_t operator () (const Lit& l) const { return Minisat::toInt(l); }
};
#endif	/* MINISAT_AUX_HH */

