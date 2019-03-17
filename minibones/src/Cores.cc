/*
 * File:  Cores.cc
 * Author:  mikolas
 * Created on:  Sat, Mar 17, 2012 5:54:18 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#include "Cores.hh"
using namespace minibones;
#define CORE_DBG(t)

Cores::Cores(ToolConfig& _tool_configuration, Var _max_id, const CNF& _clauses)
  : tool_configuration(_tool_configuration)
  , max_id(_max_id)
  , clauses(_clauses)
{
  //initialize the solver
  for (Var i=0; i<=max_id; ++i) solver.newVar();
  vec<Lit> ls;
  FOR_EACH(ci, clauses) {
    ls.clear ();
    FOR_EACH(li, *ci) {
      const Lit l = *li;
      assert (var(l) <= max_id);
      ls.push(l);
    }
    solver.addClause(ls);
  }
}

Cores::~Cores() {}

bool Cores::try_to_flip(const LitBitSet& might_be, vec<Lit>& reasons) {
  vec<Lit> assumptions(might_be.size());
  const_infinite_LitBitSetIterator it = might_be.infinite_iterator();
  size_t i = might_be.size();
  while (i) {
    --i;
    const Lit l = *it;
    assumptions[i] = ~l;
    ++it;
  }
  CORE_DBG( cerr << "assumptions: " << assumptions << endl; )
  const bool satisfiable = solver.solve(assumptions);
  reasons.clear();
  if (satisfiable) return true;
  CORE_DBG( cerr << "core: " << solver.conflict << endl; )
  solver.conflict.toVec().copyTo(reasons);
  return false;
}
