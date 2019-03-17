/*
 * File:  Rotatable.cc
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 11:34:02 AM
 * Copyright (C) 2012, Mikolas Janota
 */
#include "minisat_aux.hh"
#include "Rotatable.hh"
using Minisat::lit_Undef;
Rotatable::Rotatable(const CNF& _input_cnf)
 : input_cnf(_input_cnf)
{}

void Rotatable::get_rotatables(const vec<lbool>& model, VarSet& rotatables) {
  // first assume everything is rotatable
  for (Var index = 1; index < model.size(); ++index) rotatables.add(index);

  // remove unit literals from rotatables
  FOR_EACH (ci, input_cnf) {
    const LitSet&  clause=*ci;
    Lit last_satisfied = lit_Undef;
    size_t satisfied_count = 0;
    FOR_EACH (li, clause) {
      const Lit l = *li;
      if (is_true(l, model)) {
        ++satisfied_count;
        last_satisfied = l;
      }
    }
    assert(satisfied_count); // must be a model
    assert(last_satisfied!=lit_Undef); // must be a model
    if (satisfied_count==1)  rotatables.remove(var(last_satisfied));  // unit literal
  }
}


