/*
 * File:  UpperBound.cc
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 3:45:21 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#include "UpperBound.hh"
using namespace minibones; 

#define UPPERBOUND_DBG(t)

/*------------------------- initialization -----------------------------------*/
UpperBound::UpperBound(ToolConfig& _tool_configuration, ostream& _output,  Var _max_id, const CNF& _clauses)
  : tool_configuration(_tool_configuration)
  , chunk_size (tool_configuration.get_chunk_size() ? tool_configuration.get_chunk_size() : _max_id)
  , output(_output)
  , max_id(_max_id)
  , clauses(_clauses)
  , solver_calls(0)
  , solver_time(0)
  , lifter(clauses)
  , rotatable_computer(clauses)
  , cores(_tool_configuration,_max_id,_clauses)
  , might_be_iterator(might_be.infinite_iterator())
  , end_of_chunk(max_id)
{}

UpperBound::~UpperBound() {}

bool UpperBound::initialize() {
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

  // get the first model
  const bool result = run_solver();
  if (!result) {// the whole instance is not satisfiable
    return false;
  }

  // TODO apply initial pruning

  const vec<lbool>& model = solver.model;
  UPPERBOUND_DBG( print_model(cerr << "initial model: ", model) << endl; );

  for (Var variable = 1; variable <= max_id; ++variable) {
    if (max_id >= model.size()) continue;
    const lbool value = model[variable];
    if (value != l_Undef) {
      const Lit l = mkLit(variable, value==l_False);
      might_be.add(l);
    }
  }

  UPPERBOUND_DBG( might_be.print(cerr << "initial might be:")<< endl; );

  process_model(solver.model);

  return true;
}

/*----------------------  run of the worker ----------------------------------*/

void UpperBound::run() {
  vec<Lit> literals;
  vec<Lit> assumptions(1);
  Lit relaxation_literal = lit_Undef; 
  while (might_be.size()) { 
    UPPERBOUND_DBG( might_be.print(cerr << "might be:")<< endl; );
    if (relaxation_literal != lit_Undef) solver.addClause(relaxation_literal); // relax the last chunk clause (forever)
    // prepare new chunk clause
    literals.clear();
    relaxation_literal = make_chunk(literals);
    // call solver on the new chunk
    solver.addClause(literals);
    if (tool_configuration.get_use_variable_bumping()) {
      for (int i=0; i<literals.size(); ++i) {
        const Lit l=literals[i];
        if (l==relaxation_literal) continue;
        const Var v=var(literals[i]);
        solver.bump(v);
        solver.setPolarity(v, sign(l) ? l_True : l_False);
      }
    }
    UPPERBOUND_DBG( cerr << "chunk clause: " << literals << endl; );
    assumptions[0] = ~relaxation_literal; // disallow relaxing current chunk
    const bool is_sat = run_solver(assumptions);
    // analyze solver's output
    if (!is_sat) { // complements of the literals in the chunk are backbones
      UPPERBOUND_DBG( cerr << "bb" << endl; );
      assert(literals[0]==relaxation_literal);
      for (int index = 1; index < literals.size (); ++index) {
        const Lit backbone = ~literals[index];
        might_be.remove(backbone);
        must_be.add(backbone);
        if (tool_configuration.get_backbone_insertion()) solver.addClause(backbone);
      }
    } else {
      const vec<lbool>& model = solver.model;
      process_model(model);
    }
  }
}

Lit UpperBound::make_chunk(vec<Lit>& literals) {
  const Lit relaxation_literal = mkLit(solver.newVar());      // generate new relaxation variable
  literals.push(relaxation_literal);                          // enable the new clause to be relaxed in the future
  const int real_chunk_size = std::min(chunk_size, (int)might_be.size()) + 1; // +1 is there for relaxation literal

  if (tool_configuration.get_use_chunk_keeping()) {
    Var variable;
    for (variable=1; variable <= end_of_chunk; ++variable) {
      if (literals.size() >= real_chunk_size) break;
      const Lit pl = mkLit(variable);
      const Lit nl = ~pl;
      const bool may_pl = might_be.get(pl);
      const bool may_nl = might_be.get(nl);
      assert(!may_nl || !may_pl);
      if (may_pl) literals.push(nl);
      if (may_nl) literals.push(pl);
      if ( (1==literals.size()) && (end_of_chunk==variable) ) end_of_chunk=max_id;
    }
    end_of_chunk=variable;
  } if (tool_configuration.get_use_random_chunks()) {
    const double probability = (double)(real_chunk_size-1) / (double)might_be.size();
    const bool add_all = (size_t)chunk_size>=might_be.size();
    for (Var variable = 1; variable <= max_id; ++variable) {
      if (real_chunk_size<=literals.size()) break;
      const Lit pl = mkLit(variable);
      const Lit nl = ~pl;
      const bool may_pl = might_be.get(pl);
      const bool may_nl = might_be.get(nl);
      assert(!may_nl || !may_pl);
      if (!may_nl && !may_pl) continue; // no literal for this variable
      const double coin = (double)rand()/(double)RAND_MAX;
      const bool should_add = add_all || (coin<probability);
      if (should_add && may_pl) literals.push(nl);
      if (should_add && may_nl) literals.push(pl);
    }
  } if (tool_configuration.get_use_cores())  {
    while (true) {
      vec<Lit> reasons;
      const bool flipped = cores.try_to_flip(might_be, reasons);
      if (flipped) { // no more backbones
        might_be.clear();
        break;
      }
      assert(reasons.size());
      if (reasons.size()==1) {
        const Lit backbone = reasons[0];
        assert(might_be.get(backbone));
        might_be.remove(backbone);
        must_be.add(backbone);
      } else {
        for (int i=0; i<reasons.size(); ++i) literals.push(~reasons[i]);
        break;
      }
    }
  } else {
    while (literals.size() < real_chunk_size) {
      ++might_be_iterator;
      const Lit lit = *might_be_iterator;
      literals.push(~lit);
    }
  }
  return relaxation_literal;
}

void UpperBound::process_model(const vec<lbool>& model) {
  if (tool_configuration.get_scdc_pruning()) {
    vec<lbool> rmodel;
    vec<Lit> assumptions;
    lifter.reduce(assumptions, model, rmodel);
    assert(lifter.is_model(rmodel));
    process_pruned_model(rmodel);
  } else if (tool_configuration.get_rotatable_pruning ()) {
    VarSet pruned;
    rotatable_computer.get_rotatables(solver.model, pruned);
    FOR_EACH (vi, pruned) {
      const Lit pl = mkLit(*vi);
      const Lit nl = ~pl;
      debone(pl);
      debone(nl);
    }
    process_pruned_model(model);
  } else {
    process_pruned_model(model);
  }
}

void UpperBound::process_pruned_model(const vec<lbool>& model) {
  for (Var variable = 1; variable < model.size(); ++variable) {
    const lbool value = model[variable];
    const Lit   pos_lit = mkLit(variable);
    const Lit   neg_lit = ~mkLit(variable);
    if (value != l_True) might_be.remove(pos_lit);
    if (value != l_False) might_be.remove(neg_lit);
  }

  for (Var variable=model.size(); variable<=max_id; ++variable) {
    const Lit pos_lit = mkLit(variable);
    might_be.remove( pos_lit);
    might_be.remove(~pos_lit);
  }
}

bool UpperBound::run_solver(const vec<Lit>& assumptions) {
  const auto t0 = read_cpu_time();
  const bool retv = solver.solve(assumptions);
  const auto t1 = read_cpu_time();
  solver_time+=(t1-t0);
  ++solver_calls;
  return retv;
}
bool UpperBound::run_solver() {
  const vec<Lit> assumptions;
  return run_solver(assumptions);
}

/*-----------------------------  getters -------------------------------------*/
bool UpperBound::is_backbone(const Lit& literal) const { return must_be.get(literal); }
