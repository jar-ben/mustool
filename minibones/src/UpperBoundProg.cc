/*
 * File:  UpperBoundProg.cc
 * Author:  mikolas
 * Created on:  Thu Feb 23 14:04:59 WET 2012
 * Copyright (C) 2012, Mikolas Janota
 */
#include "UpperBoundProg.hh"
using namespace minibones; 

/*------------------------- initialization -----------------------------------*/
UpperBoundProg::UpperBoundProg(ToolConfig& _tool_configuration, ostream& _output,  Var _max_id, const CNF& _clauses)
  : tool_configuration(_tool_configuration)
  , chunk_size (tool_configuration.get_chunk_size() ? tool_configuration.get_chunk_size() : _max_id)
  , output(_output)
  , max_id(_max_id)
  , clauses(_clauses)
  , solver_calls(0)
  , solver_time(0)
  , lifter(clauses)
  , rotatable_computer(clauses)
{}

UpperBoundProg::~UpperBoundProg() {}

bool UpperBoundProg::initialize() {
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

  VarSet pruned;
  if (tool_configuration.get_rotatable_pruning()) {
    rotatable_computer.get_rotatables(model, pruned);
  }

  vec<Lit> prog_cls;
  for (Var variable = 1; variable <= max_id; ++variable) {
    if (max_id >= model.size()) continue;
    const lbool value = model[variable];
    if (value==l_Undef || pruned.get(variable)) continue;
    // insert complement of the might be into programmable clause
    const Lit mbl = mkLit(variable, value==l_False);
    might_be.add(mbl);
    const Lit cl  = ~mbl;
    const Lit cl_selector = mkLit(solver.newVar());
    assert(!CONTAINS(lit_selectors, cl));
    lit_selectors[cl]=cl_selector;
    const Lit aux_l = mkLit(solver.newVar());
    prog_cls.push(aux_l);
    solver.addClause(~aux_l, cl_selector);
    solver.addClause(~aux_l, cl);
    UPPERBOUND_DBG( 
                   cerr << "pcl: " << ~aux_l << " " <<  cl_selector << endl;
                   cerr << "pcl: " << ~aux_l << " " <<  cl << endl;);
  }
  solver.addClause(prog_cls);
  UPPERBOUND_DBG(cerr << "pcl: " << prog_cls  << endl;);

  UPPERBOUND_DBG( might_be.print(cerr << "initial might be:")<< endl; );

  process_model(solver.model);

  tool_configuration.prefix(output) << "UpperBoundProg initialized" << endl;

  return true;
}

/*----------------------  run of the worker ----------------------------------*/

void UpperBoundProg::run() {
  vec<Lit> assumptions;
  Var chunk_end = 0;
  while (might_be.size()) { 
    UPPERBOUND_DBG( might_be.print(cerr << "might be:")<< endl; );
    // prepare new chunk
    Var chunk_start=chunk_end+1;
    assert(chunk_start<=max_id);
    chunk_end = std::min(chunk_start+chunk_size-1,max_id);
    while (true) {
      assumptions.clear();
      size_t mb_chunk=make_chunk(chunk_start, chunk_end, assumptions); 
      UPPERBOUND_DBG(cerr << "cassum: " << assumptions << endl;);
      UPPERBOUND_DBG(cerr << "csz: " << mb_chunk << endl;);
      if (!mb_chunk) break; // empty chunk
      // call solver on the new chunk
      const bool is_sat = run_solver(assumptions);
      // analyze solver's output
      if (!is_sat) { // complements of the literals in the chunk are backbones
        UPPERBOUND_DBG( cerr << "bb" << endl; );
        for (Var v=chunk_start; v<=chunk_end; ++v) {
          const Lit pl = mkLit(v);
          const Lit nl = ~pl;
          const bool mbp = might_be.get(pl);
          const bool mbn = might_be.get(nl);
          assert (!mbn || !mbp);
          if (!mbn && !mbp) continue;
          const Lit backbone = mbn ? nl : pl;
          might_be.remove(backbone);
          must_be.add(backbone);
          solver.addClause(backbone);
          solver.addClause(~get_selector(~backbone));
        }
        break; // chunk emptied
      } else {
        const vec<lbool>& model = solver.model;
        UPPERBOUND_DBG( print_model(cerr << "model found: ", model) << endl; );
        process_model(model);
      }
    }
  }
}

size_t UpperBoundProg::make_chunk(const Var chunk_start, const Var chunk_end, vec<Lit>& assumptions) {
  size_t mb_chunk=0; 
  for (Var v=1; v <= max_id; ++v) {
    const Lit pl = mkLit(v);
    const Lit nl = ~pl;
    const bool mbp = might_be.get(pl);
    const bool mbn = might_be.get(nl);
    assert (!mbn || !mbp);
    const bool in_chunk = (chunk_start<=v) && (v<=chunk_end);
    if (in_chunk) {
      if (mbn || mbp) ++mb_chunk;
      if (mbn) assumptions.push(get_selector(pl));
      if (mbp) assumptions.push(get_selector(nl));
    } else {
      if (mbn) assumptions.push(~get_selector(pl));
      if (mbp) assumptions.push(~get_selector(nl));
    }
  }
  return mb_chunk;
}

void UpperBoundProg::process_model(const vec<lbool>& model) {
  if (tool_configuration.get_scdc_pruning()) {
    vec<lbool> rmodel;
    vec<Lit> assumptions;
    lifter.reduce(assumptions, model, rmodel);
    assert(lifter.is_model(rmodel));
    process_pruned_model(rmodel);
  } else if (tool_configuration.get_rotatable_pruning()) {
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

void UpperBoundProg::process_pruned_model(const vec<lbool>& model) {
  for (Var variable = 1; variable < model.size(); ++variable) {
    const lbool value = model[variable];
    const Lit   pos_lit = mkLit(variable);
    const Lit   neg_lit = ~mkLit(variable);
    if (value != l_True) debone(pos_lit);
    if (value != l_False) debone(neg_lit);
  }

  for (Var variable=model.size(); variable<=max_id; ++variable) {
    const Lit pos_lit = mkLit(variable);
    debone( pos_lit);
    debone(~pos_lit);
  }
}

bool UpperBoundProg::run_solver(const vec<Lit>& assumptions) {
  const auto t0 = read_cpu_time();
  const bool retv = solver.solve(assumptions);
  const auto t1 = read_cpu_time();
  solver_time+=(t1-t0);
  ++solver_calls;
  return retv;
}
bool UpperBoundProg::run_solver() {
  const vec<Lit> assumptions;
  return run_solver(assumptions);
}
/*-----------------------------  getters -------------------------------------*/
bool UpperBoundProg::is_backbone(const Lit& literal) const { return must_be.get(literal); }
