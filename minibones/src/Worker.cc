/* 
 * File:   Worker.cc
 * Author: mikolas
 * 
 * Created on March 19, 2011, 4:51 PM
 */
#include "Worker.hh"
#include <assert.h>
#include "minisat_aux.hh"
#include "BBInfo.hh"
#include <limits>
#include <time.h>
using namespace minibones;

/*------------------------- initialization -----------------------------------*/
Worker::Worker(ToolConfig& _tool_configuration,
        ostream& _output,        
        Var _max_id,
        const CNF& _clauses,
        const Range& _variable_range)
: tool_configuration(_tool_configuration)
, output(_output)
, max_id(_max_id)
, clauses(_clauses)
, variable_range(_variable_range)
, solver_time(0)
, solver_calls(0)
, bbInfo(max_id)
, lifter(clauses)
, rotatable_computer(clauses)
{
    assert (variable_range.first <= 1);
    assert (variable_range.first <= variable_range.second);
    assert (variable_range.second<=max_id);
}

bool Worker::initialize() {
  //initialize testing progress
  to_test.resize(max_id+1, false);
  to_test_count=0;
  for (Var variable = variable_range.first; variable <= variable_range.second; ++variable) {
    to_test[(size_t)variable]=true;
    to_test_count+=2; //two literals for each variable
  }
  //initialize the solver
  for (Var i=0; i<=max_id; ++i) solver.newVar();
  vec<Lit> ls;
  FOR_EACH(ci, clauses) {
    ls.clear ();
    FOR_EACH(li, *ci) ls.push(*li);
    solver.addClause(ls);
  }

  // get the first model
  const bool result = run_solver();
  if (!result) {// the whole instance is not satisfiable
    return false;
  }
  process_model(solver.model);
  return true;
}

/*----------------------  run of the worker ----------------------------------*/

void Worker::run() {
  output << "c starting computation"<<endl;
  for (Var variable=1; variable<=variable_range.second; ++variable) {
    if (should_stop()) break; //stop the work
    if (!to_test[variable]) continue;
    const time_t sst1 = time(NULL);
    test_backbone(mkLit(variable));
    test_backbone(~mkLit(variable));
    const time_t sst2 = time(NULL);
    const time_t d = sst2-sst1;
    if (d > 5) output << "c long sat[s] " << d << endl;
    to_test[variable]    = false;
  }
  assert(is_complete());
}

bool Worker::should_stop() const { return to_test_count==0; }
/*----------------------  backbone testing  ----------------------------------*/

bool Worker::test_backbone(const Lit& literal) {
  if (!bbInfo.might_be_bb(literal)) return false; // we already know that it is not a backbone
  if (bbInfo.is_bb(literal)) return true;         // we already know that it is a backbone
  // call the solver
  vec<Lit> assumptions(1);
  assumptions[0]=~literal;
  const bool sat = run_solver(assumptions);
  // analyze solver's result
  VarSet pruned;

  if (!sat) { // is a backbone
    mark_bone(literal);
    return true;
  } else { // not a backbone        
    assert(lifter.is_model(solver.model));
    if (tool_configuration.get_scdc_pruning()) {
      vec<lbool> rmodel;
      lifter.reduce(assumptions, solver.model, rmodel);
      assert(lifter.is_model(rmodel));
      process_model(rmodel);
    } else if (tool_configuration.get_rotatable_pruning ()) {
      pruned.clear();
      rotatable_computer.get_rotatables(solver.model, pruned);
      FOR_EACH (vi, pruned) {
        const Lit l = mkLit(*vi);
        debone( l);
        debone(~l);
      }
      process_model(solver.model);
    } else {
      process_model(solver.model);
    }
    return false;
  }
}

void Worker::process_model(const vec<lbool>& model) {
  assert (model.size()>max_id);
  for (Var index=1; index<=max_id; ++index) {
    const lbool v = model[index];
    if (v != l_True)  debone(mkLit(index));
    if (v != l_False) debone(~mkLit(index));
  }
}

/*-----------------------------  getters -------------------------------------*/

bool Worker::is_backbone(const Lit& literal) const
{return bbInfo.is_bb(literal);}

/*------------------------  debugging ----------------------------------------*/

bool Worker::is_complete() const {
  for (Var variable=1; variable<=max_id; ++variable) {
    const Lit pl = mkLit(variable);
    const Lit nl = ~mkLit(variable);
    const bool ok = !to_test[variable]
      || (   (!bbInfo.might_be_bb(pl) || bbInfo.is_bb(pl))
             && (!bbInfo.might_be_bb(nl) || bbInfo.is_bb(nl)) );
    assert (!bbInfo.is_bb(pl) || !bbInfo.is_bb(nl));
    if (!ok) return false;
  }
  return true;
}


bool Worker::run_solver(const vec<Lit>& assumptions) {
  const auto t0 = read_cpu_time();
  const bool retv = solver.solve(assumptions);
  const auto t1 = read_cpu_time();
  solver_time+=(t1-t0);
  ++solver_calls;
  return retv;
}
bool Worker::run_solver() {
  const vec<Lit> assumptions;
  return run_solver(assumptions);
}


/*----------------------------  dummy ----------------------------------------*/
//Worker::Worker(const Worker& orig)
//: output(orig.output)
//, rank (orig.rank)
//, processes_count(orig.processes_count)
//, max_id(orig.max_id)
//, clauses(orig.clauses)
//, variable_range(orig.variable_range)
//, solver_calls(orig.solver_calls)
//, lifter(orig.lifter)
//, data_sz(orig.data_sz)
//{assert(false);}

Worker::~Worker() {}
