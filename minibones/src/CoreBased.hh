#line 10 "CoreBased.nw"
#ifndef COREBASED_W_3637
#define COREBASED_W_3637
#line 116 "CoreBased.nw"
#include <vector>
#include <functional>
#include "MiniSatExt.hh"
#include "auxiliary.hh"
#include "Lifter.hh"
#include "ToolConfig.hh"
#include "Rotatable.hh"
#include "BackboneInformation.hh"
#include "LitBitSet.hh"
using MinibonesMinisat::MiniSatExt;
using MinibonesMinisat::Var;
using MinibonesMinisat::Lit;
using MinibonesMinisat::vec;

#line 18 "CoreBased.nw"
namespace minibones {
  class CoreBased : public BackboneInformation {
  public: 
#line 28 "CoreBased.nw"
CoreBased(ToolConfig& tool_configuration, 
          std::ostream& output, 
          Var max_id, const CNF& clauses);
virtual ~CoreBased();
bool initialize();// initialize, returns true iff the instance is SAT
void run();// start the worker, initialize must be called first.
virtual bool is_backbone(const Lit& literal) const;
size_t  get_solver_calls() const { return solver_calls; }
size_t  get_solver_time() const  { return solver_time; }

#line 21 "CoreBased.nw"
  private: 
#line 42 "CoreBased.nw"
LitBitSet           might_be; 
LitBitSet           must_be;  

#line 47 "CoreBased.nw"
const ToolConfig&   tool_configuration;
ostream&            output;
const Var           max_id;
const CNF&          clauses;

#line 54 "CoreBased.nw"
size_t              solver_calls;// number of solver calls, for statistical purposes
double              solver_time;// for statistical purposes


#line 60 "CoreBased.nw"
Lifter              lifter;// used to reduce models via lifting
Rotatable           rotatable_computer;// used to get rotatable variables


#line 66 "CoreBased.nw"
MiniSatExt sat_solver;
void process_model(const vec<lbool>& model);// calls pruning techniques and debones based of the model
void process_pruned_model(const vec<lbool>& model);// debones everything not in the model
size_t make_chunk(vec<Lit>& literals);
inline bool debone(Lit literal);// mark a literal as not a backbone
inline void mark_backbone(Lit literal);// mark a literal as a backbone
void test_backbones(const LitBitSet& literals);
inline bool run_sat_solver(const vec<Lit>& assumptions);
inline bool run_sat_solver();


#line 22 "CoreBased.nw"
  };
  
#line 79 "CoreBased.nw"
inline bool CoreBased::debone(Lit literal) { 
  assert(!must_be.get(literal));
  return might_be.remove(literal);
}

#line 86 "CoreBased.nw"
inline void CoreBased::mark_backbone(Lit literal) { 
  assert(might_be.get(literal));
  might_be.remove(literal);
  might_be.remove(~literal);
  if (must_be.add(literal)) {
    if (tool_configuration.get_backbone_insertion()) {
      sat_solver.addClause(literal);
    }
  }
}

#line 99 "CoreBased.nw"
bool CoreBased::run_sat_solver(const vec<Lit>& assumptions) {
  const auto t0 = read_cpu_time();
  const bool retv = sat_solver.solve(assumptions);
  const auto t1 = read_cpu_time();
  solver_time+=(t1-t0);
  ++solver_calls;
  return retv;
}

bool CoreBased::run_sat_solver() {
  const vec<Lit> assumptions;
  return run_sat_solver(assumptions);
}


#line 24 "CoreBased.nw"
}

#line 14 "CoreBased.nw"
#endif

