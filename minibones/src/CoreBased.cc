#line 329 "CoreBased.nw"
#include "CoreBased.hh"
using namespace minibones;

#line 137 "CoreBased.nw"
CoreBased::CoreBased(ToolConfig& ptool_configuration, ostream& poutput, 
   Var pmax_id, const CNF& pclauses)
: tool_configuration(ptool_configuration)
, output(poutput)
, max_id(pmax_id)
, clauses(pclauses)
, solver_calls(0)
, solver_time(0)
, lifter(clauses)
, rotatable_computer(clauses)
{ /* */ }

CoreBased::~CoreBased() { /* */ }


#line 154 "CoreBased.nw"
bool CoreBased::initialize() {
  //initialize the solver
  sat_solver.new_variables(max_id);
  vec<Lit> ls;
  FOR_EACH(ci, clauses) {
    ls.clear ();
    FOR_EACH(li, *ci) {
      const Lit l = *li;
      assert (var(l) <= max_id);
      ls.push(l);
    }
    sat_solver.addClause(ls);
  }

  // get the first model
  const bool result = run_sat_solver();
  if (!result) {// the whole instance is not satisfiable
    return false;
  }

  const vec<lbool>& model = sat_solver.model;
  const Var max_model_var = std::min(model.size()-1,max_id);
  for (Var variable = 1; variable <= max_model_var; ++variable) {
    const lbool value = model[variable];
    if (value != l_Undef) {
      const Lit l = mkLit(variable, value==l_False);
      might_be.add(l);
    }
  }
  process_model(sat_solver.model);

  return true;
}

#line 192 "CoreBased.nw"
void CoreBased::run() {
  LitBitSet relaxed; // literals removed from assumptions
  vec<Lit>  assumptions; // assumptions passed onto the SAT solver
  vec<Lit>  old_assumptions; // assumptions in prev.\ iteration of inner loop
  while (might_be.size()) {
    const size_t chunk_size = make_chunk(old_assumptions);
    
#line 206 "CoreBased.nw"
relaxed.clear();
while (true) {
  
#line 220 "CoreBased.nw"
assumptions.clear();
for (int i=0; i<old_assumptions.size(); ++i) {
  const Lit l = old_assumptions[i];
  if (relaxed.get(l)) continue;
  assumptions.push(l);
}

#line 208 "CoreBased.nw"
                         ;
  const bool flipped = run_sat_solver(assumptions);
  if (flipped) { 
#line 230 "CoreBased.nw"
// cerr<<"SAT"<<endl;
process_model(sat_solver.model);
break; // we are done with the chunk

#line 210 "CoreBased.nw"
                               ; }
  else { 
#line 238 "CoreBased.nw"
const auto& conflict = sat_solver.conflict;
// cerr<<"RX:"<<relaxed.size()<<" MBS: "<<might_be.size()<<" CS: "<<conflict.size()<<endl;
if (sat_solver.conflict.size()==1) // A unit conflict means backbone. 
  mark_backbone(conflict[0]); 
for (int i=0; i<conflict.size(); ++i) 
  relaxed.add(~conflict[i]);
assert(relaxed.size()<=chunk_size);
if (relaxed.size()==chunk_size) {
  test_backbones(relaxed);
  break; // we are done with the chunk
}

#line 211 "CoreBased.nw"
                         ; }
  old_assumptions.clear();
  assumptions.copyTo(old_assumptions);
}


#line 198 "CoreBased.nw"
                     ;
  }
}

#line 252 "CoreBased.nw"
void CoreBased::test_backbones(const LitBitSet& negliterals) {
  vec<Lit> assumptions(1);
  FOR_EACH (li, negliterals) {
    const Lit negl=*li;
    const Lit    l=~negl;
    if (!might_be.get(l)) continue;
    assumptions[0]=negl;
    const bool flipped = run_sat_solver(assumptions);
    if (flipped) process_model(sat_solver.model);
    else mark_backbone(l);
  }
}

#line 268 "CoreBased.nw"
void CoreBased::process_model(const vec<lbool>& model) {
  if (tool_configuration.get_scdc_pruning()) {
    vec<lbool> rmodel;
    vec<Lit> assumptions;
    lifter.reduce(assumptions, model, rmodel);
    assert(lifter.is_model(rmodel));
    process_pruned_model(rmodel);
  } else if (tool_configuration.get_rotatable_pruning ()) {
    VarSet pruned;
    rotatable_computer.get_rotatables(model, pruned);
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

#line 293 "CoreBased.nw"
void CoreBased::process_pruned_model(const vec<lbool>& model) {
  for (Var variable = 1; variable < model.size(); ++variable) {
    const lbool value = model[variable];
    const Lit   pos_lit = mkLit(variable);
    const Lit   neg_lit = ~mkLit(variable);
    if (value != l_True) might_be.remove(pos_lit);
    if (value != l_False) might_be.remove(neg_lit);
  }

  for (Var variable=std::max(model.size(),1); variable<=max_id; ++variable) {
    const Lit pos_lit = mkLit(variable);
    might_be.remove( pos_lit);
    might_be.remove(~pos_lit);
  }
}

#line 311 "CoreBased.nw"
size_t CoreBased::make_chunk(vec<Lit>& literals) {
  const size_t chunk_size(tool_configuration.get_chunk_size() ? 
                          tool_configuration.get_chunk_size() : max_id);
  const int real_chunk_size = (int)std::min(chunk_size, might_be.size()); 
  literals.clear();
  FOR_EACH (li, might_be) {
    if (literals.size() >= real_chunk_size) break;
    const Lit lit = *li;
    literals.push(~lit);
  }
  return literals.size();
}

#line 325 "CoreBased.nw"
bool CoreBased::is_backbone(const Lit& literal) const 
{ return must_be.get(literal); }


