/*
 * File:  UpperBoundProg.hh
 * Author:  mikolas
 * Created on:  Thu Feb 23 14:04:54 WET 2012
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef UPPERBOUNDPROG_HH_24653
#define UPPERBOUNDPROG_HH_24653
#include <vector>
#include <functional>
#include "minisat/core/Solver.h"
#include "auxiliary.hh"
#include "Lifter.hh"
#include "LitBitSet.hh"
#include "ToolConfig.hh"
#include "Rotatable.hh"
#include "BackboneInformation.hh"
#define UPPERBOUND_DBG(t)
namespace minibones {
  using Minisat::Solver;
  using Minisat::Var;
  using Minisat::Lit;
  using Minisat::vec;

  /** Class for computing backbones with the upper bound approach. */
  class UpperBoundProg : public BackboneInformation {
  public:
    UpperBoundProg(ToolConfig& tool_configuration, ostream& output, Var max_id, const CNF& clauses);
    virtual ~UpperBoundProg();
  public:
    bool initialize(); // initialize the worker, returns true iff the instance is SAT
    void run();        // start the worker, initialize must be called first.
    virtual bool is_backbone(const Lit& literal) const;
    size_t  get_solver_calls() const { return solver_calls; }
    size_t  get_solver_time() const  { return solver_time; }
  private:// initial
    const ToolConfig&   tool_configuration;
    int                 chunk_size;
    ostream&            output;
    const Var           max_id;
    const CNF&          clauses;
    const Range         variable_range;
  private:// state
    LitBitSet           might_be; // literals that might still be backbones   
    LitBitSet           must_be;  // literals that must be backbones   
  private:// stats
    size_t              solver_calls;  // number of solver calls, for statistical purposes
    double              solver_time;   // for statistical purposes
  private:// sub-objects
    Solver              solver;
    Lifter              lifter;               // used to reduce models via lifting
    Rotatable           rotatable_computer;   // used to get rotatable variables
  private:// backbone testing
    typedef std::unordered_map<Lit,Lit,Lit_hash, Lit_equal>              Lit2Lit;
    Lit2Lit             lit_selectors;
    size_t              make_chunk(const Var chunk_start, const Var chunk_end, vec<Lit>& assumptions);
    void                process_model(const vec<lbool>& model);        // calls pruning techniques and debones based of the model
    void                process_pruned_model(const vec<lbool>& model); // debones everything not in the model
    inline bool         debone(Lit literal);                    // mark a literal as not a backbone
    inline Lit          get_selector(Lit literal);
    bool                run_solver(const vec<Lit>& assumptions);
    bool                run_solver();
  };

  inline bool UpperBoundProg::debone(Lit literal) { 
    assert(!must_be.get(literal));
    const bool rr = might_be.remove(literal);
    if (rr) {
      solver.addClause(~get_selector(~literal));
      UPPERBOUND_DBG( cerr << "dbcl: " << ~get_selector(~literal) << endl;);
    }
    return rr;
  }

  inline Lit UpperBoundProg::get_selector(Lit literal) { 
    const auto i = lit_selectors.find(literal);
    assert(i!=lit_selectors.end());
    return i->second;
  }

} /* end of namespace minibones */
#endif /* UPPERBOUNDPROG_HH_24653 */
