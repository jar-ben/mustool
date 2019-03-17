/*
 * File:  UpperBound.hh
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 3:44:59 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef UPPERBOUND_HH_28341
#define UPPERBOUND_HH_28341
#include <vector>
#include <functional>
#include "MiniSatExt.hh"
#include "auxiliary.hh"
#include "Lifter.hh"
#include "LitBitSet.hh"
#include "ToolConfig.hh"
#include "Rotatable.hh"
#include "BackboneInformation.hh"
#include "Cores.hh"
namespace minibones {
  using Minisat::MiniSatExt;
  using Minisat::Var;
  using Minisat::Lit;
  using Minisat::vec;

  /** Class for computing backbones with the upper bound approach. */
  class UpperBound : public BackboneInformation {
  public:
    UpperBound(ToolConfig& tool_configuration, ostream& output, Var max_id, const CNF& clauses);
    virtual ~UpperBound();
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
    MiniSatExt          solver;
    Lifter              lifter;               // used to reduce models via lifting
    Rotatable           rotatable_computer;   // used to get rotatable variables
    Cores           cores;   // used to get rotatable variables
  private:// backbone testing
    const_infinite_LitBitSetIterator might_be_iterator; // = might_be.infinite_iterator();
    Var  end_of_chunk;
    void process_model(const vec<lbool>& model);        // calls pruning techniques and debones based of the model
    void process_pruned_model(const vec<lbool>& model); // debones everything not in the model
    Lit  make_chunk(vec<Lit>& literals);                // prepares a chunk and returns the relaxation literal for it
    inline bool debone(Lit literal);                    // mark a literal as not a backbone
    bool run_solver(const vec<Lit>& assumptions);
    bool run_solver();
  };

  inline bool UpperBound::debone(Lit literal) { 
    assert(!must_be.get(literal));
    return might_be.remove(literal);
  }
} /* end of namespace minibones */
#endif /* UPPERBOUND_HH_28341 */
