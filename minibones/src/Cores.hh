/*
 * File:  Cores.hh
 * Author:  mikolas
 * Created on:  Sat, Mar 17, 2012 5:54:13 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef CORES_HH_13827
#define CORES_HH_13827
#include <vector>
#include <functional>
#include "MiniSatExt.hh"
#include "auxiliary.hh"
#include "Lifter.hh"
#include "LitBitSet.hh"
#include "ToolConfig.hh"
#include "Rotatable.hh"
#include "BackboneInformation.hh"
namespace minibones {
  using Minisat::MiniSatExt;
  using Minisat::Var;
  using Minisat::Lit;
  using Minisat::vec;

  /** Class for computing backbones with the upper bound approach. */
  class Cores {
  public:
    Cores(ToolConfig& tool_configuration, Var max_id, const CNF& clauses);
    virtual ~Cores();
  public:
    bool try_to_flip(const LitBitSet& might_be, vec<Lit>& reasons);
  private:
    const ToolConfig&   tool_configuration;
    const Var           max_id;
    const CNF&          clauses;
  private:// sub-objects
    MiniSatExt          solver;
  };
}
#endif /* CORES_HH_13827 */
