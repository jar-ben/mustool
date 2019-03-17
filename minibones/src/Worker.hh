/*
 * File:   Worker.hh
 * Author: mikolas
 *
 * Created on March 19, 2011, 4:51 PM
 */
#ifndef WORKER_HH
#define	WORKER_HH
#include <vector>
#include <functional>

#include "auxiliary.hh"
#include "Lifter.hh"
#include "BBInfo.hh"
#include "ToolConfig.hh"
#include "Rotatable.hh"
#include "BackboneInformation.hh"
#include "minisat/core/Solver.h"

using Minisat::Solver;
using Minisat::Lit;
using Minisat::lbool;
using Minisat::Var;
using Minisat::toInt;
using Minisat::sign;
using Minisat::vec;
using std::vector;
using std::pair;

namespace minibones {
  /**
   * Class for computing backbones.
   */
  class Worker : public BackboneInformation {
  public:
    Worker(ToolConfig& tool_configuration,
           ostream& output,
           Var max_id, const CNF& clauses,
           const Range& variable_range);
    /**unsupported*/
    //Worker(const Worker& orig);
    virtual ~Worker();
  public:
    /** Initialize the worker, returns true iff the instance is SAT. */
    bool initialize();
    /**Start the worker, Run only after  {@code initialize} is called */
    void run();
    virtual bool is_backbone(const Lit& literal) const;
    size_t  get_solver_time() const  { return solver_time; }
    size_t  get_solver_calls() const { return solver_calls; }
  private:// initial
    const ToolConfig&   tool_configuration;
    ostream&            output;
    const Var           max_id;
    const CNF&          clauses;
    const Range&        variable_range;
  private:// state
    double       solver_time;   // for statistical purposes
    size_t       solver_calls;  // number of solver calls, for statistical purposes
    vector<bool> to_test;       // to be tested whether they are backbones or not
    size_t       to_test_count; // number of literals still to be tested
    BBInfo       bbInfo;
  private:// sub-objects
    Solver      solver;
    Lifter      lifter;      // used to reduce models
    Rotatable   rotatable_computer;      // used to get rotatable variables
  private:// backbone testing
    /** Tests if the given literal is a backbone or not.
     * @param literal tested literal
     * @return true iff {@code literal} is a backbone */
    bool test_backbone(const Lit& literal);
    /** Debones everything not in the model.*/
    void process_model(const vec<lbool>& model);
    inline bool debone(const Lit& literal);
    inline bool mark_bone(const Lit& literal);
    bool run_solver(const vec<Lit>& assumptions);
    bool run_solver();
  private:// running
    bool should_stop() const;
  private:// debugging
    bool is_complete () const;
  };

  bool Worker::debone(const Lit& literal) {
    bool return_value = bbInfo.debone(literal);
    if (return_value && to_test[(size_t)var(literal)])  {
      assert (to_test_count>0);
      --to_test_count;
    }
    return return_value;
  }

  bool Worker::mark_bone(const Lit& literal) {
    bool return_value =  bbInfo.mark_bone(literal);
    if (return_value && to_test[(size_t)var(literal)])  {
      assert (to_test_count>0);
      --to_test_count;
    }
    return return_value;
  }

} /* namespace minibones */
#endif	/* WORKER_HH */
