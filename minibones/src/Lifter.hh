/* 
 * File:   Lifter.hh
 * Author: mikolas
 *
 * Created on March 24, 2011, 6:41 PM
 */

#ifndef MODELREDUCE_HH
#define	MODELREDUCE_HH
#include <vector>
#include "minisat/core/SolverTypes.h"
#include "auxiliary.hh"
#include "heap.hh"
#include "ClausesInversion.hh"
using std::vector;
using Minisat::lbool;
using Minisat::sign;
using Minisat::vec;
using Minisat::Lit;
using Minisat::lit_Undef;

#define MR_LOG( t )  

struct LiteralScore;//auxiliary data structure used to keep literals in a heap, defined later on

class Lifter {
public:
  /** Initialize the object with a formula that is to be operated on.
   *@param formula CNF formula, reference is kept, so must not be modified*/
  Lifter(const CNF& formula);
  virtual ~Lifter();
public:
  /**
   * Remove literals from a given assignment while preserving satisfiability of
   * the given formula.
   * @param input variable assignment which must satisfy the given formula
   * @param output the input assignment where some literals are replaced with the
   *               press l_Undef value
   */
  void reduce(const vec<Lit>& assumptions, const vec<lbool>& input, vec <lbool>& output);
  bool is_model(const vec<lbool>& model) const;
private:
  const CNF&                formula;
  const ClausesInversion    inversion; // lits->clauses map of formula
  inline bool is_satisfied(const Lit literal, lbool value) const;
  inline bool is_satisfied(const Lit literal, const vec<lbool>& model) const;
  bool is_sat(const LitSet& clause,const vec<lbool>& model) const;
  bool is_unit(const LitSet& clause,
               const vec<lbool>& model,
               Lit& a_satisfied_literal) const;
  void satisfy_literal(const Lit literal,vec<lbool>& output, vector<bool>& satisfied_clauses,
                       LiteralScore* occurrences,
                       heap<LiteralScore*,size_t>& occurrence_heap);
public:
  Lifter(const Lifter& orig);
};

inline bool Lifter::is_satisfied(const Lit literal, const vec<lbool>& model) const {
  const Var variable = var(literal);
  MR_LOG( if (variable>=model.size()) print_model(cerr << "error in is_satisfied: " << variable << " not in ", model) << endl; )
  assert(variable<model.size());
  return is_satisfied(literal, model[variable]);
}

inline bool Lifter::is_satisfied(const Lit literal, lbool value) const {
  assert(literal!=lit_Undef);
  if (sign(literal)  && (value==l_False)) return true;
  if (!sign(literal) && (value==l_True))  return true;
  return false;
}

struct LiteralScore {
  Lit       literal;
  size_t    score;
  size_t    _index;
  inline LiteralScore(const LiteralScore& orig) : literal(orig.literal),score(orig.score),_index(orig._index) {}
  inline LiteralScore() : literal(Minisat::lit_Undef), score(0), _index(0) {}
  inline LiteralScore(Lit _literal, size_t _occurrence, size_t __index)
    : literal(_literal), score(_occurrence), _index(__index) {}
  inline size_t& key()   { return score; }
  inline size_t& index() { return _index; }

  friend ostream& operator << (ostream& output, const LiteralScore& score) {
    return output <<"[l:"<<score.literal<< ",s" << score.score << ",i" <<  score._index << "]";
  }
};
#endif	/* MODELREDUCE_HH */

