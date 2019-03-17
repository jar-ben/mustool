/* 
 * File:   Lifter.cc
 * Author: mikolas
 * 
 * Created on March 24, 2011, 6:41 PM
 */

#include "Lifter.hh"
#include "minisat_aux.hh"
#include "ClausesInversion.hh"
using Minisat::lit_Undef;
using Minisat::mkLit;

Lifter::Lifter(const CNF& _formula)
:formula(_formula) 
,inversion(formula)
{}

void Lifter::reduce(const vec<Lit>& assumptions,
        const vec<lbool>& input, vec<lbool>& output) {
    output.growTo(input.size(),l_Undef);
    // mark assumptions as fixed
    for (int index = 0; index < assumptions.size(); ++index) {
        const Lit& literal=assumptions[index];
        MR_LOG( cerr << "assumption " << literal << endl; )
        output[var(literal)] = sign(literal) ? l_False : l_True;
    }
    // find unit clauses
    FOR_EACH (clause_index,formula){         
         Lit satisfied_literal = lit_Undef;
         const bool unit = is_unit(*clause_index, input, satisfied_literal);
         assert(satisfied_literal!=lit_Undef);// all clauses must be satisfied
         if (unit) {
              MR_LOG( cerr << "unit lit " << satisfied_literal << endl; )
              const Var variable = var(satisfied_literal);
              output[variable]=input[variable];
         } 
    }
    //  things assigned in [output] until now will be assigned
    vector<bool> satisfied_clauses;
    satisfied_clauses.resize(formula.size(), false);
    // input literals go into the heap according to their occurrence
    LiteralScore* const scores = new LiteralScore[input.size()];
    vector<LiteralScore*> scores_pointers;//pointers to indivifual fields of [scores]
    scores_pointers.resize(input.size(), NULL);
    for (int index = 0; index<input.size(); ++index) {
        scores_pointers[index]=scores + index;
        const Var variable      = (Var)index;
        const Lit literal       = (input[index]==l_False) ? ~mkLit(index) : mkLit(index);
        const vector<size_t>* ixs  = inversion[literal];
        const size_t occurrences = (ixs==NULL) ? 0 : ixs->size();
        scores[variable].literal = literal;
        scores[variable].index() = 0;
        assert( occurrences < (2*formula.size()+1) );
        scores[variable].score = output[index]!=l_Undef ? 2*formula.size()+1 // push fixed literals  to the top
                                                        : occurrences;
        MR_LOG( cerr << "putting on heap " << literal << " with score " << (scores[variable].score) << endl; )
    }
    
    heap<LiteralScore*,size_t> occurrence_heap(scores_pointers);
    // satisfy literals, going from the most occurring ones
    // stop when all clauses satisfied
    while (occurrence_heap.size()>0 && occurrence_heap.maximum()->score>0) {
        LiteralScore* occurrence = occurrence_heap.extract_max();        
        MR_LOG( cerr << "satisfying literal " << occurrence->literal << " with score " << (occurrence->score) << endl; )
        if (occurrence->score==0) break;
        satisfy_literal(occurrence->literal, output, satisfied_clauses, scores, occurrence_heap);
    }
    delete[] scores;
}

Lifter::~Lifter() {}

bool Lifter::is_sat(const LitSet& clause, const vec<lbool>& model) const {
  FOR_EACH(literal_index, clause)
    if (is_satisfied(*literal_index,model)) return true;
  return false;
}

bool Lifter::is_unit(const LitSet& clause,
                     const vec<lbool>& model,
                     Lit& a_satisfied_literal) const {
  size_t counter = 0;
  a_satisfied_literal = lit_Undef;
  FOR_EACH(literal_index, clause) {
    const Lit literal = *literal_index;
    if (is_satisfied(literal,model)) {
      ++counter;
      a_satisfied_literal = literal;
    }
    if (counter > 1) break;
  }
  return counter==1;
}

void Lifter::satisfy_literal(const Lit literal,
                             vec<lbool>& output, vector<bool>& satisfied_clauses,
                             LiteralScore* occurrences,
                             heap<LiteralScore*,size_t>& occurrence_heap) {

  MR_LOG( cerr << "satisfying " << literal << endl; )
    const Var variable = var(literal);
  const lbool lit_val = sign(literal) ? l_False : l_True;
  assert(output[variable]==l_Undef || output[variable]==lit_val);
  output[variable] = lit_val;//copy to output
  const vector<size_t>& occurrence_indexes = *inversion[literal];//get literal's occurrences
  if (inversion[literal]==NULL) {// literal not occurring, nothing to do
    return;
  }

  //make sure that all clauses containing the literal are satisfied
  FOR_EACH(clause_index,occurrence_indexes) {
    if (satisfied_clauses[*clause_index]) continue;
    satisfied_clauses[*clause_index]=true;
    // decrease occurrences for all the contained literals
    FOR_EACH(literal_index, formula[*clause_index]) {
      const Lit lit = *literal_index;
      const Var vid = var(lit);
      LiteralScore* lo = occurrences + vid;
      assert(var(lo->literal)==vid);
      if ((lo->literal!=lit) || !occurrence_heap.contains(lo)) continue;
      assert(lo->score > 0);
      occurrence_heap.update_key(lo, (lo->score)-1);
    }
  }
}

bool Lifter::is_model(const vec<lbool>& model) const {
  FOR_EACH (clause_index,formula) 
    if (!is_sat(*clause_index, model)) return false;
  return true;
}
