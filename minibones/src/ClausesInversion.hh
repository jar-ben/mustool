/* 
 * File:   ClausesInversion.hh
 * Author: mikolas
 *
 * Created on April 1, 2011, 5:21 PM
 */
#ifndef CLAUSESINVERSION_HH
#define	CLAUSESINVERSION_HH
#include "minisat_aux.hh"
#include "types.hh"
using std::vector;
/** Represents a given formula as a map from literals to
    occurrences.*/
class ClausesInversion {
public:
  typedef vector<const std::vector<size_t>*>::const_iterator const_iterator;
  ClausesInversion(const CNF& source);
  ClausesInversion(const ClausesInversion& orig);
  virtual ~ClausesInversion();
  void invert(const CNF& source);
  inline size_t lit_index(const Lit& l) const { return sign(l)? var(l)<<1 : (var(l)<<1) + 1; }
  const vector<size_t>* operator []  (Lit l) const { return get(l); };
  const vector<size_t>* get(Lit l) const;
  void clear();
private:
  std::vector<std::vector< size_t >*> inversion;
  inline vector<size_t>* _get(Lit l) const;
  inline bool put(Lit l, vector<size_t>* c);
};

inline const vector<size_t>* ClausesInversion::get(Lit l) const { 
  const size_t i = lit_index(l);        
  return i<inversion.size() ? inversion[i] : NULL;
}

inline vector<size_t>* ClausesInversion::_get(Lit l) const {
  const size_t i = lit_index(l);
  return i<inversion.size() ? inversion[lit_index(l)] : NULL;
}

inline bool ClausesInversion::put(Lit l, vector<size_t>* c) {
  const size_t i = lit_index(l);
  if (inversion.size() <= i) inversion.resize(i+1, NULL);
  const vector<size_t>* const o=inversion[i];
  inversion[i]=c;
  return o==NULL;
}
#endif	/* CLAUSESINVERSION_HH */

