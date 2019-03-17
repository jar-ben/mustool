/* 
 * File:   ClausesInversion.cc
 * Author: mikolas
 * 
 * Created on April 1, 2011, 5:21 PM
 */

#include "ClausesInversion.hh"

ClausesInversion::ClausesInversion(const CNF& source) {invert(source);}

ClausesInversion::ClausesInversion(const ClausesInversion&) {assert(0);}

ClausesInversion::~ClausesInversion() { clear(); }

void ClausesInversion::clear() {
  for (size_t i=0; i<inversion.size(); ++i) {
    if (inversion[i]==NULL) continue;
    delete inversion[i];
  }
  inversion.clear();
}

void ClausesInversion::invert(const CNF& source) {
  for (size_t clause_index=0; clause_index<source.size();++clause_index) {
    const LitSet& clause = source[clause_index];
    FOR_EACH(literal_index, clause) {
      const Lit& literal = *literal_index;
      vector<size_t>* index_list = _get(literal);
      if (index_list == NULL) {
        index_list = new vector<size_t>();
        put(literal,index_list);
      }
      index_list->push_back(clause_index);
    }
  }
}

