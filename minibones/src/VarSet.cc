#include "auxiliary.hh"
#include "VarSet.hh"

VarSet::VarSet() 
  : _size(0)
{}

VarSet::VarSet(const VariableVector& variables)
: _size(0) 
{ FOR_EACH(vi, variables) add(*vi); }

const_VarIterator::const_VarIterator(const VarSet& ls, size_t x) : ls(ls), i(x) {
  while ((i<ls.physical_size()) && !ls.get(i)) ++i;
}

