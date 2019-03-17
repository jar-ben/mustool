/*
 * File:  LitBitSet.cc
 * Author:  mikolas
 * Created on:  Mon, Feb 20, 2012 11:47:06 AM
 * Copyright (C) 2012, Mikolas Janota
 */
#include "LitBitSet.hh"

ostream& LitBitSet::print(ostream& out) {
  for (size_t index = 0; index < set.size(); ++index) {
    if (!set[index]) continue;
    if (index) out << " ";
    out << index2literal(index);
  }
  return out;
}

const_infinite_LitBitSetIterator::const_infinite_LitBitSetIterator(const LitBitSet& ls, size_t x) : ls(ls), i(x) {
  if (!ls.size()) return;
  while (i<2 || !ls.get(index2literal(i))) i=(i+1) % ls.physical_size();
  assert (ls.get(index2literal(i)));
}

const_LitBitSetIterator::const_LitBitSetIterator(const LitBitSet& ls, size_t x) : ls(ls), i(x) {
  while ((i<ls.physical_size()) && !ls.get(index2literal(i))) ++i;
}
