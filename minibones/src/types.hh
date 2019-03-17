/*
 * File:  types.hh
 * Author:  mikolas
 * Created on:  Fri, Feb 17, 2012 10:59:33 AM
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef TYPES_HH_5674
#define TYPES_HH_5674
#include <iostream>
#include <vector>
#include <utility>
#include "LitSet.hh"
using std::ostream;
using std::pair;
typedef vector<LitSet>                        CNF;
typedef pair<Var,Var>                         Range;
typedef vector<Var>                           VariableVector;

ostream & operator << (ostream& outs, const CNF& f);
#endif /* TYPES_HH_5674 */
