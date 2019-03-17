/*
 * File:  Rotatable.hh
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 11:33:58 AM
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef ROTATABLE_HH_32718
#define ROTATABLE_HH_32718
#include "types.hh"
#include "VarSet.hh"
class Rotatable {
public:
  Rotatable(const CNF& input_cnf);
  void get_rotatables(const vec<lbool>& model, VarSet& rotatables);
private:
  const CNF& input_cnf;
};
#endif /* ROTATABLE_HH_32718 */
