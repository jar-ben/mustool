/*
 * File:  BackboneInformation.cc
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 4:59:15 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#include "BackboneInformation.hh"
using namespace minibones;

bool BackboneInformation::is_backbone(Var var) const  {return is_backbone(~mkLit(var)) || is_backbone(mkLit(var));}

bool BackboneInformation::backbone_sign(Var var) const {
  assert(is_backbone(var));
  return is_backbone(mkLit(var));
}
