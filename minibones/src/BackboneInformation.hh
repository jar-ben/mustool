/*
 * File:  BackboneInformation.hh
 * Author:  mikolas
 * Created on:  Tue, Feb 21, 2012 4:38:23 PM
 * Copyright (C) 2012, Mikolas Janota
 */
#ifndef BACKBONEINFORMATION_HH_13096
#define BACKBONEINFORMATION_HH_13096
#include "minisat_aux.hh"
namespace minibones {
  class BackboneInformation {
  public:
    virtual bool is_backbone(const Lit& literal) const =0;
    virtual bool is_backbone(Var var) const;
    virtual bool backbone_sign(Var var) const;
  public:
    BackboneInformation() {}
    virtual ~BackboneInformation() {}
  };
} /* end of namespace minibones */
#endif /* BACKBONEINFORMATION_HH_13096 */
