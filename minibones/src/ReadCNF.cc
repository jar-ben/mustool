/* 
 * File:   ReadCNF.cc
 * Author: mikolas
 * 
 * Created on November 11, 2010, 1:19 PM
 */

#include "ReadCNF.hh"
#include "Reader.hh"
#include "minisat_aux.hh"
#include <vector>
#include <string>
#include <assert.h>
using Minisat::lit_Undef;

ReadCNF::ReadCNF(Reader& input_file)
: in(input_file), mxid(0)
{}

ReadCNF::ReadCNF(const ReadCNF& orig)
: in(orig.in), mxid(orig.mxid)
{assert(false);}

ReadCNF::~ReadCNF() {}

void ReadCNF::read_cnf_clause(vector<Lit>& lits) {
  Lit parsed_lit;
  lits.clear();
  for (;;){
    in.skip_whitespace();
    parsed_lit = parse_lit(in);
    if (parsed_lit == lit_Undef) break;
    const Var v = var(parsed_lit);
    if (v > mxid) { mxid = v; }
    lits.push_back(parsed_lit);
  }
}

void ReadCNF::read() {
  vector<Lit> literals;
  for (;;){
    in.skip_whitespace();
    if (*in == EOF)
      break;
    else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else {
      literals.clear();
      read_cnf_clause(literals);
      clause_vector.push_back(LitSet(literals));
    }
  }
}

Lit ReadCNF::parse_lit(Reader& in) {
    Var  v = 0;
    bool neg = false;
    if (*in=='-') { neg=true; ++in; }
    else if (*in == '+') ++in;
    if ((*in < '0') || (*in > '9')) {
        string s("unexpected char in place of a literal: ");  s+=*in;
        throw ReadCNFException(s);
    }
    while (*in >= '0' && *in <= '9') {
        v = v*10 + (*in - '0');
        ++in;
    }
    return v ? mkLit(v, neg) : lit_Undef;
}
