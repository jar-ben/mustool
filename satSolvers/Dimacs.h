#ifndef SatSolver_Dimacs_h
#define SatSolver_Dimacs_h

#include <stdio.h>
#include <vector>
#include <iostream>
#include "core/misc.h"
#include "satSolvers/ParseUtils.h"
#include <string>

using namespace ReMUS; 

void UreadClause(StreamBuffer& in, std::vector<int>& lits);
void UreadGroupClause(StreamBuffer& in, std::vector<int>& lits);
void Uparse_DIMACS_main(StreamBuffer& in, std::vector<std::vector<int>>& cls, bool gcnf);
void Uparse_DIMACS(std::string filename, std::vector<std::vector<int>>& cls);
#endif

