/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef ReMUS_Dimacs_h
#define ReMUS_Dimacs_h

#include <stdio.h>
#include <vector>
#include <iostream>
#include "satSolvers/ParseUtils.h"
#include <string>
#include "core/misc.h"

namespace ReMUS{

//=================================================================================================
// DIMACS Parser:

bool ends_with(std::string const & value, std::string const & ending) {
    	if (ending.size() > value.size()) return false;
	return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}


static void readClause(StreamBuffer& in, std::vector<int>& lits) {
	int parsed_lit, var;
	lits.clear();
    	for (;;){
		parsed_lit = parseInt(in);
		if (parsed_lit == 0) break;
		var = abs(parsed_lit)-1;
		lits.push_back(parsed_lit);
	}
}

void readGroupClause(StreamBuffer& in, std::vector<int>& lits)
{
    	int parsed_lit, var;
    lits.clear();
    if (!eagerMatch(in, "{")) {
        printf("PARSE ERROR! Expected {groupid}, got char: %c\n", *in);
        exit(3);
    }
    parseInt(in);
    if (!eagerMatch(in, "}")) {
        printf("PARSE ERROR! Unexpected char: %c\n", *in);
        exit(3);
    }
    for (;;) {
        parsed_lit = parseInt(in);
        if (parsed_lit == 0)
            break;
        var = abs(parsed_lit) - 1;
	lits.push_back(parsed_lit);
    }
    return;
}

static void parse_DIMACS_main(StreamBuffer& in, std::vector<std::vector<int>>& cls, bool gcnf) {
	int vars    = 0;
   	int clauses = 0;
    	int cnt     = 0;
	std::vector<int> clause;
    	for (;;){
        	skipWhitespace(in);
		if (*in == EOF) break;
		else if (*in == 'p'){
		    if ((gcnf &&  eagerMatch(in, "p gcnf")) || (!gcnf &&  eagerMatch(in, "p cnf"))){
			vars    = parseInt(in);
			clauses = parseInt(in);
		   	if(gcnf) parseInt(in); //the number of groups
		    }else{
			printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
		    }
		} else if (*in == 'c' || *in == 'p')
		    skipLine(in);
		else{
		    cnt++;
		    if(gcnf)
			readGroupClause(in, clause);
		    else
		    	readClause(in, clause);
		    cls.push_back(clause);
		    }
	    }
    if (cnt  != cls.size())
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");
}

// Inserts problem into solver.
//
static void parse_DIMACS(std::string filename, std::vector<std::vector<int>>& cls) {
     	gzFile in = gzopen(filename.c_str(), "rb");

	bool group = ends_with(filename, "gcnf");
	if(in == NULL)
		print_err("could not open the file " + filename);
    	StreamBuffer ins(in);
       	parse_DIMACS_main(ins, cls, group); 
	gzclose_r(in);
}

//=================================================================================================
}
#endif

