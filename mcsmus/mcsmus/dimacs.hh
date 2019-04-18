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

#ifndef MCSMUS_DIMACS_HH
#define MCSMUS_DIMACS_HH

#include <stdio.h>
#include <vector>

#include "mcsmus/parseutils.hh"

namespace mcsmus {

//=================================================================================================
// DIMACS Parser:

template<class B>
void readClause(B& in, std::vector<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        // GK: newVar() will be called by addClause, if needed
        //while (var >= S.nVars()) S.newVar();
        lits.push_back( mkLit(var, parsed_lit<0) );
    }
}

template<class B, class Cnf>
void parse_DIMACS_main(B& in, Cnf& S, bool strictp = false) {
    std::vector<Lit> lits;
    int clauses = 0;
    int cnt     = 0;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
            if (eagerMatch(in, "p cnf")){
                /*vars=*/ parseInt(in);
                clauses = parseInt(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
            }else{
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else{
            cnt++;
            readClause(in, lits);
            S.addClause(lits, cnt); // everything is soft in a cnf
        }
    }
    if (strictp && cnt != clauses)
        printf("PARSE ERROR! DIMACS header mismatch: wrong number of clauses\n");
}

// Inserts problem into solver.
//
void parse_DIMACS(gzFile input_stream, Wcnf& S, bool strictp = false) {
    StreamBuffer in(input_stream);
    parse_DIMACS_main(in, S, strictp); }

//=================================================================================================
// group CNF Parser:

template <class B> int readGroupClause(B& in, std::vector<Lit>& lits)
{
    int parsed_lit, var;
    int group;
    lits.clear();
    if (!eagerMatch(in, "{")) {
        printf("PARSE ERROR! Expected {groupid}, got char: %c\n", *in);
        exit(3);
    }
    group = parseInt(in);
    if (!eagerMatch(in, "}")) {
        printf("PARSE ERROR! Unexpected char: %c\n", *in);
        exit(3);
    }
    for (;;) {
        parsed_lit = parseInt(in);
        if (parsed_lit == 0)
            break;
        var = abs(parsed_lit) - 1;
        lits.push_back(mkLit(var, parsed_lit < 0));
    }
    return group;
}

template <class B, class Cnf>
void parse_GCNF_main(B& in, Cnf& S, bool strictp = false)
{
    std::vector<Lit> lits;
    int clauses = 0;
    int groups = 0;
    int cnt = 0;
    for (;;) {
        skipWhitespace(in);
        if (*in == EOF)
            break;
        else if (*in == 'p') {
            if (eagerMatch(in, "p gcnf")) {
                /*vars=*/parseInt(in);
                clauses = parseInt(in);
                groups = parseInt(in);
            } else {
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else {
            cnt++;
            int group = readGroupClause(in, lits);
            if (strictp && group > groups)
                printf("PARSE ERROR! GCNF header mismatch: clause %d has group "
                       "%d, %d groups declared\n",
                    cnt, group, groups);
            S.addClause(lits, group);
        }
    }
    if (strictp && cnt != clauses)
        printf(
            "PARSE ERROR! DIMACS header mismatch: wrong number of clauses\n");
}

// Read group-CNF file
void parse_GCNF(gzFile input_stream, Wcnf& S, bool strictp)
{
    StreamBuffer in(input_stream);
    parse_GCNF_main(in, S, strictp);
}

//=================================================================================================
// auto-detect type of input and read the correct thing
enum class input_type { CNF_FILE, GCNF_FILE, UNKNOWN };

template <class B> input_type determineFileType(B& in)
{
    for (;;) {
        skipWhitespace(in);
        if (*in == EOF)
            break;
        else if (*in == 'p') {
            if( lookaheadMatch(in, "p cnf") )
                return input_type::CNF_FILE;
            else if( lookaheadMatch(in, "p gcnf") )
                return input_type::GCNF_FILE;
            else
                return input_type::UNKNOWN;
        } else
            return input_type::UNKNOWN;
    }
    return input_type::UNKNOWN;
}

input_type parse_auto_detect(gzFile input_stream, Wcnf& S, bool strictp)
{
    StreamBuffer in(input_stream);
    auto filetype = determineFileType(in);
    if( filetype == input_type::CNF_FILE )
        parse_DIMACS_main(in, S, strictp);
    else if( filetype == input_type::GCNF_FILE )
        parse_GCNF_main(in, S, strictp);
    else {
        printf("PARSE ERROR! Could not detect input file type\n");
        exit(3);
    }
    return filetype;
}

//=================================================================================================
}

#endif
