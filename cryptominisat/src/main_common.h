/*************************************************************
MiniSat       --- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
CryptoMiniSat --- Copyright (c) 2014, Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***************************************************************/

#ifndef __MAIN_COMMON_H__
#define __MAIN_COMMON_H__

#include "cryptominisat5/cryptominisat.h"
#include "solverconf.h"
#include <iostream>
#include <cmath>

class MainCommon
{
public:
    uint32_t print_model(CMSat::SATSolver* solver,
                         std::ostream* os,
                         std::vector<uint32_t>* only = NULL);
    void handle_drat_option();

    string dratfilname;
    bool dratDebug = false;
    std::ostream* dratf = NULL;
    bool zero_exit_status = false;
    CMSat::SolverConf conf;
    unsigned num_threads = 1;
};

#endif //__MAIN_COMMON_H__
