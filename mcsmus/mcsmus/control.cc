/***********[control.cc]
Copyright (c) 2015 George Katsirelos

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********/

#include "mcsmus/control.hh"
#include "minisat/utils/mcsmus_System.h"
#include <functional>
#include <cstdlib>
#include <cmath>
#include <signal.h>

namespace {
    std::function<void()> set_interrupt_function;
    std::function<void()> interrupt_function;

    void status_alarm(int)
    {
        mcsmus::getGlobalControl()->status_timeout = true;
    }
}

// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int) { set_interrupt_function(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int) {
    interrupt_function();
    std::exit(1);
}

mcsmus::Control::Control()
{
    signal(SIGALRM, status_alarm);

    set_interrupt_function = [&]() { interrupt(); };
    interrupt_function = [&]() {
        printf("\n"); printf("*** INTERRUPTED ***\n");
    };

    Minisat::sigTerm(SIGINT_interrupt);
}

mcsmus::Control::~Control()
{
    Minisat::sigTerm(SIG_DFL);
}

void mcsmus::Control::force_exit_on_interrupt()
{
    Minisat::sigTerm(SIGINT_exit);
}

void mcsmus::Control::notify_on_interrupt()
{
    Minisat::sigTerm(SIGINT_interrupt);
}

void mcsmus::Control::interrupt()
{
    asynch_interrupt = true;
}

void mcsmus::Control::checkPrintStatusLine()
{
    if( status_timeout ) {
        if( fn_printStatusLine && verbosity >= 1 )
            fn_printStatusLine();
        status_timeout = false;
        if( std::round(Minisat::cpuTime()) < status_interval )
            alarm(ceil(status_interval/10.0));
        else
            alarm(status_interval);
    }
}

mcsmus::Control* mcsmus::getGlobalControl()
{
    static Control control;
    return &control;
}
