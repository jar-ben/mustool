/***********[control.hh]
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

#ifndef MCSMUS_CONTROL_H
#define MCSMUS_CONTROL_H

#include <functional>

namespace mcsmus {

struct Control;
Control *getGlobalControl();

/* There should be a single instance of this, shared among all the
   solvers/musers/whatevers */
struct Control {
    std::function<void()> fn_printStats;
    std::function<void(double)> fn_printHeader;
    std::function<void()> fn_printFooter;
    std::function<void()> fn_printStatusLine;

    bool status_timeout{true};
    int status_interval{10};
    bool asynch_interrupt{false};

    int verbosity{0};

    void interrupt();


    // Use signal handlers that forcibly quit. Used in cases where the
    // solver is unable to respond (mostly parsing)
    void force_exit_on_interrupt();

    // Use signal handlers that notify the solver to quit.
    void notify_on_interrupt();

    void checkPrintStatusLine();
protected:
    Control();
    virtual ~Control();
    friend Control *getGlobalControl();
};

}

#endif
