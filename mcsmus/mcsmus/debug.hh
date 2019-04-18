/***********[debug.hh]
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

#ifndef MCSMUS_DEBUG_HH
#define MCSMUS_DEBUG_HH

/* macros for debugging

 we can write

 DOUT << "hello\n";
 DBG {
   // do some expensive debugging stuff here
 }

 and this code will be removed if MCSMUS_DEBUG is not defined (compiled
 but 0 runtime cost)
*/

#if defined(MCSMUS_DEBUG)
#define DOUT cout
#define DBG
#else
#define DOUT if(1); else std::cout
#define DBG if(1); else
#endif

#define PV(x) "("#x")=" << (x)

#endif
