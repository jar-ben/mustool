This is a so far unnamed tool for online enumeration of minimal unsatisfiable subsets (MUSes) of a given unsatisfiable set of constraints. The tool currently implements three online MUS enumeration algorithms: MARCO[1], TOME[2], and ReMUS[3], and supports MUS enumeration in three constraint domains: SAT, SMT, and LTL.


We distribute this source code under the MIT licence. See the file ./COPYRIGHTS for mode details.

In case of any troubles, do not hesitate to contact me, Jaroslav Bendik, the developer of the tool, at xbendik=at=fi.muni.cz.

## Third-party Tools

Besides our original code, the tool uses several external tools that we in some form distribute with our code. Therefore, by using our tool, you are also using the third parties tools and it is your responsibility to obey the lincenses of these tools. In particular, our tool builds on the following third-party tools:
* miniSAT[4] - miniSAT is a well-known SAT CNF solver. We use it in our tool as a satisfiability solver int the SAT domain and we also use it as a solver for maintaining the symbolic representation of the unexplored subsets during the MUS enumeration.
* muser2[5] - muser2 is a SAT CNF single MUS extraction tool and we use it for an implementation of the shrinking procedure in the three MUS enumeration algorithms. We distribute only a linux complied binary of muser2 (./muser2-para) alongside our tool. The source code of muser2 is available at https://bitbucket.org/anton_belov/muser2. We recommend you, the user of our tool, download the source code of muser2 and build it on your own to make sure that you have a build that is suitable for your architecture.
* z3[6] - we use the SMT satisfiability solver z3 as a satisfiability solver in the SMT domain. However, we do not distribute the code of z3 (and nor a binary); we include it as a C++ library and you have to download and install it on your own.
* SPOT[7] - SPOT is a C++ library for LTL, omega-automata manipulation, and LTL model checking. We use SPOT as a satisfiability solver in the LTL domain. You have to download and install the library on your own.
* MiniBones[8] - MiniBones is a tool for computing complete bakcbones of a given satisfiable SAT CNF formula. We use it as a subroutine in our heuristic called Backbone heuristic; the heuristic can be optionally used to speed up a MUS enumeration in the SAT domain. Unfortunately,
the license of MiniBones is currently not precisely stated and thus we do not distribute it with our tool.
You can download its binary from its author's webpage http://sat.inesc-id.pt/~mikolas/sw/minibones/. Name it minibones_binary and put it into this directory. Even a better option is to contact the author of MiniBones and ask him for a source code so that you can use MiniBones directly (it is faster). We provide a wrapper, Bones.h and Bones.cpp, for running MiniBones from source code. See Bones.h for the folder structure that needs to be set up to incorporate MiniBones from source code. Also, you have to change the MINIBONESTYPE flag in Makefile.

## Installation
If you have already installed the above mentioned tools, just run "make"

## Running the Tool
To run the tool use: ./mvc -v 1 <path_to_input_file>

The -v argument determines the algorithm that should be used. In particular, 1 - ReMUS, 2 - TOME, 3 - Marco. The name of the input file should end either as .cnf, .smt, or .ltl. See the example inputs in ./examples/

For more options, invoke ./mvc -h.

For trying our new heuristics, Backbone and Matchmaker, for in the SAT domain, use flag "-w 1" or "-w 2" for Backbone and Matchmaker, respectively.

## Final Note
The tool is still under development and there is a lot of things that need to be done (or improved). In case of any troubles, please contact me at xbendik=at=fi.muni.cz.

* [1] Mark H. Liffiton, Alessandro Previti, Ammar Malik, João Marques-Silva: Fast, flexible MUS enumeration. Constraints 21(2), 2016.
* [2] Jaroslav Bendík, Nikola Beneš, Ivana Černá, Jiří Barnat: Tunable Online MUS/MSS Enumeration. FSTTCS 2016.
* [3] Jaroslav Bendík, Ivana Černá, Nikola Beneš: Recursive Online Enumeration of All Minimal Unsatisfiable Subsets. ATVA 2018.
* [4] http://minisat.se/
* [5] https://bitbucket.org/anton_belov/muser2
* [6] https://github.com/Z3Prover/z3
* [7] https://spot.lrde.epita.fr/
* [8] http://sat.inesc-id.pt/~mikolas/sw/minibones/
