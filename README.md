# MUST
MUST is a  tool for online enumeration of minimal unsatisfiable subsets (MUSes) of a given unsatisfiable set of constraints. The tool currently implements three online MUS enumeration algorithms: MARCO [1], TOME [2], and ReMUS [3], and supports MUS enumeration in three constraint domains: SAT, SMT, and LTL.

We distribute this source code under the MIT licence. See ./LICENSE for mode details.

## Installation

To be able to deal with the three constraint domains, our tool employs several external libraries. In particular, we use:
- miniSAT [4] and zlib library for dealing with the SAT domain
- z3 [6] for dealing with the SMT domain
- and SPOT [7] for dealing with the LTL domain.

MiniSAT is packed directly with the source code of our tool; you do not install it separately. Zlib has to be installed. Finally, z3 and/or SPOT need to be installed only if you want to use our tool in the SMT and/or LTL domains. Note that installation of z3 and SPOT might take several hours.


### Installation of zlib
Zlib is a part of the package zlib1g-dev, you should be able to install it with:
```
sudo apt install zlib1g-dev
```

### Installation of z3
z3 can be downloaded at [https://github.com/Z3Prover/z3](https://github.com/Z3Prover/z3)
Please follow its README file for installation instructions. Basically, the following should do the trick:
```
python scripts/mk_make.py
cd build; make
sudo make install
```


### Installation of SPOT
SPOT can be downloaded from [https://spot.lrde.epita.fr/](https://spot.lrde.epita.fr/)
Again, follow installation instructions that are provided by its authors. 
Basically, the following should do the trick:
```
./configure --disable-python
make
sudo make install
sudo ldconfig
```
Note the "sudo ldconfig", that is important. Do not forget to run it, otherwise SPOT might not be visible for our tool.


### Building our tool
Building our tool requires at least zlib to be installed.
To build our tool, go to its main directory (the one where is this README.md file) and run:
```
make
```
This will make our tool only with support for the SAT domain. If you want to support also the other domains, you have to install z3 and/or SPOT first. Then run one of the following:
```
make USESMT=YES
make USELTL=YES
make USESMT=YES USELTL=YES
```
If you run make repeatedly, e.g. if you have decided to allow another domain, run first "make clean" and then appropriate make.


## Running our tool
in the main directory, run "./must _file_" where _file_ is an input file of constraints. You can use one of our examples, e.g.:
```
./must examples/test.cnf
./must examples/test.smt2
./must examples/test.ltl
./must examples/bf1355-228.cnf
```
To run the tool with a time limit (always recommended), use e.g.:
```
timeout 20 examples/bf1355-228.cnf
```
To save the identified MUSes into a file, run:
```
timeout 10 ./must -o output_file examples/test.cnf
```
To see all the available parameters, run:
```
./must -h
```

## Other Third-Party Tools
Besides the above mentioned tools, we also use two single MUS extraction algorithms: muser2 [5] and mcsmus [8]. You do not have to install these. Muser2 is presented in our repo in a binary form, and mcsmus is compiled with our code. 

## References
* [1] Mark H. Liffiton, Alessandro Previti, Ammar Malik, João Marques-Silva: Fast, flexible MUS enumeration. Constraints 21(2), 2016.
* [2] Jaroslav Bendík, Nikola Beneš, Ivana Černá, Jiří Barnat: Tunable Online MUS/MSS Enumeration. FSTTCS 2016.
* [3] Jaroslav Bendík, Ivana Černá, Nikola Beneš: Recursive Online Enumeration of All Minimal Unsatisfiable Subsets. ATVA 2018.
* [4] http://minisat.se/
* [5] https://bitbucket.org/anton_belov/muser2
* [6] https://github.com/Z3Prover/z3
* [7] https://spot.lrde.epita.fr/
* [8] https://bitbucket.org/gkatsi/mcsmus

## Contact
In case of any troubles, do not hesitate to contact me, Jaroslav Bendik, the developer of the tool, at xbendik=at=fi.muni.cz.
