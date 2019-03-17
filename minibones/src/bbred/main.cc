#include <stdio.h>
#include <time.h> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iosfwd>
#include <sstream>
#include <signal.h>
#include <unistd.h>
#include <algorithm>
#include "ReadCNF.hh"
#include "auxiliary.hh"
#include "BBRed.hh"

using std::cout;
using std::cin;
using std::setprecision;

Reader* make_reader(string flafile);
void print_usage(const char*);

int main(int argc, char** argv) {
  if (argc!=2) {
    cerr<<"Error parsing options."<<endl;
    print_usage(argc ? argv[0] : "bbred");
    exit(100);
  }
  //read input
  Reader* fr = make_reader(string(argv[1]));  
  ReadCNF reader(*fr);
  reader.read();
  const CNF& cnf = reader.get_clause_vector();
  BBRed r(reader.get_max_id(),reader.get_clause_vector());
  r.run();

  assert(r.get_redundant_count()>=0);
  size_t countdown=(cnf.size()-r.get_redundant_count());
  cout<<"c number of removed clauses: "<<r.get_redundant_count()<<std::endl;
  cout<<"p cnf "<<reader.get_max_id()<<" "<<countdown<<std::endl;
  for (size_t i=0; countdown; ++i) {
    if (r.is_redundant(i)) continue;
    --countdown;
    std::cout<<cnf[i]<<" 0"<<endl;
  }
  assert(countdown==0);
  delete fr;
}

Reader* make_reader(string flafile) {
  gzFile ff=Z_NULL;
  if (flafile.size()==1 && flafile[0]=='-') {
    return new Reader(std::cin);
  } else {
    ff = gzopen(flafile.c_str(), "rb");
    if (ff==Z_NULL) {
      cerr << "Unable to open file: " << flafile << endl;
      exit(100);
    }
    return new Reader(ff);
  }
  assert(0);
  return NULL;
}

void print_usage(const char* program_name) {
  cout<<"USAGE: "<<program_name<<" <qdimacs_filename>"<<endl;
  cout<<"  Note: If filename is -, then the formula is read from standard input."<<endl;
}
