/* 
 * File:   ReadCNF.hh
 * Author: mikolas
 *
 * Created on November 11, 2010, 1:19 PM
 */

#ifndef READCNF_HH
#define	READCNF_HH
#include <zlib.h>
#include "auxiliary.hh"
#include "Reader.hh"
#include "ReadCNFException.hh"
#include "types.hh"
using std::vector;
using std::string;

class ReadCNF {
public:
  ReadCNF(Reader& input_file);
  ReadCNF(const ReadCNF& orig);
  virtual ~ReadCNF();
  void                       read();
  Var                        get_max_id() const         {return mxid;}
  const CNF&                 get_clause_vector() const  {return clause_vector;}
private:
  Reader&           in;
  Var               mxid;
  CNF               clause_vector;
  Lit               parse_lit(Reader& in);
  void              read_cnf_clause(vector<Lit>& lits);
};
#endif	/* READCNF_HH */

