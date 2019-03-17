/* 
 * File:   ReadCNFException.hh
 * Author: mikolas
 *
 * Created on March 20, 2011, 10:03 PM
 */

#ifndef READCNFEXCEPTION_HH
#define	READCNFEXCEPTION_HH
#include <string>
using std::string;
using std::exception;
class ReadCNFException : public exception {
public:
    ReadCNFException(const string& message);
//    ReadQException(const ReadQException& other);
    const char* what() const throw();
    ~ReadCNFException()  throw() { delete[] s; }
private:
    char* s;
};

#endif	/* READCNFEXCEPTION_HH */

