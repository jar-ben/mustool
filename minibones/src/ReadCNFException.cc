/* 
 * File:   ReadCNFException.cc
 * Author: mikolas
 * 
 * Created on March 20, 2011, 10:03 PM
 */

#include "ReadCNFException.hh"

ReadCNFException::ReadCNFException(const string& message) {
    s=new char[message.size()+1];
    message.copy(s,message.size(),0);
}

const char* ReadCNFException::what() const throw() {return s;}



