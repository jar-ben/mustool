/* 
 * File:   BBInfo.cc
 * Author: mikolas
 * 
 * Created on May 30, 2011, 6:23 PM
 */

#include <stddef.h>

#include "BBInfo.hh"
using namespace minibones;

BBInfo::BBInfo(Var _max_id)
: max_id(_max_id)
, work_buf_sz((max_id+9)/8)
, done(new uint8_t[work_buf_sz])
{
    // initialize information about backbones
    const size_t literal_sz = (1+max_id)*2;
    might_be.resize(literal_sz,true); // anything can be backbone
    is.resize(literal_sz,false);      // nothing is backbone so far
}

BBInfo::~BBInfo() {}

