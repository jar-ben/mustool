/* 
 * File:   Reader.hh
 * Author: mikolas
 *
 * Created on January 12, 2011, 4:19 PM
 */

#ifndef READER_HH
#define	READER_HH
#include <iostream>
#include <zlib.h>
#include <utility>
#include <stdio.h>
#include "parse_utils.hh"
using std::istream;
class Reader {
public:
    Reader(gzFile& zf);
    Reader(StreamBuffer& zipStream);
    Reader(istream& stream);
    Reader(const Reader& orig);
    virtual ~Reader();
    int  operator*();
    void operator++();
    void skip_whitespace();
    inline size_t get_line_number();
private:
    size_t        lnn;
    StreamBuffer* zip;
    istream* s;
    int c;
};

inline size_t Reader::get_line_number() { return lnn; }
#endif	/* READER_HH */

