/* 
 * File:   Reader.cc
 * Author: mikolas
 * 
 * Created on January 12, 2011, 4:19 PM
 */

#include "Reader.hh"
Reader::Reader(gzFile& zf) :  lnn(0), zip(new StreamBuffer(zf)), s(NULL) {}
Reader::Reader(StreamBuffer& zipStream) : lnn(0), zip(&zipStream), s(NULL) {}
Reader::Reader(istream& stream) :  lnn(0), zip(NULL), s(&stream) { c=s->get(); }
Reader::Reader(const Reader& orig) :  lnn(orig.lnn), zip(orig.zip), s(orig.s), c(orig.c) {}
Reader::~Reader() {}

int  Reader::operator *  () {
  const int r = s==NULL ? **zip :  c;
  return r; }

void Reader::operator ++ () { 
    if (s==NULL) {
       ++(*zip);
       if ((**zip)=='\n') ++lnn;
    }
    else {
       if (s->eof()) c=EOF;
       else c = s->get();
       if (c=='\n') ++lnn;
    }
}
void Reader::skip_whitespace()  {
      while (((**this) >= 9 && (**this) <= 13) || (**this) == 32) ++(*this);
}
