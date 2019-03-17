/* 
 * File:   BBInfo.hh
 * Author: mikolas
 *
 * Created on May 30, 2011, 6:23 PM
 */

#ifndef BBINFO_HH
#define	BBINFO_HH
#include "minisat_aux.hh"

namespace minibones {
#include <vector>
using std::vector;


class BBInfo {
public:
    BBInfo(Var max_id);
    virtual ~BBInfo();
    /** Marks the given literal as certainly not a backbone.
     * @param literal tested literal
     * @return true iff {@code literal} could have been a backbone before */
    inline bool debone(const Lit& literal);
    inline bool mark_bone(const Lit& literal);
    inline bool might_be_bb(const Lit& literal) const;
    inline bool is_done(size_t variable) const;
    inline bool set_done(size_t variable);
    bool backbone_sign(Var var) const;
    inline bool is_bb(const Lit& literal) const;
    inline void* cpy_done(void* dest) const;

private:
    const Var    max_id;
    const size_t work_buf_sz;
    vector<bool> might_be;    // literals that might be backbone
    vector<bool> is;          // literals that are backbone
    uint8_t* const done;
    inline bool update_done(const Lit& literal);
    inline bool _debone(const Lit& literal);
    inline bool _mark_bone(const Lit& literal);
};

bool BBInfo::set_done(size_t variable) {
    const size_t by = variable / 8;
    const size_t bi = variable % 8;
    const bool   r  = done[by] & (1<<bi);
    done[by] |= (1<<bi);
    return r;
}

bool BBInfo::is_done(size_t variable) const {
    const size_t by = variable / 8;
    const size_t bi = variable % 8;
    return done[by] & (1<<bi);
}

bool BBInfo::might_be_bb(const Lit& literal) const {
   return might_be[literal_index(literal)]; }

bool BBInfo::is_bb(const Lit& literal) const {
    return is[literal_index(literal)]; }

bool BBInfo::mark_bone(const Lit& literal) {
    const bool return_value = _mark_bone(literal);
    if (return_value) update_done(literal);
    return return_value;
}

bool BBInfo::debone(const Lit& literal) {
    const bool return_value = _debone(literal);
    if (return_value) update_done(literal);
    return return_value;
}

bool BBInfo::_mark_bone(const Lit& literal) {
    const size_t index = literal_index(literal);
    if (is[index]) return false; // we already know it's a backbone
    is[index] = true;
    return true;
}

bool BBInfo::_debone(const Lit& literal) {
    const size_t index = literal_index(literal);
    if (!might_be[index]) return false; // we already know it's not a backbone
    might_be[index] = false;
    return true;
}

bool BBInfo::update_done(const Lit& literal) {
    const size_t index =  literal_index(literal);
    const size_t indexo = literal_index(~literal);
    const size_t v= (size_t)( var(literal) );
    const size_t by = v / 8;
    const size_t bi = v % 8;
    if (done[by] & (1<<bi)) return false;
    if ( (is[index] || is [indexo] )
       || (!might_be[index] && !might_be[indexo]) ) {
        done[by] |= (1<<bi);
    }
    return true;
}

void* BBInfo::cpy_done(void* dest) const {
   return memcpy(dest,done,work_buf_sz);
}

}
#endif	/* BBINFO_HH */

