#ifndef MCSMUS_ALGORITHM_HH
#define MCSMUS_ALGORITHM_HH

#include <vector>
#include <set>
#include <unordered_set>
#include <iostream>
#include <cassert>

namespace mcsmus {

//--------------------------------------------------
// Find min element with value given by a predicate
template<typename Iter, typename ValuePred>
auto min_element_value(Iter begin, Iter end, ValuePred p) -> decltype(p(*begin))
{
    assert( begin != end );
    auto rv = p(*begin);
    ++begin;
    while( begin != end ) {
        auto bv = p(*begin);
        if( bv < rv )
            rv = bv;
        ++begin;
    }
    return rv;
}

//--------------------------------------------------
// erase elements from a vector/vec, subject to a Pred
template<typename T, typename A, typename Pred>
void erase_if(std::vector<T, A>& v, Pred p)
{
    v.erase( remove_if(begin(v), end(v), p), end(v) );
}

//--------------------------------------------------
// move elements from one vector to another, subject to a Pred
template <typename T, typename A, typename Pred>
void move_if(std::vector<T, A>& from, std::vector<T, A>& to, Pred p)
{
    erase_if(from, [&](auto&& e) {
        if (p(e)) {
            to.emplace_back(e);
            return true;
        }
        return false;
    });
}

//--------------------------------------------------
// an abbreviation for find(x) != end for set,map, unordered_set,
// unordered_map
template<typename T>
bool contains(std::set<T> const& s, T const& t)
{
    return s.find(t) != s.end();
}

template<typename T, typename H>
bool contains(std::unordered_set<T, H> const& s, T const& t)
{
    return s.find(t) != s.end();
}

//--------------------------------------------------
// transform_if
//
// like std::transform, but also takes a unary predicate p, so that it
// only performs the transformation for elements evaluate p to
// true. Returns an iterator past the last element in the output
// sequence
template<typename InIter, typename OutIter, typename Op, typename Pred>
OutIter transform_if(InIter begin, InIter end, OutIter obegin,
                     Pred pred, Op op)
{
    for(; begin != end; ++begin)
        if( pred(*begin) )
            *obegin++ = op(*begin);
    return obegin;
}

}

//--------------------------------------------------
// output operators

namespace std {

template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& v)
{
    os << "v(sz=" << v.size() << ")[";
    bool first = true;
    for( auto&& t: v ) {
        if( first )
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, std::set<T> const& v)
{
    os << "s[";
    bool first = true;
    for( auto&& t: v ) {
        if( first )
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template<typename T, typename H>
std::ostream& operator<<(std::ostream& os, std::unordered_set<T, H> const& v)
{
    os << "s[";
    bool first = true;
    for( auto&& t: v ) {
        if( first )
            first = false;
        else
            os << ",";
        os << t;
    }
    os << "]";
    return os;
}

template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, std::pair<T, U> const& p)
{
    os << "p<" << p.first << "," << p.second << ">";
    return os;
}

}

namespace mcsmus {

//--------------------------------------------------
template <typename F>
struct on_scope_exit_
{
    F f;
    on_scope_exit_(F pf) : f(pf) {}
    ~on_scope_exit_() { f(); }
};

template<typename F>
on_scope_exit_<F> on_scope_exit(F f)
{
    return on_scope_exit_<F>(f);
}

//--------------------------------------------------

template <typename Int, typename Diff> struct int_iterator {
    Int i;

    explicit int_iterator(Int ii)
        : i(ii)
    {
    }

    int_iterator<Int, Diff> operator++()
    {
        ++i;
        return *this;
    }
    Int operator*() const { return i; }
    bool operator==(int_iterator o) const { return i == o.i; }
    bool operator!=(int_iterator o) const { return i != o.i; }

    typedef std::input_iterator_tag iterator_category;
    typedef Int value_type;
    typedef Diff difference_type;
    typedef Int* pointer;
    typedef Int& reference;
};

template <typename Int> struct int_range {
    static_assert(std::is_integral<Int>::value,
        "int_range<Int> requires Int to be an integral type");

    typedef typename std::make_signed<Int>::type difference_type;

    Int b, e;

    int_range(Int bb, Int ee)
        : b(bb)
        , e(ee)
    {
        assert(e >= b);
    }

    int_iterator<Int, difference_type> begin() const
    {
        return int_iterator<Int, difference_type>(b);
    }
    int_iterator<Int, difference_type> end() const
    {
        return int_iterator<Int, difference_type>(e);
    }
    size_t size() const { return e - b; }
};

template <typename Int> inline auto begin(int_range<Int> ir)
{
    return ir.begin();
}
template <typename Int> inline auto end(int_range<Int> ir)
{
    return ir.end();
}

//--------------------------------------------------

template <typename T> auto indices(std::vector<T> const& v)
{
    typedef typename std::vector<T>::size_type Int;
    return int_range<Int>(Int(0), v.size());
}

//=================================================================================================
}

#endif
