// Copyright (c) 1998-2005,2007  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Sylvain Pion, Michael Hemmer

#ifndef CGAL_INTERVAL_NT_H
#define CGAL_INTERVAL_NT_H

// This file contains the description of the following classes:
// - Interval_nt<false>  It's a number type that needs the FPU rounding mode
//                       to be set to +inf.  It is also typedef'd to
//                       Interval_nt_advanced for backward compatibility.
// - Interval_nt<true>   Same but it does the rounding mode itself so you
//                       don't have to worry about it.  But it's slower.
//
// Note: When rounding is towards +infinity, to make an operation rounded
// towards -infinity, it's enough to take the opposite of some of the operand,
// and the opposite of the result (see operator+, operator*,...).

// TODO : 
// - test whether stopping constant propagation only in functions taking
//   double as arguments, improves performance.

#include <x86intrin.h>
#include <utility> // for std::pair
#include <CGAL/number_type_config.h>
#include <CGAL/number_utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/number_utils.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/double.h>
#include <CGAL/FPU.h>
#include <CGAL/IO/io.h>
#include <iostream>

namespace CGAL {

template <bool Protected = true>
class Interval_nt
{
  typedef Interval_nt<Protected>     IA;
  typedef std::pair<double, double>  Pair;

public:

  typedef double      value_type;

  typedef Uncertain_conversion_exception            unsafe_comparison;
  typedef Checked_protect_FPU_rounding<Protected>   Internal_protector;
  typedef Protect_FPU_rounding<!Protected>          Protector;

  Interval_nt()
#ifndef CGAL_NO_ASSERTIONS
#ifdef CGAL_USE_SSE2
      : val(CGAL_OPACIFY_M128D_CST(_mm_setr_pd(-1, 0)))
#else
      : _inf(-1), _sup(0)
#endif
             // to early and deterministically detect use of uninitialized
#endif
    {}

  Interval_nt(int i)
  {
    *this = Interval_nt(static_cast<double>(i));
  }

  Interval_nt(unsigned i)
  {
    *this = Interval_nt(static_cast<double>(i));
  }

  Interval_nt(long long i)
  {
    // Is this safe against excess precision? -- Marc Glisse, Dec 2012
    double d = static_cast<double>(i);
    *this = Interval_nt(d);
#ifdef __GNUC__
    long long safe = 1LL << 52; // Use numeric_limits?
    bool exact = ((long long)d == i) || (i <= safe && i >= -safe);
    if (!(__builtin_constant_p(exact) && exact))
#endif
    // gcc ignores -frounding-math when converting integers to floats.
      *this += smallest();
  }

  Interval_nt(unsigned long long i)
  {
    double d = static_cast<double>(i);
    *this = Interval_nt(d);
#ifdef __GNUC__
    unsigned long long safe = 1ULL << 52; // Use numeric_limits?
    bool exact = ((unsigned long long)d == i) || (i <= safe);
    if (!(__builtin_constant_p(exact) && exact))
#endif
      *this += smallest();
  }

  Interval_nt(long i)
  {
    *this = (sizeof(int)==sizeof(long)) ?
      Interval_nt((int)i) :
      Interval_nt((long long)i);
  }

  Interval_nt(unsigned long i)
  {
    *this = (sizeof(int)==sizeof(long)) ?
      Interval_nt((unsigned)i) :
      Interval_nt((unsigned long long)i);
  }

  Interval_nt(double d)
  {
    CGAL_assertion(is_finite(d));
    *this = Interval_nt(d, d);
  }

// The Intel compiler on Linux is aggressive with constant propagation and
// it seems there is no flag to stop it, so disable this check for it.
#if !defined(CGAL_DISABLE_ROUNDING_MATH_CHECK) && \
    defined(__INTEL_COMPILER) && defined(__linux)
#  define CGAL_DISABLE_ROUNDING_MATH_CHECK
#endif

#ifdef CGAL_USE_SSE2
  // This constructor should really be private, like the simd() function, but
  // that would mean a lot of new friends, so they are only undocumented.
  explicit Interval_nt(__m128d v) : val(CGAL_OPACIFY_M128D_CST2(v)) {}
  // Opacifying here instead of in the operations experimentally works.
#endif

  // WARNING: the version of clang I tested not only propagates constants, but
  // moves fesetenv across operations and optimizes a^x+a^y to a^(x+y).
  Interval_nt(double i, double s)
#ifdef CGAL_USE_SSE2
    : val(CGAL_OPACIFY_M128D_CST(_mm_setr_pd(-i, s)))
#else
    : _inf(-i), _sup(s)
#endif
  {
      // VC++ should use instead : (i<=s) || !is_valid(i) || !is_valid(s)
      // Or should I use is_valid() ? or is_valid_or_nan() ?
    CGAL_assertion_msg(!(i>s),
	      " Variable used before being initialized (or CGAL bug)");
#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
    CGAL_assertion_code((void) tester;) // Necessary to trigger a runtime test of rounding modes.
#endif
  }

  Interval_nt(const Pair & p)
  {
    *this = Interval_nt(p.first, p.second);
  }

  IA operator-() const {
#ifdef CGAL_USE_SSE2
    return IA (_mm_shuffle_pd(val, val, 1));
#else
    return IA (-sup(), -inf());
#endif
  }

  IA & operator+= (const IA &d) { return *this = *this + d; }
  IA & operator-= (const IA &d) { return *this = *this - d; }
  IA & operator*= (const IA &d) { return *this = *this * d; }
  IA & operator/= (const IA &d) { return *this = *this / d; }

  bool is_point() const
  {
    // TODO: sup()-inf()==0 might help a compiler generate the SSE3 haddpd.
    return sup() == inf();
  }

  bool is_same (const IA & d) const
  {
#ifdef CGAL_USE_SSE2
    // Faster to answer yes, but slower to answer no.
    return _mm_movemask_pd (_mm_cmpneq_pd (val, d.val)) == 0;
#else
    return inf() == d.inf() && sup() == d.sup();
#endif
  }

  bool do_overlap (const IA & d) const
  {
#ifdef CGAL_USE_SSE2
    __m128d m = _mm_set1_pd (-0.);
    __m128d y = _mm_xor_pd ((-d).val, m); // {-ds,di}
    __m128d c = _mm_cmplt_pd (val, y); // {i>ds,s<di}
    return _mm_movemask_pd (c) == 0;
#else
    return !(d.inf() > sup() || d.sup() < inf());
#endif
  }

  double inf() const {
#ifdef CGAL_USE_SSE2
    return -_mm_cvtsd_f64(val);
#else
    return -_inf;
#endif
  }
  double sup() const {
#ifdef CGAL_USE_SSE2
    // Should we use shufpd because it is already used by operator- ?
    return _mm_cvtsd_f64(_mm_unpackhi_pd(val, val));
#else
    return _sup;
#endif
  }
#ifdef CGAL_USE_SSE2
  __m128d simd() const { return val; }
#endif

  std::pair<double, double> pair() const
  {
    return std::pair<double, double>(inf(), sup());
  }

  static IA largest()
  {
    return IA(-internal::infinity, internal::infinity);
  }

  static IA smallest()
  {
    return IA(-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE);
  }

#if 0 // def CGAL_HISTOGRAM_PROFILER  // not yet ready
  ~Interval_nt()
  {
    CGAL_HISTOGRAM_PROFILER("[Interval_nt relative precision in log2 scale]",
                             (unsigned) ( ::log(relative_precision(*this))) / ::log(2.0) )  );
  }
#endif

private:
  // Pair inf_sup;
  // The value stored in _inf is the negated lower bound.
  // TODO: experiment with different orders of the values in the SSE2 register,
  // for instance {sup, -inf}, or {inf, -sup}, and adapt users to query the low
  // value in priority. {-inf, sup} has the drawback that neither inf nor sup
  // is free to access.
#ifdef CGAL_USE_SSE2
  __m128d val;
#else
  double _inf, _sup;
#endif

  struct Test_runtime_rounding_modes {
    Test_runtime_rounding_modes()
    {
      // We test whether GCC's -frounding-math option has been forgotten.
      // The macros CGAL_IA_MUL and CGAL_IA_DIV stop constant propagation only
      // on the second argument, so if -fno-rounding-math, the compiler optimizes
      // the 2 negations and we get wrong rounding.
      typename Interval_nt<>::Internal_protector P;
      CGAL_assertion_msg(-CGAL_IA_MUL(-1.1, 10.1) != CGAL_IA_MUL(1.1, 10.1),
                         "Wrong rounding: did you forget the  -frounding-math  option if you use GCC (or  -fp-model strict  for Intel)?");
      CGAL_assertion_msg(-CGAL_IA_DIV(-1., 10) != CGAL_IA_DIV(1., 10),
                         "Wrong rounding: did you forget the  -frounding-math  option if you use GCC (or  -fp-model strict  for Intel)?");
    }
  };

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
  static Test_runtime_rounding_modes tester;
#endif
};

#ifndef CGAL_DISABLE_ROUNDING_MATH_CHECK
template <bool Protected>
typename Interval_nt<Protected>::Test_runtime_rounding_modes
Interval_nt<Protected>::tester;
#endif

template <bool Protected>
inline
Uncertain<bool>
operator<(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (a.sup()  < b.inf()) return true;
  if (a.inf() >= b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (a.sup() <= b.inf()) return true;
  if (a.inf() >  b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (b.inf() >  a.sup() || b.sup() <  a.inf()) return false;
  if (b.inf() == a.sup() && b.sup() == a.inf()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return ! (a == b); }


// Mixed operators with double.

template <bool Protected>
inline
Uncertain<bool>
operator<(double a, const Interval_nt<Protected> &b)
{
  if (a < b.inf()) return true;
  if (a >= b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(double a, const Interval_nt<Protected> &b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(double a, const Interval_nt<Protected> &b)
{
  if (a <= b.inf()) return true;
  if (a >  b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(double a, const Interval_nt<Protected> &b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(double a, const Interval_nt<Protected> &b)
{
  if (b.inf() >  a || b.sup() <  a) return false;
  if (b.is_point()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(double a, const Interval_nt<Protected> &b)
{ return ! (a == b); }

template <bool Protected>
inline
Uncertain<bool>
operator<(const Interval_nt<Protected> &a, double b)
{
  if (a.inf() >= b) return false;
  if (a.sup()  < b) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(const Interval_nt<Protected> &a, double b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(const Interval_nt<Protected> &a, double b)
{
  if (a.inf() >  b) return false;
  if (a.sup() <= b) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(const Interval_nt<Protected> &a, double b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(const Interval_nt<Protected> &a, double b)
{ return b == a; }

template <bool Protected>
inline
Uncertain<bool>
operator!=(const Interval_nt<Protected> &a, double b)
{ return ! (b == a); }


// Mixed operators with int.

template <bool Protected>
inline
Uncertain<bool>
operator<(int a, const Interval_nt<Protected> &b)
{
  return static_cast<double>(a) < b;
}

template <bool Protected>
inline
Uncertain<bool>
operator>(int a, const Interval_nt<Protected> &b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(int a, const Interval_nt<Protected> &b)
{
  return static_cast<double>(a) <= b;
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(int a, const Interval_nt<Protected> &b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(int a, const Interval_nt<Protected> &b)
{
  return static_cast<double>(a) == b;
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(int a, const Interval_nt<Protected> &b)
{ return ! (a == b); }

template <bool Protected>
inline
Uncertain<bool>
operator<(const Interval_nt<Protected> &a, int b)
{
  return a < static_cast<double>(b);
}

template <bool Protected>
inline
Uncertain<bool>
operator>(const Interval_nt<Protected> &a, int b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(const Interval_nt<Protected> &a, int b)
{
  return a <= static_cast<double>(b);
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(const Interval_nt<Protected> &a, int b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(const Interval_nt<Protected> &a, int b)
{
  return a == static_cast<double>(b);
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(const Interval_nt<Protected> &a, int b)
{ return ! (a == b); }



// Non-documented
// Returns true if the interval is a unique representable double.
template <bool Protected>
inline
bool
fit_in_double (const Interval_nt<Protected> & d, double &r)
{
  bool b = d.is_point();
  if (b)
    r = d.inf();
  return b;
}

// Non-documented
template <bool Protected>
inline
bool
is_singleton (const Interval_nt<Protected> & d)
{
  return d.is_point();
}

// Non-documented
template <bool Protected>
inline
double
magnitude (const Interval_nt<Protected> & d)
{
#ifdef CGAL_USE_SSE2
#if 0
  //FIXME: Intel's compiler seems to be missing _mm_set1_epi64x ???
  const __m128d m = _mm_castsi128_pd (_mm_set1_epi64x (0x7fffffffffffffff));
#else
  union { long long l; double d; } b;
  b.l = 0x7fffffffffffffff;
  const __m128d m = _mm_set1_pd(b.d);
#endif
  __m128d x = _mm_and_pd (d.simd(), m); // { abs(inf), abs(sup) }
  __m128d y = _mm_unpackhi_pd (x, x);
  return _mm_cvtsd_f64 (_mm_max_sd (x, y));
#else
  return (std::max)(CGAL::abs(d.inf()), CGAL::abs(d.sup()));
#endif
}

// Non-documented
template <bool Protected>
inline
double
width (const Interval_nt<Protected> & d)
{
  return d.sup() - d.inf();
}

// Non-documented
template <bool Protected>
inline
double
radius (const Interval_nt<Protected> & d)
{
  return width(d)/2; // This could be improved to avoid overflow.
}

// Non-documented
// This is the relative precision of to_double() (the center of the interval),
// hence we use radius() instead of width().
template <bool Protected>
inline
bool
has_smaller_relative_precision(const Interval_nt<Protected> & d, double prec)
{
  return magnitude(d) == 0 || radius(d) < prec * magnitude(d);
}

// Non-documented
template <bool Protected>
double
relative_precision(const Interval_nt<Protected> & d)
{
  if (magnitude(d) == 0.0)
    return 0.0;
  return radius(d) / magnitude(d);
}


template< bool Protected >
class Is_valid< Interval_nt<Protected> >
  : public std::unary_function< Interval_nt<Protected>, bool > {
  public :
    bool operator()( const Interval_nt<Protected>& x ) const {
      return is_valid(-x.inf()) &&
             is_valid(x.sup()) &&
             x.inf() <= x.sup();
    }
};

template <bool Protected>
std::ostream & operator<< (std::ostream &os, const Interval_nt<Protected> & I )
{
  return os << "[" << I.inf() << ";" << I.sup() << "]";
}

#define CGAL_SWALLOW(IS,CHAR)                           \
    {                                                   \
        char c;                                         \
        do c = is.get(); while (isspace(c));            \
        if (c != CHAR) {                                \
            is.setstate(std::ios_base::failbit);        \
        }                                               \
    }                                                   \

template <bool Protected>
std::istream & operator>> (std::istream &is, Interval_nt<Protected> & I)
{
    char c;
    do c = is.get(); while (isspace(c));
    is.putback(c);
    if(c == '['){ // read original output from operator <<
        double inf,sup;
        CGAL_SWALLOW(is, '[');// read the "["
        is >> iformat(inf);
        CGAL_SWALLOW(is, ';');// read the ";"
        is >> iformat(sup);
        CGAL_SWALLOW(is, ']');// read the "]"
        I = Interval_nt<Protected>(inf,sup);
    }else{ //read double (backward compatibility)
        double d;
        is >> d;
        I = d;
    }
    return is;
}
#undef CGAL_SWALLOW

typedef Interval_nt<false> Interval_nt_advanced;  // for backward-compatibility


template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  typename Interval_nt<Protected>::Internal_protector P;
#ifdef CGAL_USE_SSE2
  return Interval_nt<Protected>(CGAL_IA_M128D_ADD(a.simd(), b.simd()));
#else
  return Interval_nt<Protected> (-CGAL_IA_ADD(-a.inf(), -b.inf()),
                                  CGAL_IA_ADD(a.sup(), b.sup()));
#endif
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> & a, double b)
{
  return a+Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (double a, const Interval_nt<Protected> & b)
{
  return b+a;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (int a, const Interval_nt<Protected> & b)
{
  return static_cast<double>(a)+b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> & a, int b)
{
  return a+static_cast<double>(b);
}

template< bool Protected >
inline
Interval_nt<Protected>
operator+( const Interval_nt<Protected>& a ) {
  return a;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  return a+(-b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (double a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)-b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> & a, double b)
{
  return a+Interval_nt<Protected>(-b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (int a, const Interval_nt<Protected> & b)
{
  return static_cast<double>(a)-b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> & a, int b)
{
  return a-static_cast<double>(b);
}

#if 0
static __m128d _mm_blendv_pd(__m128d n, __m128d p, __m128d m){
  //FIXME: needs m=_mm_cmplt_pd(m,_mm_setzero_pd());
  __m128d pp=_mm_and_pd(m,p);
  __m128d nn=_mm_andnot_pd(m,n);
  return _mm_or_pd(pp,nn);
}
#endif

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  typedef Interval_nt<Protected> IA;
  typename Interval_nt<Protected>::Internal_protector P;
//FIXME: use protection for all variants
#ifdef CGAL_USE_SSE2
# if 0
  // With branches, better than scalar but still not so fast.
  __m128d aa = a.simd();
  double t = -b.inf(); // -bi
  double u = b.sup();  //  bs
  if (t <= 0) { // b >= 0
    __m128d res1 = CGAL_IA_M128D_MUL (_mm_set1_pd (-t), aa); // {-ai*bi,as*bi}
    __m128d res2 = CGAL_IA_M128D_MUL (_mm_set1_pd ( u), aa); // {-ai*bs,as*bs}
    return IA (_mm_max_pd (res1, res2));
  }
  __m128d     ap = _mm_shuffle_pd (aa, aa, 1); // {as,-ai}
  __m128d      x = CGAL_IA_M128D_MUL (_mm_set1_pd ( t), ap); // {-as*bi,ai*bi}

  __m128d y;
  if (u < 0)   y = CGAL_IA_M128D_MUL (_mm_set1_pd (-u), ap); // {-as*bs,ai*bs}
  else         y = CGAL_IA_M128D_MUL (_mm_set1_pd ( u), aa); // {-ai*bs,as*bs}
  return IA (_mm_max_pd (x, y));
# elif 0
  // Another branchy version with less _mm_set1_pd
  // barely better
  // swapping a and b sometimes improves the timing...
  __m128d aa = a.simd();				// {-ai,as}
  __m128d bb = b.simd();				// {-bi,bs}
  double nai = _mm_cvtsd_f64(aa);			// -ai
  __m128d ap = _mm_shuffle_pd (aa, aa, 1);		// {as,-ai}
  if(nai<=0){
    __m128d mi = _mm_set_sd(-0.);			// {-0,+0}
    __m128d c = _mm_xor_pd(aa, mi);			// {ai,as}
    __m128d cp = _mm_shuffle_pd (c, c, 1);		// {as,ai}
    __m128d x = _mm_mul_pd(c,bb);			// {-ai*bi,as*bs}
    __m128d y = _mm_mul_pd(cp,bb);			// {-as*bi,ai*bs}
    return IA(_mm_max_pd(x,y));
  }else if(_mm_cvtsd_f64(ap)<=0){
    __m128d ms = _mm_setr_pd(0.,-0.);			// {+0,-0}
    __m128d c = _mm_xor_pd(aa, ms);			// {-ai,-as}
    __m128d cp = _mm_shuffle_pd (c, c, 1);		// {-as,-ai}
    __m128d bp = _mm_shuffle_pd (bb, bb, 1);		// {bs,-bi}
    __m128d x = _mm_mul_pd(c,bp);			// {-ai*bs,as*bi}
    __m128d y = _mm_mul_pd(cp,bp);			// {-as*bs,ai*bi}
    return IA(_mm_max_pd(x,y));
  }else{
    __m128d c = _mm_unpackhi_pd(bb, bb);		// {bs,bs}
    __m128d d = _mm_set1_pd(_mm_cvtsd_f64(bb));		// {-bi,-bi}
    __m128d x = _mm_mul_pd(c, aa);			// {-ai*bs,as*bs}
    __m128d y = _mm_mul_pd(d, ap);			// {-as*bi,ai*bi}
    return IA(_mm_max_pd(x,y));
  }
# elif 0
// we want to multiply -ai,as with {ai<0?bs:bi,as<0?bi:bs}
// we want to multiply -as,ai with {as<0?bs:bi,ai<0?bi:bs}
// too many xor
  __m128d m = _mm_set_sd(-0.);				// {-0,+0}
  __m128d m1 = _mm_set1_pd(-0.);			// {-0,-0}
  __m128d aa = a.simd();				// {-ai,as}
  __m128d ax = _mm_shuffle_pd (aa, aa, 1);		// {as,-ai}
  __m128d ap = _mm_xor_pd (ax, m1);			// {-as,ai}
  __m128d bb = _mm_xor_pd(b.simd(), m);			// {bi,bs}
  __m128d bp = _mm_shuffle_pd (bb, bb, 1);		// {bs,bi}
  __m128d az = _mm_xor_pd(aa, m);			// {ai,as}
  __m128d neg = _mm_cmplt_pd (az, _mm_setzero_pd());	// {ai<0,as<0}
  __m128d x = _mm_blendv_pd (bb, bp, neg);		// {ai<0?bs:bi,as<0?bi:bs}
  __m128d negp = _mm_shuffle_pd (neg, neg, 1);		// {as<0,ai<0}
  __m128d y = _mm_blendv_pd (bb, bp, negp);		// {as<0?bs:bi,ai<0?bi:bs}
  __m128d p1 = _mm_mul_pd (aa, x);
  __m128d p2 = _mm_mul_pd (ap, y);
  return IA (_mm_max_pd (p1, p2));
# elif 1
// we want to multiply ai,as with {ai<0?-bs:-bi,as<0?bi:bs}
// we want to multiply as,ai with {as<0?-bs:-bi,ai<0?bi:bs}
// best?
  __m128d m = _mm_set_sd(-0.);				// {-0,+0}
  __m128d m1 = _mm_set1_pd(-0.);			// {-0,-0}
  __m128d aa = a.simd();				// {-ai,as}
  __m128d az = _mm_xor_pd(aa, m);			// {ai,as}
  __m128d azp = _mm_shuffle_pd (az, az, 1);		// {as,ai}
  __m128d bb = b.simd();				// {-bi,bs}
  __m128d bx = _mm_shuffle_pd (bb, bb, 1);		// {bs,-bi}
  __m128d bp = _mm_xor_pd(bx, m1);			// {-bs,bi}
  __m128d x = _mm_blendv_pd (bb, bp, az);		// {ai<0?-bs:-bi,as<0?bi:bs}
  __m128d y = _mm_blendv_pd (bb, bp, azp);		// {as<0?-bs:-bi,ai<0?bi:bs}
  __m128d p1 = _mm_mul_pd (az, x);
  __m128d p2 = _mm_mul_pd (azp, y);
  return IA (_mm_max_pd (p1, p2));
# elif 0
// we want to multiply -ai,as with {ai>0?bi:bs,as<0?bi:bs}
// we want to multiply -as,ai with {as<0?bs:bi,ai>0?bs:bi}
// comparable
  __m128d m1 = _mm_set1_pd(-0.);			// {-0,-0}
  __m128d aa = a.simd();				// {-ai,as}
  __m128d ax = _mm_shuffle_pd (aa, aa, 1);		// {as,-ai}
  __m128d ap = _mm_xor_pd (ax, m1);			// {-as,ai}
  __m128d bb = b.simd();				// {-bi,bs}
  double bi = -_mm_cvtsd_f64(bb);
  double bs = _mm_cvtsd_f64(_mm_unpackhi_pd(bb,bb));
  __m128d bbi = _mm_set1_pd(bi);			// {bi,bi}
  __m128d bbs = _mm_set1_pd(bs);			// {bs,bs}
  __m128d x = _mm_blendv_pd (bbs, bbi, aa);		// {ai>0?bi:bs,as<0?bi:bs}
  __m128d y = _mm_blendv_pd (bbi, bbs, ax);		// {as<0?bs:bi,ai>0?bs:bi}
  __m128d p1 = _mm_mul_pd (aa, x);
  __m128d p2 = _mm_mul_pd (ap, y);
  return IA (_mm_max_pd (p1, p2));
# elif 0
  // Brutal, compute all products in all directions.
  // The actual winner (by a hair) on recent hardware
  __m128d aa = a.simd();				// {-ai,as}
  __m128d bb = b.simd();				// {-bi,bs}
  __m128d m = _mm_set_sd(-0.);				// {-0,+0}
  __m128d m1 = _mm_set1_pd(-0.);			// {-0,-0}
  __m128d ax = _mm_shuffle_pd (aa, aa, 1);		// {as,-ai}
  __m128d ap = _mm_xor_pd (ax, m1);			// {-as,ai}
  __m128d bz = _mm_xor_pd(bb, m);			// {bi,bs}
  __m128d c = _mm_shuffle_pd (bz, bz, 1);		// {bs,bi}
  __m128d x1 = _mm_mul_pd(aa,bz);			// {-ai*bi,as*bs}
  __m128d x2 = _mm_mul_pd(aa,c);			// {-ai*bs,as*bi}
  __m128d x3 = _mm_mul_pd(ap,bz);			// {-as*bi,ai*bs}
  __m128d x4 = _mm_mul_pd(ap,c);			// {-as*bs,ai*bi}
  __m128d y1 = _mm_max_pd(x1,x2);
  __m128d y2 = _mm_max_pd(x3,x4);
  return IA (_mm_max_pd (y1, y2));
# else
  // AVX version of the brutal method, same running time or slower
  __m128d aa = a.simd();				// {-ai,as}
  __m128d bb = b.simd();				// {-bi,bs}
  __m128d m = _mm_set_sd(-0.);				// {-0,+0}
  __m128d m1 = _mm_set1_pd(-0.);			// {-0,-0}
  __m128d ax = _mm_shuffle_pd (aa, aa, 1);		// {as,-ai}
  __m128d ap = _mm_xor_pd (ax, m1);			// {-as,ai}
  __m128d bz = _mm_xor_pd(bb, m);			// {bi,bs}
  __m256d X = _mm256_set_m128d(ap,aa);			// {-ai,as,-as,ai}
  __m256d Y1 = _mm256_set_m128d(bz,bz);			// {bi,bs,bi,bs}
  __m256d Y2 = _mm256_permute_pd(Y1,5);			// {bs,bi,bs,bi}
  __m256d Z1 = _mm256_mul_pd(X,Y1);
  __m256d Z2 = _mm256_mul_pd(X,Y2);
  __m256d Z = _mm256_max_pd(Z1,Z2);
  __m128d z1 = _mm256_castpd256_pd128(Z);
  __m128d z2 = _mm256_extractf128_pd(Z,1);
  return IA (_mm_max_pd (z1, z2));
# endif
#else

  if (a.inf() >= 0.0)					// a>=0
  {
    // b>=0     [a.inf()*b.inf(); a.sup()*b.sup()]
    // b<=0     [a.sup()*b.inf(); a.inf()*b.sup()]
    // b~=0     [a.sup()*b.inf(); a.sup()*b.sup()]
    double aa = a.inf(), bb = a.sup();
    if (b.inf() < 0.0)
    {
	aa = bb;
	if (b.sup() < 0.0)
	    bb = a.inf();
    }
    return IA(-CGAL_IA_MUL(aa, -b.inf()), CGAL_IA_MUL(bb, b.sup()));
  }
  else if (a.sup()<=0.0)				// a<=0
  {
    // b>=0     [a.inf()*b.sup(); a.sup()*b.inf()]
    // b<=0     [a.sup()*b.sup(); a.inf()*b.inf()]
    // b~=0     [a.inf()*b.sup(); a.inf()*b.inf()]
    double aa = a.sup(), bb = a.inf();
    if (b.inf() < 0.0)
    {
	aa=bb;
	if (b.sup() < 0.0)
	    bb=a.sup();
    }
    return IA(-CGAL_IA_MUL(-bb, b.sup()), CGAL_IA_MUL(-aa, -b.inf()));
  }
  else						// 0 \in a
  {
    if (b.inf()>=0.0)				// b>=0
      return IA(-CGAL_IA_MUL(-a.inf(), b.sup()),
                 CGAL_IA_MUL( a.sup(), b.sup()));
    if (b.sup()<=0.0)				// b<=0
      return IA(-CGAL_IA_MUL( a.sup(), -b.inf()),
                 CGAL_IA_MUL(-a.inf(), -b.inf()));
        					// 0 \in b
    double tmp1 = CGAL_IA_MUL(-a.inf(),  b.sup());
    double tmp2 = CGAL_IA_MUL( a.sup(), -b.inf());
    double tmp3 = CGAL_IA_MUL(-a.inf(), -b.inf());
    double tmp4 = CGAL_IA_MUL( a.sup(),  b.sup());
    return IA(-(std::max)(tmp1,tmp2), (std::max)(tmp3,tmp4));
  }
#endif
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (double a, Interval_nt<Protected> b)
{
  typedef Interval_nt<Protected> IA;
  typename IA::Internal_protector P;
  if (a < 0) { a = -a; b = -b; }
  // Now a >= 0
#ifdef CGAL_USE_SSE2
  return IA(CGAL_IA_M128D_MUL (_mm_set1_pd(a), b.simd()));
#else
  return IA(-CGAL_IA_MUL(a, -b.inf()), CGAL_IA_MUL(a, b.sup()));
#endif
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> & a, double b)
{
  return b*a;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (int a, const Interval_nt<Protected> & b)
{
  return static_cast<double>(a)*b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> & a, int b)
{
  return a*static_cast<double>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  typedef Interval_nt<Protected> IA;
  typename Interval_nt<Protected>::Internal_protector P;
  // TODO: add an sse2 version similar to the multiplication.

#ifdef CGAL_USE_SSE2
  __m128d aa = a.simd();
  double t = -b.inf(); // -bi
  double u = b.sup();  //  bs
  if (t < 0) { // b > 0
    __m128d res1 = CGAL_IA_M128D_DIV (aa, _mm_set1_pd (-t)); // {-ai/bi,as/bi}
    __m128d res2 = CGAL_IA_M128D_DIV (aa, _mm_set1_pd ( u)); // {-ai/bs,as/bs}
    return IA (_mm_max_pd (res1, res2));
  }
  else if (u < 0) { // b < 0
    aa = _mm_shuffle_pd (aa, aa, 1); // {as,-ai}
    __m128d res1 = CGAL_IA_M128D_DIV (aa, _mm_set1_pd ( t)); // {-as/bi,ai/bi}
    __m128d res2 = CGAL_IA_M128D_DIV (aa, _mm_set1_pd (-u)); // {-as/bs,ai/bs}
    return IA (_mm_max_pd (res1, res2));
  }

#else
  if (b.inf() > 0.0)				// b>0
  {
    // e>=0	[a.inf()/b.sup(); a.sup()/b.inf()]
    // e<=0	[a.inf()/b.inf(); a.sup()/b.sup()]
    // e~=0	[a.inf()/b.inf(); a.sup()/b.inf()]
    double aa = b.sup(), bb = b.inf();
    if (a.inf() < 0.0)
    {
	aa = bb;
	if (a.sup() < 0.0)
	    bb = b.sup();
    }
    return IA(-CGAL_IA_DIV(-a.inf(), aa), CGAL_IA_DIV(a.sup(), bb));
  }
  else if (b.sup()<0.0)			// b<0
  {
    // e>=0	[a.sup()/b.sup(); a.inf()/b.inf()]
    // e<=0	[a.sup()/b.inf(); a.inf()/b.sup()]
    // e~=0	[a.sup()/b.sup(); a.inf()/b.sup()]
    double aa = b.sup(), bb = b.inf();
    if (a.inf() < 0.0)
    {
	bb = aa;
	if (a.sup() < 0.0)
	    aa = b.inf();
    }
    return IA(-CGAL_IA_DIV(a.sup(), -aa), CGAL_IA_DIV(a.inf(), bb));
  }
#endif
  else					// b~0
    return IA::largest();
	   // We could do slightly better -> [0;infinity] when b.sup()==0,
	   // but is this worth ?
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (double a, const Interval_nt<Protected> & b)
{
  // TODO: specialize
  return Interval_nt<Protected>(a)/b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (Interval_nt<Protected> a, double b)
{
  typedef Interval_nt<Protected> IA;
  typename IA::Internal_protector P;
  if (b < 0) { a = -a; b = -b; }
  else if (b == 0) return IA::largest();
  // Now b > 0
#ifdef CGAL_USE_SSE2
  return IA(CGAL_IA_M128D_DIV (a.simd(), _mm_set1_pd(b)));
#else
  return IA(-CGAL_IA_DIV(-a.inf(), b), CGAL_IA_DIV(a.sup(), b));
#endif
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (int a, const Interval_nt<Protected> & b)
{
  return static_cast<double>(a)/b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (const Interval_nt<Protected> & a, int b)
{
  return a/static_cast<double>(b);
}

// TODO: What about these two guys? Where do they belong to?
template <bool Protected>
struct Min <Interval_nt<Protected> >
    : public std::binary_function<Interval_nt<Protected>,
                             Interval_nt<Protected>,
                             Interval_nt<Protected> >
{
    Interval_nt<Protected> operator()( const Interval_nt<Protected>& d,
                                       const Interval_nt<Protected>& e) const
    {
#ifdef CGAL_USE_SSE2
        __m128d x = _mm_min_pd (d.simd(), e.simd());
        // Use _mm_max_sd instead?
        __m128d y = _mm_max_pd (d.simd(), e.simd());
        return Interval_nt<Protected> (_mm_move_sd (x, y));
#else
        return Interval_nt<Protected>(
                -(std::max)(-d.inf(), -e.inf()),
                 (std::min)( d.sup(),  e.sup()));
#endif
    }
};

template <bool Protected>
struct Max <Interval_nt<Protected> >
    : public std::binary_function<Interval_nt<Protected>,
                             Interval_nt<Protected>,
                             Interval_nt<Protected> >
{
    Interval_nt<Protected> operator()( const Interval_nt<Protected>& d,
                                       const Interval_nt<Protected>& e) const
    {
#ifdef CGAL_USE_SSE2
        // Use _mm_min_sd instead?
        __m128d x = _mm_min_pd (d.simd(), e.simd());
        __m128d y = _mm_max_pd (d.simd(), e.simd());
        return Interval_nt<Protected> (_mm_move_sd (y, x));
#else
        return Interval_nt<Protected>(
                -(std::min)(-d.inf(), -e.inf()),
                 (std::max)( d.sup(),  e.sup()));
#endif
    }
};

template<bool Protected> inline 
Interval_nt<Protected> min BOOST_PREVENT_MACRO_SUBSTITUTION(
const Interval_nt<Protected> & x,
const Interval_nt<Protected> & y){
  return CGAL::Min<Interval_nt<Protected> > ()(x,y);
}
template<bool Protected> inline 
Interval_nt<Protected> max BOOST_PREVENT_MACRO_SUBSTITUTION(
const Interval_nt<Protected> & x,
const Interval_nt<Protected> & y){
  return CGAL::Max<Interval_nt<Protected> > ()(x,y);
}



// TODO : document, when we are OK with the interface.
// - should it allow other number types for the exponent ?
template < bool b >
Interval_nt<b>
ldexp(const Interval_nt<b> &i, int e)
{
  double scale = std::ldexp(1.0, e);
  Interval_nt<b> scale_interval (
                      CGAL_NTS is_finite(scale) ? scale : CGAL_IA_MAX_DOUBLE,
                      scale == 0 ? CGAL_IA_MIN_DOUBLE : scale);
  return i * scale_interval;
}


// We also specialize some corresponding functors returning Uncertain<>.

// TODO: To which concept do these functors belong? Can we remove them??
template < bool b >
struct Equal_to < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x == y; }
};

template < bool b >
struct Not_equal_to < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x != y; }
};

template < bool b >
struct Greater < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x > y; }
};

template < bool b >
struct Less < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x < y; }
};

template < bool b >
struct Greater_equal < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x >= y; }
};

template < bool b >
struct Less_equal < Interval_nt<b>, Interval_nt<b> >
  : public std::binary_function< Interval_nt<b>, Interval_nt<b>, Uncertain<bool> >
{
  Uncertain<bool> operator()( const Interval_nt<b>& x,
                              const Interval_nt<b>& y) const
  { return x <= y; }
};


// As in MP_float.h, the namespace INTERN_INTERVAL_NT contains (now) global
// functions like square or sqrt which would have collided with the new
// global functions from AST/RET
//
// TODO: IMHO, a better solution would be to put the INTERN_MP_FLOAT-functions
//       into the MP_Float-class... But there is surely a reason why this is not
//       the case..?


namespace INTERN_INTERVAL_NT {

  template <bool Protected>
  inline
  double
  to_double (const Interval_nt<Protected> & d)
  {
    return (d.sup() + d.inf()) * 0.5;
    // This may overflow...
  }

  template <bool Protected>
  inline
  std::pair<double, double>
  to_interval (const Interval_nt<Protected> & d)
  {
    return d.pair();
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  sqrt (const Interval_nt<Protected> & d)
  {
    typename Interval_nt<Protected>::Internal_protector P;  // not optimal here.
    // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
    // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
    // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
    FPU_set_cw(CGAL_FE_DOWNWARD);
    double i = (d.inf() > 0.0) ? CGAL_IA_SQRT(d.inf()) : 0.0;
    FPU_set_cw(CGAL_FE_UPWARD);
    return Interval_nt<Protected>(i, CGAL_IA_SQRT(d.sup()));
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  square (Interval_nt<Protected> d)
  {
    typename Interval_nt<Protected>::Internal_protector P;
    if (d.inf()>=0.0)
#ifdef CGAL_USE_SSE2
      // TODO: try a branchless version
      {
	__m128d x = d.simd();
	__m128d r = CGAL_IA_M128D_MUL (_mm_xor_pd (x, _mm_set_sd (-0.)), x);
	return Interval_nt<Protected> (r);
      }
#else
        return Interval_nt<Protected>(-CGAL_IA_MUL(d.inf(), -d.inf()),
                                 CGAL_IA_MUL(d.sup(), d.sup()));
#endif
    if (d.sup()<=0.0)
#ifdef CGAL_USE_SSE2
      {
	__m128d x = (-d).simd();
	__m128d r = CGAL_IA_M128D_MUL (_mm_xor_pd (x, _mm_set_sd (-0.)), x);
	return Interval_nt<Protected> (r);
      }
#else
        return Interval_nt<Protected>(-CGAL_IA_MUL(d.sup(), -d.sup()),
                               CGAL_IA_MUL(-d.inf(), -d.inf()));
#endif
    return Interval_nt<Protected>(0.0, CGAL_IA_SQUARE((std::max)(-d.inf(),
                     d.sup())));
  }

  template <bool Protected>
  inline
  Interval_nt<Protected>
  abs (const Interval_nt<Protected> & d)
  {
    if (d.inf() >= 0.0) return d;
    if (d.sup() <= 0.0) return -d;
    return Interval_nt<Protected>(0.0, (std::max)(-d.inf(), d.sup()));
  }

  template <bool Protected>
  inline
  Uncertain<Sign>
  sign (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 0.0) return POSITIVE;
    if (d.sup() < 0.0) return NEGATIVE;
    if (d.inf() == d.sup()) return ZERO;
    return Uncertain<Sign>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<Comparison_result>
  compare (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
  {
    if (d.inf() > e.sup()) return LARGER;
    if (e.inf() > d.sup()) return SMALLER;
    if (e.inf() == d.sup() && d.inf() == e.sup()) return EQUAL;
    return Uncertain<Comparison_result>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_zero (const Interval_nt<Protected> & d)
  {
    // TODO: try using _mm_movemask_pd+_mm_cmp*_pd and is_point
    if (d.inf() > 0.0) return false;
    if (d.sup() < 0.0) return false;
    if (d.inf() == d.sup()) return true;
    return Uncertain<bool>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_one (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 1) return false;
    if (d.sup() < 1) return false;
    if (d.inf() == d.sup()) return true;
    return Uncertain<bool>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_positive (const Interval_nt<Protected> & d)
  {
    if (d.inf() > 0.0) return true;
    if (d.sup() <= 0.0) return false;
    return Uncertain<bool>::indeterminate();
  }

  template <bool Protected>
  inline
  Uncertain<bool>
  is_negative (const Interval_nt<Protected> & d)
  {
    if (d.inf() >= 0.0) return false;
    if (d.sup() < 0.0) return true;
    return Uncertain<bool>::indeterminate();
  }

 // TODO: Whats this for? Why is this in this file??
  inline
  std::pair<double, double>
  to_interval (const long & l)
  {
    if (sizeof(double) > sizeof(long)) {
      // On 64bit platforms, a long doesn't fit exactly in a double.
      // Well, a perfect fix would be to use std::numeric_limits<>, but...
      Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
      Interval_nt<false> approx ((double) l);
      FPU_set_cw(CGAL_FE_UPWARD);
      approx += Interval_nt<false>::smallest();
      return approx.pair();
    }
    else
      return std::pair<double,double>(l,l);
  }
} // namespace INTERN_INTERVAL_NT


template< bool B > class Real_embeddable_traits< Interval_nt<B> >
  : public INTERN_RET::Real_embeddable_traits_base< Interval_nt<B> , CGAL::Tag_true> {
  public:
    typedef Interval_nt<B>  Type;
  typedef Uncertain<CGAL::Sign> Sign;
  typedef Uncertain<bool> Boolean;
  typedef Uncertain<CGAL::Comparison_result> Comparison_result; 

    class Abs
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::abs( x );
        }
    };

    class Sgn
        : public std::unary_function< Type, Uncertain< ::CGAL::Sign > > {
      public:
        Uncertain< ::CGAL::Sign > operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::sign( x );
        }
    };

    class Is_positive
      : public std::unary_function< Type, Uncertain<bool> > {
      public:
        Uncertain<bool> operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_positive( x );
        }
    };

    class Is_negative
      : public std::unary_function< Type, Uncertain<bool> > {
      public:
        Uncertain<bool> operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_negative( x );
        }
    };

    class Compare
      : public std::binary_function< Type, Type, Comparison_result > {
      public:
      Comparison_result operator()( const Type& x, const Type& y ) const {
        return INTERN_INTERVAL_NT::compare( x, y );
      }
      CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
          Comparison_result )
    };

    class To_double
      : public std::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::to_double( x );
        }
    };

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
            return INTERN_INTERVAL_NT::to_interval( x );
        }
    };

    class Is_finite
      : public std::unary_function< Type, Boolean > {
      public :
        Boolean operator()( const Type& x ) const {
          return CGAL_NTS is_finite(-x.inf() ) && CGAL_NTS is_finite( x.sup() );
        }
    };

};

// Algebraic structure traits
template< bool B >
class Algebraic_structure_traits< Interval_nt<B> >
  : public Algebraic_structure_traits_base< Interval_nt<B>,
                                            Field_with_sqrt_tag >  {
  public:
    typedef Interval_nt<B>      Type;
    typedef Tag_false           Is_exact;
    typedef Tag_true            Is_numerical_sensitive;
    typedef Uncertain<bool>     Boolean; 

    class Is_zero
      : public std::unary_function< Type, Boolean > {
      public:
        Boolean operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_zero( x );
        }
    };

    class Is_one
      : public std::unary_function< Type, Boolean > {
      public:
        Boolean operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::is_one( x );
        }
    };

    class Square
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::square( x );
        }
    };

    class Sqrt
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return INTERN_INTERVAL_NT::sqrt( x );
        }
    };

    struct Is_square
        :public std::binary_function<Interval_nt<B>,Interval_nt<B>&,Boolean >
    {
        Boolean operator()(const Interval_nt<B>& x) const {
            return INTERN_INTERVAL_NT::is_positive( x );
        }

        Boolean operator()(
                const Interval_nt<B>& x,
                Interval_nt<B>      & result) const {
            Boolean is_positive = INTERN_INTERVAL_NT::is_positive( x );
            if ( is_positive.inf() == true ){
                typename Algebraic_structure_traits<Interval_nt<B> >::Sqrt sqrt;
                result = sqrt(x);
            }else{
                typename Real_embeddable_traits<Interval_nt<B> >::Abs  abs;
                typename Algebraic_structure_traits<Interval_nt<B> >::Sqrt sqrt;
                result = sqrt(abs(x));
            }
            return is_positive;
        }
    };

  class Divides
    : public std::binary_function< Type, Type, Boolean > { 
  public:
    Boolean operator()( const Type& x, const Type&) const {
      return ! Is_zero()(x);
    } 
    // second operator computing q
    Boolean operator()( const Type& x, const Type& y, Type& q) const {
      if (! Is_zero()(x) )
        q  = y/x ;
      return Boolean(true);
    }
  };
  
};


// COERCION_TRAITS BEGIN
template < class A, class B , int > struct Coercion_traits_for_level;
template < class A, class B, class C> struct Coercion_traits_interval_nt;

template<class A ,bool P >
struct Coercion_traits_for_level<A,Interval_nt<P>,CTL_INTERVAL>
    :public Coercion_traits_interval_nt<A,Interval_nt<P>,
            typename Real_embeddable_traits<A>::Is_real_embeddable>{};

template<class A , bool P>
struct Coercion_traits_for_level<Interval_nt<P>,A,CTL_INTERVAL>
    :public Coercion_traits_for_level<A,Interval_nt<P>, CTL_INTERVAL>{};

template<class A , bool P >
struct Coercion_traits_interval_nt<A, Interval_nt<P>,Tag_false>
    :public Coercion_traits_for_level<A,Interval_nt<P>,0>{};

template<class A , bool P>
struct Coercion_traits_interval_nt<A, Interval_nt<P>, Tag_true>{
    typedef Tag_true Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef Interval_nt<P> Type;
    struct Cast {
        typedef Interval_nt<P> result_type;
        Interval_nt<P> inline operator()(const Interval_nt<P>& x ) const {
            return x;
        }
        Interval_nt<P> inline operator()(const A& x ) const {
            return typename Real_embeddable_traits<A>::To_interval()(x);
        }
    };
};

// COERCION_TRAITS END

template< bool B >
class Interval_traits< Interval_nt<B> >
  : public internal::Interval_traits_base< Interval_nt<B> >  {
public: 
  typedef Interval_traits<Interval_nt<B> > Self; 
  typedef Interval_nt<B> Interval; 
  typedef double Bound; 
  typedef CGAL::Tag_false With_empty_interval; 
  typedef CGAL::Tag_true  Is_interval; 

 struct Construct :public std::binary_function<Bound,Bound,Interval>{
    Interval operator()( const Bound& l,const Bound& r) const {
      CGAL_precondition( l < r ); 
      return Interval(l,r);
    }
  };

  struct Lower :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.inf();
    }
  };

  struct Upper :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.sup();
    }
  };

  struct Width :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return width(a); 
    }
  };

  struct Median :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return (Lower()(a)+Upper()(a))/2.0;
    }
  };
    
  struct Norm :public std::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return magnitude(a);
    }
  };

  struct Singleton :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return Lower()(a) == Upper()(a);
    }
  };

  struct Zero_in :public std::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return Lower()(a) <= 0  &&  0 <= Upper()(a);
    }
  };

  struct In :public std::binary_function<Bound,Interval,bool>{
    bool operator()( Bound x, const Interval& a ) const {
      return Lower()(a) <= x && x <= Upper()(a);
    }
  };

  struct Equal :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.is_same(b);
    }
  };
    
  struct Overlap :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.do_overlap(b);
    }
  };
    
  struct Subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Lower()(b) <= Lower()(a) && Upper()(a) <= Upper()(b) ;  
    }
  };
    
  struct Proper_subset :public std::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Subset()(a,b) && ! Equal()(a,b); 
    }
  };
    
  struct Hull :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      return Interval( 
          std::make_pair(
              min BOOST_PREVENT_MACRO_SUBSTITUTION (Lower()(a),b.inf()), 
              max BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(a),b.sup())));
    }
  };
    
  
//  struct Empty is Null_functor 
  
  struct Intersection :public std::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      Bound l(max BOOST_PREVENT_MACRO_SUBSTITUTION (Lower()(a),Lower()(b)));
      Bound u(min BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(a),Upper()(b)));
      if(u < l ) throw Exception_intersection_is_empty();
      return Construct()(l,u);
    }
  };
};

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<bool b> struct NumTraits<CGAL::Interval_nt<b> >
  {
    typedef CGAL::Interval_nt<b> Real;
    typedef CGAL::Interval_nt<b> NonInteger;
    typedef CGAL::Interval_nt<b> Nested;

    static inline Real epsilon() { return 0; }

    // Costs could depend on b.
    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 0,
      ReadCost = 2,
      AddCost = 2,
      MulCost = 10
    };
  };
}

#endif // CGAL_INTERVAL_NT_H
