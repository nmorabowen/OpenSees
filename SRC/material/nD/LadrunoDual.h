/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 11/2026  (WP-262)
//
// LadrunoDual — forward-mode automatic differentiation by dual numbers.
//
// A Dual<N> carries a value v and N directional derivatives d[0..N-1]:
//
//     x = v + sum_i d_i eps_i ,   eps_i eps_j = 0
//
// so every arithmetic operation propagates the chain rule exactly (no step
// size, no truncation, no round-off amplification):  f(x) = f(v) + f'(v) d eps.
// Seeding one input per direction (Dual<6>::variable(strain_i, i)) and pushing
// it through a kernel written ONCE as `template <class T>` returns the full
// Jacobian of the kernel's outputs in one pass.
//
// Design (WP-262; research bench gitAPE/dual-numbers-opensees, survey 01 §5):
//   * fixed N, plain array storage, no heap (the Sacado SFad model);
//   * nestable: Dual<N, Dual<M>> gives second derivatives (energy Hessians);
//   * C++17, header-only, OpenSees-free; compiles under MSVC and g++ (the
//     kernel checks build with g++, the binaries with MSVC);
//   * NO implicit conversion to double: a Dual can never silently lose its
//     derivative. Branches read the value explicitly through value(x);
//   * the value part follows EXACTLY the operations of the double code path,
//     so value(f(Dual)) == f(double) bit for bit for the operations below.
//
// Math functions are found by argument-dependent lookup: kernel code calls
// unqualified sqrt/exp/... and gets ::sqrt for double, ladruno_dual::sqrt for
// Dual. Kernels must therefore NOT qualify them as std::sqrt.
//
// Known limits (research/05_downsides.md D4-D6):
//   * sqrt at 0, pow(0, a<1) and similar: the value is finite but the
//     derivative is infinite. Kernels keep their existing guards (normFloor).
//   * fabs/fmax/fmin and every value-based branch return the derivative of the
//     active branch: one-sided at the kink itself.
//   * Iterative solves: converge, then take ONE Newton step in Dual at the
//     converged point (implicit-function theorem) so the derivative does not
//     lag the value by an iterate.

#ifndef LadrunoDual_h
#define LadrunoDual_h

#include <math.h>
#include <type_traits>

namespace ladruno_dual {

template <int N, class S = double>
struct Dual {
  static_assert(N >= 1, "Dual<N> needs at least one derivative direction");

  S v;      // value
  S d[N];   // derivatives along the N seeded directions

  Dual() : v(0.0) { for (int i = 0; i < N; i++) d[i] = S(0.0); }
  Dual(const S& x) : v(x) { for (int i = 0; i < N; i++) d[i] = S(0.0); }
  template <class A, typename std::enable_if<std::is_arithmetic<A>::value, int>::type = 0>
  Dual(A x) : v(S(x)) { for (int i = 0; i < N; i++) d[i] = S(0.0); }

  // an independent variable: value x, unit derivative along direction i
  static Dual variable(const S& x, int i)
  {
    Dual r(x);
    r.d[i] = S(1.0);
    return r;
  }

  Dual& operator+=(const Dual& b) { v += b.v; for (int i = 0; i < N; i++) d[i] += b.d[i]; return *this; }
  Dual& operator-=(const Dual& b) { v -= b.v; for (int i = 0; i < N; i++) d[i] -= b.d[i]; return *this; }
  Dual& operator*=(const Dual& b)
  {
    for (int i = 0; i < N; i++) d[i] = d[i]*b.v + v*b.d[i];
    v *= b.v;
    return *this;
  }
  Dual& operator/=(const Dual& b)
  {
    const S inv = S(1.0)/b.v;
    v = v/b.v;                                     // the value op of the double path
    for (int i = 0; i < N; i++) d[i] = (d[i] - v*b.d[i])*inv;
    return *this;
  }
};

// ---------------------------------------------------------------- traits ---
template <class T> struct is_dual : std::false_type {};
template <int N, class S> struct is_dual<Dual<N, S> > : std::true_type {};

#define LADRUNO_DUAL_ARITH(A) typename std::enable_if<std::is_arithmetic<A>::value, int>::type = 0

// value(): the plain double underneath any (possibly nested) Dual; identity on double
inline double value(double x) { return x; }
template <int N, class S> inline double value(const Dual<N, S>& x) { return value(x.v); }

// zero_value(): set the innermost scalar to 0 and keep every derivative part.
// Used for the implicit-function-theorem step at a converged local solve:
//     corr = R/dR;  zero_value(corr);  x -= corr;
// corrects x's derivative to -R_eps/R_x while leaving its value bit-unchanged.
inline void zero_value(double& x) { x = 0.0; }
template <int N, class S> inline void zero_value(Dual<N, S>& x) { zero_value(x.v); }

// ---------------------------------------------------------- arithmetic ----
template <int N, class S> inline Dual<N, S> operator+(const Dual<N, S>& a) { return a; }
template <int N, class S> inline Dual<N, S> operator-(const Dual<N, S>& a)
{
  Dual<N, S> r(-a.v);
  for (int i = 0; i < N; i++) r.d[i] = -a.d[i];
  return r;
}

template <int N, class S>
inline Dual<N, S> operator+(const Dual<N, S>& a, const Dual<N, S>& b)
{
  Dual<N, S> r(a.v + b.v);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i] + b.d[i];
  return r;
}
template <int N, class S>
inline Dual<N, S> operator-(const Dual<N, S>& a, const Dual<N, S>& b)
{
  Dual<N, S> r(a.v - b.v);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i] - b.d[i];
  return r;
}
template <int N, class S>
inline Dual<N, S> operator*(const Dual<N, S>& a, const Dual<N, S>& b)
{
  Dual<N, S> r(a.v * b.v);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*b.v + a.v*b.d[i];
  return r;
}
template <int N, class S>
inline Dual<N, S> operator/(const Dual<N, S>& a, const Dual<N, S>& b)
{
  Dual<N, S> r(a.v / b.v);
  const S inv = S(1.0)/b.v;
  for (int i = 0; i < N; i++) r.d[i] = (a.d[i] - r.v*b.d[i])*inv;
  return r;
}

// Mixed Dual / arithmetic scalar. Spelled out because template argument
// deduction never applies the double -> Dual converting constructor.
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator+(const Dual<N, S>& a, A b) { Dual<N, S> r(a); r.v = a.v + b; return r; }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator+(A a, const Dual<N, S>& b) { Dual<N, S> r(b); r.v = a + b.v; return r; }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator-(const Dual<N, S>& a, A b) { Dual<N, S> r(a); r.v = a.v - b; return r; }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator-(A a, const Dual<N, S>& b)
{
  Dual<N, S> r(a - b.v);
  for (int i = 0; i < N; i++) r.d[i] = -b.d[i];
  return r;
}
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator*(const Dual<N, S>& a, A b)
{
  Dual<N, S> r(a.v * b);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i] * b;
  return r;
}
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator*(A a, const Dual<N, S>& b)
{
  Dual<N, S> r(a * b.v);
  for (int i = 0; i < N; i++) r.d[i] = a * b.d[i];
  return r;
}
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator/(const Dual<N, S>& a, A b)
{
  Dual<N, S> r(a.v / b);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i] / b;
  return r;
}
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> operator/(A a, const Dual<N, S>& b)
{
  Dual<N, S> r(a / b.v);                         // d(a/b) = -(a/b) db / b
  const S inv = S(1.0)/b.v;
  for (int i = 0; i < N; i++) r.d[i] = -(r.v*b.d[i])*inv;
  return r;
}

// ------------------------------------------------ comparisons (on value) ---
#define LADRUNO_DUAL_CMP(OP)                                                          \
  template <int N, class S>                                                           \
  inline bool operator OP(const Dual<N, S>& a, const Dual<N, S>& b) { return a.v OP b.v; } \
  template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>                           \
  inline bool operator OP(const Dual<N, S>& a, A b) { return a.v OP b; }              \
  template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>                           \
  inline bool operator OP(A a, const Dual<N, S>& b) { return a OP b.v; }
LADRUNO_DUAL_CMP(<)
LADRUNO_DUAL_CMP(>)
LADRUNO_DUAL_CMP(<=)
LADRUNO_DUAL_CMP(>=)
LADRUNO_DUAL_CMP(==)
LADRUNO_DUAL_CMP(!=)
#undef LADRUNO_DUAL_CMP

// -------------------------------------------------------------- math ------
// Each one: value through the inner type's own function (so nesting works and
// the double value path is the plain <math.h> call), derivative by chain rule.
template <int N, class S>
inline Dual<N, S> sqrt(const Dual<N, S>& a)
{
  using ::sqrt;
  Dual<N, S> r(sqrt(a.v));
  const S k = S(0.5)/r.v;                        // infinite at a.v == 0 (see header)
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*k;
  return r;
}
template <int N, class S>
inline Dual<N, S> exp(const Dual<N, S>& a)
{
  using ::exp;
  Dual<N, S> r(exp(a.v));
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*r.v;
  return r;
}
template <int N, class S>
inline Dual<N, S> log(const Dual<N, S>& a)
{
  using ::log;
  Dual<N, S> r(log(a.v));
  const S inv = S(1.0)/a.v;
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*inv;
  return r;
}
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> pow(const Dual<N, S>& a, A b)
{
  using ::pow;
  Dual<N, S> r(pow(a.v, b));
  const S k = S(b)*pow(a.v, b - 1.0);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*k;
  return r;
}
template <int N, class S>
inline Dual<N, S> sin(const Dual<N, S>& a)
{
  using ::sin; using ::cos;
  Dual<N, S> r(sin(a.v));
  const S c = cos(a.v);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*c;
  return r;
}
template <int N, class S>
inline Dual<N, S> cos(const Dual<N, S>& a)
{
  using ::sin; using ::cos;
  Dual<N, S> r(cos(a.v));
  const S s = -sin(a.v);
  for (int i = 0; i < N; i++) r.d[i] = a.d[i]*s;
  return r;
}
// fabs / fmax / fmin: derivative of the active branch (one-sided at the kink).
template <int N, class S>
inline Dual<N, S> fabs(const Dual<N, S>& a) { return (a.v < 0.0) ? -a : a; }
template <int N, class S>
inline Dual<N, S> fmax(const Dual<N, S>& a, const Dual<N, S>& b) { return (a.v >= b.v) ? a : b; }
template <int N, class S>
inline Dual<N, S> fmin(const Dual<N, S>& a, const Dual<N, S>& b) { return (a.v <= b.v) ? a : b; }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> fmax(const Dual<N, S>& a, A b) { return (a.v >= b) ? a : Dual<N, S>(b); }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> fmax(A a, const Dual<N, S>& b) { return (b.v >= a) ? b : Dual<N, S>(a); }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> fmin(const Dual<N, S>& a, A b) { return (a.v <= b) ? a : Dual<N, S>(b); }
template <int N, class S, class A, LADRUNO_DUAL_ARITH(A)>
inline Dual<N, S> fmin(A a, const Dual<N, S>& b) { return (b.v <= a) ? b : Dual<N, S>(a); }

#undef LADRUNO_DUAL_ARITH

} // namespace ladruno_dual

#endif
