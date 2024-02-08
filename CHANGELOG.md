# Changelog

## Version 0.2.1

* New feature: Cantor-Zassenhaus algorithm for factorization of polynomials in
  big finite fields. An algorithm for factorization can be specified passing the
  third argument to CL-POLYNOMIAL/FPX:FACTOR (can be :BERLEKAMP or
  :CANTOR-ZASSENHAUS).
* New feature: Computation of cyclotomic polynomials in ℤ[x].
* New feature: Basic number-theoretic functions MOEBIUS and TOTIENT in a new
  package CL-POLYNOMIAL/Z. They are rather slow (I use them for testing).
* Optimization: BERLEKAMP-FACTOR performs tests with reducing polynomials much
  faster.
* Optimization: Linear Hensel lifting is replaced with quadratic Hensel lifting
  which seems a bit faster, especially in cases where a polynomial does not have
  a non-trivial factorization in ℤ[x].
* Optimization: Addition of polynomials is slightly improved.
* Improvement: A new constant CL-POLYNOMIAL/POLYNOMIAL:+VARIABLE+ in which the
  polynomial X is stored.
  

## Version 0.2

Optimization and bugfix release.

* Optimization: FACTOR-SQUARE-FREE performs lesser amount of Hensel liftings and
  is much faster in some cases.
* Optimization: DESTRUCTURING-BIND was replaced with faster BIND-MONOMIAL where
  it is important for performance.
* Bug fix: Unreachable code is removed from FACTOR-SQUARE-FREE.
* Optimization: MOD-SYM is now a bit faster.
* Bug fix: SQUARE-FREE (in both packages) correctly handles all corner cases.

## Version 0.1

* Factorization of univariate polynomials over finite field (Berlekamp
  algorithm) and over integers (Zassenhaus algorithm).
