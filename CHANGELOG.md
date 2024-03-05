# Changelog

## Version 0.2.3

* Improvement: Do not use floating point functions when computing Mignotte
  bound.
* Improvement: When factoring in ℤ[x], try different primes if the first prime
  produces to much factors in 𝔽_p[x].

## Version 0.2.2

* New feature: POLYNOMIAL/= function which behaves like a composition of NOT and
  POLYNOMIAL= (saves typing)
* Improvement: use a new source of primes, suitable for parallel applications.
* Optimization: Fewer call to the Hensel lifting. Now all lifting is performed
  before the recombination step.
* Optimization: Faster POLYNOMIAL= and EXPT.

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
