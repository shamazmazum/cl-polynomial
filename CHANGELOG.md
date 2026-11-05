# Changelog

## Version 0.4.1

* Improvement: Fermat primality test, CL-POLYNOMIAL/Z:MAYBE-PRIME-P was added.
* Improvement: Cantor-Zassenhaus algorithm for factorization in a finite field
  with big number of elements is now used in CL-POLYNOMIAL/ZX:FACTOR when
  appropriate.
* Improvement: Factor recombination from a finite ring to ℤ[x] is somewhat
  faster.
* Bug fix: CL-POLYNOMIAL/FPX:IRREDUCIBLEP returns T when called with constant
  polynomials.
* Bug fix: MAP-POLY can be used with functions which return zero (strictly
  speaking, POLYNOMIAL is not a functor).

## Version 0.4.0

* Incompatible change: FFT is now calculated for arrays of fixnums rather than
  integers.
* FFT: the second and the third arguments are flipped (i.e. you must write `(fft
  array ω p)` and not `(fft array p ω)`. This is for consistency with other
  functions in this library.
* Add FFT! and IFFT! functions which modify their argument ("true" inplace FFT).

## Version 0.3.0

New features:

* RECIPROCAL and PALYNDROMEP functions were added.
* Evaluation (any ring) and interpolation (over finite fields) of polynomials.
* FFT as a special case of above.

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
