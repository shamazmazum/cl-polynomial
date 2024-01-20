# Changelog

## Version 0.1

* Factorization of univariate polynomials over finite field (Berlekamp
  algorithm) and over integers (Zassenhaus algorithm).

## Version 0.2

Optimization and bugfix release.

* Optimization: FACTOR-SQUARE-FREE performs lesser amount of Hensel liftings and
  is much faster in some cases.
* Optimization: DESTRUCTURING-BIND was replaced with faster BIND-MONOMIAL where
  it is important for performance.
* Bug fix: Unreachable code is removed from FACTOR-SQUARE-FREE.
* Optimization: MOD-SYM is now a bit faster.
* Bug fix: SQUARE-FREE (in both packages) correctly handles all corner cases.
