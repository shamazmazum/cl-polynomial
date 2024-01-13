(defpackage polynomial
  (:use #:cl)
  (:shadow #:gcd #:constantp)
  (:local-nicknames (#:sera   #:serapeum)
                    (#:alex   #:alexandria)
                    (#:si     #:stateless-iterators)
                    (#:primes #:cl-prime-maker))
  (:export #:polynomial
           #:polynomial=
           #:+zero+
           #:+one+
           #:degree
           #:leading-coeff
           #:list->polynomial
           #:sequence->polynomial
           #:polynomial->list
           #:map-poly
           #:negate
           #:add
           #:subtract
           #:multiply
           #:scale
           #:modulo
           #:invert-integer
           #:mod-sym
           #:divide
           #:remainder
           #:monic-polynomial
           #:monicp
           #:constantp
           #:gcd
           #:gcdex
           #:derivative
           #:square-free
           #:reducing-polynomials
           #:irreduciblep
           #:factor

           #:lift-factors
           #:remove-content
           #:suitable-primes
           #:suitable-bound))
