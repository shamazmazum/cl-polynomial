(defpackage polynomial
  (:use #:cl)
  (:shadow #:gcd #:constantp)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:si   #:stateless-iterators))
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
           #:divide
           #:remainder
           #:monic-polynomial
           #:monicp
           #:constantp
           #:gcd
           #:derivative
           #:square-free
           #:reducing-polynomials
           #:irreduciblep
           #:factor))
