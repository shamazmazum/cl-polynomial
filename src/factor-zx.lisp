(defpackage cl-polynomial/factor-zx
  (:use #:cl)
  (:local-nicknames (#:sera   #:serapeum)
                    (#:alex   #:alexandria)
                    (#:si     #:stateless-iterators)
                    (#:primes #:cl-prime-maker)
                    (#:u      #:cl-polynomial/util)
                    (#:p      #:cl-polynomial/polynomial)
                    (#:ff     #:cl-polynomial/factor))
  (:export #:lift-factors
           #:remove-content
           #:suitable-primes
           #:suitable-bound))
(in-package :cl-polynomial/factor-zx)

(sera:-> remove-content (p:polynomial)
         (values p:polynomial fixnum &optional))
(defun remove-content (f)
  "For polynomial \\(f \\in \\mathbb{Z}[x]\\) return its primitive
part and content."
  (let ((content (reduce
                  (lambda (acc m)
                    (declare (type u:monomial m))
                    (gcd acc (cdr m)))
                  (p:polynomial-coeffs f)
                  :initial-value 0)))
    (values
     (p:map-poly (lambda (x) (/ x content)) f)
     content)))

(sera:-> replace-lc (p:polynomial fixnum)
         (values p:polynomial &optional))
(defun replace-lc (f c)
  "Replace leading coefficient of non-zero @c(f) with @c(c)."
  (if (p:polynomial= f p:+zero+)
      (p:scale p:+one+ c)
      (p:polynomial
       (destructuring-bind ((d . %c) &rest rest)
           (p:polynomial-coeffs f)
         (declare (ignore %c))
         (cons (cons d c) rest)))))

(sera:-> lift-factors
         (p:polynomial p:polynomial p:polynomial u:prime alex:non-negative-fixnum)
         (values p:polynomial p:polynomial boolean alex:non-negative-fixnum &optional))
(defun lift-factors (f f1 f2 p d)
  "For a primitive polynomial \\(f \\in \\mathbb{Z}[x]\\) with leading
coefficient > 0 and \\(f_1, f_2 \\in \\mathbb{Z}_p[x]\\), \\(p\\)
being prime, find polynomials \\(\\hat{f_1}, \\hat{f_2} \\in
\\mathbb{Z}[x]\\) such that \\(f = \\hat{f_1} \\hat{f_2}\\) if \\(f =
f_1 f_2 \\mod p\\). \\(d\\) is a maximal absolute value of
coefficients in \\(f\\) or its factors.

The first two values returned are the desired factors. The third value
is a boolean being equal to @c(T) if the algorithm has successfully
found a solution in \\(\\mathbb{Z}[x]\\). If this value is @c(NIL)
then the algorithm is not successful and \\(f = \\hat{f_1} \\hat{f_2}
\\mod p^N\\) where \\(N\\) is the forth returned value."
  ;; This algorithm is called the linear Hensel lifting. Having the
  ;; factorization in ℤ_{p^k}[x] it finds the factorization in
  ;; ℤ_{p^{k+1}}[x] until the factorization in ℤ[x] is found or a
  ;; limit of iterations is reached.

  ;; KLUDGE: Check that leading coefficient is > 0. Otherwise
  ;; factorization will be incorrect by -1 multiple.
  (unless (> (p:leading-coeff f) 0)
    (error "Leading coefficient < 0: ~a" f))
  (let* ((lc (p:leading-coeff f))
         ;; p^max-steps > 2D*lc must hold
         (max-steps (1+ (ceiling (log (abs (* 2 d lc)) p))))
         ;; Adjust f, f1 and f2 to correctly solve non-monic case
         ;; f ∈ ℤ[x], f₁, f₂ ∈ ℤ_p[x]
         (f  (identity (p:scale f  (* lc))))
         (f1 (p:modulo (p:scale f1 (* lc (u:invert-integer (p:leading-coeff f1) p))) p))
         (f2 (p:modulo (p:scale f2 (* lc (u:invert-integer (p:leading-coeff f2) p))) p)))
    (multiple-value-bind (gcd s d)
        (p:gcdex f1 f2 p)
      (declare (ignore gcd))
      (labels ((%lift-factors (%f1 %f2 q step)
                 (let* ((%f1 (replace-lc %f1 lc))
                        (%f2 (replace-lc %f2 lc))
                        (prod (p:multiply %f1 %f2))
                        (diff (p:subtract f prod)))
                   (if (or (p:polynomial= diff p:+zero+)
                           (= step max-steps))
                       (values (remove-content %f1)
                               (remove-content %f2)
                               (p:polynomial= diff p:+zero+)
                               step)
                       (let* ((rhs (p:map-poly (lambda (x) (/ x q)) diff))
                              (%δf2 (p:modulo (p:multiply s rhs) p))
                              (%δf1 (p:modulo (p:multiply d rhs) p)))
                         (multiple-value-bind (quo δf2)
                             (p:divide %δf2 %f2 p)
                           (let ((δf1 (p:modulo (p:add (p:multiply %f1 quo) %δf1) p)))
                             (%lift-factors (p:add %f1 (p:scale δf1 q))
                                            (p:add %f2 (p:scale δf2 q))
                                            (* p q)
                                            (1+ step)))))))))
        (%lift-factors f1 f2 p 0)))))

(sera:-> suitable-primes (p:polynomial)
         (values si:iterator &optional))
(defun suitable-primes (polynomial)
  "Return an iterator for suitable primes for reducing a factorization
in \\(\\mathbb{Z}[x]\\) to a factorization in \\(\\mathbb{F}_p[x]\\)."
  (assert (not (p:polynomial= polynomial p:+zero+)))
  ;; A suitable prime is the first prime which does not divide LC
  (let ((lc (p:leading-coeff polynomial)))
    (si:drop-while
     (lambda (prime) (zerop (rem lc prime)))
     (si:imap
      (lambda (n) (primes:get-nth-prime n))
      (si:count-from 2)))))

(sera:-> suitable-bound (p:polynomial)
         (values alex:positive-fixnum &optional))
(defun suitable-bound (polynomial)
  (let ((degree (p:degree polynomial)))
    (* (floor (sqrt (1+ degree)))
       (expt 2 degree)
       (reduce #'max (mapcar (alex:compose #'abs #'cdr)
                             (p:polynomial-coeffs polynomial))))))

(sera:-> possible-factors (p:polynomial)
         (values list &optional))
(defun possible-factors (polynomial)
  ;; Polynomial must be content-free and square-free
  (labels ((find-prime (primes-source)
             (multiple-value-bind (prime primes-source)
                 (si:consume-one primes-source)
               (let* ((f (p:monic-polynomial (p:modulo polynomial prime) prime))
                      (sf-factors (ff:square-free f prime)))
                 (if (and (= (length sf-factors) 1)
                          (= (caar sf-factors) 1))
                     (values prime f)
                     (find-prime primes-source))))))
    ;; Find a prime which gets a square-free factorization in the
    ;; corresponding finite field and get POLYNOMIAL modulo that
    ;; prime.
    (multiple-value-bind (p f)
        (find-prime (suitable-primes polynomial))
      (let ((factors (ff:berlekamp-factor f p))
            (bound (suitable-bound polynomial)))
        (mapcar
         (lambda (f1)
           ;; KLUDGE: This seems to be inefficient. I just loose one of
           ;; lifted factors in this procedure.
           (let ((f2 (p:modulo (apply #'p:multiply (remove f1 factors)) p)))
             (lift-factors polynomial f1 f2 p bound)))
         factors)))))
