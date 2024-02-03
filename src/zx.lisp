;; Operations and factorization in ℤ[x]

(defpackage cl-polynomial/zx
  (:use #:cl)
  (:shadow #:gcd)
  (:local-nicknames (#:sera   #:serapeum)
                    (#:si     #:stateless-iterators)
                    (#:primes #:cl-prime-maker)
                    (#:u      #:cl-polynomial/util)
                    (#:p      #:cl-polynomial/polynomial)
                    (#:fpx    #:cl-polynomial/fpx))
  (:export #:lift-factors
           #:remove-content
           #:suitable-primes
           #:suitable-bound
           #:divide
           #:remainder
           #:gcd
           #:square-free
           #:factor
           #:irreduciblep))
(in-package :cl-polynomial/zx)

(sera:-> scale-divide (p:polynomial integer)
         (values p:polynomial &optional))
(declaim (inline scale-divide))
(defun scale-divide (poly a)
  (p:map-poly (lambda (x) (/ x a)) poly))

(sera:-> remove-content (p:polynomial)
         (values p:polynomial integer &optional))
(defun remove-content (f)
  "For polynomial \\(f \\in \\mathbb{Z}[x]\\) return its primitive
part and content."
  (let ((content (reduce
                  (lambda (acc m)
                    (declare (type u:monomial m))
                    (cl:gcd acc (cdr m)))
                  (p:polynomial-coeffs f)
                  :initial-value 0)))
    (values (scale-divide f content) content)))

(sera:-> replace-lc (p:polynomial integer)
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

(sera:-> fix-unit (p:polynomial integer u:prime-power)
         (values p:polynomial &optional))
(declaim (inline fix-unit))
(defun fix-unit (f lc q)
  (fpx:modulo (p:scale f (* lc (u:invert-integer (p:leading-coeff f) q))) q))

;; ========================================
;; Quadratic lifting
;; https://www.csd.uwo.ca/~mmorenom/CS874/Lectures/Newton2Hensel.html/node17.html
(sera:-> lifting-step
         (p:polynomial p:polynomial p:polynomial p:polynomial p:polynomial u:prime-power)
         (values p:polynomial p:polynomial p:polynomial p:polynomial &optional))
(defun lifting-step (f g h s p m)
  (let ((e (fpx:modulo (p:subtract f (p:multiply g h)) m)))
    (multiple-value-bind (q r)
        (fpx:divide (p:multiply s e) h m)
      (let* ((%g (fpx:modulo (p:add g (p:multiply p e) (p:multiply q g)) m))
             (%h (fpx:modulo (p:add h r) m))
             (b (fpx:modulo (p:subtract (p:add (p:multiply s %g) (p:multiply p %h)) p:+one+)
                            m)))
              (multiple-value-bind (c d)
                  (fpx:divide (p:multiply s b) %h m)
                (values
                 %g %h
                 (fpx:modulo (p:subtract s d) m)
                 (fpx:modulo (p:subtract p (p:multiply p b) (p:multiply c %g)) m)))))))

(sera:-> lift-factors (p:polynomial p:polynomial p:polynomial u:prime-power)
         (values p:polynomial p:polynomial boolean &optional))
(defun lift-factors (f g h q)
  "For a primitive polynomial \\(f \\in \\mathbb{Z}[x]\\) with leading
coefficient > 0 and \\(g, h \\in \\mathbb{Z}_q[x]\\) find polynomials
\\(\\hat{g}, \\hat{h} \\in \\mathbb{Z}[x]\\) such that \\(f = \\hat{g}
\\hat{h}\\) if \\(f = g h \\mod q\\).

The first two values returned are the desired factors. The third value
is a boolean being equal to @c(T) if the algorithm has successfully
found a solution in \\(\\mathbb{Z}[x]\\). If this value is @c(NIL)
the first two values can be ignored."
  ;; KLUDGE: Check that leading coefficient is > 0. Otherwise
  ;; factorization will be incorrect by -1 multiple.
  ;; FIXME: Is it still true?
  (unless (> (p:leading-coeff f) 0)
    (error "Leading coefficient < 0: ~a" f))
  (let* ((lc (p:leading-coeff f))
         (f (p:scale f lc))
         (g (fix-unit g lc q))
         (h (fix-unit h lc q)))
    (multiple-value-bind (gcd s p)
        (fpx:gcdex g h q)
      (declare (ignore gcd))
      (let ((bound (* 2 (suitable-bound f))))
        (labels ((%lf (g h s p q)
                   (let ((g (replace-lc g lc))
                         (h (replace-lc h lc))
                         (convp (p:polynomial= f (p:multiply g h))))
                     (if (or (> q bound) convp)
                         (values (remove-content g)
                                 (remove-content h)
                                 convp)
                         (let ((q (expt q 2)))
                           (multiple-value-bind (g h s p)
                               (lifting-step f g h s p q)
                             (%lf g h s p q)))))))
          (%lf g h s p q))))))

(sera:-> suitable-primes (p:polynomial)
         (values si:iterator &optional))
(defun suitable-primes (polynomial)
  "Return an iterator for suitable primes for reducing a factorization
in \\(\\mathbb{Z}[x]\\) to a factorization in \\(\\mathbb{F}_p[x]\\)."
  (assert (not (p:polynomial= polynomial p:+zero+)))
  ;; A suitable prime is the first prime which does not divide LC
  (let ((lc (p:leading-coeff polynomial)))
    (si:filter
     (lambda (prime)
       (not (zerop (rem lc prime))))
     (si:imap
      (lambda (n) (primes:get-nth-prime n))
      (si:count-from 2)))))

(sera:-> suitable-bound (p:polynomial)
         (values (integer 1) &optional))
(defun suitable-bound (polynomial)
  (let ((degree (p:degree polynomial)))
    (* (floor (sqrt (1+ degree)))
       (expt 2 degree)
       (reduce #'max (p:polynomial-coeffs polynomial)
               :key (lambda (x) (abs (cdr x)))))))

(sera:-> combinations (list (integer 0))
         (values list &optional))
(defun combinations (list n)
  (cond ((zerop n) (list nil))
        ((null list) nil)
        (t (append
            (mapcar (lambda (comb)
                      (cons (car list) comb))
                    (combinations (cdr list) (1- n)))
            (combinations (cdr list) n)))))

(sera:-> factor-square-free (p:polynomial)
         (values list &optional))
(defun factor-square-free (polynomial)
  ;; Polynomial must be content-free and square-free
  (labels ((find-prime (primes-source)
             ;; Find a suitable prime which does not divide the
             ;; leading coefficient and provides a square-free
             ;; factorization in a finite field (this factorization
             ;; exists because POLYNOMIAL itself is square-free).
             (multiple-value-bind (prime primes-source)
                 (si:consume-one primes-source)
               (let ((sf-factors (fpx:square-free (fpx:modulo polynomial prime) prime)))
                 (if (and (= (length sf-factors) 1)
                          (= (caar sf-factors) 1))
                     (values prime (cdar sf-factors))
                     (find-prime primes-source))))))
    (multiple-value-bind (p f)
        (find-prime (suitable-primes polynomial))
      (labels ((lift-factor (f f1 f2)
                 (lift-factors f f1 f2 p))
               (try-combinations (f combs acc all-factors)
                 ;; F ∈ ℤ[x]. Here ALL-FACTORS is a list of factors of F mod P and COMB
                 ;; contain all combinations of N elements from ALL-FACTORS. For example,
                 ;; F mod P factors like `abcd` and COMB contains combinations of 2
                 ;; elements. Then COMB is equal to `((ab) (ac) (ad) (bc) (bd)
                 ;; (cd))`. Each combination is lifted to ℤ[x] and tested if it divides
                 ;; F. If it does, it is removed from ALL-FACTORS and RECOMBINE is called
                 ;; again with F/∏COMB and ALL-FACTORS \ COMB.  This function does not
                 ;; return anything useful because RECOMBINE jumps directly to the end of
                 ;; FACTOR-SQUARE-FREE.
                 (when combs
                   (let* ((c1 (car combs))
                          (c2 (set-difference all-factors c1 :test #'p:polynomial=))
                          (f1 (apply #'p:multiply c1))
                          (f2 (apply #'p:multiply c2)))
                     (multiple-value-bind (f1l f2l convp)
                         (lift-factor f f1 f2)
                       (if convp
                           (recombine f2l c2 (cons f1l acc) (length c1))
                           (try-combinations f (cdr combs) acc all-factors))))))
               (recombine (f factors acc i)
                 ;; F ∈ ℤ[x]. Generate combinations from factors of F mod P and try to
                 ;; lift these combinations to factors of F.  For example, suppose F mod P
                 ;; factors as abcd, where a,b,c and d ∈ 𝔽_p[x].  Firstly combinations (a)
                 ;; (b) (c) and (d) are constructed and we try to lift them to factors of
                 ;; F. If we are not successful, combinations (ab) (ac) ... (cd) are
                 ;; constructed and one of these combinations can be lifted to a factor of
                 ;; F.
                 (try-combinations f (combinations factors i) acc factors)
                 (if (or (p:polynomial= f p:+one+)
                         (= i (length factors)))
                     (return-from factor-square-free acc)
                     (recombine f factors acc (1+ i)))))
        (recombine polynomial (fpx:berlekamp-factor f p) nil 1)))))

;; DIVIDE and GCD have their own versions for polynomials over
;; integers. DIVIDE can be used only for exact division externally.
(sera:-> divide (p:polynomial p:polynomial)
         (values p:polynomial p:polynomial &optional))
(defun divide (poly1 poly2)
  (let ((degree (p:degree poly2))
        (lc (p:leading-coeff poly2)))
    (if (zerop degree)
        ;; Division by a constant is a special case
        (values
         (scale-divide poly1 lc)
         p:+zero+)
        (labels ((division-step (quotient-coeffs remainder)
                   (declare (type list quotient-coeffs)
                            (type p:polynomial remainder))
                   (let* ((remainder-coeffs (p:polynomial-coeffs remainder))
                          (d (caar remainder-coeffs))
                          (c (cdar remainder-coeffs)))
                     (if (or (< (p:degree remainder) degree)
                             (not (zerop (rem c lc))))
                         (values (p:polynomial (reverse quotient-coeffs))
                                 remainder)
                         (let* ((quotient-degree (- d degree))
                                (quotient-coeff (/ c lc))
                                (monomial (cons quotient-degree quotient-coeff)))
                           (division-step
                            (cons monomial quotient-coeffs)
                            (p:subtract remainder
                                        (p::multiply-monomial monomial poly2))))))))
          (division-step nil poly1)))))

(sera:-> remainder (p:polynomial p:polynomial)
         (values p:polynomial &optional))
(declaim (inline remainder))
(defun remainder (poly1 poly2)
  (nth-value 1 (divide poly1 poly2)))

(sera:-> gcd-primitive (p:polynomial p:polynomial)
         (values p:polynomial &optional))
(defun gcd-primitive (poly1 poly2)
  (p:positive-lc
   (remove-content
    (labels ((%gcd (f1 f2 &optional %d %γ %ψ contp)
               (if (p:polynomial= f2 p:+zero+) f1
                   (let* ((d (- (p:degree f1) (p:degree f2)))
                          (γ (p:leading-coeff f2))
                          (ψ (if contp (/ (expt (- %γ) %d)
                                          (expt %ψ (1- %d)))
                                 -1))
                          (β (if contp (- (* %γ (expt ψ d)))
                                 (expt -1 (1+ d)))))
                     (%gcd
                      f2 (scale-divide (remainder (p:scale f1 (expt γ (1+ d))) f2) β)
                      d γ ψ t)))))
      (if (> (p:degree poly1)
             (p:degree poly2))
          (%gcd poly1 poly2)
          (%gcd poly2 poly1))))))

(sera:-> gcd (p:polynomial p:polynomial)
         (values p:polynomial &optional))
(defun gcd (poly1 poly2)
  "Calculate GCD of two polynomials in \\(\\mathbb{Z}[x]\\)."
  (multiple-value-bind (p1 c1)
      (remove-content poly1)
    (multiple-value-bind (p2 c2)
        (remove-content poly2)
      (p:scale
       (gcd-primitive p1 p2)
       (cl:gcd c1 c2)))))

(sera:-> square-free (p:polynomial)
         (values list integer &optional))
(defun square-free (polynomial)
  "Perform square-free factorization for a polynomial \\(f \\in
\\mathbb{Z}[x]\\). The first returned value is a list of tuples \\((d_i
. f_i)\\), so the polynomial is equal to \\(\\sum_i f_i^{d_i}\\)
multiplied by a second returned value."
  (multiple-value-bind (f cont)
      (remove-content (p:positive-lc polynomial))
    (labels ((%square-free (b d acc deg)
               (if (p:polynomial= b p:+one+) acc
                   (let* ((gcd (gcd b d))
                          (%b (divide b gcd))
                          (c  (divide d gcd))
                          (%d (p:subtract c (p:derivative %b))))
                     (%square-free %b %d
                                   (if (p:polynomial= gcd p:+one+) acc
                                       (cons (cons deg gcd) acc))
                                   (1+ deg))))))
      (cond
        ((p:polynomial= polynomial p:+zero+)
         (error "Cannot factor zero polynomial"))
        ((= (p:degree polynomial) 0)
         (values (list (cons 1 p:+one+))
                 (p:leading-coeff polynomial)))
        (t
         (let ((factors (%square-free f (p:derivative f) nil 0)))
           (values (remove 0 factors :key #'car)
                   (* cont (signum (p:leading-coeff polynomial))))))))))

(sera:-> factor (p:polynomial)
         (values list integer &optional))
(defun factor (polynomial)
  "Factor a polynomial in \\(\\mathbb{Z}[x]\\)."
  (if (p:polynomial= polynomial p:+zero+)
      (error "Cannot factor zero polynomial")
      (multiple-value-bind (f c)
          (square-free polynomial)
        (values
         (reduce #'append f
                 :key
                 (lambda (factor)
                   (destructuring-bind (m . f) factor
                     (mapcar
                      (lambda (factor)
                        (cons m factor))
                      (factor-square-free f)))))
         c))))

(sera:-> irreducible-p (p:polynomial)
         (values boolean &optional))
(defun irreduciblep (polynomial)
  "Test if polynomial is irreducible in \\(\\mathbb{Z}[x]\\)."
  (multiple-value-bind (p c)
      (remove-content polynomial)
    (if (= c 1)
        (destructuring-bind ((m . p) &rest rest)
            (square-free p)
          (and (null rest)
               (= m 1)
               (= (length (factor p)) 1))))))
