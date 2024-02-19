;; Operations and factorization in ℤ[x]

(defpackage cl-polynomial/zx
  (:use #:cl)
  (:shadow #:gcd)
  (:local-nicknames (#:sera   #:serapeum)
                    (#:alex   #:alexandria)
                    (#:si     #:stateless-iterators)
                    (#:u      #:cl-polynomial/util)
                    (#:p      #:cl-polynomial/polynomial)
                    (#:z      #:cl-polynomial/z)
                    (#:fpx    #:cl-polynomial/fpx))
  (:export #:lift-factors
           #:lifting-steps
           #:remove-content
           #:suitable-prime
           #:suitable-bound
           #:divide
           #:remainder
           #:gcd
           #:square-free
           #:factor
           #:irreduciblep
           #:cyclotomic))
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
                    (declare (type p:monomial m))
                    (cl:gcd acc (cdr m)))
                  (p:polynomial-coeffs f)
                  :initial-value 0)))
    (values (scale-divide f content) content)))

;; ========================================
;; Quadratic lifting
;; https://www.csd.uwo.ca/~mmorenom/CS874/Lectures/Newton2Hensel.html/node17.html
(sera:-> lifting-step
         (p:polynomial p:polynomial p:polynomial p:polynomial p:polynomial u:prime-power)
         (values p:polynomial p:polynomial p:polynomial p:polynomial &optional))
(defun lifting-step (f g h s p m)
  "Having \\(f \\in \\mathbb{Z}\\), \\(g, h \\in \\mathbb{Z}_m\\),
\\(\\gcd(g, h) = 1\\), \\(sg + ph = 1 \\mod m\\) and \\(f = gh \\mod
m\\), find \\(g^*, h^*, s^*, p^*\\), so that \\(s^*g^* + p^*h^* = 1
\\mod m^2\\) and \\(f = g^*h^* \\mod m^2\\)."
  (unless (= (p:leading-coeff f) 1)
    (error "F is not monic: ~a" f))
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

(sera:-> lifting-steps ((integer 1) u:prime)
         (values (integer 1) (integer 0) &optional))
(defun lifting-steps (b p)
  (labels ((%count (n q)
             (if (>= q b) (values q n)
                 (%count (1+ n) (* q q)))))
    (%count 0 p)))

(sera:-> lift-two-factors (p:polynomial p:polynomial p:polynomial u:prime (integer 0))
         (values p:polynomial p:polynomial &optional))
(defun lift-two-factors (f g h p n)
  "For a primitive polynomial \\(f \\in \\mathbb{Z}[x]\\) with leading
coefficient > 0 and two relatively prime polynomials \\(g, h \\in
\\mathbb{F}_p[x]\\) find \\(g^*, h^* \\in \\mathbb{Z}_{p^{2^n}}[x]\\)
such that \\(f = g^* h^* \\mod p^{2^n}\\) if \\(f = g h \\mod p\\)."
  (unless (= (p:leading-coeff f) 1)
    (error "Leading coefficient ≠ 1: ~a" f))
  (multiple-value-bind (gcd s d)
      (fpx:gcdex g h p)
    (declare (ignore gcd))
    (labels ((%lf (g h s d q m)
               (if (= m n)
                   (values g h)
                   (let ((q (expt q 2)))
                     (multiple-value-bind (g h s d)
                         (lifting-step f g h s d q)
                       (%lf g h s d q (1+ m)))))))
      (%lf g h s d p 0))))

(sera:-> lift-factors (p:polynomial list u:prime (integer 0))
         (values list &optional))
(defun lift-factors (f gs p n)
  "For a monic polynomial \\(f \\in \\mathbb{Z}[x]\\) and a list of
its factors in \\(\\mathbb{F}_p[x]\\), lift these factors to
\\(\\mathbb{Z}_{p^{2^n}}[x]\\)."
  (unless gs
    (error "The list of factors is empty"))
  (labels ((%lift (f gs acc)
             (if (null gs) acc
                 (let ((g (car gs))
                       (h (fpx:modulo (apply #'p:multiply (cdr gs)) p)))
                   (multiple-value-bind (g h)
                       (lift-two-factors f g h p n)
                     (%lift h (cdr gs) (cons g acc)))))))
    (%lift f gs nil)))

(sera:-> suitable-prime (p:polynomial)
         (values u:prime &optional))
(defun suitable-prime (polynomial)
  "Return a suitable prime for reducing a factorization in
\\(\\mathbb{Z}[x]\\) to a factorization in \\(\\mathbb{F}_p[x]\\)."
  (assert (not (p:polynomial= polynomial p:+zero+)))
  ;; Find a suitable prime which does not divide the leading
  ;; coefficient and provides a square-free factorization in a finite
  ;; field (this factorization exists because POLYNOMIAL itself is
  ;; square-free).
  (let ((lc (p:leading-coeff polynomial)))
    (nth-value
     0 (si:consume-one
        (si:filter
         (lambda (prime)
           (and (not (zerop (rem lc prime)))
                (let ((sf-factors (fpx:square-free (fpx:modulo polynomial prime) prime)))
                  (and (= (length sf-factors) 1)
                       (= (caar sf-factors) 1)))))
         z:*prime-source*)))))

;; https://en.wikipedia.org/wiki/Landau-Mignotte_bound
(sera:-> suitable-bound (p:polynomial)
         (values (integer 1) &optional))
(defun suitable-bound (f)
  "Return a suitable bound for absolute values of coefficients of
factors of \\(f\\)."
  (let ((degree (p:degree f)))
    (* 2 (u:root (1+ degree) 2)
       (expt 2 degree)
       (reduce #'max (p:polynomial-coeffs f)
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

(defun factor-square-free (polynomial)
  (declare (optimize (speed 3)))
  ;; Polynomial must be content-free and square-free
  (labels ((restore-factor (comb f q)
             (remove-content
              (fpx:modulo (p:scale (apply #'p:multiply comb) (p:leading-coeff f)) q)))
           (try-combinations (f combs all-factors comb-length q)
             (declare (type alex:positive-fixnum comb-length))
             ;; F ∈ ℤ[x]. Here ALL-FACTORS is a list of factors of F mod Q in ℤ_q[x] and
             ;; COMB contain all combinations of COMB-LENGTH elements from
             ;; ALL-FACTORS. For example, F mod P factors like `abcd` and COMB contains
             ;; combinations of 2 elements. Then COMB is equal to `((ab) (ac) (ad) (bc)
             ;; (bd) (cd))`. Each combination is tested if it divides F. If it does, it is
             ;; removed from ALL-FACTORS and a newly found "true" factor is returned.
             (if (null combs)
                 (values f all-factors (1+ comb-length))
                 (let* ((c1 (car combs))
                        (c2 (set-difference all-factors c1 :test #'p:polynomial=))
                        (f1 (restore-factor c1 f q))
                        (f2 (restore-factor c2 f q)))
                   (if (p:polynomial= f (p:multiply f1 f2))
                       (values f2 c2 comb-length f1)
                       (try-combinations f (cdr combs) all-factors comb-length q)))))
           (recombine (f factors acc comb-length q)
             ;; F ∈ ℤ[x]. Generate combinations from factors of F mod Q in ℤ_q[x] and try
             ;; to detect products of which combinations are the factors of F in ℤ[x]. For
             ;; example, suppose F mod Q factors as abcd, where a,b,c and d ∈ ℤ_q[x].
             ;; Firstly combinations (a) (b) (c) and (d) are constructed and we check if
             ;; they are factors of F. If we are not successful, combinations (ab) (ac)
             ;; ... (cd) are constructed and checked for being the factors.
             (multiple-value-bind (f factors comb-length maybe-factor)
                 (try-combinations
                  f (combinations factors comb-length) factors comb-length q)
               (let ((acc (if maybe-factor (cons maybe-factor acc) acc)))
                 (if (p:polynomial= f p:+one+) acc
                     (recombine f factors acc comb-length q))))))
    (let ((p (suitable-prime polynomial))
          (b (suitable-bound polynomial)))
      (multiple-value-bind (q n)
          (lifting-steps b p)
        (let ((monic-p (fpx:monic-polynomial (fpx:modulo polynomial p) p))
              (monic-q (fpx:monic-polynomial polynomial q)))
          (recombine polynomial (lift-factors monic-q (fpx:berlekamp-factor monic-p p) p n)
                     nil 1 q))))))

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

;; https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Subresultant_pseudo-remainder_sequence

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

(sera:-> cyclotomic ((integer 1))
         (values p:polynomial &optional))
(defun cyclotomic (n)
  "Get n-th cyclotomic polynomial in \\(\\mathbb{Z}[x]\\)."
  (labels ((%cyclotomic (acc-m acc-d m)
             (if (> m n) (values acc-m acc-d)
                 (multiple-value-bind (q r)
                     (floor n m)
                   (if (zerop r)
                       (let ((p (p:polynomial (list (cons m 1) '(0 . -1)))))
                         (case (z:moebius q)
                           ( 0 (%cyclotomic acc-m acc-d (1+ m)))
                           ( 1 (%cyclotomic (p:multiply acc-m p) acc-d (1+ m)))
                           (-1 (%cyclotomic acc-m (p:multiply acc-d p) (1+ m)))))
                       (%cyclotomic acc-m acc-d (1+ m)))))))
    (multiple-value-bind (m d)
        (%cyclotomic p:+one+ p:+one+ 1)
      (nth-value 0 (divide m d)))))
