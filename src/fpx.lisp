;; Operations and factorization in 𝔽_p[x]

(defpackage cl-polynomial/fpx
  (:use #:cl)
  (:shadow #:gcd)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:p    #:cl-polynomial/polynomial)
                    (#:la   #:cl-polynomial/linalg))
  (:export #:modulo ; Operations mod p
           #:divide
           #:remainder
           #:monic-polynomial
           #:gcd
           #:gcdex

           #:square-free ; Factorization
           #:berlekamp-factor
           #:reducing-polynomials
           #:irreduciblep
           #:factor))
(in-package :cl-polynomial/fpx)

;; OPERATIONS

(sera:-> modulo (p:polynomial (integer 1))
         (values p:polynomial &optional))
(defun modulo (polynomial n)
  "Return a polynomial with every coefficient of @c(polynomial) being
taken modulo @c(n)."
  (p:polynomial
   (reduce
    (lambda (monomial acc)
      (declare (type u:monomial monomial))
      (destructuring-bind (d . c) monomial
        (let ((c (u:mod-sym c n)))
          (if (zerop c) acc (cons (cons d c) acc)))))
    (p:polynomial-coeffs polynomial)
    :from-end t
    :initial-value nil)))

(sera:-> divide (p:polynomial p:polynomial u:prime)
         (values p:polynomial p:polynomial &optional))
(defun divide (poly1 poly2 p)
  "Calculate \\(p_1 / p_2\\) where \\(p_1, p_2 \\in
\\mathbb{F}_p[x]\\), \\(p\\) being prime. A quotient and a remainder
are returned as 2 values."
  (let ((degree (p:degree poly2))
        (i (u:invert-integer (p:leading-coeff poly2) p)))
    (if (zerop degree)
        ;; Division by a constant is a special case
        (values
         (modulo (p:scale poly1 i) p)
         p:+zero+)
        (labels ((division-step (quotient-coeffs remainder)
                   (declare (type list quotient-coeffs)
                            (type p:polynomial remainder))
                   (if (< (p:degree remainder) degree)
                       (values (p:polynomial (reverse quotient-coeffs))
                               remainder)
                       (let ((remainder-coeffs (p:polynomial-coeffs remainder)))
                         (destructuring-bind (d . c) (car remainder-coeffs)
                           (let* ((quotient-degree (- d degree))
                                  (quotient-coeff (u:mod-sym (* c i) p))
                                  (monomial (cons quotient-degree quotient-coeff)))
                             (division-step
                              (cons monomial quotient-coeffs)
                              (modulo
                               (p:subtract remainder
                                         (p::multiply-monomial monomial poly2))
                               p))))))))
          (division-step nil poly1)))))

(sera:-> remainder (p:polynomial p:polynomial u:prime)
         (values p:polynomial &optional))
(declaim (inline remainder))
(defun remainder (poly1 poly2 p)
  "Calculate a remainder of \\(p_1 / p_2\\) where \\(p_1, p_2 \\in
\\mathbb{F}_p[x]\\), \\(p\\) being prime.

This function returns the second value of @c(divide)."
  (nth-value 1 (divide poly1 poly2 p)))

(sera:-> monic-polynomial (p:polynomial u:prime)
         (values p:polynomial integer &optional))
(defun monic-polynomial (polynomial p)
  "Factor an arbitrary non-zero polynomial in \\(\\mathbb{F}_p[x]\\),
\\(p\\) being prime, into a monic polynomial and a constant factor."
  (let ((c (p:leading-coeff polynomial)))
    (values
     (if (zerop c)
         ;; This polynomial equals to 0
         p:+zero+
         (modulo
          (p:scale polynomial (u:invert-integer c p)) p))
     c)))

(sera:-> gcd (p:polynomial p:polynomial u:prime)
         (values p:polynomial &optional))
(defun gcd (poly1 poly2 p)
  "Calculate the greatest common divisor of two polynomials in
\\(\\mathbb{F}_p[x]\\), p being prime."
  ;; Two first cases are special cases
  ;; gcd(a, 0) = positive-lc(a) and
  ;; gcd(0, a) = positive-lc(a) and
  ;; gcd(0, 0) = 0
  (p:positive-lc
   (cond
     ((p:polynomial= poly1 p:+zero+)
      poly2)
     ((p:polynomial= poly2 p:+zero+)
      poly1)
     (t
      ;; The rest is the Euclidean algorithm
      (let ((degree1 (p:degree poly1))
            (degree2 (p:degree poly2)))
        (labels ((%gcd (p1 p2)
                   (let ((r (remainder p1 p2 p)))
                     (if (p:polynomial= r p:+zero+)
                         (nth-value 0 (monic-polynomial p2 p))
                         (%gcd p2 r)))))
          (if (> degree1 degree2)
              (%gcd poly1 poly2)
              (%gcd poly2 poly1))))))))

(sera:-> gcdex (p:polynomial p:polynomial u:prime)
         (values p:polynomial p:polynomial p:polynomial &optional))
(defun gcdex (poly1 poly2 p)
  "Find \\(\\gcd(p_1, p_2)\\), \\(p_1, p_2 \\in \\mathbb{F}_p[x]\\)
and also find a solution of Bezout's equation \\(a p_1 + b p_2 =
\\gcd(p_1, p_2)\\) with minimal possible degree of \\(a\\) and
\\(b\\)."
  ;; POLY1 == 0 and/or POLY2 == 0 are special cases covered by COND
  (cond
    ((p:polynomial= poly1 p:+zero+)
     (values (p:positive-lc poly2) p:+zero+
             (p:scale p:+one+ (signum (p:leading-coeff poly2)))))
    ((p:polynomial= poly2 p:+zero+)
     (values (p:positive-lc poly1)
             (p:scale p:+one+ (signum (p:leading-coeff poly1))) p:+zero+))
    (t
     ;; This is an extended Euclidean algorithm
     (labels ((%gcd (p1 p2 s0 s1 d0 d1)
                (multiple-value-bind (q r)
                    (divide p1 p2 p)
                  (let ((s (modulo (p:subtract s0 (p:multiply q s1)) p))
                        (d (modulo (p:subtract d0 (p:multiply q d1)) p)))
                    (if (p:polynomial= r p:+zero+)
                        ;; Convert GCD to monic polynomial
                        (multiple-value-bind (m s)
                            (monic-polynomial p2 p)
                          (let ((s^-1 (u:invert-integer s p)))
                            (values
                             (p:positive-lc m)
                             (modulo (p:scale s1 s^-1) p)
                             (modulo (p:scale d1 s^-1) p))))
                        (%gcd p2 r s1 s d1 d))))))
       (if (> (p:degree poly1)
              (p:degree poly2))
           (%gcd poly1 poly2 p:+one+ p:+zero+ p:+zero+ p:+one+)
           (multiple-value-bind (gcd s d)
               (%gcd poly2 poly1 p:+one+ p:+zero+ p:+zero+ p:+one+)
             (values gcd d s)))))))

;; FACTORIZATION

(sera:-> x^pk-case (p:polynomial u:prime)
         (values p:polynomial &optional))
(defun x^pk-case (poly p)
  "Replace a polynomial in the form \\(\\sum_k a_k x^{b_k p}\\) with
\\(\\sum_k a_k x^{b_k}\\)."
  (p:polynomial
   (mapcar
    (lambda (m)
      (declare (type u:monomial m))
      (destructuring-bind (d . c) m
        (assert (zerop (rem d p)))
        (cons (/ d p) c)))
    (p:polynomial-coeffs poly))))

(sera:-> square-free (p:polynomial u:prime)
         (values list integer &optional))
(defun square-free (polynomial p)
  "Perform square-free factorization of a polynomial in
\\(\\mathbb{F}_p[x]\\), \\(p\\) being prime, with \\(\\deg f > 0\\). A list of
tuples \\((d_i . f_i)\\) is returned, so the supplied polynomial is equal to
\\(\\prod_i f_i^{d_i}\\) multiplied by the second returned value."
  (labels ((%%collect (p1 p2 acc n multiplicity)
             (if (p:polynomial= p2 p:+one+)
                 (values acc p1)
                 (let* ((gcd (gcd p1 p2 p))
                        (%p1 (divide p1 gcd p))
                        (%p2 (divide p2 gcd p)))
                   (%%collect %p1 gcd
                              (if (p:polynomial= %p2 p:+one+) acc
                                  (cons (cons (* n multiplicity) %p2) acc))
                              (1+ n) multiplicity))))
           (%collect (poly acc multiplicity)
             (let* ((gcd (gcd poly (modulo (p:derivative poly) p) p))
                    (w (divide poly gcd p)))
               (multiple-value-bind (%acc rest)
                   (%%collect gcd w acc 1 multiplicity)
                 (if (p:polynomial= rest p:+one+)
                     %acc
                     (%collect (x^pk-case rest p)
                               %acc (* multiplicity p)))))))
    (multiple-value-bind (monic m)
        (monic-polynomial polynomial p)
      (values (reverse (%collect monic nil 1)) m))))

;; I really don't know how to name it
(sera:-> berlekamp-matrix (p:polynomial u:prime)
         (values u:matrix &optional))
(defun berlekamp-matrix (poly p)
  (let* ((degree (p:degree poly))
         (matrix (make-array (list degree degree)
                             :element-type '(signed-byte 32))))
    (dotimes (i degree)
      (let ((x^n-mod-poly (remainder (p:polynomial (list (cons (* i p) 1))) poly p)))
        (dotimes (j degree)
          (setf (aref matrix i j)
                (mod
                 (- (or (cdr (assoc j (p:polynomial-coeffs x^n-mod-poly))) 0)
                    (if (= i j) 1 0))
                 p)))))
    ;; The Berlekamp matrix + I has coefficients of polynomials x^{i
    ;; prime} mod polynomial for 0≤i<deg polynomial in its columns.
    matrix))

(sera:-> reducing-polynomials (p:polynomial u:prime)
         (values list &optional))
(defun reducing-polynomials (f p)
  "Return a list of all f-reducing polynomials for a monic non-constant
square-free polynomial \\(f \\in \\mathbb{F}_p[x]\\)."
  (let ((nullspace (la:nullspace (berlekamp-matrix f p) p)))
    (mapcar #'p:sequence->polynomial nullspace)))

;; TODO: Generalize to all polynomials
(sera:-> irreduciblep (p:polynomial u:prime)
         (values boolean &optional))
(defun irreduciblep (f p)
  "Test if a monic non-constant square-free polynomial \\(f \\in
\\mathbb{F}_p[x]\\) is irreducible."
  (= (length (la:nullspace (berlekamp-matrix f p) p)) 1))

;; TODO: This implementation does a lot of redundant work. Optimize it.
(sera:-> berlekamp-factor (p:polynomial u:prime)
         (values list &optional))
(defun berlekamp-factor (f p)
  "Given a monic square-free non-constant polynomial, return a list of its
factors."
  (let ((rps (remove p:+one+ (reducing-polynomials f p)
                     :test #'p:polynomial=)))
    (if (null rps)
        (list f)
        (labels ((collect-factors (g)
                   ;; Return non-trivial factors of g. They are not necessarily
                   ;; irreducible and not all irreducible factors are contained
                   ;; in this factorization.
                   (reduce
                    (lambda (acc rp)
                      (labels ((add-factors (n acc)
                                 (if (= n p) acc
                                     (add-factors
                                      (1+ n)
                                      (let ((gcd
                                             (gcd g (modulo (p:add rp (p:scale p:+one+ n)) p)
                                                  p)))
                                        (if (p:polynomial= gcd p:+one+) acc
                                            (cons gcd acc)))))))
                        (add-factors 0 acc)))
                    rps
                    :initial-value nil))
                 (more-factors (factors)
                   ;; (MORE-FACTORS (COLLECT-FACTORS G)) return *all*
                   ;; non-trivial factors.
                   (let* ((old-factors (remove-duplicates factors :test #'p:polynomial=))
                          (new-factors (remove-duplicates
                                        (alex:flatten
                                         ;; Here we calculate some redundant GCDs because
                                         ;; gcd(f, gcd(f, g)) = gcd(f, g)
                                         (mapcar #'collect-factors old-factors))
                                        :test #'p:polynomial=)))
                     ;; If there are nothing new, stop
                     (if (null (set-difference new-factors old-factors :test #'p:polynomial=))
                         new-factors
                         (more-factors new-factors)))))
          ;; Now we need to filter irreducible factors from this mess.
          (let ((factors (more-factors (collect-factors f))))
            (remove-if
             (lambda (g)
               (some
                (lambda (h)
                  (let ((gcd (gcd h g p)))
                    (not
                     (or (p:polynomial= gcd p:+one+)
                         (p:polynomial= gcd g)))))
                factors))
             factors))))))

;; TODO: Generalize to constant polynomials?
(sera:-> factor (p:polynomial u:prime)
         (values list integer &optional))
(defun factor (polynomial p)
  "Factor a non-constant polynomial in \\(\\mathbb{F}_p[x]\\) into irreducible factors."
  (multiple-value-bind (sqfr-factors m)
      (square-free polynomial p)
    (values
     (reduce #'append
             (mapcar
              (lambda (factor)
                (destructuring-bind (m . f) factor
                  (mapcar
                   (lambda (irreducible-factor)
                     (cons m irreducible-factor))
                   (berlekamp-factor f p))))
              sqfr-factors))
     m)))
