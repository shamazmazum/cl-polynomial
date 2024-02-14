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
           #:distinct-degree
           #:berlekamp-factor
           #:reducing-polynomials
           #:irreduciblep
           #:factor))
(in-package :cl-polynomial/fpx)

;; OPERATIONS

(sera:-> modulo (p:polynomial u:prime-power)
         (values p:polynomial &optional))
(defun modulo (polynomial q)
  "Return a polynomial with every coefficient of @c(polynomial) being
taken modulo @c(q)."
  (p:polynomial
   (reduce
    (lambda (monomial acc)
      (p:bind-monomial (d c) monomial
        (let ((c (u:mod-sym c q)))
          (if (zerop c) acc (cons (cons d c) acc)))))
    (p:polynomial-coeffs polynomial)
    :from-end t
    :initial-value nil)))

(sera:-> divide (p:polynomial p:polynomial u:prime-power)
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
                         (p:bind-monomial (d c) (car remainder-coeffs)
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

(sera:-> remainder (p:polynomial p:polynomial u:prime-power)
         (values p:polynomial &optional))
(declaim (inline remainder))
(defun remainder (poly1 poly2 p)
  "Calculate a remainder of \\(p_1 / p_2\\) where \\(p_1, p_2 \\in
\\mathbb{F}_p[x]\\), \\(p\\) being prime.

This function returns the second value of @c(divide)."
  (nth-value 1 (divide poly1 poly2 p)))

(sera:-> monic-polynomial (p:polynomial u:prime-power)
         (values p:polynomial integer &optional))
(defun monic-polynomial (polynomial p)
  "Factor an arbitrary non-zero polynomial in \\(\\mathbb{F}_p[x]\\),
\\(p\\) being prime, into a monic polynomial and a constant factor."
  (let ((c (p:leading-coeff polynomial)))
    (values
     (if (zerop c) p:+zero+
         (modulo
          (p:scale polynomial (u:invert-integer c p)) p))
     c)))

(sera:-> gcd (p:polynomial p:polynomial u:prime-power)
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
     ((p:polynomial= poly1 p:+zero+) poly2)
     ((p:polynomial= poly2 p:+zero+) poly1)
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

(sera:-> gcdex (p:polynomial p:polynomial u:prime-power)
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
      (declare (type p:monomial m))
      (p:bind-monomial (d c) m
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
                 (if (p:polynomial= rest p:+one+) %acc
                     (%collect (x^pk-case rest p)
                               %acc (* multiplicity p)))))))
    (cond
      ((p:polynomial= polynomial p:+zero+)
       (error "Cannot factor zero polynomial"))
      ((= (p:degree polynomial) 0)
       (values (list (cons 1 p:+one+))
               (p:leading-coeff polynomial)))
      (t
       (multiple-value-bind (monic m)
           (monic-polynomial polynomial p)
         (values (%collect monic nil 1) m))))))

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

(sera:-> irreduciblep (p:polynomial u:prime)
         (values boolean &optional))
(defun irreduciblep (f p)
  "Test if a polynomial \\(f \\in \\mathbb{F}_p[x]\\) is irreducible."
  (multiple-value-bind (factors c)
      (square-free f p)
    (declare (ignore c))
    (if (= (length factors) 1)
        (destructuring-bind (m . f) (car factors)
          (and (= m 1)
               (= (length (la:nullspace (berlekamp-matrix f p) p)) 1))))))

(sera:-> berlekamp-factor (p:polynomial u:prime)
         (values list &optional))
(defun berlekamp-factor (f p)
  "Given a monic square-free non-constant polynomial, return a list of its
factors."
  (let* ((rps (remove p:+one+ (reducing-polynomials f p)
                      :test #'p:polynomial=))
         (nfactors (1+ (length rps))))
    (labels ((constant-poly (c)
               (p:polynomial (list (cons 0 c))))
             ;; Given a polynomial FACTOR and F-reducing polynomial RP, try to find
             ;; non-trivial factors of FACTOR. If there are such factors, remove
             ;; FACTOR from ACC, add those factors to ACC and try again.
             (%collect-factors (factor rp element acc)
               (if (or (= element p)
                       (= (length acc) nfactors))
                   acc
                   (let ((gcd (gcd factor (p:add rp (constant-poly element)) p))
                         (next-element (1+ element)))
                     (if (or (p:polynomial= gcd p:+one+)
                             (p:polynomial= gcd factor))
                         ;; GCD is a trivial factor of FACTOR, try the next element
                         ;; of a finite field.
                         (%collect-factors factor rp next-element acc)
                         ;; GCD is a non-trivial factor of FACTOR, continue
                         ;; recursively with factors GCD and FACTOR/GCD.
                         (let ((new-factor (divide factor gcd p)))
                           (%collect-factors
                            gcd rp next-element
                            (%collect-factors
                             new-factor rp next-element
                             (append (list new-factor gcd) (remove factor acc)))))))))
             ;; Try to split already found factors of F in ACC using F-reducing
             ;; polynomial RP.
             (collect-factors (rp element acc)
               (reduce
                (lambda (acc factor)
                  (if (= (length acc) nfactors) acc
                      (%collect-factors factor rp element acc)))
                acc :initial-value acc)))
      ;; Starting with (F) as a list of known factors of F, try to use F-reducing
      ;; polynomials to find NFACTORS distinct non-trivial factors of F.
      (reduce
       (lambda (acc rp)
         (if (= (length acc) nfactors) acc
             (collect-factors rp 0 acc)))
       rps :initial-value (list f)))))

(sera:-> expt-rem (p:polynomial integer p:polynomial u:prime)
         (values p:polynomial &optional))
(defun expt-rem (f n g p)
  "Calculate \\(f^n(x) \\mod g(x)\\) for a non-negative integer
\\(n\\) and \\(f,g \\in \\mathbb{F}_p[x]\\)."
  (labels ((mul-rem (f1 f2)
             (remainder (modulo (p:multiply f1 f2) p) g p))
           (%expt-rem (f n acc)
             (cond
               ((zerop n) acc)
               ((evenp n)
                (%expt-rem (mul-rem f f) (floor n 2) acc))
               (t
                (%expt-rem f (1- n) (mul-rem f acc))))))
    (%expt-rem f n p:+one+)))

;; f(x)^q mod q can be done really fast, no need to call slow EXPT-REM.
(sera:-> expt-q (p:polynomial u:prime)
         (values p:polynomial &optional))
(defun expt-q (f q)
  "Calculate \\(f(x)^q mod q\\)."
  (p:polynomial
   (mapcar
    (lambda (monomial)
      (p:bind-monomial (d c)
          monomial
        (cons (* d q) c)))
    (p:polynomial-coeffs f))))

(sera:-> cantor-zassenhaus (p:polynomial u:prime u:degree)
         (values list &optional))
(defun cantor-zassenhaus (f p deg)
  "Given a monic square-free non-constant polynomial, return a list of
its factors in \\(\\mathbb{F}_p[x]\\) for a prime \\(p > 2\\). Each
factor must be of degree @c(deg)."
  (assert (> p 2))
  (let ((nfactors (/ (p:degree f) deg)))
    (labels ((random-poly ()
               (p:add
                (p:list->polynomial
                 (loop repeat (1- (* 2 deg))
                       collect (- (random p) (floor p 2))))
                (p:polynomial (list (cons (1- (* 2 deg)) 1)))))
             (w (f)
               (modulo
                (p:add (expt-rem (random-poly) (floor (expt p deg) 2) f p) p:+one+) p))
             (collect-factors (factor acc)
               (if (= (length acc) nfactors) acc
                   (let ((gcd (gcd factor (w f) p)))
                     (if (or (p:polynomial= gcd p:+one+)
                             (p:polynomial= gcd factor))
                         ;; GCD is a trivial factor of FACTOR, try another
                         ;; random polynomial.
                         (collect-factors factor acc)
                         ;; GCD is a non-trivial factor of FACTOR, continue
                         ;; recursively with factors GCD and FACTOR/GCD.
                         (let* ((new-factor (divide factor gcd p))
                                (%acc (append (list new-factor gcd)
                                              (remove factor acc)))
                                ;; Factor GCD if needed
                                (%acc (if (= (p:degree gcd) deg) %acc
                                          (collect-factors gcd %acc)))
                                ;; Factor NEW-FACTOR if needed
                                (%acc (if (= (p:degree new-factor) deg) %acc
                                          (collect-factors new-factor %acc))))
                           %acc))))))
      (collect-factors f (list f)))))

(sera:-> distinct-degree (p:polynomial u:prime)
         (values list &optional))
(defun distinct-degree (f p)
  "Perform distinct degree factorization of a non-constant square-free
polynomial \\(f \\in \\mathbb{F}_p[x]\\). Return a list of pairs
\\((d_i . p_i)\\) where \\(p_i\\) is a product of all factors of
\\(f\\) of degree \\(d_i\\)."
  (let ((x (p:polynomial '((1 . 1)))))
    (labels ((collect (f w deg acc)
               (cond
                 ((p:polynomial= f p:+one+) acc)
               ((> deg (floor (p:degree f) 2))
                (cons (cons (p:degree f) f) acc))
               (t
                (let* ((w (remainder (expt-q w p) f p))
                       (gcd (gcd f (p:subtract w x) p)))
                  (if (p:polynomial= gcd p:+one+)
                      (collect f w (1+ deg) acc)
                      (collect (divide f gcd p)
                        w (1+ deg)
                        (cons (cons deg gcd) acc))))))))
      (collect f x 1 nil))))

(sera:-> big-field-factor (p:polynomial u:prime)
         (values list &optional))
(defun big-field-factor (f p)
  "Perform factorization of a non-constant square-free polynomial in a
big finite field by firstly apply DISTINCT-DEGREE and then
CANTOR-ZASSENHAUS to F."
  (let ((factors (distinct-degree f p)))
    (reduce #'append factors
            :from-end t ; Less consing
            :key (lambda (factor)
                   (cantor-zassenhaus (cdr factor) p (car factor))))))

(deftype factorization-algorithm ()
  '(member :berlekamp :cantor-zassenhaus))

;; TODO: Generalize to constant polynomials?
(sera:-> factor (p:polynomial u:prime &optional factorization-algorithm)
         (values list integer &optional))
(defun factor (polynomial p &optional (algorithm :berlekamp))
  "Factor a polynomial in \\(\\mathbb{F}_p[x]\\) into irreducible
factors. An optional parameter @c(algorithm) can be either
@c(:berlekamp) (default) or @c(:cantor-zassenhaus)."
  (multiple-value-bind (sqfr-factors m)
      (square-free polynomial p)
    (values
     (reduce #'append sqfr-factors
             :key
             (lambda (factor)
               (destructuring-bind (m . f) factor
                 (mapcar
                  (lambda (irreducible-factor)
                    (cons m irreducible-factor))
                  (cond
                    ((or (eq algorithm :berlekamp) (= p 2))
                     (berlekamp-factor f p))
                    ((eq algorithm :cantor-zassenhaus)
                     (big-field-factor f p))
                    (t
                     (error "ALGORITHM must be either :BERLEKAMP or :CANTOR-ZASSENHAUS")))))))
     m)))
