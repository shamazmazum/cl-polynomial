(in-package :polynomial)

(defstruct polynomial
  (coeffs nil :type list :read-only t))

(sera:-> polynomial= (polynomial polynomial)
         (values boolean &optional))
(defun polynomial= (poly1 poly2)
  "Test if two polynomials are equal."
  (equalp (polynomial-coeffs poly1)
          (polynomial-coeffs poly2)))

;; Convenient constructor
(sera:-> polynomial (list)
         (values polynomial &optional))
(declaim (inline polynomial))
(defun polynomial (list)
  (make-polynomial :coeffs list))

(sera:-> degree (polynomial)
         (values alex:non-negative-fixnum &optional))
(defun degree (polynomial)
  "Return degree of a polynomial"
  (let ((first (first (polynomial-coeffs polynomial))))
    (if first (car first) 0)))

(sera:-> leading-coeff (polynomial)
         (values fixnum &optional))
(defun leading-coeff (polynomial)
  "Return the leading coefficient of a polynomial"
  (let ((first (first (polynomial-coeffs polynomial))))
    (if first (cdr first) 0)))

(sera:-> monicp (polynomial)
         (values boolean &optional))
(defun monicp (polynomial)
  "Test if a polynomial is monic"
  (= (leading-coeff polynomial) 1))

(sera:-> constantp (polynomial)
         (values boolean &optional))
(defun constantp (polynomial)
  "Test if a polynomial has degree 0."
  (zerop (degree polynomial)))

(sera:-> list->polynomial (list)
         (values polynomial &optional))
(defun list->polynomial (coeffs)
  "Given a list of integer coefficients convert it to a polynomial. An
index into the list corresponds with a power for that coefficient, e.g.:

@begin[lang=lisp](code)
CL-USER> (polynomial:list->polynomial '(1 2))
2X^1 + X^0
@end(code)"
  (polynomial
   (si:foldl
    (lambda (acc monomial)
      (declare (type monomial monomial))
      (destructuring-bind (degree . coeff) monomial
        (declare (ignore degree))
        (if (zerop coeff) acc (cons monomial acc))))
    nil
    (si:enumerate (si:list->iterator coeffs)))))

(sera:-> sequence->polynomial (sequence)
         (values polynomial &optional))
(defun sequence->polynomial (seq)
  "Like list->polynomial, but for all types of sequences.

@begin[lang=lisp](code)
CL-USER> (polynomial:sequence->polynomial #*101)
X^2 + X^0
@end(code)"
  (list->polynomial (coerce seq 'list)))

(sera:-> polynomial->list (polynomial)
         (values list &optional))
(defun polynomial->list (polynomial)
  "Convert a polynomial to a list of coefficients. This is an inverse
of @c(list->polynomial)."
  (let ((degree (degree polynomial))
        (coeffs (polynomial-coeffs polynomial)))
    (loop for d from 0 to degree collect
          (or (cdr (assoc d coeffs)) 0))))

(defparameter +zero+ (list->polynomial '(0))
  "Additive identity")

(defparameter +one+ (list->polynomial '(1))
  "Multiplicative identity")

(sera:-> print-monomial (monomial boolean stream)
         (values &optional))
(defun print-monomial (monomial lastp stream)
  (declare (ignore lastp))
  (destructuring-bind (degree . coeff) monomial
    (unless (zerop coeff)
      (format stream "~:[~d~;~]" (= coeff 1) coeff)
      (format stream "X^~d" degree)))
  (values))

(defmethod print-object ((polynomial polynomial) stream)
  (let* ((coeffs (polynomial-coeffs polynomial))
         (last (car (last coeffs))))
    (if (polynomial= +zero+ polynomial)
        (princ "0X^0" stream)
        (dolist (monomial coeffs)
          (let ((lastp (eq monomial last)))
            (print-monomial monomial lastp stream)
            (unless lastp
              (princ " + " stream)))))))

(sera:-> %add (polynomial polynomial)
         (values polynomial &optional))
(defun %add (poly1 poly2)
  (declare (optimize (speed 3)))
  ;; Polynomials coefficients are stored in a assoc list sorted by
  ;; power of X, so we need to compare heads of two lists when adding
  ;; two polynomials.
  (labels ((collect-coeffs (acc ms1 ms2)
             (cond
               ((null ms1)
                (append (reverse ms2) acc))
               ((null ms2)
                (append (reverse ms1) acc))
               (t
                (destructuring-bind ((d1 . c1) &rest rest1) ms1
                  (destructuring-bind ((d2 . c2) &rest rest2) ms2
                    (declare (type alex:non-negative-fixnum d1 d2)
                             (type fixnum c1 c2))
                    (cond
                      ((> d1 d2)
                       (collect-coeffs (cons (cons d1 c1) acc)
                                       rest1 ms2))
                      ((> d2 d1)
                       (collect-coeffs (cons (cons d2 c2) acc)
                                       rest2 ms1))
                      (t
                       (let ((sum (+ c1 c2)))
                         (collect-coeffs (if (zerop sum)
                                             acc (cons (cons d1 (+ c1 c2)) acc))
                                         rest1 rest2))))))))))
    (polynomial
     (reverse
      (collect-coeffs
       nil
       (polynomial-coeffs poly1)
       (polynomial-coeffs poly2))))))
                    
(sera:-> add (&rest polynomial)
         (values polynomial &optional))
(declaim (inline add))
(defun add (&rest polynomials)
  "Add polynomials together"
  (reduce #'%add polynomials :initial-value +zero+))

(sera:-> map-poly ((sera:-> (fixnum) (values fixnum &optional)) polynomial)
         (values polynomial &optional))
(defun map-poly (fn polynomial)
  "Given a polynomial \\(\\sum_n a_n x^n\\) return a polynomial
\\(\\sum_n fn(a_n)x^n\\).

This is @c(fmap) for a type @c(data Poly a b = Poly [(a, b)])."
  (polynomial
   (mapcar
    (lambda (m)
      (declare (type monomial m))
      (destructuring-bind (d . c) m
        (cons d (funcall fn c))))
    (polynomial-coeffs polynomial))))

(sera:-> negate (polynomial)
         (values polynomial &optional))
(declaim (inline negate))
(defun negate (polynomial)
  "Given a polynomial \\(\\sum_n a_n x^n\\) return a polynomial
\\(\\sum_n -a_n x^n\\).

This is function is equivalent to
@begin[lang=lisp](code)
(alexandria:curry #'map-polynomial #'-)
@end(code)"
  (map-poly #'- polynomial))

(sera:-> subtract (polynomial &rest polynomial)
         (values polynomial &optional))
(declaim (inline subtract))
(defun subtract (polynomial &rest polynomials)
  "Subtract @c(polynomials) from @c(polynomial)."
  (reduce #'%add (mapcar #'negate polynomials)
          :initial-value polynomial))

(sera:-> multiply-monomial (monomial polynomial)
         (values polynomial &optional))
(defun multiply-monomial (monomial polynomial)
  (destructuring-bind (d1 . c1) monomial
    (polynomial
     (mapcar
      (lambda (m)
        (declare (type monomial m))
        (destructuring-bind (d2 . c2) m
          (cons (+ d1 d2) (* c1 c2))))
      (polynomial-coeffs polynomial)))))

(sera:-> %multiply (polynomial polynomial)
         (values polynomial &optional))
(defun %multiply (poly1 poly2)
  (let ((polys (mapcar
                (lambda (monomial)
                  (multiply-monomial monomial poly2))
                (polynomial-coeffs poly1))))
    (reduce #'%add polys :initial-value +zero+)))

(sera:-> multiply (&rest polynomial)
         (values polynomial &optional))
(declaim (inline multiply))
(defun multiply (&rest polynomials)
  "Multiply polynomials"
  (reduce #'%multiply polynomials :initial-value +one+))

(sera:-> scale (polynomial fixnum)
         (values polynomial &optional))
(declaim (inline scale))
(defun scale (polynomial c)
  "Multiply a polynomial by a constant."
  (map-poly (lambda (a) (* a c)) polynomial))

(sera:-> modulo (polynomial alex:positive-fixnum)
         (values polynomial &optional))
(defun modulo (polynomial n)
  "Return a polynomial with every coefficient of @c(polynomial) being
taken modulo @c(n)."
  (polynomial
   (reduce
    (lambda (monomial acc)
      (declare (type monomial monomial))
      (destructuring-bind (d . c) monomial
        (let ((c (mod c n)))
          (if (zerop c) acc (cons (cons d c) acc)))))
    (polynomial-coeffs polynomial)
    :from-end t
    :initial-value nil)))

(sera:-> invert-integer (alex:positive-fixnum prime)
         (values alex:positive-fixnum &optional))
(defun invert-integer (n p)
  "Find a multiplicative inverse of \\(n\\) in \\(\\mathbb{F}_p\\), p
being prime, i.e. find \\(x\\) such that \\(xn = nx = 1\\)."
  ;; Remember that n^p = n
  (declare (optimize (speed 3)))
  (mod (expt n (- p 2)) p))

(sera:-> divide (polynomial polynomial prime)
         (values polynomial polynomial &optional))
(defun divide (poly1 poly2 p)
  "Calculate \\(p_1 / p_2\\) where \\(p_1, p_2 \\in
\\mathbb{F}_p[x]\\), \\(p\\) being prime. A quotient and a remainder
are returned as 2 values."
  (let ((degree (degree poly2))
        (i (invert-integer (leading-coeff poly2) p)))
    (if (zerop degree)
        ;; Division by a constant is a special case
        (values
         (modulo (scale poly1 i) p)
         +zero+)
        (labels ((division-step (quotient-coeffs remainder)
                   (declare (type list quotient-coeffs)
                            (type polynomial remainder))
                   (if (< (degree remainder) degree)
                       (values (polynomial (reverse quotient-coeffs))
                               remainder)
                       (let ((remainder-coeffs (polynomial-coeffs remainder)))
                         (destructuring-bind (d . c) (car remainder-coeffs)
                           (let* ((quotient-degree (- d degree))
                                  (quotient-coeff (mod (* c i) p))
                                  (monomial (cons quotient-degree quotient-coeff)))
                             (division-step
                              (cons monomial quotient-coeffs)
                              (modulo
                               (subtract remainder
                                         (multiply-monomial monomial poly2))
                               p))))))))
          (division-step nil poly1)))))

(sera:-> remainder (polynomial polynomial prime)
         (values polynomial &optional))
(declaim (inline remainder))
(defun remainder (poly1 poly2 p)
  "Calculate a remainder of \\(p_1 / p_2\\) where \\(p_1, p_2 \\in
\\mathbb{F}_p[x]\\), \\(p\\) being prime.

This function returns the second value of @c(divide)."
  (nth-value 1 (divide poly1 poly2 p)))

(sera:-> monic-polynomial (polynomial prime)
         (values polynomial fixnum &optional))
(defun monic-polynomial (polynomial p)
  "Factor an arbitrary non-zero polynomial in \\(\\mathbb{F}_p[x]\\),
\\(p\\) being prime, into a monic polynomial and a constant factor."
  (let ((c (leading-coeff polynomial)))
    (values
     (if (zerop c)
         ;; This polynomial equals to 0
         +zero+
         (modulo
          (scale polynomial (invert-integer c p)) p))
     c)))

(sera:-> gcd (polynomial polynomial prime)
         (values polynomial &optional))
(defun gcd (poly1 poly2 p)
  "Calculate thr greatest common divisor of two polynomials in
\\(\\mathbb{F}_p[x]\\), p being prime."
  (nth-value
   0 (monic-polynomial
      ;; Two first cases are special cases
      ;; gcd(a, 0) = monic(a) and
      ;; gcd(0, a) = monic(a) and
      ;; gcd(0, 0) = 0
      (cond
        ((polynomial= poly1 +zero+)
         poly2)
        ((polynomial= poly2 +zero+)
         poly1)
        (t
         ;; The rest is the Euclidean algorithm
         (let ((degree1 (degree poly1))
               (degree2 (degree poly2)))
           (labels ((%gcd (p1 p2)
                      (let ((r (remainder p1 p2 p)))
                        (if (polynomial= r +zero+) p2
                            (%gcd p2 r)))))
             (if (> degree1 degree2)
                 (%gcd poly1 poly2)
                 (%gcd poly2 poly1))))))
      p)))

(sera:-> derivative (polynomial)
         (values polynomial &optional))
(defun derivative (polynomial)
  "Return a formal derivative of a polynomial."
  (polynomial
   (reduce
    (lambda (m acc)
      (declare (type monomial m))
      (destructuring-bind (d . c) m
        (if (zerop d) acc
            (cons
             (cons (1- d) (* d c))
             acc))))
    (polynomial-coeffs polynomial)
    :from-end t
    :initial-value nil)))
