;; Polynomials and opeartions in R[x], R being a ring such as ℤ.

(defpackage cl-polynomial/polynomial
  (:use #:cl)
  (:shadow #:constantp #:expt)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:si   #:stateless-iterators))
  (:export #:monomial
           #:bind-monomial
           #:polynomial
           #:polynomial=
           #:polynomial/=
           #:+zero+
           #:+one+
           #:+variable+
           #:degree
           #:polynomial-coeffs
           #:leading-coeff
           #:list->polynomial
           #:monomials->polynomial
           #:sequence->polynomial
           #:polynomial->list
           #:positive-lc
           #:map-poly
           #:negate
           #:add
           #:subtract
           #:multiply
           #:scale
           #:monicp
           #:constantp
           #:derivative
           #:expt
           #:reciprocal
           #:palindromep
           #:evaluate))
(in-package :cl-polynomial/polynomial)

;; FIXME: let monomial be just an alias for a tuple?
(deftype monomial () '(cons u:degree integer))

;; Like destructuring-bind, but without additional checks
(defmacro bind-monomial ((deg coeff) monomial &body body)
  (let ((m (gensym)))
    `(let ((,m ,monomial))
       (declare (type monomial ,m))
       (let ((,deg   (car ,m))
             (,coeff (cdr ,m)))
         ,@body))))

(defstruct (polynomial
             (:constructor polynomial (coeffs)))
  (coeffs nil :type list :read-only t))

(sera:-> monomials->polynomial (list)
         (values polynomial &optional))
(defun monomials->polynomial (list)
  "Convert a list of monomials in the form @c((degree . coefficient))
to a polynomial. A list must not contain two monomials of the same
degree."
  (when (/= (length list)
            (length (remove-duplicates list :key #'car :test #'=)))
    (error "A list contains monomials of the same degree"))
  (polynomial
   (sort (copy-list list) #'> :key #'car)))

(sera:-> polynomial= (polynomial polynomial)
         (values boolean &optional))
(defun polynomial= (poly1 poly2)
  "Test if two polynomials are equal."
  (declare (optimize (speed 3)))
  (or (eq poly1 poly2)
      (labels ((%= (ms1 ms2)
                 (cond
                   ((and (null ms1) (null ms2)) t)
                   ((or  (null ms1) (null ms2)) nil)
                   (t
                    (bind-monomial (d1 c1) (car ms1)
                      (bind-monomial (d2 c2) (car ms2)
                        (if (and (= d1 d2) (= c1 c2))
                            (%= (cdr ms1) (cdr ms2)))))))))
        (%=
         (polynomial-coeffs poly1)
         (polynomial-coeffs poly2)))))

(sera:-> polynomial/= (polynomial polynomial)
         (values boolean &optional))
(defun polynomial/= (poly1 poly2)
  "Test if two polynomials are not equal."
  (declare (optimize (speed 3)))
  (and (not (eq poly1 poly2))
       (labels ((%/= (ms1 ms2)
                  (cond
                    ((and (null ms1) (null ms2)) nil)
                    ((or (null ms2) (null ms1)) t)
                    (t
                     (bind-monomial (d1 c1) (car ms1)
                       (bind-monomial (d2 c2) (car ms2)
                         (if (or (/= d1 d2) (/= c1 c2)) t
                             (%/= (cdr ms1) (cdr ms2)))))))))
         (%/=
          (polynomial-coeffs poly1)
          (polynomial-coeffs poly2)))))

(sera:-> degree (polynomial)
         (values alex:non-negative-integer &optional))
(defun degree (polynomial)
  "Return degree of a polynomial"
  (let ((first (first (polynomial-coeffs polynomial))))
    (if first (car first) 0)))

(sera:-> leading-coeff (polynomial)
         (values integer &optional))
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
#<POLYNOMIAL:POLYNOMIAL 2X + 1>
@end(code)"
  (polynomial
   (si:foldl
    (lambda (acc monomial)
      (declare (type monomial monomial))
      (bind-monomial (degree coeff) monomial
        (declare (ignore degree))
        (if (zerop coeff) acc (cons monomial acc))))
    nil
    (si:enumerate (si:list->iterator coeffs)))))

(sera:-> sequence->polynomial (sequence)
         (values polynomial &optional))
(defun sequence->polynomial (seq)
  "Like @c(list->polynomial), but for all types of sequences.

@begin[lang=lisp](code)
CL-USER> (polynomial:sequence->polynomial #*101)
#<POLYNOMIAL:POLYNOMIAL X^2 + 1>
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

(defparameter +variable+ (list->polynomial '(0 1))
  "The polynomial \\(x\\)")

(sera:-> print-monomial (monomial boolean stream)
         (values &optional))
(defun print-monomial (monomial firstp stream)
  (bind-monomial (degree coeff) monomial
    (unless (zerop coeff)
      (if firstp
          (when (< coeff 0) (princ "- " stream))
          (format stream " ~:[+~;-~] " (< coeff 0)))
      (let ((abs (abs coeff)))
        (format stream "~:[~d~;~]"
                (and (= abs 1) (not (zerop degree)))
                abs))
      (format stream "~[~;X~:;X^~d~]" degree degree)))
  (values))

(defmethod print-object ((polynomial polynomial) stream)
  (let* ((coeffs (polynomial-coeffs polynomial))
         (first (first coeffs)))
    (print-unreadable-object (polynomial stream :type t)
      (if (polynomial= +zero+ polynomial)
        (princ 0 stream)
        (dolist (monomial coeffs)
          (print-monomial monomial (eq monomial first) stream))))))

(sera:-> %add (polynomial polynomial)
         (values polynomial &optional))
(defun %add (poly1 poly2)
  (declare (optimize (speed 3)))
  ;; Polynomials coefficients are stored in a assoc list sorted by
  ;; power of X, so we need to compare heads of two lists when adding
  ;; two polynomials.
  (labels ((%do (ms1 ms2)
             (cond
               ((null ms1) ms2)
               ((null ms2) ms1)
               (t
                (bind-monomial (d1 c1) (car ms1)
                  (bind-monomial (d2 c2) (car ms2)
                    (cond
                      ((> d1 d2)
                       (cons (car ms1) (%do (cdr ms1) ms2)))
                      ((> d2 d1)
                       (cons (car ms2) (%do ms1 (cdr ms2))))
                      (t
                       (let ((sum (+ c1 c2)))
                         (if (zerop sum)
                             (%do (cdr ms1) (cdr ms2))
                             (cons (cons d1 sum) (%do (cdr ms1) (cdr ms2)))))))))))))
    (polynomial (%do (polynomial-coeffs poly1) (polynomial-coeffs poly2)))))
                    
(sera:-> add (&rest polynomial)
         (values polynomial &optional))
(declaim (inline add))
(defun add (&rest polynomials)
  "Add polynomials together"
  (reduce #'%add polynomials :initial-value +zero+))
(define-compiler-macro add (&whole whole &rest polynomials)
  (case (length polynomials)
    (0 +zero+)
    (1 (first polynomials))
    (2 `(%add ,(first polynomials)
              ,(second polynomials)))
    (t whole)))

(sera:-> map-poly ((sera:-> (integer) (values integer &optional)) polynomial)
         (values polynomial &optional))
(defun map-poly (fn polynomial)
  "Given a polynomial \\(\\sum_n a_n x^n\\) return a polynomial
\\(\\sum_n fn(a_n)x^n\\).

This is @c(fmap) for a type @c(data Poly a b = Poly [(a, b)])."
  (declare (optimize (speed 3)))
  (polynomial
   (mapcar
    (lambda (m)
      (bind-monomial (d c) m
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
(alexandria:curry #'map-poly #'-)
@end(code)"
  (map-poly #'- polynomial))

(sera:-> positive-lc (polynomial)
         (values polynomial &optional))
(declaim (inline positive-lc))
(defun positive-lc (polynomial)
  "Return @c(polynomial) or its negation, so that the leading
coefficient is positive."
  (let ((lc (leading-coeff polynomial)))
    (if (> lc 0) polynomial
        (negate polynomial))))

(sera:-> subtract (polynomial &rest polynomial)
         (values polynomial &optional))
(declaim (inline subtract))
(defun subtract (polynomial &rest polynomials)
  "Subtract @c(polynomials) from @c(polynomial)."
  (reduce #'%add polynomials
          :key #'negate
          :initial-value polynomial))
(define-compiler-macro subtract (&whole whole polynomial &rest polynomials)
  (case (length polynomials)
    (0 polynomial)
    (1 `(%add ,polynomial (negate ,(first polynomials))))
    (t whole)))

(sera:-> multiply-monomial (monomial polynomial)
         (values polynomial &optional))
(defun multiply-monomial (monomial polynomial)
  (bind-monomial (d1 c1) monomial
    (polynomial
     (mapcar
      (lambda (m)
        (bind-monomial (d2 c2) m
          (cons (+ d1 d2) (* c1 c2))))
      (polynomial-coeffs polynomial)))))

(sera:-> %multiply (polynomial polynomial)
         (values polynomial &optional))
(defun %multiply (poly1 poly2)
  (reduce #'%add (polynomial-coeffs poly1)
          :initial-value +zero+
          :key (lambda (monomial)
                 (multiply-monomial monomial poly2))))

(sera:-> multiply (&rest polynomial)
         (values polynomial &optional))
(declaim (inline multiply))
(defun multiply (&rest polynomials)
  "Multiply polynomials"
  (reduce #'%multiply polynomials :initial-value +one+))
(define-compiler-macro multiply (&whole whole &rest polynomials)
  (case (length polynomials)
    (0 +one+)
    (1 (first polynomials))
    (2 `(%multiply ,(first polynomials)
                   ,(second polynomials)))
    (t whole)))

(sera:-> scale (polynomial integer)
         (values polynomial &optional))
(declaim (inline scale))
(defun scale (polynomial c)
  "Multiply a polynomial by a constant."
  (map-poly (lambda (a) (* a c)) polynomial))

(sera:-> derivative (polynomial)
         (values polynomial &optional))
(defun derivative (polynomial)
  "Return a formal derivative of a polynomial."
  (polynomial
   (reduce
    (lambda (m acc)
      (bind-monomial (d  c) m
        (if (zerop d) acc
            (cons
             (cons (1- d) (* d c))
             acc))))
    (polynomial-coeffs polynomial)
    :from-end t
    :initial-value nil)))

(sera:-> expt (polynomial unsigned-byte)
         (values polynomial &optional))
(defun expt (f n)
  "For polynomial \\(f(x)\\) and a non-negative integer number \\(n\\)
calculate \\(f^n(x)\\)."
  (labels ((%expt (f n acc)
             (cond
               ((zerop n) acc)
               ((evenp n)
                (%expt (multiply f f) (floor n 2) acc))
               (t
                (%expt f (1- n) (multiply f acc))))))
    (%expt f n +one+)))

(sera:-> reciprocal (polynomial)
         (values polynomial &optional))
(defun reciprocal (f)
  "Return a reciprocal polynomial of \\(f\\), i.e. given \\(f(x)\\),
return \\(x^{\\deg f} f(x^{-1})\\)."
  (let ((degree (degree f)))
    (polynomial
     (reduce
      (lambda (acc m)
        (declare (type monomial m))
        (bind-monomial (d c) m
          (cons (cons (- degree d) c) acc)))
      (polynomial-coeffs f)
      :initial-value nil))))

(sera:-> palindromep (polynomial)
         (values boolean &optional))
(defun palindromep (f)
  "Return @c(T) if \\(f\\) is a palindrome, i.e. it is equal to its
reciprocal polynomial."
  (polynomial= f (reciprocal f)))


(sera:-> evaluate (polynomial integer)
         (values integer &optional))
(defun evaluate (p x)
  "Evaluate a polynomial in \\(\\mathbb{Z}[x]\\) at point \\(x\\)."
  (labels ((%go (acc coeffs)
             (if (null coeffs) acc
                 (bind-monomial (deg coeff)
                     (car coeffs)
                   (let ((next-deg (or (car (second coeffs)) 0)))
                     (%go (* (+ acc coeff)
                             (cl:expt x (- deg next-deg)))
                          (cdr coeffs)))))))
    (%go 0 (polynomial-coeffs p))))
