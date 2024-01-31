(defpackage cl-polynomial/util
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria))
  (:export #:matrix
           #:row
           #:prime
           #:degree
           #:monomial

           #:mod-sym
           #:invert-integer
           #:bind-monomial))
(in-package :cl-polynomial/util)

;; A type for matrices
(deftype matrix () '(simple-array (signed-byte 32) (* *)))

;; A type for rows
(deftype row () '(simple-array (signed-byte 32) (*)))

;; A bit more narrow range for primes
(deftype prime () '(integer 2 #.most-positive-fixnum))

;; A type for degree of a polynomial
(deftype degree () 'alex:non-negative-fixnum)

;; FIXME: let monomial be just an alias for a tuple?
(deftype monomial () '(cons degree integer))

;; Symmetric modulo operation widely used in cl-polynomial instead of
;; MOD.
(sera:-> mod-sym (fixnum prime)
         (values fixnum &optional))
(defun mod-sym (x n)
  "Compute \\(x \\mod n\\). The result is in a range \\(0 \\dots 1\\)
if \\(n = 2\\) or \\(-(n-1)/2 \\dots (n-1)/2\\) if \\(n > 2\\)."
  (declare (optimize (speed 3)))
  (let ((mod (mod x n))
        (half (floor n 2)))
    (if (or (= n 2)
            (<= mod half))
        mod (- mod n))))

(sera:-> expt-rem (fixnum alex:non-negative-fixnum prime)
         (values fixnum &optional))
(defun expt-rem (x n p)
  "Calculate \\(x^n \\mod p\\) with \\(n \\ge 0\\)."
  (declare (optimize (speed 3)))
  (labels ((%expt-rem (x n acc)
             (declare (type alex:non-negative-fixnum n)
                      (type fixnum x acc))
             (cond
               ((zerop n) acc)
               ((evenp n)
                (%expt-rem (mod-sym (* x x) p) (floor n 2) acc))
               (t
                (%expt-rem x (1- n) (mod-sym (* acc x) p))))))
    (%expt-rem x n 1)))

(sera:-> invert-integer (fixnum prime)
         (values fixnum &optional))
(declaim (inline invert-integer))
(defun invert-integer (n p)
  "Find a multiplicative inverse of \\(n\\) in \\(\\mathbb{F}_p\\), p
being prime, i.e. find \\(x\\) such that \\(xn = nx = 1\\)."
  (declare (optimize (speed 3)))
  ;; Remember that n^p = n
  (expt-rem n (- p 2) p))

;; Like destructuring-bind, but without additional checks
(defmacro bind-monomial ((deg coeff) monomial &body body)
  (let ((m (gensym)))
    `(let ((,m ,monomial))
       (declare (type monomial ,m))
       (let ((,deg   (car ,m))
             (,coeff (cdr ,m)))
         ,@body))))
