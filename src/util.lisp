(defpackage cl-polynomial/util
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria))
  (:export #:matrix
           #:row
           #:prime
           #:prime-power
           #:degree
           #:monomial

           #:gcdex
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
;; A power of a prime. In the Hensel lifting this often goes above the
;; most positive fixnum.
(deftype prime-power () '(integer 2))

;; A type for degree of a polynomial
(deftype degree () 'alex:non-negative-fixnum)

;; FIXME: let monomial be just an alias for a tuple?
(deftype monomial () '(cons degree integer))

;; Symmetric modulo operation widely used in cl-polynomial instead of
;; MOD.
;; NB: INTEGER may be replaced with FIXNUM only in the case of finite
;; fields. ZX:LIFT-FACTORS *REQUIRES* integers here.
(sera:-> mod-sym (integer prime-power)
         (values integer &optional))
(defun mod-sym (n q)
  "Compute \\(n \\mod q\\). The result is in a range \\(0 \\dots 1\\)
if \\(q = 2\\) or \\(-(q-1)/2 \\dots (q-1)/2\\) if \\(q > 2\\)."
  (let ((mod (mod n q))
        (half (floor q 2)))
    (if (or (= q 2)
            (<= mod half))
        mod (- mod q))))

(sera:-> gcdex (integer integer)
         (values integer integer integer &optional))
(defun gcdex (u v)
  "For \\(u, v \\in \\mathbb{Z}\\) find \\(\\gcd(u,v)\\) and solutions
of an equation \\(au + bv = \\gcd(u,v)\\) with minimal absolute
values."
  (cond
    ((zerop u)
     (values (abs v) 0 (signum v)))
    ((zerop v)
     (values (abs u) (signum u) 0))
    (t
     (labels ((%gcd (u v s0 s1 d0 d1)
                (multiple-value-bind (q r)
                    (round u v)
                  (let ((s (- s0 (* q s1)))
                        (d (- d0 (* q d1))))
                    (if (zerop r)
                        (values (abs v)
                                (* (signum v) s1)
                                (* (signum v) d1))
                        (%gcd v r s1 s d1 d))))))
       (if (> (abs u) (abs v))
           (%gcd u v 1 0 0 1)
           (multiple-value-bind (gcd s d)
               (%gcd v u 1 0 0 1)
             (values gcd d s)))))))

;; TODO: Update docs and types: works for any ring if the inverse exists.
(sera:-> invert-integer (integer prime-power)
         (values integer &optional))
(declaim (inline invert-integer))
(defun invert-integer (n q)
  "Find a multiplicative inverse of \\(n\\) in \\(\\mathbb{Z}_q\\), q
being a power of prime, i.e. find \\(x\\) such that \\(xn = nx =
1\\). Signal an error, if there is no such inverse."
  (when (zerop n)
    (error "Zero does not have a multiplicative inverse"))
  (multiple-value-bind (gcd a b)
      (gcdex n q)
    (declare (ignore b))
    (when (/= gcd 1)
      (error "N does not have a multiplicative inverse: N is not coprime with P"))
    (mod-sym a q)))

;; Like destructuring-bind, but without additional checks
(defmacro bind-monomial ((deg coeff) monomial &body body)
  (let ((m (gensym)))
    `(let ((,m ,monomial))
       (declare (type monomial ,m))
       (let ((,deg   (car ,m))
             (,coeff (cdr ,m)))
         ,@body))))
