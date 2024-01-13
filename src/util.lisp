(defpackage cl-polynomial/util
  (:use #:cl)
  (:local-nicknames (#:sera   #:serapeum)
                    (#:alex   #:alexandria))
  (:export #:matrix
           #:row
           #:prime
           #:monomial

           #:mod-sym
           #:invert-integer))
(in-package :cl-polynomial/util)

;; A type for matrices
(deftype matrix () '(simple-array (signed-byte 32) (* *)))

;; A type for rows
(deftype row () '(simple-array (signed-byte 32) (*)))

;; A bit more narrow range for primes
(deftype prime () '(integer 2 #.most-positive-fixnum))

;; FIXME: let monomial be just an alias for a tuple?
(deftype monomial () '(cons alex:non-negative-fixnum fixnum))

;; Symmetric modulo operation widely used in cl-polynomial instead of
;; MOD.
(sera:-> mod-sym (fixnum alex:positive-fixnum)
         (values fixnum &optional))
(defun mod-sym (x n)
  "Compute \\(x \\mod n\\). The result is in a range \\(0 \\dots 1\\)
if \\(n = 2\\) or \\(-(n-1)/2 \\dots (n-1)/2\\) if \\(n > 2\\)."
  (declare (optimize (speed 3)))
  (let ((mod (mod x n))
        (half (floor n 2)))
    (if (or (= n 2)
            (<= mod  half))
        mod (- mod n))))

(sera:-> invert-integer (fixnum prime)
         (values fixnum &optional))
(defun invert-integer (n p)
  "Find a multiplicative inverse of \\(n\\) in \\(\\mathbb{F}_p\\), p
being prime, i.e. find \\(x\\) such that \\(xn = nx = 1\\)."
  ;; Remember that n^p = n
  (declare (optimize (speed 3)))
  (mod-sym (expt n (- p 2)) p))
