;; Some very basic and slow(!) number-theoretic functions

(defpackage cl-polynomial/z
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:si   #:stateless-iterators)
                    (#:u    #:cl-polynomial/util))
  (:export #:factor
           #:totient
           #:moebius
           #:*prime-source*))
(in-package :cl-polynomial/z)

(sera:-> factor ((integer 0))
         (values list &optional))
(defun factor (n)
  "Factor a non-negative integer @c(n) into prime numbers."
  (declare (optimize (speed 3)))
  (if (< n 2) (list n)
      (labels ((stop (n) (floor (sqrt n)))
               (%factor (acc n x stop)
                 (declare (type integer n x stop))
                 (multiple-value-bind (q r)
                     (floor n x)
                   (cond
                     ((> x stop) (cons n acc))
                     ((zerop r)
                      (%factor (cons x acc) q 2 (stop q)))
                     (t
                      (%factor acc n (1+ x) stop))))))
        (%factor nil n 2 (stop n)))))

(sera:-> totient ((integer 0))
         (values (integer 0) &optional))
(defun totient (n)
  "Euler's totient function."
  (declare (optimize (speed 3)))
  (if (< n 2) n
      (labels ((%totient (acc m)
                 (declare (type (integer 0) acc m))
                 (if (<= m n)
                     (%totient (if (= (gcd m n) 1) (1+ acc) acc) (1+ m))
                     acc)))
        (%totient 0 1))))

(sera:-> moebius ((integer 1))
         (values (integer -1 1) &optional))
(defun moebius (n)
  "Moebius function."
  (let* ((factors (factor n))
         (length (length factors)))
    (cond
      ((= n 1) 1)
      ((/= length (length (remove-duplicates factors :test #'=))) 0)
      ((evenp length) 1)
      (t -1))))

;; https://en.wikipedia.org/wiki/Formula_for_primes
(sera:-> get-prime ((cons unsigned-byte alex:positive-fixnum))
         (values u:prime (cons unsigned-byte alex:positive-fixnum) &optional))
(defun get-prime (state)
  (let* ((n! (car state))
         (n  (cdr state))
         (%n (1+ n)))
    (values
     (+ (* (1- n) (floor (mod n! %n) n)) 2)
     (cons (* n! %n) %n))))

(declaim (type si:iterator *prime-source*))
(defparameter *prime-source*
  (si:concat
   (si:singleton 2)
   (si:filter
    (lambda (x) (/= x 2))
    (si:unfold #'get-prime '(2 . 2))))
  "Iterator which returns primes.

@begin[lang=lisp](code)
(stateless-iterators:collect
  (stateless-iterators:take 100 cl-polynomial/z:*prime-source*)) =>
(2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103
 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199
 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313
 317 331 337 347 349 353 359 367 373 379 383 389 397 401 409 419 421 431 433
 439 443 449 457 461 463 467 479 487 491 499 503 509 521 523 541)
@end(code)")
