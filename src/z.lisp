;; Some very basic and slow(!) number-theoretic functions

(defpackage cl-polynomial/z
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum))
  (:export #:factor #:totient #:moebius))
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
