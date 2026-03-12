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
           #:maybe-prime-p
           #:*prime-source*))
(in-package :cl-polynomial/z)

(sera:-> factor ((integer 0))
         (values list &optional))
(defun factor (n)
  "Factor a non-negative integer @c(n) into prime numbers."
  (declare (optimize (speed 3)))
  (if (< n 2) (list n)
      (labels ((%factor (acc n x stop)
                 (declare (type integer n x stop))
                 (multiple-value-bind (q r)
                     (floor n x)
                   (cond
                     ((> x stop) (cons n acc))
                     ((zerop r)
                      (%factor (cons x acc) q 2 (isqrt q)))
                     (t
                      (%factor acc n (1+ x) stop))))))
        (%factor nil n 2 (isqrt n)))))

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

;; The sieve of eratosthenes
(sera:-> sieve (alex:positive-fixnum)
         (values (simple-array bit (*)) &optional))
(defun sieve (n)
  (declare (optimize (speed 3)))
  (let ((sieve (make-array n :element-type 'bit :initial-element 0)))
    (labels ((%go (p)
               (loop for i from (expt p 2) below n by p do
                     (setf (aref sieve i) 1))
               (let ((%p (position 0 sieve :start (1+ p))))
                 (if %p (%go %p) sieve))))
      (%go 2))))

(sera:-> eratosthenes-primes (alex:positive-fixnum)
         (values si:iterator &optional))
(defun eratosthenes-primes (n)
  "Make a fast iterator over primes less than or equal to @c(n)"
  (let ((sieve (sieve (1+ n))))
    (si:filter
     (lambda (n) (zerop (aref sieve n)))
     (si:range 2 (1+ n)))))

(declaim (type si:iterator *prime-source*))
(defparameter *prime-source*
  (si:concat
   (eratosthenes-primes 100000000)
   (si:undefined "Cannot compute primes beyond this limit ;("))
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


;; Additional thanks to nakrakiiya!
;; https://github.com/nakrakiiya/cl-prime-maker
(serapeum:-> expt-mod (unsigned-byte unsigned-byte unsigned-byte)
             (values unsigned-byte &optional))
(defun expt-mod (a n q)
  "Compute \\(a^n \\mod q\\)."
  (labels ((mul-rem (a b)
             (mod (* a b) q))
           (%expt-mod (a n acc)
             (cond
               ((zerop n) acc)
               ((evenp n)
                (%expt-mod (mul-rem a a) (floor n 2) acc))
               (t
                (%expt-mod a (1- n) (mul-rem a acc))))))
    (%expt-mod a n 1)))

(declaim (inline %random))
(defun %random (n)
  (+ (random (- n 3)) 2))

;; Fermat's little theorem states that if N is prime then A^N mod N = A. So
;; to test if N is prime we choose some random A which is less than N and
;; compute A^N mod N. If this is not equal to A then N is definitely not a
;; prime. If the test succeeds then A might be a prime (certain composite
;; numbers pass the Fermat test, these are called pseudo-primes), if we
;; perform the test over and over again then the probability of mis-classifying
;; the number reduces by roughly one half each time we perform the test. After
;; (say) one hundred iterations the probability of mis-classifying a number
;; is approximately 2^-100. So we can be fairly sure that the classification
;; is correct.

(declaim (inline test-prime-p))
(defun test-prime-p (n a)
  (= (expt-mod a (1- n) n) 1))

(serapeum:-> maybe-prime-p ((integer 1) &optional alex:positive-fixnum)
             (values boolean &optional))
(defun maybe-prime-p (a &optional (ntests 100))
  "Return @c(T) if @c(a) is probably prime.

This is the Fermat primality test which can have false positives."
  (labels ((%go (n)
             (declare (type fixnum n))
             (cond
               ((zerop n) t)
               ((test-prime-p a (%random a))
                (%go (1- n))))))
    (cond
      ((< a 2) nil)
      ((< a 4) t)
      (t (%go ntests)))))
