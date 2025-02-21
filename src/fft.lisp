(defpackage cl-polynomial/fft
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:si   #:stateless-iterators)
                    (#:z    #:cl-polynomial/z))
  (:export #:fourier-primes
           #:*-group-generator
           #:primitive-root-of-unity
           #:pad-array
           #:fft
           #:ifft))
(in-package :cl-polynomial/fft)

;; TODO: requires a faster prime source
(sera:-> fourier-primes (unsigned-byte)
         (values si:iterator &optional))
(defun fourier-primes (n)
  "Return primes in the form \\(k n + 1\\) where \\(k\\) is odd and \\(n = 2^l, l > 0\\)."
  (unless (zerop (logand n (1- n)))
    (error "N must be a power of 2."))
  (si:filter
   (lambda (p)
     (let ((m (/ (1- p) n)))
       (and (integerp m) (oddp m))))
   z:*prime-source*))

(sera:-> unique-factors ((integer 0))
         (values list &optional))
(defun unique-factors (n)
  "Return unique prime factors of \\(n\\) in ascending order."
  (declare (optimize (speed 3)))
  (reduce
   (lambda (acc x)
     (let ((last (car acc)))
       (declare (type integer x)
                (type (or null integer) last))
       (if (and last (= x last)) acc
           (cons x acc))))
   (z:factor n)
   :initial-value nil))

(sera:-> *-group-generator (u:prime)
         (values integer &optional))
(defun *-group-generator (p)
  "Get a generator of a multiplicative group of a field \\(\\mathbb{F}_p\\)."
  (declare (optimize (speed 3)))
  (if (= p 2) 1
      (let ((factors (unique-factors (1- p))))
        (flet ((generatorp (n)
                 (every (lambda (m)
                          (declare (type u:prime m))
                          (let ((a (/ (1- p) m)))
                            (declare (type integer a))
                            (/= (u:mod-sym (expt n a) p) 1)))
                        factors)))
          (loop for gen = (+ (random (- p 2)) 2)
                until (generatorp gen)
                finally (return gen))))))

(sera:-> primitive-root-of-unity (u:prime unsigned-byte)
         (values integer &optional))
(defun primitive-root-of-unity (p n)
  (declare (optimize (speed 3)))
  (unless (zerop (rem (1- p) n))
    (error "There is no ~d-th roots of unity in this field" n))
  (let ((gen (*-group-generator p))
        (m (/ (1- p) n)))
    (declare (type integer m))
    (u:mod-sym (expt gen m) p)))

(sera:-> padded-length (alex:positive-fixnum)
         (values alex:positive-fixnum &optional))
(declaim (inline padded-length))
(defun padded-length (n)
  (ash 1 (integer-length (1- n))))

(sera:-> pad-array ((simple-array integer (*)) &optional alex:positive-fixnum)
         (values (simple-array integer (*)) &optional))
(defun pad-array (array &optional (n (length array)))
  "Pad array with zeros to the smallest power of two which is equal or
greater than \\(n\\)."
  (declare (optimize (speed 3)))
  (let ((%array (make-array (padded-length n)
                            :element-type 'integer
                            :initial-element 0)))
    (map-into %array #'identity array)))

(sera:-> renormalize ((simple-array integer (*)) u:prime)
         (values (simple-array integer (*)) &optional))
(defun renormalize (array p)
  (let ((m (u:invert-integer (length array) p)))
    (map '(vector integer)
         (lambda (x)
           (u:mod-sym (* x m) p))
         array)))

;; Requirement: Array length is a power of 2
(sera:-> reverse-bits ((integer 1) unsigned-byte)
         (values unsigned-byte &optional))
(declaim (inline reverse-bits))
(defun reverse-bits (length n)
  (let ((length (integer-length (1- length))))
    (si:foldl
     (lambda (acc i)
       (logior acc (ash (ldb (byte 1 i) n) (- length i 1))))
     0 (si:range 0 length))))

;; Requirement: Array length is a power of 2
(sera:-> reorder-input ((simple-array integer (*)))
         (values (simple-array integer (*)) &optional))
(defun reorder-input (array)
  (declare (optimize (speed 3)))
  (let* ((length (length array))
         (result (make-array length :element-type 'integer)))
    (loop for i below length do
          (setf (aref result i)
                (aref array (reverse-bits length i))))
    result))

;; Requirement: Array length is a power of 2
(sera:-> %fft! ((simple-array integer (*)) u:prime integer)
         (values (simple-array integer (*)) &optional))
(defun %fft! (array p ω)
  (declare (optimize (speed 3)))
  (let* ((length (length array))
         (steps  (integer-length (1- length))))
    (loop for s below steps
          for pair-offset  = (ash 1 s)
          for ngroups = (ash 1 (- steps s 1))
          for group-size = pair-offset ;; number of evens
          for m = (u:mod-sym (expt ω ngroups) p) do
          (loop for i below ngroups
                for group-offset fixnum from 0 by (* group-size 2) do
                (loop for j below group-size
                      with %m = 1
                      for even-idx = (+ group-offset j)
                      for odd-idx  = (+ even-idx pair-offset)
                      for even = (aref array even-idx)
                      for odd  = (aref array odd-idx)
                      for %odd = (* odd %m) do
                      (setf (aref array even-idx)
                            (u:mod-sym (+ even %odd) p)
                            (aref array odd-idx)
                            (u:mod-sym (- even %odd) p)
                            %m (u:mod-sym (* %m m) p)))))
    array))

(declaim (inline sanity-checks))
(defun sanity-checks (array p ω)
  (let ((n (length array)))
    (unless (zerop (logand n (1- n)))
      (error "Length of the input array is not a power of 2"))
    (unless (= (u:mod-sym (expt ω n) p) 1)
      (error "ω is not a root of unity we need"))))

(sera:-> fft ((simple-array integer (*)) u:prime integer)
         (values (simple-array integer (*)) &optional))
(defun fft (array p ω)
  "Perform forward FFT of ARRAY in a field
\\(\\mathbb{F}_p\\). \\(\\omega\\) is an \\(n\\)-th primitive root of
unity in that field, where \\(n\\) is the length of ARRAY. The length
\\(n\\) must be a positive integer power of two."
  (sanity-checks array p ω)
  (%fft! (reorder-input array) p ω))

(sera:-> ifft ((simple-array integer (*)) u:prime integer)
         (values (simple-array integer (*)) &optional))
(defun ifft (array p ω)
  "Invert FFT.

@begin[lang=lisp](code)
(every #'= a (ifft (fft a p ω) p ω)) ; Evaluates to T
@end(code)"
  (sanity-checks array p ω)
  (renormalize
   (%fft! (reorder-input array) p (u:invert-integer ω p)) p))
