(defpackage cl-polynomial/fft
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:p    #:cl-polynomial/polynomial)
                    (#:si   #:stateless-iterators)
                    (#:z    #:cl-polynomial/z))
  (:export #:fourier-primes
           #:*-group-generator
           #:primitive-root-of-unity
           #:pad-array
           #:fft
           #:ifft
           #:polynomial->vector))
(in-package :cl-polynomial/fft)

;; TODO: requires a faster prime source
(sera:-> fourier-primes (fixnum)
         (values si:iterator &optional))
(defun fourier-primes (n)
  "Return primes in the form \\(2^l k + 1\\) where \\(k\\) is odd and
\\(l\\) is an integer greater than zero."
  (unless (zerop (logand n (1- n)))
    (error "N must be a power of 2."))
  (si:filter
   (lambda (p)
     (let ((m (/ (1- p) n)))
       (and (integerp m) (oddp m))))
   z:*prime-source*))

(sera:-> unique-factors (alex:positive-fixnum)
         (values list &optional))
(defun unique-factors (n)
  "Return unique prime factors of \\(n\\) in ascending order."
  (declare (optimize (speed 3)))
  (reduce
   (lambda (acc x)
     (declare (type fixnum x))
     (let ((last (car acc)))
       (declare (type (or null fixnum) last))
       (if (and last (= x last)) acc
           (cons x acc))))
   (z:factor n)
   :initial-value nil))

(sera:-> *-group-generator (u:prime)
         (values fixnum &optional))
(defun *-group-generator (p)
  "Get a generator of a multiplicative group of a field \\(\\mathbb{F}_p\\)."
  (declare (optimize (speed 3)))
  (if (= p 2) 1
      (let ((factors (unique-factors (1- p))))
        (flet ((generatorp (n)
                 (every (lambda (m)
                          (declare (type u:prime m))
                          (let ((a (/ (1- p) m)))
                            (/= (u:expt-mod n a p) 1)))
                        factors)))
          (loop for gen = (+ (random (- p 2)) 2)
                until (generatorp gen)
                finally (return gen))))))

(sera:-> primitive-root-of-unity (alex:positive-fixnum u:prime)
         (values fixnum &optional))
(defun primitive-root-of-unity (n p)
  "Get an \\(n\\)-th primitive root of unity in \\(\\mathbb{F}_p\\)."
  (declare (optimize (speed 3)))
  (unless (zerop (rem (1- p) n))
    (error "There is no ~d-th roots of unity in this field" n))
  (let ((gen (*-group-generator p))
        (m (/ (1- p) n)))
    (u:expt-mod gen m p)))

(sera:-> padded-length (alex:positive-fixnum)
         (values alex:positive-fixnum &optional))
(declaim (inline padded-length))
(defun padded-length (n)
  (ash 1 (integer-length (1- n))))

(sera:-> pad-array ((simple-array fixnum (*)) &optional alex:positive-fixnum)
         (values (simple-array fixnum (*)) &optional))
(defun pad-array (array &optional (n (length array)))
  "Pad array with zeros to the smallest power of two which is equal or
greater than \\(n\\)."
  (declare (optimize (speed 3)))
  (let ((%array (make-array (padded-length n)
                            :element-type 'fixnum
                            :initial-element 0)))
    (map-into %array #'identity array)))

(sera:-> renormalize ((simple-array fixnum (*)) u:prime)
         (values (simple-array fixnum (*)) &optional))
(defun renormalize (array p)
  (let ((m (u:invert-integer (length array) p)))
    (map '(vector fixnum)
         (lambda (x)
           (u:mod-sym (* x m) p))
         array)))

;; Requirement: Array length is a power of 2
(sera:-> reverse-bits (alex:positive-fixnum alex:non-negative-fixnum)
         (values alex:non-negative-fixnum &optional))
(declaim (inline reverse-bits))
(defun reverse-bits (length n)
  (loop with length = (integer-length (1- length))
        for i below length sum
        (ash (ldb (byte 1 i) n) (- length i 1))
        fixnum))

;; Requirement: Array length is a power of 2
(sera:-> reorder-input ((simple-array fixnum (*)))
         (values (simple-array fixnum (*)) &optional))
(defun reorder-input (array)
  (declare (optimize (speed 3)))
  (let* ((length (length array))
         (result (make-array length :element-type 'fixnum)))
    (loop for i below length do
          (setf (aref result i)
                (aref array (reverse-bits length i))))
    result))

;; Collect a list of "group ωs", that is ω^{2^{steps-k-1}} for k = 0,1,…,steps-1
(sera:-> ωs (fixnum alex:positive-fixnum u:prime)
         (values list &optional))
(defun ωs (ω n p)
  (declare (optimize (speed 3)))
  (labels ((%go (acc k)
             (declare (type fixnum k))
             (if (zerop k) acc
                 (let ((ω (car acc)))
                   (%go (cons (u:mod-sym (expt ω 2) p) acc)
                        (1- k))))))
    (%go (list ω) (1- n))))

;; Requirement: Array length is a power of 2
(sera:-> %fft! ((simple-array fixnum (*)) fixnum u:prime)
         (values (simple-array fixnum (*)) &optional))
(defun %fft! (array ω p)
  (declare (optimize (speed 3)))
  (let* ((length (length array))
         (steps  (integer-length (1- length)))
         (ωs (ωs ω steps p)))
    (loop for s below steps
          for pair-offset  = (ash 1 s)
          for ngroups = (ash 1 (- steps s 1))
          for group-size = pair-offset ;; number of evens
          for m fixnum in ωs do
          (loop for i below ngroups
                for group-offset fixnum from 0 by (* group-size 2) do
                (loop for j below group-size
                      with %m = 1
                      for even-idx = (+ group-offset j)
                      for odd-idx  = (+ even-idx pair-offset)
                      for even = (aref array even-idx)
                      for odd  = (aref array odd-idx)
                      for %odd fixnum = (* odd %m) do
                      (setf (aref array even-idx)
                            (u:mod-sym (+ even %odd) p)
                            (aref array odd-idx)
                            (u:mod-sym (- even %odd) p)
                            %m (u:mod-sym (the fixnum (* %m m)) p)))))
    array))

(declaim (inline sanity-checks))
(defun sanity-checks (array ω p)
  (let ((n (length array)))
    (unless (zerop (logand n (1- n)))
      (error "Length of the input array is not a power of 2"))
    (unless (= (u:expt-mod ω n p) 1)
      (error "ω is not a root of unity we need"))))

(sera:-> fft ((simple-array fixnum (*)) fixnum u:prime)
         (values (simple-array fixnum (*)) &optional))
(defun fft (array ω p)
  "Perform forward FFT of @c(array) in a field
\\(\\mathbb{F}_p\\). \\(\\omega\\) is an \\(n\\)-th primitive root of
unity in that field, where \\(n\\) is the length of @c(array). The
length \\(n\\) must be a positive integral power of two."
  (sanity-checks array ω p)
  (%fft! (reorder-input array) ω p))

(sera:-> ifft ((simple-array fixnum (*)) fixnum u:prime)
         (values (simple-array fixnum (*)) &optional))
(defun ifft (array ω p)
  "Invert FFT.

@begin[lang=lisp](code)
(every #'= a (ifft (fft a p ω) p ω)) ; Evaluates to T
@end(code)"
  (sanity-checks array ω p)
  (renormalize
   (%fft! (reorder-input array) (u:invert-integer ω p) p) p))

(sera:-> polynomial->vector (p:polynomial &optional alex:positive-fixnum)
         (values (simple-array fixnum (*)) &optional))
(defun polynomial->vector (p &optional (n (1+ (p:degree p))))
  "Convert a polynomial to a vector suitable for use in @c(fft)
function. The resulting vector has a size which is the smallest
positive integral power of two which is greater or equal to \\(n\\). If
\\(n\\) is supplied, it must be greater than the degree of \\(p\\)."
  (unless (> n (p:degree p))
    (error "N must be greater than the degree of P."))
  (let ((array (make-array (padded-length n) :element-type 'fixnum :initial-element 0)))
    (loop for (deg . c) in (p:polynomial-coeffs p) do
          (setf (aref array deg) c))
    array))
