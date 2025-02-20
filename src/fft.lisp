(defpackage cl-polynomial/fft
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:z    #:cl-polynomial/z))
  (:export #:*-group-generator
           #:primitive-root-of-unity))
(in-package :cl-polynomial/fft)

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

(sera:-> pad-array ((simple-array integer (*)))
         (values (simple-array integer (*)) &optional))
(defun pad-array (array)
  (declare (optimize (speed 3)))
  (let ((%array (make-array (padded-length (length array))
                            :element-type 'integer
                            :initial-element 0)))
    (map-into %array #'identity array)))

(sera:-> phase-split ((simple-array integer (*)))
         (values (simple-array integer (*))
                 (simple-array integer (*))
                 &optional))
(defun phase-split (array)
  (declare (optimize (speed 3)))
  (let* ((length (length array))
         (phase1 (make-array (/ length 2) :element-type 'integer))
         (phase2 (make-array (/ length 2) :element-type 'integer)))
    (loop for i below length by 2
          for j from 0 by 1 do
          (setf (aref phase1 j) (aref array (+ i 0))
                (aref phase2 j) (aref array (+ i 1))))
    (values phase1 phase2)))

;; Output size = input size!
;; ω = n-th primitive root of unity in F_p, n being the input's length
(sera:-> %fft ((simple-array integer (*)) u:prime integer)
         (values (simple-array integer (*)) &optional))
(defun %fft (array p ω)
  (declare (optimize (speed 3)))
  (let* ((length   (length array))
         (length/2 (/ length 2)))
    (if (= length 1) array
        (multiple-value-bind (phase1 phase2)
            (phase-split array)
          (let ((fft1 (%fft phase1 p (u:mod-sym (expt ω 2) p)))
                (fft2 (%fft phase2 p (u:mod-sym (expt ω 2) p)))
                (result (make-array length :element-type 'integer)))
            (loop with ω^i = 1
                  for i below length/2
                  for even = (aref fft1 i)
                  for odd  = (aref fft2 i)
                  for %odd = (* odd ω^i) do
                  (setf (aref result i)
                        (u:mod-sym (+ even %odd) p)
                        (aref result (+ i length/2))
                        (u:mod-sym (- even %odd) p)
                        ω^i (u:mod-sym (* ω^i ω) p)))
            result)))))

(sera:-> renormalize ((simple-array integer (*)) u:prime)
         (values (simple-array integer (*)) &optional))
(defun renormalize (array p)
  (let ((m (u:invert-integer (length array) p)))
    (map '(vector integer)
         (lambda (x)
           (u:mod-sym (* x m) p))
         array)))
