(defpackage cl-polynomial/fp
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:z    #:cl-polynomial/z))
  (:export #:*-group-generator
           #:primitive-root-of-unity))
(in-package :cl-polynomial/fp)

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
