(defpackage cl-polynomial/factor
  (:use #:cl)
  (:local-nicknames (#:sera #:serapeum)
                    (#:alex #:alexandria)
                    (#:u    #:cl-polynomial/util)
                    (#:p    #:cl-polynomial/polynomial)
                    (#:la   #:cl-polynomial/linalg))
  (:export #:square-free
           #:berlekamp-factor
           #:reducing-polynomials
           #:irreduciblep
           #:factor))
(in-package :cl-polynomial/factor)

(sera:-> x^pk-case (p:polynomial u:prime)
         (values p:polynomial &optional))
(defun x^pk-case (poly p)
  "Replace a polynomial in the form \\(\\sum_k a_k x^{b_k p}\\) with
\\(\\sum_k a_k x^{b_k}\\)."
  (p:polynomial
   (mapcar
    (lambda (m)
      (declare (type u:monomial m))
      (destructuring-bind (d . c) m
        (assert (zerop (rem d p)))
        (cons (/ d p) c)))
    (p:polynomial-coeffs poly))))

(sera:-> square-free (p:polynomial u:prime)
         (values list &optional))
(defun square-free (polynomial p)
  "Perform square-free factorization of a monic polynomial in
\\(\\mathbb{F}_p[x]\\), \\(p\\) being prime, with \\(\\deg f > 0\\). A list of
tuples \\((d_i . f_i)\\) is returned, so the supplied polynomial is equal to
\\(\\prod_i f_i^{d_i}\\)."
  (labels ((%%collect (p1 p2 acc n multiplicity)
             (if (p:polynomial= p2 p:+one+)
                 (values acc p1)
                 (let* ((gcd (p:gcd p1 p2 p))
                        (%p1 (p:divide p1 gcd p))
                        (%p2 (p:divide p2 gcd p)))
                   (%%collect %p1 gcd
                              (if (p:polynomial= %p2 p:+one+) acc
                                  (cons (cons (* n multiplicity) %p2) acc))
                              (1+ n) multiplicity))))
           (%collect (poly acc multiplicity)
             (let* ((gcd (p:gcd poly (p:modulo (p:derivative poly) p) p))
                    (w (p:divide poly gcd p)))
               (multiple-value-bind (%acc rest)
                   (%%collect gcd w acc 1 multiplicity)
                 (if (p:polynomial= rest p:+one+)
                     %acc
                     (%collect (x^pk-case rest p)
                               %acc (* multiplicity p)))))))
    (reverse
     (%collect polynomial nil 1))))

;; I really don't know how to name it
(sera:-> berlekamp-matrix (p:polynomial u:prime)
         (values u:matrix &optional))
(defun berlekamp-matrix (poly p)
  (let* ((degree (p:degree poly))
         (matrix (make-array (list degree degree)
                             :element-type '(signed-byte 32))))
    (dotimes (i degree)
      (let ((x^n-mod-poly (p:remainder (p:polynomial (list (cons (* i p) 1))) poly p)))
        (dotimes (j degree)
          (setf (aref matrix i j)
                (mod
                 (- (or (cdr (assoc j (p:polynomial-coeffs x^n-mod-poly))) 0)
                    (if (= i j) 1 0))
                 p)))))
    ;; The Berlekamp matrix + I has coefficients of polynomials x^{i
    ;; prime} mod polynomial for 0≤i<deg polynomial in its columns.
    matrix))

(sera:-> reducing-polynomials (p:polynomial u:prime)
         (values list &optional))
(defun reducing-polynomials (f p)
  "Return a list of all f-reducing polynomials for a monic non-constant
square-free polynomial \\(f \\in \\mathbb{F}_p[x]\\)."
  (let ((nullspace (la:nullspace (berlekamp-matrix f p) p)))
    (mapcar #'p:sequence->polynomial nullspace)))

;; TODO: Generalize to all polynomials
(sera:-> irreduciblep (p:polynomial u:prime)
         (values boolean &optional))
(defun irreduciblep (f p)
  "Test if a monic non-constant square-free polynomial \\(f \\in
\\mathbb{F}_p[x]\\) is irreducible."
  (= (length (la:nullspace (berlekamp-matrix f p) p)) 1))

;; TODO: This implementation does a lot of redundant work. Optimize it.
(sera:-> berlekamp-factor (p:polynomial u:prime)
         (values list &optional))
(defun berlekamp-factor (f p)
  "Given a monic square-free non-constant polynomial, return a list of its
factors."
  (let ((rps (remove p:+one+ (reducing-polynomials f p)
                     :test #'p:polynomial=)))
    (if (null rps)
        (list f)
        (labels ((collect-factors (g)
                   ;; Return non-trivial factors of g. They are not necessarily
                   ;; irreducible and not all irreducible factors are contained
                   ;; in this factorization.
                   (reduce
                    (lambda (acc rp)
                      (labels ((add-factors (n acc)
                                 (if (= n p) acc
                                     (add-factors
                                      (1+ n)
                                      (let ((gcd (p:gcd
                                                  g (p:modulo
                                                     (p:add rp (p:scale p:+one+ n)) p)
                                                  p)))
                                        (if (p:polynomial= gcd p:+one+) acc
                                            (cons gcd acc)))))))
                        (add-factors 0 acc)))
                    rps
                    :initial-value nil))
                 (more-factors (factors)
                   ;; (MORE-FACTORS (COLLECT-FACTORS G)) return *all*
                   ;; non-trivial factors.
                   (let* ((old-factors (remove-duplicates factors :test #'p:polynomial=))
                          (new-factors (remove-duplicates
                                        (alex:flatten
                                         ;; Here we calculate some redundant GCDs because
                                         ;; gcd(f, gcd(f, g)) = gcd(f, g)
                                         (mapcar #'collect-factors old-factors))
                                        :test #'p:polynomial=)))
                     ;; If there are nothing new, stop
                     (if (null (set-difference new-factors old-factors :test #'p:polynomial=))
                         new-factors
                         (more-factors new-factors)))))
          ;; Now we need to filter irreducible factors from this mess.
          (let ((factors (more-factors (collect-factors f))))
            (remove-if
             (lambda (g)
               (some
                (lambda (h)
                  (let ((gcd (p:gcd h g p)))
                    (not
                     (or (p:polynomial= gcd p:+one+)
                         (p:polynomial= gcd g)))))
                factors))
             factors))))))

;; TODO: Generalize to constant polynomials?
(sera:-> factor (p:polynomial u:prime)
         (values list fixnum &optional))
(defun factor (polynomial p)
  "Factor a non-constant polynomial in \\(\\mathbb{F}_p[x]\\) into irreducible factors."
  (multiple-value-bind (monic m)
      (p:monic-polynomial polynomial p)
    (values
     (reduce #'append
             (mapcar
              (lambda (factor)
                (destructuring-bind (m . f) factor
                  (mapcar
                   (lambda (irreducible-factor)
                     (cons m irreducible-factor))
                   (berlekamp-factor f p))))
              (square-free monic p)))
     m)))
