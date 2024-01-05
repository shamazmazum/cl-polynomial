(in-package :polynomial)

(sera:-> root (polynomial prime)
         (values polynomial &optional))
(defun root (poly p)
  (polynomial
   (mapcar
    (lambda (m)
      (declare (type monomial m))
      (destructuring-bind (d . c) m
        (assert (zerop (rem d p)))
        (cons (/ d p) c)))
    (polynomial-coeffs poly))))

(sera:-> square-free (polynomial prime)
         (values list &optional))
(defun square-free (polynomial p)
  "Perform square-free factorization of a monic polynomial in
\\(\\mathbb{F}_p[x]\\), \\(p\\) being prime, with \\(\\deg f > 0\\). A list of
tuples \\((d_i . f_i)\\) is returned, so the supplied polynomial is equal to
\\(\\prod_i f_i^{d_i}\\)."
  (labels ((%%collect (p1 p2 acc n multiplicity)
             (if (polynomial= p2 +one+)
                 (values acc p1)
                 (let* ((gcd (gcd p1 p2 p))
                        (%p1 (divide p1 gcd p))
                        (%p2 (divide p2 gcd p)))
                   (%%collect %p1 gcd
                              (if (polynomial= %p2 +one+) acc
                                  (cons (cons (* n multiplicity) %p2) acc))
                              (1+ n) multiplicity))))
           (%collect (poly acc multiplicity)
             (let* ((gcd (gcd poly (modulo (derivative poly) p) p))
                    (w (divide poly gcd p)))
               (multiple-value-bind (%acc rest)
                   (%%collect gcd w acc 1 multiplicity)
                 (if (polynomial= rest +one+)
                     %acc
                     (%collect (root rest p) %acc (* multiplicity p)))))))
    (reverse
     (%collect polynomial nil 1))))

;; I really don't know how to name it
(sera:-> berlekamp-matrix (polynomial prime)
         (values matrix &optional))
(defun berlekamp-matrix (poly p)
  (let* ((degree (degree poly))
         (matrix (make-array (list degree degree)
                             :element-type '(signed-byte 32))))
    (dotimes (i degree)
      (let ((x^n-mod-poly (remainder (polynomial (list (cons (* i p) 1))) poly p)))
        (dotimes (j degree)
          (setf (aref matrix i j)
                (mod
                 (- (or (cdr (assoc j (polynomial-coeffs x^n-mod-poly))) 0)
                    (if (= i j) 1 0))
                 p)))))
    matrix))

(sera:-> reducing-polynomials (polynomial prime)
         (values list &optional))
(defun reducing-polynomials (f p)
  "Return a list of all f-reducing polynomials for a monic non-constant
square-free polynomial \\(f \\in \\mathbb{F}_p[x]\\)."
  (let ((nullspace (nullspace (berlekamp-matrix f p) p)))
    (mapcar #'sequence->polynomial nullspace)))

;; TODO: Generalize to all polynomials
(sera:-> irreduciblep (polynomial prime)
         (values boolean &optional))
(defun irreduciblep (f p)
  "Test if a monic non-constant square-free polynomial \\(f \\in
\\mathbb{F}_p[x]\\) is irreducible."
  (= (length (nullspace (berlekamp-matrix f p) p)) 1))

;; TODO: This implementation does a lot of redundant work. Optimize it.
(sera:-> berlekamp-factor (polynomial prime)
         (values list &optional))
(defun berlekamp-factor (f p)
  "Given a monic square-free non-constant polynomial, return a list of its
factors."
  (let ((rps (remove +one+ (reducing-polynomials f p)
                     :test #'polynomial=)))
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
                                      (let ((gcd
                                             (gcd g (modulo (add rp (scale +one+ n)) p) p)))
                                        (if (polynomial= gcd +one+) acc (cons gcd acc)))))))
                        (add-factors 0 acc)))
                    rps
                    :initial-value nil))
                 (more-factors (factors)
                   ;; (MORE-FACTORS (COLLECT-FACTORS G)) return *all*
                   ;; non-trivial factors.
                   (let* ((old-factors (remove-duplicates factors :test #'polynomial=))
                          (new-factors (remove-duplicates
                                        (alex:flatten
                                         ;; Here we calculate some redundant GCDs because
                                         ;; gcd(f, gcd(f, g)) = gcd(f, g)
                                         (mapcar #'collect-factors old-factors))
                                        :test #'polynomial=)))
                     ;; If there are nothing new, stop
                     (if (null (set-difference new-factors old-factors :test #'polynomial=))
                         new-factors
                         (more-factors new-factors)))))
          ;; Now we need to filter irreducible factors from this mess.
          (let ((factors (more-factors (collect-factors f))))
            (remove-if
             (lambda (g)
               (some
                (lambda (h)
                  (let ((gcd (gcd h g p)))
                    (not
                     (or (polynomial= gcd +one+)
                         (polynomial= gcd g)))))
                factors))
             factors))))))

;; TODO: Generalize to constant polynomials?
(sera:-> factor (polynomial prime)
         (values list fixnum &optional))
(defun factor (polynomial p)
  "Factor a non-constant polynomial in \\(\\mathbb{F}_p[x]\\) into irreducible factors."
  (multiple-value-bind (monic m)
      (monic-polynomial polynomial p)
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
