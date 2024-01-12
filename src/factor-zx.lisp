(in-package :polynomial)

(sera:-> content-free (polynomial)
         (values polynomial fixnum &optional))
(defun content-free (f)
  "For polynomial \\(f \\in \\mathbb{Z}[x]\\) return its primitive
part and content."
  (let ((content (reduce
                  (lambda (acc m)
                    (declare (type monomial m))
                    (cl:gcd acc (cdr m)))
                  (polynomial-coeffs f)
                  :initial-value 0)))
    (values
     (map-poly (lambda (x) (/ x content)) f)
     content)))

(defun replace-lc (f c)
  "Replace leading coefficient of non-zero @c(f) with @c(c)."
  (if (polynomial= f +zero+)
      (scale +one+ c)
      (polynomial
       (destructuring-bind ((d . %c) &rest rest)
           (polynomial-coeffs f)
         (declare (ignore %c))
         (cons (cons d c) rest)))))

;; TODO: Rewrite this
(sera:-> modulo-symmetric (polynomial prime)
         (values polynomial &optional))
(defun modulo-symmetric (polynomial prime)
  (if (= prime 2)
      ;; 2 is the only even prime, so there is no symmetric
      ;; representation.
      (modulo polynomial prime)
      (let ((half (/ (1- prime) 2)))
        (map-poly (lambda (x)
                    (if (<= x half) x (- x prime)))
                  (modulo polynomial prime)))))

(sera:-> lift-factors (polynomial polynomial polynomial prime alex:non-negative-fixnum)
         (values polynomial polynomial boolean alex:non-negative-fixnum &optional))
(defun lift-factors (f f1 f2 p d)
  "For a polynomial \\(f \\in \\mathbb{Z}[x]\\) and \\(f_1, f_2 \\in
\\mathbb{Z}_p[x]\\), \\(p\\) being prime, find polynomials
\\(\\hat{f_1}, \\hat{f_2} \\in \\mathbb{Z}[x]\\) such that \\(f =
\\hat{f_1} \\hat{f_2}\\) if \\(f = f_1 f_2 \\mod p\\). \\(d\\) is a
maximal absolute value of coefficients in \\(f\\) or its factors.

The first two values returned are the desired factors. The third value
is a boolean being equal to @c(T) if the algorithm has successfully
found a solution in \\(\\mathbb{Z}[x]\\). If this value is @c(NIL)
then the algorithm is not successful and \\(f = \\hat{f_1} \\hat{f_2}
\\mod p^N\\) where \\(N\\) is the forth returned value."
  ;; Adjust f, f1 and f2 to correctly solve non-monic case
  ;; f ∈ ℤ[x], f₁, f₂ ∈ ℤ_p[x]
  (let* ((lc (leading-coeff f))
         ;; p^max-steps > 2D*lc must hold
         (max-steps (1+ (ceiling (log (* 2 d lc) p))))
         (f  (identity (scale f  (* lc))))
         (f1 (modulo   (scale f1 (* lc (invert-integer (leading-coeff f1) p))) p))
         (f2 (modulo   (scale f2 (* lc (invert-integer (leading-coeff f2) p))) p)))
    (multiple-value-bind (gcd s d)
        (gcdex f1 f2 p)
      (declare (ignore gcd))
      (labels ((%lift-factors (%f1 %f2 q step)
                 (let* ((%f1 (replace-lc %f1 lc))
                        (%f2 (replace-lc %f2 lc))
                        (prod (multiply %f1 %f2))
                        (diff (subtract f prod)))
                   (if (or (polynomial= diff +zero+)
                           (= step max-steps))
                       (values (content-free %f1)
                               (content-free %f2)
                               (polynomial= diff +zero+)
                               step)
                       (let* ((rhs (map-poly (lambda (x) (/ x q)) diff))
                              (%δf2 (modulo-symmetric (multiply s rhs) p))
                              (%δf1 (modulo-symmetric (multiply d rhs) p)))
                         (multiple-value-bind (quo r)
                             (divide %δf2 %f2 p)
                           (let ((δf2 (modulo-symmetric r p))
                                 (δf1 (modulo-symmetric (add (multiply %f1 quo) %δf1) p)))
                             #+nil
                             (assert
                              (polynomial=
                               (modulo-symmetric rhs p)
                               (modulo-symmetric
                                (add
                                 (multiply %f1 δf2)
                                 (multiply %f2 δf1))
                                p)))
                             (%lift-factors (add %f1 (scale δf1 q))
                                            (add %f2 (scale δf2 q))
                                            (* p q)
                                            (1+ step)))))))))
        (%lift-factors f1 f2 p 0)))))

(sera:-> suitable-prime (polynomial)
         (values prime &optional))
(defun suitable-prime (polynomial)
  "Return a suitable prime for reducing a factorization in
\\(\\mathbb{Z}[x]\\) to a factorization in \\(mathbb{F}_p[x]\\)."
  (assert (not (polynomial= polynomial +zero+)))
  (let ((lc (leading-coeff polynomial)))
    (nth-value
     0 (si:consume-one
        (si:drop-while
         (lambda (prime) (zerop (rem lc prime)))
         (si:imap
          (lambda (n) (primes:get-nth-prime n))
          (si:count-from 1)))))))
