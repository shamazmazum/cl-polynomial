(in-package :cl-polynomial/tests)

(alex:define-constant +some-primes+ '(2 3 5 7 11 17 19)
  :test #'equalp)

(defun run-tests ()
  (every #'identity
         (mapcar (lambda (suite)
                   (let ((status (run suite)))
                     (explain! status)
                     (results-status status)))
                 '(algebra factor))))

(defun coeffs-sorted-p (polynomial)
  (let* ((coeffs (p::polynomial-coeffs polynomial))
         (sorted (sort (copy-seq coeffs) #'> :key #'car)))
    (equalp coeffs sorted)))

(defun random-poly (max-coeff state &optional (max-degree 100))
  (p:list->polynomial
   (loop repeat (1+ (random max-degree state))
         collect (- (random max-coeff state)
                    (if (= max-coeff 2) 0 (floor max-coeff 2))))))

(defun random-prime (state)
  (let ((length (length +some-primes+)))
    (nth (random length state) +some-primes+)))

(defun ratsimp (factors)
  (apply
   #'p:multiply
   (mapcar
    (lambda (factor)
      (destructuring-bind (m . f) factor
        (apply #'p:multiply
               (loop repeat m collect f))))
    factors)))

(def-suite algebra :description "Generic algebraic operations on polynomials")
(def-suite factor  :description "Factorization of polynomials over finite fields")

(in-suite algebra)

(test addition/subtraction
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((poly1 (random-poly 100 state))
               (poly2 (random-poly 100 state))
               (sum (p:add poly1 poly2)))
          (is-true (coeffs-sorted-p sum))
          (is (p:polynomial= sum (p:add poly2 poly1)))
          (is (p:polynomial= (p:add poly1 p:+zero+) poly1))
          (is (<= (p:degree sum)
                  (max (p:degree poly1)
                       (p:degree poly2))))
          (is (p:polynomial= poly1 (p:subtract sum poly2)))
          (is (p:polynomial= poly2 (p:subtract sum poly1))))))
          
(test multiplication
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state))
               (product (p:modulo (p:multiply poly1 poly2) prime)))
          (is-true (coeffs-sorted-p product))
          (is (p:polynomial= (p:multiply poly1 p:+zero+) p:+zero+))
          (is (p:polynomial= (p:multiply poly1 p:+one+) poly1))
          (is (p:polynomial= product (p:modulo (p:multiply poly2 poly1) prime)))
          (unless (or (p:polynomial= poly1 p:+zero+)
                      (p:polynomial= poly2 p:+zero+))
            (is (= (p:degree product)
                   (+ (p:degree poly1)
                      (p:degree poly2)))))
          (unless (p:polynomial= poly1 p:+zero+)
            (multiple-value-bind (q r)
                (p:divide product poly1 prime)
              (is (p:polynomial= q poly2))
              (is (p:polynomial= r p:+zero+)))))))

(test division
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state)))
          (unless (p:polynomial= poly2 p:+zero+)
            (multiple-value-bind (q r)
                (p:divide poly1 poly2 prime)
              (is-true (coeffs-sorted-p q))
              (is-true (coeffs-sorted-p r))
              (is-true (or (< (p:degree r) (p:degree poly2))
                           (= 0 (p:degree r) (p:degree poly2))))
              (is (p:polynomial=
                   poly1
                   (p:modulo (p:add (p:multiply q poly2) r) prime))))))))

(test monic
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly (random-poly prime state)))
          (unless (p:polynomial= poly p:+zero+)
            (multiple-value-bind (monic c)
                (p:monic-polynomial poly prime)
              (is-true (coeffs-sorted-p monic))
              (is (= (p:leading-coeff monic) 1))
              (is (p:polynomial=
                   poly
                   (p:modulo (p:scale monic c) prime))))))))

(test gcd
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state)))
          (unless (or (p:polynomial= poly1 p:+zero+)
                      (p:polynomial= poly2 p:+zero+))
            (let* ((gcd (p:gcd poly1 poly2 prime))
                   (lcm (p:divide (p:modulo (p:multiply poly1 poly2) prime) gcd prime)))
              (is (p:polynomial= gcd (p:gcd poly2 poly1 prime)))
              (is-true (p:monicp gcd))
              (unless (p:constantp gcd)
                (is (p:polynomial= (p:remainder poly1 gcd prime) p:+zero+))
                (is (p:polynomial= (p:remainder poly2 gcd prime) p:+zero+))
                (is (p:polynomial= (p:remainder lcm poly1 prime) p:+zero+))
                (is (p:polynomial= (p:remainder lcm poly2 prime) p:+zero+))))))))

(test gcdex
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state)))
          (multiple-value-bind (gcd a b)
              (p:gcdex poly1 poly2 prime)
            (is (p:polynomial= gcd (p:gcd poly1 poly2 prime)))
            (is (p:polynomial= (p:modulo
                                (p:add (p:multiply poly1 a)
                                       (p:multiply poly2 b))
                                prime)
                               gcd))))))

(in-suite factor)

(test square-free
  (loop with state = (make-random-state t)
        repeat 10000 do
        ;; Generate polynomials with low degree in 𝔽_p[x] with p equal
        ;; to 2 or 3, so the result may have non-trivial factorization
        ;; with high amount of probability.
        (let* ((prime (+ 2 (random 2 state)))
               (polynomial (p:monic-polynomial (random-poly prime state 20) prime)))
          (unless (p:polynomial= polynomial p:+zero+)
            (is (p:polynomial=
                 (p:modulo (ratsimp (p:square-free polynomial prime)) prime)
                 polynomial))))))

(test reducing-polys
  (loop with state = (make-random-state t)
        repeat 10000 do
        ;; Generate polynomials with low degree in 𝔽_p[x] with p equal
        ;; to 2 or 3, so the result may have non-trivial factorization
        ;; with high amount of probability.
        (let* ((prime (+ 2 (random 2 state)))
               (polynomial (p:monic-polynomial (random-poly prime state 20) prime)))
          (unless (zerop (p:degree polynomial))
            (dolist (rp (p:reducing-polynomials polynomial prime))
              (is (p:polynomial=
                   (p:remainder (p:modulo
                                 (apply #'p:multiply
                                        (loop repeat prime collect rp))
                                 prime)
                                polynomial prime)
                   rp)))))))

(test factor-finite
  (loop with state = (make-random-state t)
        repeat 10000 do
        ;; The same comment as above
        (let* ((prime (+ 2 (random 2 state)))
               (polynomial (p:monic-polynomial (random-poly prime state 20) prime)))
          (unless (zerop (p:degree polynomial))
            (multiple-value-bind (factors c)
                (p:factor polynomial prime)
              (is (p:polynomial= polynomial
                                 (p:modulo (p:scale (ratsimp factors) c) prime)))
              (dolist (factor factors)
                (is-true (p:irreduciblep (cdr factor) prime))))))))

(test lifting
  (loop with state = (make-random-state t)
        repeat 400000 do
        (let* (;; Generate a random primitive polynomial
               (poly (p:remove-content (random-poly 20 state 10)))
               ;; Make sure that the leading coefficient is > 0
               (poly (if (< (p:leading-coeff poly) 0) (p:negate poly) poly)))
          (unless (p:polynomial= poly p:+zero+)
            (let* ((prime (si:consume-one (p:suitable-primes poly)))
                   ;; Constant multiplier in the factorization in
                   ;; 𝔽_p[x] can be ignored.
                   (factors (p:factor (p:modulo poly prime) prime)))
              ;; For simplicity choose a situation with only 2 factors
              (when (and (= (length factors) 2))
                (destructuring-bind ((m1 . f1) (m2 . f2)) factors
                  (when (= m1 m2 1)
                    (multiple-value-bind (f1zx f2zx convp steps)
                        (p:lift-factors poly f1 f2 prime (p:suitable-bound poly))
                      (declare (ignore steps))
                      ;; When a factorization in ℤ[x] exists...
                      (when convp
                        (is (p:polynomial= poly (p:multiply f1zx f2zx)))))))))))))
