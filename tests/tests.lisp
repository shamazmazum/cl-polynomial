(in-package :cl-polynomial/tests)

(alex:define-constant +some-primes+ '(2 3 5 7 11 17 19)
  :test #'equalp)

(defun run-tests ()
  (every #'identity
         (mapcar (lambda (suite)
                   (let ((status (run suite)))
                     (explain! status)
                     (results-status status)))
                 '(number algebra factor))))

(defun coeffs-sorted-p (polynomial)
  (let* ((coeffs (p:polynomial-coeffs polynomial))
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

(defun make-square-free (polynomial prime)
  (ratsimp
   (mapcar
    (lambda (factor)
      (destructuring-bind (m . f) factor
        (declare (ignore m))
        (cons 1 f)))
    (fpx:square-free polynomial prime))))

(defun make-square-free-zx (polynomial)
  (ratsimp
   (mapcar
    (lambda (factor)
      (destructuring-bind (m . f) factor
        (declare (ignore m))
        (cons 1 f)))
    (zx:square-free polynomial))))

(defun set-equal-p (s1 s2)
  (and (subsetp s1 s2 :test #'equalp)
       (subsetp s2 s1 :test #'equalp)))

(defun factor-subset-p (small big)
  (every
   (lambda (small-factor)
     (let ((big-factor (find (cdr small-factor) big :key #'cdr :test #'p:polynomial=)))
       (and big-factor (<= (car small-factor) (car big-factor)))))
   small))

(defun factor-x^n-1 (n)
  (si:foldl
   (lambda (acc m)
     (if (zerop (rem n m)) (cons (zx:cyclotomic m) acc) acc))
   nil (si:range 1 (1+ n))))

(def-suite number  :description "Basic number-theoretic functions")
(def-suite algebra :description "Generic algebraic operations on polynomials")
(def-suite factor  :description "Factorization of polynomials over finite fields")

(in-suite number)

(test factor-numbers
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((number (random 100000))
               (factors (z:factor number)))
          (is (= number (reduce #'* factors)))
          (is-true
           (every
            (lambda (number)
              (= (length (z:factor number)) 1))
            factors)))))

(test totient
  (loop with state = (make-random-state t)
        repeat 1000 do
        (let ((n1 (random 500))
              (n2 (random 500)))
          (when (= (gcd n1 n2) 1)
            (is (= (* (z:totient n1) (z:totient n2))
                   (z:totient (* n1 n2)))))
          (is (= (loop for i from 1 to n1 sum
                       (if (zerop (rem n1 i)) (z:totient i) 0))
                 n1)))))

(test moebius
  (loop with state = (make-random-state t)
        repeat 1000 do
        (let ((n1 (random 500))
              (n2 (random 500)))
          (when (= (gcd n1 n2) 1)
            (is (= (* (z:moebius n1) (z:moebius n2))
                   (z:moebius (* n1 n2)))))
          (when (> n1 1)
            (is (= (loop for i from 1 to n1 sum
                         (if (zerop (rem n1 i)) (z:moebius i) 0))
                   0))))))

(in-suite algebra)

(test mod-sym-homomorphism
  (loop for p in '(2 3 5 7) do
        (loop for n in '(1 2 3 4) do
              (let ((x1 (1+ (random 1000)))
                    (x2 (1+ (random 1000)))
                    (q (expt p n)))
                (is (= (u:mod-sym 0 q) 0))
                (is (= (u:mod-sym 1 q) 1))
                (is (= (u:mod-sym (* x1 x2) q)
                       (u:mod-sym
                        (* (u:mod-sym x1 q)
                           (u:mod-sym x2 q))
                        q)))
                (is (= (u:mod-sym (+ x1 x2) q)
                       (u:mod-sym
                        (+ (u:mod-sym x1 q)
                           (u:mod-sym x2 q))
                        q)))))))

(test equality
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let ((poly1 (random-poly 100 state (1+ (random 100))))
              (poly2 (random-poly 100 state (1+ (random 100)))))
          (is-true  (p:polynomial=  poly1 poly1))
          (is-false (p:polynomial/= poly1 poly1))
          (is-true  (p:polynomial= poly1 (p:polynomial (p:polynomial-coeffs poly1))))
          (is-false (p:polynomial/= poly1 (p:polynomial (p:polynomial-coeffs poly1))))
          (when (equalp poly1 poly2)
            (is-true  (p:polynomial= poly1 poly2))
            (is-false (p:polynomial/= poly1 poly2)))
          (when (not (equalp poly1 poly2))
            (is-true  (p:polynomial/= poly1 poly2))
            (is-false (p:polynomial= poly1 poly2))))))

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
               (product (fpx:modulo (p:multiply poly1 poly2) prime)))
          (is-true (coeffs-sorted-p product))
          (is (p:polynomial= (p:multiply poly1 p:+zero+) p:+zero+))
          (is (p:polynomial= (p:multiply poly1 p:+one+) poly1))
          (is (p:polynomial= product (fpx:modulo (p:multiply poly2 poly1) prime)))
          (unless (or (p:polynomial= poly1 p:+zero+)
                      (p:polynomial= poly2 p:+zero+))
            (is (= (p:degree product)
                   (+ (p:degree poly1)
                      (p:degree poly2)))))
          (unless (p:polynomial= poly1 p:+zero+)
            (multiple-value-bind (q r)
                (fpx:divide product poly1 prime)
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
                (fpx:divide poly1 poly2 prime)
              (is-true (coeffs-sorted-p q))
              (is-true (coeffs-sorted-p r))
              (is-true (or (< (p:degree r) (p:degree poly2))
                           (= 0 (p:degree r) (p:degree poly2))))
              (is (p:polynomial=
                   poly1
                   (fpx:modulo (p:add (p:multiply q poly2) r) prime))))))))

(test monic
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly (random-poly prime state)))
          (unless (p:polynomial= poly p:+zero+)
            (multiple-value-bind (monic c)
                (fpx:monic-polynomial poly prime)
              (is-true (coeffs-sorted-p monic))
              (is (= (p:leading-coeff monic) 1))
              (is (p:polynomial=
                   poly
                   (fpx:modulo (p:scale monic c) prime))))))))

(test gcd
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state)))
          (unless (or (p:polynomial= poly1 p:+zero+)
                      (p:polynomial= poly2 p:+zero+))
            (let* ((gcd (fpx:gcd poly1 poly2 prime))
                   (lcm (fpx:divide (fpx:modulo (p:multiply poly1 poly2) prime) gcd prime)))
              (is (p:polynomial= gcd (fpx:gcd poly2 poly1 prime)))
              (is-true (p:monicp gcd))
              (unless (p:constantp gcd)
                (is (p:polynomial= (fpx:remainder poly1 gcd prime) p:+zero+))
                (is (p:polynomial= (fpx:remainder poly2 gcd prime) p:+zero+))
                (is (p:polynomial= (fpx:remainder lcm poly1 prime) p:+zero+))
                (is (p:polynomial= (fpx:remainder lcm poly2 prime) p:+zero+))))))))

(test gcd2
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state))
               (poly3 (random-poly prime state)))
          (unless (p:polynomial= poly3 p:+zero+)
            (let* ((poly1 (fpx:modulo (p:multiply poly1 poly3) prime))
                   (poly2 (fpx:modulo (p:multiply poly2 poly3) prime))
                   (gcd (fpx:gcd poly1 poly2 prime)))
              (is (p:polynomial= (fpx:remainder gcd poly3 prime) p:+zero+))
              (is (p:polynomial= (fpx:remainder gcd poly3 prime) p:+zero+)))))))

(test gcdex
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state))
               (poly2 (random-poly prime state)))
          (multiple-value-bind (gcd a b)
              (fpx:gcdex poly1 poly2 prime)
            (is (p:polynomial= gcd (fpx:gcd poly1 poly2 prime)))
            (is (p:polynomial= (fpx:modulo
                                (p:add (p:multiply poly1 a)
                                       (p:multiply poly2 b))
                                prime)
                               gcd))))))

(test gcd-zx
  (loop with state = (make-random-state t)
        repeat 300000 do
        (let ((poly1 (random-poly 20 state 20))
              (poly2 (random-poly 20 state 20)))
          (unless (or (p:polynomial= poly1 p:+zero+)
                      (p:polynomial= poly2 p:+zero+))
            (let* ((gcd (zx:gcd poly1 poly2))
                   (lcm (zx:divide (p:multiply poly1 poly2) gcd)))
              (is (p:polynomial= gcd (zx:gcd poly2 poly1)))
              (unless (p:constantp gcd)
                (is (p:polynomial= (zx:remainder poly1 gcd) p:+zero+))
                (is (p:polynomial= (zx:remainder poly2 gcd) p:+zero+))
                (is (p:polynomial= (zx:remainder lcm poly1) p:+zero+))
                (is (p:polynomial= (zx:remainder lcm poly2) p:+zero+))))))))

(test gcd-zx2
  (loop with state = (make-random-state t)
        repeat 50000 do
        (let ((poly1 (random-poly 20 state 20))
              (poly2 (random-poly 20 state 20))
              (poly3 (random-poly 20 state 20)))
          (unless (p:polynomial= poly3 p:+zero+)
            (let* ((poly1 (p:multiply poly1 poly3))
                   (poly2 (p:multiply poly2 poly3))
                   (gcd (zx:gcd poly1 poly2)))
              (is (p:polynomial= (zx:remainder gcd poly3) p:+zero+))
              (is (p:polynomial= (zx:remainder gcd poly3) p:+zero+)))))))

(test expt
  (loop with state = (make-random-state t)
        repeat 10000 do
        (let ((p (random-poly 20 state 10))
              (n (random 10)))
          (unless (p:constantp p)
            (is-true (every (lambda (x) (>= (car x) n))
                            (zx:factor (p:expt p n))))))))

(in-suite factor)

(test square-free
  (loop with state = (make-random-state t)
        repeat 40000 do
        (let* ((prime (random-prime state))
               (polynomial (fpx:monic-polynomial (random-poly prime state 20) prime)))
          (unless (p:polynomial= polynomial p:+zero+)
            (is (p:polynomial=
                 (fpx:modulo (ratsimp (fpx:square-free polynomial prime)) prime)
                 polynomial))))))

(test distinct-degree
  (loop with state = (make-random-state t)
        repeat 40000 do
        (let* ((prime (random-prime state))
               (poly (random-poly prime state 20)))
          (unless (p:constantp poly)
            (let* ((poly (make-square-free poly prime))
                   (factors (fpx:distinct-degree poly prime))
                   (rps (fpx:reducing-polynomials poly prime)))
              (dolist (factor factors)
                (destructuring-bind (deg . f) factor
                  (is (not (zerop (p:degree f))))
                  (is (zerop (rem (p:degree f) deg)))))
              (is (= (length rps)
                     (reduce #'+ factors
                             :key (lambda (factor)
                                    (/ (p:degree (cdr factor))
                                       (car factor)))))))))))

(test reducing-polys
  (loop with state = (make-random-state t)
        repeat 40000 do
        (let* ((prime (random-prime state))
               (polynomial (fpx:monic-polynomial (random-poly prime state 20) prime)))
          (unless (zerop (p:degree polynomial))
            (dolist (rp (fpx:reducing-polynomials polynomial prime))
              (is (p:polynomial=
                   ;; rp^prime mod polynomial
                   (fpx:modulo (fpx::expt-rem rp prime polynomial prime) prime)
                   rp)))))))

(test factor-finite
  (loop with state = (make-random-state t)
        repeat 40000 do
        (let* ((prime (random-prime state))
               (polynomial (fpx:monic-polynomial (random-poly prime state 20) prime)))
          (unless (zerop (p:degree polynomial))
            (multiple-value-bind (factors c)
                (fpx:factor polynomial prime :berlekamp)
              (is (p:polynomial= polynomial (fpx:modulo (p:scale (ratsimp factors) c) prime)))
              (dolist (factor factors)
                (is-true (fpx:irreduciblep (cdr factor) prime)))
              (when (/= prime 2)
                (let ((%factors (fpx:factor polynomial prime :cantor-zassenhaus)))
                  (is (set-equal-p factors %factors)))))))))

(test factor-finite2
  (loop with state = (make-random-state t)
        repeat 40000 do
        (let* ((prime (random-prime state))
               (poly1 (random-poly prime state 20))
               (poly2 (random-poly prime state 20))
               (poly3 (fpx:modulo (p:multiply poly1 poly2) prime)))
          (unless (or (p:polynomial= poly3 p:+zero+)
                      (= (p:degree poly1) 0)
                      (= (p:degree poly2) 0))
            (dolist (method '(:berlekamp :cantor-zassenhaus))
              (let ((factors1 (fpx:factor poly1 prime method))
                    (factors2 (fpx:factor poly2 prime method))
                    (factors3 (fpx:factor poly3 prime method)))
                (is-true (factor-subset-p factors1 factors3))
                (is-true (factor-subset-p factors2 factors3))
                (is (> (reduce #'+ factors3 :key #'car) 1))))))))

(test lifting
  (loop with state = (make-random-state t)
        repeat 100000 do
        ;; Generate a random square-free monic polynomial
        (let* ((poly (make-square-free-zx
                      (p:add (random-poly 20 state 10) (p:polynomial '((10 . 1))))))
               (prime (zx:suitable-prime poly))
               (fpx-factors (mapcar #'cdr (fpx:factor (fpx:modulo poly prime) prime)))
               (lifting-steps (nth-value 1 (zx:suitable-bound poly prime)))
               (lifted-factors (zx:lift-factors poly fpx-factors prime lifting-steps)))
          (is (p:polynomial= poly (fpx:modulo (apply #'p:multiply lifted-factors)
                                              (expt prime (expt 2 lifting-steps))))))))

(test square-free-zx
  (loop with state = (make-random-state t)
        repeat 100000 do
        ;; Generate a polynomial of relatively low degree to increase
        ;; the probability of not being irreducible.
        (let ((polynomial (random-poly 10 state 10)))
          (unless (p:polynomial= polynomial p:+zero+)
            (multiple-value-bind (factors c)
                (zx:square-free polynomial)
              (is (p:polynomial=
                   (p:scale (ratsimp factors) c)
                   polynomial)))))))

(test factor-zx
  (loop with state = (make-random-state t)
        repeat 100000 do
        ;; Generate a polynomial of relatively low degree to increase
        ;; the probability of not being irreducible.
        (let ((polynomial (random-poly 10 state 10)))
          (unless (p:polynomial= polynomial p:+zero+)
            (multiple-value-bind (factors c)
                (zx:factor polynomial)
              (is (p:polynomial=
                   (p:scale (ratsimp factors) c)
                   polynomial))
              (is-true (every (alex:compose #'zx:irreduciblep #'cdr) factors)))))))

(test factor-zx2
  (loop with state = (make-random-state t)
        repeat 100000 do
        (let* ((poly1 (random-poly 20 state 20))
               (poly2 (random-poly 20 state 20))
               (poly3 (p:multiply poly1 poly2)))
          (unless (or (p:polynomial= poly3 p:+zero+)
                      (= (p:degree poly1) 0)
                      (= (p:degree poly2) 0))
            (let ((factors1 (zx:factor poly1))
                  (factors2 (zx:factor poly2))
                  (factors3 (zx:factor poly3)))
              (is-true (factor-subset-p factors1 factors3))
              (is-true (factor-subset-p factors2 factors3))
              (is (> (reduce #'+ factors3 :key #'car) 1)))))))

(test factor-cyclotomic
  (loop for n from 1 to 100 do
        (let ((factors1 (factor-x^n-1 n))
              (factors2 (zx:factor (p:polynomial (list (cons n 1) '(0 . -1))))))
          (is-true (every (lambda (f) (= (car f) 1)) factors2))
          (is (set-equal-p factors1 (mapcar #'cdr factors2))))))
