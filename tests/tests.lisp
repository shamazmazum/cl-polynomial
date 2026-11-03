(in-package :cl-polynomial/tests)

(defun run-tests ()
  (let ((*random-state* (make-random-state t))
        (*num-trials*   5000)
        (*max-trials*   50000))
    (every #'identity
           (mapcar (lambda (suite)
                     (let ((status (run suite)))
                       (explain! status)
                       (results-status status)))
                   '(number algebra factor fft)))))

(defun coeffs-sorted-p (polynomial)
  (let* ((coeffs (p:polynomial-coeffs polynomial))
         (sorted (sort (copy-seq coeffs) #'> :key #'car)))
    (equalp coeffs sorted)))

(defun gen-poly (&key coeff degree)
  (lambda ()
    (let ((degree (funcall degree)))
      (p:list->polynomial
       (loop repeat (1+ degree)
             collect (funcall coeff))))))

(defun gen-coeff/ff (prime)
  (if (= prime 2)
      (gen-integer :min 0 :max 1)
      (let ((p/2 (floor prime 2)))
        (gen-integer :min (- p/2) :max (+ p/2)))))

(defun gen-prime ()
  (gen-one-element 2 3 5 7 11 17 19))

(defun gen-points (&key n x y)
  (lambda ()
    (fpx:list->points
     (remove-duplicates
      (loop repeat (funcall n) collect (cons (funcall x) (funcall y)))
      :key #'car :test #'=))))

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
(def-suite fft     :description "Fast Fourier transform")

(in-suite number)

(test factor-numbers
  (for-all ((number (gen-integer :min 0 :max 100000)))
  (let ((factors (z:factor number)))
    (is (= number (reduce #'* factors)))
    (is-true
     (every
      (lambda (number)
        (let ((factors (z:factor number)))
          (and factors (null (cdr factors)))))
      factors)))))

(test totient
  (for-all ((n1 (gen-integer :min 0 :max 500))
            (n2 (gen-integer :min 0 :max 500)))
    (when (= (gcd n1 n2) 1)
      (is (= (* (z:totient n1)
                (z:totient n2))
             (z:totient (* n1 n2)))))
    (is (= (loop for i from 1 to n1
                 sum (if (zerop (rem n1 i)) (z:totient i) 0))
           n1))))

(test moebius
  (for-all ((n1 (gen-integer :min 1 :max 500))
            (n2 (gen-integer :min 1 :max 500)))
    (when (= (gcd n1 n2) 1)
      (is (= (* (z:moebius n1)
                (z:moebius n2))
             (z:moebius (* n1 n2)))))
    (when (> n1 1)
      (is (= (loop for i from 1 to n1
                   sum (if (zerop (rem n1 i)) (z:moebius i) 0))
             0)))))

(test prime
  (for-all ((n (gen-integer :min 0 :max 100000000000)))
    (let ((factors (z:factor n)))
      (if (z:maybe-prime-p n)
          (is (null (cdr factors)))
          (is       (cdr factors))))))

(in-suite algebra)

(test evaluation
  (flet ((simple-eval (p x)
           (loop for (deg . coeff) in (p:polynomial-coeffs p)
                 sum (* coeff (expt x deg)))))
    (for-all ((poly (gen-poly :coeff  (gen-integer :min -10 :max 10)
                              :degree (gen-integer :min   0 :max 50)))
              (x    (gen-integer :min 0 :max 100)))
      (is (= (p:evaluate poly x) (simple-eval poly x))))))

(test interpolation
  (for-all ((points (gen-points  :n (gen-integer :min   1 :max 20)
                                 :x (gen-integer :min -20 :max 20)
                                 :y (gen-integer :min -20 :max 20))))
    (let ((f (fpx:interpolate points 211)))
      (is-true
       (every
        (lambda (point)
          (fpx:bind-point (x y) point
            (= y (fpx:evaluate f x 211))))
        (fpx:points-list points))))))

(test reciprocal
  (for-all ((p (gen-poly :coeff  (gen-integer :min -100 :max 100)
                         :degree (gen-integer :min    0 :max 100))))
    (let* ((p1 (p:add (p:multiply p p:+variable+) p:+one+))
           (p2 (p:reciprocal p1)))
      (is-true (coeffs-sorted-p p2))
      (is (= (p:degree p1) (p:degree p2)))
      (is (= (p:evaluate p1 1)
             (p:evaluate p2 1))))))

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
  (let ((gen (gen-poly :coeff  (gen-integer :min -100 :max 100)
                       :degree (gen-integer :min    0 :max 100))))
    (for-all ((p1 gen)
              (p2 gen))
      (is-true  (p:polynomial=  p1 p1))
      (is-false (p:polynomial/= p1 p1))
      (is-true  (p:polynomial=  p1 (p:polynomial (p:polynomial-coeffs p1))))
      (is-false (p:polynomial/= p1 (p:polynomial (p:polynomial-coeffs p1))))
      (when (equalp p1 p2)
        (is-true  (p:polynomial=  p1 p2))
        (is-false (p:polynomial/= p1 p2)))
      (when (not (equalp p1 p2))
        (is-true  (p:polynomial/= p1 p2))
        (is-false (p:polynomial=  p1 p2))))))

(test addition/subtraction
  (let ((gen (gen-poly :coeff  (gen-integer :min -100 :max 100)
                       :degree (gen-integer :min    0 :max 100))))
    (for-all ((p1 gen)
              (p2 gen))
      (let ((sum (p:add p1 p2)))
        (is-true (coeffs-sorted-p sum))
        (is (p:polynomial= sum (p:add p2 p1)))
        (is (p:polynomial= (p:add p1 p:+zero+) p1))
        (is (<= (p:degree sum)
                (max (p:degree p1)
                     (p:degree p2))))
        (is (p:polynomial= p1 (p:subtract sum p2)))
        (is (p:polynomial= p2 (p:subtract sum p1)))))))

(test multiplication
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100)))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))))
    (let ((f (fpx:modulo (p:multiply f1 f2) p)))
      (is-true (coeffs-sorted-p f))
      (is (p:polynomial= (p:multiply f1 p:+zero+) p:+zero+))
      (is (p:polynomial= (p:multiply f1 p:+one+)  f1))
      (is (p:polynomial= f (fpx:modulo (p:multiply f2 f1) p))))))

(test multiplication/nz
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f1 p:+zero+))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f2 p:+zero+)))
    (let ((f (fpx:modulo (p:multiply f1 f2) p)))
      (is (= (p:degree f)
             (+ (p:degree f1)
                (p:degree f2))))
      (multiple-value-bind (q r)
          (fpx:divide f f1 p)
        (is (p:polynomial= q f2))
        (is (p:polynomial= r p:+zero+))))))

(test division
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100)))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f2 p:+zero+)))
    (multiple-value-bind (q r)
        (fpx:divide f1 f2 p)
      (is-true (coeffs-sorted-p q))
      (is-true (coeffs-sorted-p r))
      (is-true (or (< (p:degree r) (p:degree f2))
                   (= 0 (p:degree r) (p:degree f2))))
      (is (p:polynomial= f1 (fpx:modulo (p:add (p:multiply q f2) r) p))))))

(test monic
  (for-all* ((p (gen-prime))
             (f (gen-poly :coeff  (gen-coeff/ff p)
                          :degree (gen-integer :min 0 :max 100))
                (p:polynomial/= f p:+zero+)))
    (multiple-value-bind (m c)
        (fpx:monic-polynomial f p)
      (is-true (coeffs-sorted-p m))
      (is (= (p:leading-coeff m) 1))
      (is (p:polynomial= f (fpx:modulo (p:scale m c) p))))))

(test gcd
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f1 p:+zero+))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f2 p:+zero+)))
    (let* ((gcd (fpx:gcd f1 f2 p))
           (lcm (fpx:divide (fpx:modulo (p:multiply f1 f2) p) gcd p)))
      (is (p:polynomial= gcd (fpx:gcd f2 f1 p)))
      (is-true (p:monicp gcd))
      (unless (p:constantp gcd)
        (is (p:polynomial= (fpx:remainder f1 gcd p) p:+zero+))
        (is (p:polynomial= (fpx:remainder f2 gcd p) p:+zero+))
        (is (p:polynomial= (fpx:remainder lcm f1 p) p:+zero+))
        (is (p:polynomial= (fpx:remainder lcm f2 p) p:+zero+))))))

(test gcd2
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100)))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100)))
             (f3 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))
                 (p:polynomial/= f3 p:+zero+)))
    (let* ((f1 (fpx:modulo (p:multiply f1 f3) p))
           (f2 (fpx:modulo (p:multiply f2 f3) p))
           (gcd (fpx:gcd f1 f2 p)))
      (is (p:polynomial= (fpx:remainder gcd f3 p) p:+zero+)))))

(test gcdex
  (for-all* ((p  (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100)))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 100))))
    (multiple-value-bind (gcd a b)
        (fpx:gcdex f1 f2 p)
      (is (p:polynomial= gcd (fpx:gcd f1 f2 p)))
      (is (p:polynomial= (fpx:modulo (p:add (p:multiply f1 a) (p:multiply f2 b)) p) gcd)))))

(test gcd-zx
  (for-all ((f1 (gen-poly :coeff  (gen-integer :min -5 :max 5)
                          :degree (gen-integer :min 0 :max 20))
                (p:polynomial/= f1 p:+zero+))
            (f2 (gen-poly :coeff  (gen-integer :min -5 :max 5)
                          :degree (gen-integer :min 0 :max 20))
                (p:polynomial/= f2 p:+zero+)))
    (let* ((gcd (zx:gcd f1 f2))
           (lcm (zx:divide (p:multiply f1 f2) gcd)))
      (is (p:polynomial= gcd (zx:gcd f2 f1)))
      (unless (p:constantp gcd)
        (is (p:polynomial= (zx:remainder f1 gcd) p:+zero+))
        (is (p:polynomial= (zx:remainder f2 gcd) p:+zero+))
        (is (p:polynomial= (zx:remainder lcm f1) p:+zero+))
        (is (p:polynomial= (zx:remainder lcm f2) p:+zero+))))))

(test gcd-zx2
  (for-all ((f1 (gen-poly :coeff  (gen-integer :min -100 :max 100)
                          :degree (gen-integer :min    0 :max 100)))
            (f2 (gen-poly :coeff  (gen-integer :min -100 :max 100)
                          :degree (gen-integer :min    0 :max 100)))
            (f3 (gen-poly :coeff  (gen-integer :min -100 :max 100)
                          :degree (gen-integer :min    0 :max 100))
                (p:polynomial/= f3 p:+zero+)))
    (let* ((f1 (p:multiply f1 f3))
           (f2 (p:multiply f2 f3))
           (gcd (zx:gcd f1 f2)))
      (is (p:polynomial= (zx:remainder gcd f3) p:+zero+)))))

(test expt
  (for-all ((f (gen-poly :coeff  (gen-integer :min -100 :max 100)
                         :degree (gen-integer :min    0 :max 10)))
            (n (gen-integer :min 0 :max 10)))
    (is (p:polynomial= (p:expt f n)
                       (labels ((%go (acc m)
                                  (if (zerop m) acc
                                      (%go (p:multiply acc f)
                                           (1- m)))))
                         (%go p:+one+ n))))))

(in-suite factor)

(test square-free
  (for-all* ((p (gen-prime))
             (f (gen-poly :coeff  (gen-coeff/ff p)
                          :degree (gen-integer :min 0 :max 20))
                (p:polynomial/= f p:+zero+)))
    (let ((f (fpx:monic-polynomial f p)))
      (is (p:polynomial= (fpx:modulo (ratsimp (fpx:square-free f p)) p) f)))))

(test distinct-degree
  (for-all* ((p (gen-prime))
             (f (gen-poly :coeff  (gen-coeff/ff p)
                          :degree (gen-integer :min 0 :max 20))
                (> (p:degree f) 0)))
    (let* ((f (make-square-free f p))
           (factors (fpx:distinct-degree f p))
           (rps (fpx:reducing-polynomials f p)))
      (dolist (factor factors)
        (destructuring-bind (deg . f) factor
          (is (not (zerop (p:degree f))))
          (is (zerop (rem (p:degree f) deg)))))
      (is (= (length rps)
             (reduce #'+ factors
                     :key (lambda (factor)
                            (/ (p:degree (cdr factor))
                               (car factor)))))))))

(test reducing-polys
  (for-all* ((p (gen-prime))
             (f (gen-poly :coeff  (gen-coeff/ff p)
                          :degree (gen-integer :min 0 :max 20))
                (> (p:degree f) 0)))
    (let ((f (fpx:monic-polynomial f p)))
      (dolist (rp (fpx:reducing-polynomials f p))
        ;; rp^prime mod polynomial
        (is (p:polynomial= (fpx::expt-mod rp p f p) rp))))))

(test factor-finite
  (for-all* ((p (gen-prime))
             (f (gen-poly :coeff  (gen-coeff/ff p)
                          :degree (gen-integer :min 0 :max 20))
                (> (p:degree f) 0)))
    (multiple-value-bind (factors c)
        (fpx:factor f p :berlekamp)
      (is (p:polynomial= f (fpx:modulo (p:scale (ratsimp factors) c) p)))
      (dolist (factor factors)
        (is-true (fpx:irreduciblep (cdr factor) p)))
      (when (/= p 2)
        (let ((%factors (fpx:factor f p :cantor-zassenhaus)))
          (is (set-equal-p factors %factors)))))))

(test factor-finite2
  (for-all* ((p (gen-prime))
             (f1 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 20))
                 (> (p:degree f1) 0))
             (f2 (gen-poly :coeff  (gen-coeff/ff p)
                           :degree (gen-integer :min 0 :max 20))
                 (> (p:degree f2) 0)))
    (let ((f3 (fpx:modulo (p:multiply f1 f2) p)))
      (dolist (method '(:berlekamp :cantor-zassenhaus))
        (let ((factors1 (fpx:factor f1 p method))
              (factors2 (fpx:factor f2 p method))
              (factors3 (fpx:factor f3 p method)))
          (is-true (factor-subset-p factors1 factors3))
          (is-true (factor-subset-p factors2 factors3))
          (is (> (reduce #'+ factors3 :key #'car) 1)))))))

(test lifting
  (for-all ((f (gen-poly :coeff  (gen-integer :min -20 :max 20)
                         :degree (gen-integer :min   0 :max 9))))
    ;; Generate a random square-free monic polynomial
    (let* ((f (make-square-free-zx (p:add f (p:polynomial '((10 . 1))))))
           (p (zx:suitable-prime f))
           (fpx-factors (mapcar #'cdr (fpx:factor (fpx:modulo f p) p)))
           (lifting-steps (nth-value 1 (zx:lifting-steps (zx:suitable-bound f) p)))
           (lifted-factors (zx:lift-factors f fpx-factors p lifting-steps)))
      (is (p:polynomial= f (fpx:modulo (apply #'p:multiply lifted-factors)
                                       (expt p (expt 2 lifting-steps))))))))


(test square-free-zx
  ;; Generate a polynomial of relatively low degree to increase the
  ;; probability of not being irreducible.
  (for-all ((f (gen-poly :coeff  (gen-integer :min -5 :max 5)
                         :degree (gen-integer :min  0 :max 20))
               (p:polynomial/= f p:+zero+)))
    (multiple-value-bind (factors c)
        (zx:square-free f)
      (is (p:polynomial= (p:scale (ratsimp factors) c) f)))))

(test factor-zx
  ;; Generate a polynomial of relatively low degree to increase
  ;; the probability of not being irreducible.
  (for-all ((f (gen-poly :coeff  (gen-integer :min -5 :max 5)
                         :degree (gen-integer :min  0 :max 20))
               (p:polynomial/= f p:+zero+)))
    (multiple-value-bind (factors c)
        (zx:factor f)
      (is (p:polynomial= (p:scale (ratsimp factors) c) f))
      (is-true (every (alex:compose #'zx:irreduciblep #'cdr) factors)))))

(test factor-zx2
  (for-all ((f1 (gen-poly :coeff  (gen-integer :min -20 :max 20)
                          :degree (gen-integer :min   0 :max 20))
                (> (p:degree f1) 0))
            (f2 (gen-poly :coeff  (gen-integer :min -20 :max 20)
                          :degree (gen-integer :min    0 :max 20))
                (> (p:degree f2) 0)))
    (let* ((f3 (p:multiply f1 f2))
           (factors1 (zx:factor f1))
           (factors2 (zx:factor f2))
           (factors3 (zx:factor f3)))
      (is-true (factor-subset-p factors1 factors3))
      (is-true (factor-subset-p factors2 factors3))
      (is (> (reduce #'+ factors3 :key #'car) 1)))))

(test factor-cyclotomic
  (loop for n from 1 to 100 do
        (let ((factors1 (factor-x^n-1 n))
              (factors2 (zx:factor (p:polynomial (list (cons n 1) '(0 . -1))))))
          (is-true (every (lambda (f) (= (car f) 1)) factors2))
          (is (set-equal-p factors1 (mapcar #'cdr factors2))))))

(in-suite fft)

(test multiplicative-group
  (si:do-iterator (p (si:take 1000 z:*prime-source*))
    (let ((g (fft:*-group-generator p)))
      (is (= (length
              (remove-duplicates
               (si:collect
                   (si:take (1- p) (si:iterate (lambda (x) (u:mod-sym (* x g) p)) 1)))
               :test #'=))
             (1- p))))))

(test primitive-root-of-unity
  (si:do-iterator (r (si:range 1 12))
    (let* ((n (expt 2 r))
           (p (si:consume-one (fft:fourier-primes n)))
           (g (fft:primitive-root-of-unity n p)))
      (is (= (u:mod-sym (expt g n) p) 1))
      (is-false (si:some (lambda (m) (= (u:mod-sym (expt g m) p) 1))
                         (si:range 1 n))))))


(test ifft-fft-identity
  (let* ((n 65536)
         (p 31916033)
         (b (floor p 2))
         (ω (fft:primitive-root-of-unity n p)))
    (for-all ((a (gen-buffer :length       (gen-integer :min n     :max n)
                             :elements     (gen-integer :min (- b) :max (+ b))
                             :element-type 'fixnum)))
      (is-true (every #'= a (fft:ifft (fft:fft a ω p) ω p))))))

(test fft-polynomial-multiplication
  (let* ((p 51713)
         (n 512)
         (ω (fft:primitive-root-of-unity n p)))
    (for-all ((a (gen-buffer :length       (gen-integer :min (/ n 2) :max (/ n 2))
                             :elements     (gen-integer :min -10     :max 10)
                             :element-type 'fixnum))
              (b (gen-buffer :length       (gen-integer :min (/ n 2) :max (/ n 2))
                             :elements     (gen-integer :min -10     :max 10)
                             :element-type 'fixnum)))
      (let* ((ft1 (fft:fft (fft:pad-array a n) ω p))
             (ft2 (fft:fft (fft:pad-array b n) ω p))
             (mul (map '(vector fixnum) #'* ft1 ft2))
             (c (fft:ifft mul ω p)))
        (is (p:polynomial= (p:sequence->polynomial c)
                           (p:multiply
                            (p:sequence->polynomial a)
                            (p:sequence->polynomial b))))))))

(test fft-values
  (let* ((p 1009)
         (n 16)
         (b (floor p 2))
         (ω (fft:primitive-root-of-unity n p)))
    (for-all ((a (gen-buffer :length       (gen-integer :min n     :max n)
                             :elements     (gen-integer :min (- b) :max (+ b))
                             :element-type 'fixnum)))
        (let* ((f (p:sequence->polynomial a))
               (fft (fft:fft a ω p))
               (slow-dft (map '(vector integer)
                              (lambda (k)
                                (u:mod-sym (p:evaluate f (expt ω k)) p))
                              (loop for k below n collect k))))
          (is-true (every #'= fft slow-dft))))))

(test inplace-fft
  (let* ((p 31916033)
         (n 65536)
         (b (floor p 2))
         (ω (fft:primitive-root-of-unity n p)))
  (for-all ((a (gen-buffer :length       (gen-integer :min n     :max n)
                           :elements     (gen-integer :min (- b) :max (+ b))
                           :element-type 'fixnum)))
    (is-true (every #'= (fft:fft a ω p) (fft:fft! (copy-seq a) ω p)))
    (is-true (every #'= (fft:ifft a ω p) (fft:ifft! (copy-seq a) ω p))))))
