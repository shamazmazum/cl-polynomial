(in-package :polynomial)

(sera:-> triangularize! (matrix prime)
         (values matrix &optional))
(defun triangularize! (matrix p)
  "Bring a matrix to an echelon form using row operations"
  (destructuring-bind (n m)
      (array-dimensions matrix)
    (loop for i below n do
          ;; Find a row with the left-most non-zero leading element in
          ;; rows from i to n.
          (multiple-value-bind (leading-row leading-column)
              (loop with row-idx = n
                    with col-idx = m
                    for j from i below n
                    for row = (select:select matrix j (select:range 0 m))
                    for pos = (position-if
                               (lambda (x) (not (zerop x)))
                               row)
                    when (and pos (< pos col-idx)) do
                    (setq col-idx pos
                          row-idx j)
                    finally (return (values row-idx col-idx)))
            ;; KLUDGE:: If there is no such row, then exit
            (when (= leading-row n)
              (return-from triangularize! matrix))

            ;; Swap this row with i-th row
            (rotatef
             (select:select matrix i (select:range 0 m))
             (select:select matrix leading-row (select:range 0 m)))

            ;; Now cancel the left-most leading element in all rows from i+1
            (loop with rowi = (select:select matrix i (select:range 0 m))
                  with elt1 = (aref rowi leading-column)
                  for j from (1+ i) below n
                  for rowj = (select:select matrix j (select:range 0 m))
                  for elt2 = (aref rowj leading-column)
                  when (not (zerop elt2)) do
                  (let* ((mul (mod (* (invert-integer elt2 p) elt1) p))
                         (new-row (map '(vector (signed-byte 32))
                                       (lambda (x1 x2)
                                         (mod (- x1 (* mul x2)) p))
                                       rowi rowj)))
                    (setf (select:select matrix j (select:range 0 m)) new-row))))))
  matrix)

(sera:-> identity-matrix (alex:positive-fixnum)
         (values matrix &optional))
(defun identity-matrix (n)
  "Return NxN identity matrix"
  (let ((matrix (make-array (list n n)
                            :element-type '(signed-byte 32)
                            :initial-element 0)))
    (loop for i below n do
          (setf (aref matrix i i) 1))
    matrix))

(sera:-> attach-identity (matrix)
         (values matrix &optional))
(defun attach-identity (matrix)
  "Attach an identity matrix to the right of a square matrix"
  (destructuring-bind (n m)
      (array-dimensions matrix)
    (assert (= n m))
    (let ((id (identity-matrix n))
          (result (make-array (list n (* n 2))
                              :element-type '(signed-byte 32))))
      (setf (select:select result
              (select:range 0 n)
              (select:range 0 n))
            matrix
            (select:select result
              (select:range 0 n)
              (select:range n (* n 2)))
            id)
      result)))

(sera:-> %nullspace (matrix prime)
         (values matrix &optional))
(defun %nullspace (matrix p)
  (triangularize! (attach-identity matrix) p))

(sera:-> nullspace (matrix prime)
         (values list &optional))
(defun nullspace (m p)
  "Return a list of vectors which span ker(M^T) with elements in
\(\mathbb{F}_p\)."
  (let ((nullspace (%nullspace m p))
        (n (array-dimension m 0)))
    (loop for i below n
          for r1 = (select:select nullspace (- n i 1) (select:range 0 n))
          for r2 = (select:select nullspace (- n i 1) (select:range n (* n 2)))
          while (every #'zerop r1) collect r2)))
