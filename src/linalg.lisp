(in-package :polynomial)

(sera:-> row (matrix alex:non-negative-fixnum
                     &optional
                     alex:non-negative-fixnum
                     alex:non-negative-fixnum)
         (values row &optional))
(defun row (matrix k &optional (start 0) (end (array-dimension matrix 1)))
  (declare (optimize (speed 3)))
  (loop with row = (make-array (- end start) :element-type '(signed-byte 32))
        for i from start below end
        for j from 0 by 1 do
        (setf (aref row j)
              (aref matrix k i))
        finally (return row)))

(sera:-> (setf row) (row matrix alex:non-negative-fixnum)
         (values row &optional))
(defun (setf row) (row matrix k)
  (declare (optimize (speed 3)))
  (dotimes (i (array-dimension matrix 1))
    (setf (aref matrix k i) (aref row i)))
  row)

(sera:-> triangularize! (matrix prime)
         (values matrix &optional))
(defun triangularize! (matrix p)
  "Bring a matrix to an echelon form using row operations"
  (declare (optimize (speed 3)))
  (destructuring-bind (n m)
      (array-dimensions matrix)
    (loop for i below n do
          ;; Find a row with the left-most non-zero leading element in
          ;; rows from i to n.
          (multiple-value-bind (leading-row leading-column)
              (loop with row-idx fixnum = n
                    with col-idx fixnum = m
                    for j from i below n
                    for row = (row matrix j)
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
             (row matrix i)
             (row matrix leading-row))

            ;; Now cancel the left-most leading element in all rows from i+1
            (loop with rowi = (row matrix i)
                  with elt1 = (aref rowi leading-column)
                  for j from (1+ i) below n
                  for rowj = (row matrix j)
                  for elt2 = (aref rowj leading-column)
                  when (not (zerop elt2)) do
                  (let* ((mul (mod (* (invert-integer elt2 p) elt1) p))
                         (new-row (map-into
                                   rowj
                                   (lambda (x1 x2)
                                     (declare (type (signed-byte 32) x1 x2))
                                     (mod (- x1 (* mul x2)) p))
                                   rowi rowj)))
                    (setf (row matrix j) new-row))))))
  matrix)

(sera:-> attach-identity (matrix)
         (values matrix &optional))
(defun attach-identity (matrix)
  "Attach an identity matrix to the right of a square matrix"
  (declare (optimize (speed 3)))
  (destructuring-bind (n m)
      (array-dimensions matrix)
    (declare (type fixnum n m))
    (assert (= n m))
    (loop with result = (make-array (list n (* n 2)) :element-type '(signed-byte 32))
          for i below n do
          (loop for j below (* n 2) do
                (setf (aref result i j)
                      (cond
                        ((< j n) (aref matrix i j))
                        ((= (- j n) i) 1)
                        (t 0))))
          finally (return result))))

(sera:-> %nullspace (matrix prime)
         (values matrix &optional))
(declaim (inline %nullspace))
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
          for r1 = (row nullspace (- n i 1) 0 n)
          for r2 = (row nullspace (- n i 1) n (* n 2))
          while (every #'zerop r1) collect r2)))
