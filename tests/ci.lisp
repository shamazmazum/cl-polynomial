(defun do-all()
  (ql:quickload :cl-polynomial/tests)
  (uiop:quit
   (if (uiop:call-function "cl-polynomial/tests:run-tests")
       0 1)))

(do-all)
