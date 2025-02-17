(defsystem :cl-polynomial
    :name :cl-polynomial
    :version "0.2.2"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :description "Factorization of polynomials over integers and finite fields"
    :license "2-clause BSD"
    :serial t
    :pathname "src"
    :class :package-inferred-system
    :depends-on (:alexandria
                 :serapeum
                 :stateless-iterators
                 :cl-polynomial/util
                 :cl-polynomial/linalg
                 :cl-polynomial/polynomial
                 :cl-polynomial/z
                 :cl-polynomial/fp
                 :cl-polynomial/fpx
                 :cl-polynomial/zx)
    :in-order-to ((test-op (load-op "cl-polynomial/tests")))
    :perform (test-op (op system)
                      (declare (ignore op system))
                      (uiop:symbol-call :cl-polynomial/tests '#:run-tests)))

(defsystem :cl-polynomial/tests
    :name :cl-polynomial/tests
    :version "0.2.2"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :license "2-clause BSD"
    :serial t
    :pathname "tests"
    :components ((:file "package")
                 (:file "tests"))
    :depends-on (:fiveam :cl-polynomial))
