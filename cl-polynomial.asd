(defsystem :cl-polynomial
    :name :cl-polynomial
    :version "0.1"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :description "Factorization of polynomials over integers and finite fields"
    :license "2-clause BSD"
    :serial t
    :pathname "src/"
    :components ((:file "package")
                 (:file "types")
                 (:file "linalg")
                 (:file "polynomial")
                 (:file "factor")
                 (:file "factor-zx"))
    :depends-on (:alexandria :serapeum :stateless-iterators :cl-prime-maker)
    :in-order-to ((test-op (load-op "cl-polynomial/tests")))
    :perform (test-op (op system)
                      (declare (ignore op system))
                      (uiop:symbol-call :cl-polynomial/tests '#:run-tests)))

(defsystem :cl-polynomial/tests
    :name :cl-polynomial/tests
    :version "0.1"
    :author "Vasily Postnicov <shamaz.mazum@gmail.com>"
    :license "2-clause BSD"
    :pathname "tests/"
    :components ((:file "package")
                 (:file "tests" :depends-on ("package")))
    :depends-on (:cl-polynomial :fiveam))
