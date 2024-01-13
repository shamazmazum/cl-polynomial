(defpackage cl-polynomial/tests
  (:use #:cl #:fiveam)
  (:local-nicknames (#:p    #:cl-polynomial/polynomial)
                    (#:ff   #:cl-polynomial/factor)
                    (#:fi   #:cl-polynomial/factor-zx)
                    (#:alex #:alexandria)
                    (#:si   #:stateless-iterators))
  (:export #:run-tests))
