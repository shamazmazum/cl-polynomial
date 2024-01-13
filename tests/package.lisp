(defpackage cl-polynomial/tests
  (:use #:cl #:fiveam)
  (:local-nicknames (#:p    #:cl-polynomial/polynomial)
                    (#:fpx  #:cl-polynomial/fpx)
                    (#:zx   #:cl-polynomial/zx)
                    (#:alex #:alexandria)
                    (#:si   #:stateless-iterators))
  (:export #:run-tests))
