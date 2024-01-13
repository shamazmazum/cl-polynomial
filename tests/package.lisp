(defpackage cl-polynomial/tests
  (:use #:cl #:fiveam)
  (:local-nicknames (#:p    #:polynomial)
                    (#:alex #:alexandria)
                    (#:si   #:stateless-iterators))
  (:export #:run-tests))
