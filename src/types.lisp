(in-package :polynomial)

;; A type for matrices
(deftype matrix () '(simple-array (signed-byte 32) (* *)))

;; A bit more narrow range for primes
(deftype prime () '(integer 2 #.most-positive-fixnum))

;; FIXME: let monomial be just an alias for a tuple?
(deftype monomial () '(cons alex:non-negative-fixnum fixnum))
