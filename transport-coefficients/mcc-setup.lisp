(in-package :tc)

;;; we define helper functions and symbol macros that will simplify
;;; the notation in the later code

(defmacro mcc-tc-coeff-denom (A B nu-ab nu-ba)
  "Expands into (A * B) - (nu-ab * nu-ba)

This is the denominator in mu-alpha, k-alpha, D-alpha, alpha-T-alpha
equations for transport coefficients Sharipov&Kalempa, 2002 (74-77)"
  `(- (* ,A ,B)
      (* ,nu-ab ,nu-ba)))

#|(define-symbol-macro omega-11% (omega-11 c-alpha-beta temp))
(define-symbol-macro omega-12% (omega-12 c-alpha-beta temp))
(define-symbol-macro omega-13% (omega-13 c-alpha-beta temp))
(define-symbol-macro omega-22% (omega-22 c-alpha-beta temp))|#