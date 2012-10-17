;;;; package.lisp

(defpackage #:transport-coefficients-user
  (:nicknames #:tcu)
  (:use #:cl #:lisp-unit
	;;
	#+probably-unnecessary #:gas-kinetics-physics-constants
	#:transport-coefficients
	#:molecular-potentials
	#:collision-integrals
	;;
	#:gnuplot-interface
	#:mv-gnuplot
	;;
	#:mv-grid)
  (:documentation
"This package provides an environement for computing and using
transport coefficients.  It can be used as a template for other
pakcages that use transport coefficients."))



