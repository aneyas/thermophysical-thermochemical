;;;; transport-coefficients.asd

(asdf:defsystem #:transport-coefficients
  :serial t
  :depends-on (#:alexandria
               #:lisp-unit
               #:physics-constants
	       #:molecular-potentials
	       #:collision-integrals
	       #:defgeneric+default
	       #:gas-kinetic-formulary
	       #:index-sugar)
  :components
  ((:module "package-setup"
	    :pathname #P"./"
	    :components ((:file "transport-coefficients-package-def")))
   (:module "setup"
	    :pathname #P"./"
	    :components ((:file "fundamental-constants")
			 (:file "generic-functions")))
   (:module "transport-coefficients"
	    :pathname #P"./"
	    :components ((:file "cck-transport-coefficients")
			 (:file "mcc-setup")
			 (:file "mcc-diffusivity")
			 (:file "mcc-viscosity")
			 (:file "mcc-thermal-conductivity")
			 (:file "mcc-thermal-diffusivity")))))


