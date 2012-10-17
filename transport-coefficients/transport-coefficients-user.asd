;;;; transport-coefficients.asd

(asdf:defsystem #:transport-coefficients-user
  :serial t
  :depends-on (#:transport-coefficients
	       #+probably-unnecessary #:gas-kinetics-physics-constants
	       #:gnuplot-interface
	       #:mv-gnuplot
	       #:mv-grid-utils)
  :components
  ((:module "package-setup"
	    :pathname #P"./"
	    :components ((:file "transport-coefficients-user-package-def")))
   (:module "setup"
	    :pathname #P"./"
	    :components ((:file "initialize-species")
			 (:file "initialize-gnuplot")
			 (:file "test-potentials")))
   (:module "example-plots"
	    :pathname #P"./"
	    :components ((:file "binary-diffusion-plots")
			 (:file "binary-diffusion-comp-with-ilt")
			 (:file "thermal-diffusivity-plots")))))


