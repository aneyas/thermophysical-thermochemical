;;;; package.lisp

(defpackage #:transport-coefficients
  (:nicknames #:tc)
  (:use #:cl #:lisp-unit
	#:defgeneric+default
	#:molecular-potentials
	#:omega-xx
	#:gas-kinetic-formulary
	#:index-sugar)
  (:shadow #:gas-kinetic-formulary :D12
	   :alpha-t
	   :delta-1)
  (:export :cck :mcc
	   :*d12-default-method*
	   :d12 :d12-1
	   :d12b :d12b-1
	   :delta :delta-1 ;; this will cause conflicts with gas-kinetic formulary
	   :*mu1-default-method*
	   :*mu2-default-method*
	   :mu1 :mu1-1
	   :mu2 :mu2-1
	   :*lambda1-default-method*
	   :*lambda-ig1-default-method*
	   :*lambda-ig2-default-method*
	   :lambda-ig1 :lambda-ig1-1
	   :lambda-ig2 :lambda-ig2-1
	   :lambda1 :lambda1-1
	   :lambda2 :lambda2-1
	   :alpha-t :alpha-t-1))

