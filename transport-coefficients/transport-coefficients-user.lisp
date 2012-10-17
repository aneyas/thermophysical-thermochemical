

(in-package :tcu)





(defun mu-cck/mcc (s1 p1 p2 temp)
  "Comparison of D12 calculated via CCK and MCC"
  (let ((lj1 (make-species-lennard-jones-6/12-potential s1)))
    (let ((mu-cck (mu1-1 'tc::cck lj1 temp))
	  (mu-mcc (multiple-value-bind (mu1 mu2)
		      (mu2-1 'tc::mcc lj1 lj1 p1 p2 temp)
		    (+ mu1 mu2))))
      (values mu-cck mu-mcc))))

(define-test mu-cck/mcc
  (let ((p1 50d0)
	(p2 5d-3)
	(temp 300d0))
    (mu-cck/mcc :Ar p1 p2 temp)))

(defun lambda2-ig-cck/mcc (s1 p1 p2 temp)
  "Comparison of D12 calculated via CCK and MCC"
  (let ((lj1 (make-species-lennard-jones-6/12-potential s1)))
    (let ((lambda-cck (lambda-ig1-1 'cck lj1 temp))
	  (lambda-mcc (lambda-ig2-1 'mcc lj1 lj1 p1 p2 temp)))
      (values lambda-cck lambda-mcc))))

(define-test lambda-ig1/2-cck/mcc
  (let ((p1 50d0)
	(p2 5d0)
	(temp 300.0d0))
    (lambda2-ig-cck/mcc :Ar p1 p2 temp)))



(defun alpha-t-ig-cck/mcc (s1 s2 p1 p2 temp)
  "Comparison of D12 calculated via CCK and MCC"
  (let ((pot-1 (make-species-lennard-jones-6/12-potential s1))
	(pot-2 (make-species-lennard-jones-6/12-potential s2)))
    (let ((alpha-t-cck (alpha-t-1 'cck pot-1 pot-2 p1 p2 temp))
	  (alpha-t-mcc (alpha-t-1 'mcc pot-1 pot-2 p1 p2 temp)))
      (values alpha-t-cck alpha-t-mcc))))

(define-test alpha-t-ig-cck/mcc
  (let ((p1 50d0)
	(p2 5d0)
	(temp 300.0d0))
    (alpha-t-ig-cck/mcc :Ar :He p1 p2 temp)))