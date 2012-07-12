;; Mirko Vukovic
;; Time-stamp: <2012-06-13 11:01:52 mcc-viscosity.lisp>
;; 
;; Copyright 2011 Mirko Vukovic
;; Distributed under the terms of the GNU General Public License
;; 
;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;; 
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; 
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

(in-package :tc)

(export '(#|mu-alpha1-hs|# mu-alpha1 delta-alpha-beta mu-binary))
;; Hard sphere calculations are in conflict with the rest of the code


(defun mu-alpha (P-alpha Psi-alpha Psi-beta
		      nu4-alpha-beta nu4-beta-alpha)
  "Viscosity of species alpha
Sharipov&Kalempa, 2002 (74)"
  (/ (* P-alpha (+ Psi-beta nu4-alpha-beta))
     (mcc-tc-coeff-denom Psi-alpha Psi-beta nu4-alpha-beta nu4-beta-alpha)))
;; &vars&eval-calls (1 2) (2 1)

(define-test mu-alpha
  (assert-numerical-equal
   1 (mu-alpha 1 1 1 0 0))
  (assert-numerical-equal
   2 (mu-alpha 2 1 1 0 0) "P-alpha")
  ;; PSI-beta and nu4-alpha-beta are related through the numerator.
  ;; Tests of their interaction
  (assert-numerical-equal
   1
   (mu-alpha 1 1 2 0 1) "PSI-beta")
  (assert-numerical-equal
   -1
   (mu-alpha 1 1 0 2 1) "nu4-alpha-beta")
  (assert-numerical-equal
   (/ (+ 1 2)
      (- 1 2))
   (mu-alpha 1 1 1 2 1) "PSI-beta/nu4-alpha-beta")
  (assert-numerical-equal
   (/ (+ 2 1)
      (- 2 1))
   (mu-alpha 1 1 2 1 1) "PSI-beta/nu4-alpha-beta")
  (assert-numerical-equal
   1/2 (mu-alpha 1 2 1 0 0) "PSI-alpha")
  (assert-numerical-equal
   -1/2 (mu-alpha 1 0 0 1 2) "nu4-beta-alpha"))


#|(defun mu-alpha1-hs (n-1 n-2 m-1 m-2 sigma-1 sigma-2 temp)
  "mu-alpha

Sharipov&Kalempa, 2002 (74)"
  (bind-p-alpha
    (bind-m*-alpha-beta
      (bind-omega11-hs-alpha-beta
	(bind-omega22-hs-alpha-beta
	  (bind-nu3-alpha-beta
	    (bind-nu4-alpha-beta
	      (bind-psi-alpha
		(bind-mu-alpha
		  (values mu-1 mu-2 ))))))))))|#


(defmethod mu2-1 ((method (eql 'mcc))
		  (s1 lennard-jones-6/12-potential)
		  (s2 lennard-jones-6/12-potential) p-1 p-2 temp
		  &key)
  "Calculate mu-alpha1,2 for `species-1' and `species-2' with respective
densities of `n-1' and `n-2' at temperature `temp'

It returns mu-alpha1,2 as a list

As a second value, the quantities (/ P-1,2 mu-1,2) are also returned
as a list"
  (let ((m-1 (mass s1))
	(m-2 (mass s2))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			(m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2 1 2) (1 1 2 2)
	     (nu3 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta omega-22-alpha-beta)
	     (nu4 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta omega-22-alpha-beta))
	  (with-f-alpha-beta*
	      ((1 2) (2 1)
	       (psi nu3-alpha-alpha nu4-alpha-alpha nu3-alpha-beta)
	       (mu-alpha p-alpha psi-alpha-beta psi-beta-alpha
			 nu4-alpha-beta nu4-beta-alpha))
	    (values (+ mu-alpha-1-2 mu-alpha-2-1)
		    (list mu-alpha-1-2 mu-alpha-2-1)
		    (list (/ p-1 mu-alpha-1-2)
			  (/ p-2 mu-alpha-2-1)))))))))

(defmethod mu2-1 ((method (eql 'mcc))
		  (s1 hard-sphere-potential)
		  (s2 hard-sphere-potential) p-1 p-2 temp &key)
  (let ((m-1 (mass s1))
	(m-2 (mass s2))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			(m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2 1 2) (1 1 2 2)
	     (nu3 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta omega-22-alpha-beta)
	     (nu4 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta omega-22-alpha-beta))
	  (with-f-alpha-beta*
	      ((1 2) (2 1)
	       (psi nu3-alpha-alpha nu4-alpha-alpha nu3-alpha-beta)
	       (mu-alpha p-alpha psi-alpha-beta psi-beta-alpha
			 nu4-alpha-beta nu4-beta-alpha))
	    (values (+ mu-alpha-1-2 mu-alpha-2-1)
		    (list mu-alpha-1-2 mu-alpha-2-1)
		    (list (/ p-1 mu-alpha-1-2)
			  (/ p-2 mu-alpha-2-1)))))))))



(defun delta-alpha-beta (R spec-alpha spec-beta p-alpha p-beta temp)
  (multiple-value-bind (mu-alpha mu-beta gamma-s)
      (mu2-1 'MCC  spec-alpha spec-beta p-alpha p-beta temp)
    (declare (ignore mu-alpha mu-beta))
    (destructuring-bind (gamma-alpha gamma-beta) gamma-s
      (list (delta1 R gamma-alpha (mass spec-alpha) temp)
	    (delta1 R gamma-beta (mass spec-beta) temp)))))



(define-test MU2-tc-vs-gkf/lj
  "Comparison between the code in TRANSPORT-COEFFICIENTS and
  GAS-KINETIC-FORMULARY"
  (let ((n1 0.9e22)
	(n2 0.1e22)
	(s1 (make-species-lennard-jones-6/12-potential :Ar))
	(temp 300d0)
	(lisp-unit:*epsilon* 1e-6))
    (let ((s2 (make-species-lennard-jones-6/12-potential :Ar))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:mu-alpha1 s1 s2 n1 n2 temp)
				(mu2-1 'MCC s1 s2 p1 p2 temp ))))
    (let ((s2 (make-species-lennard-jones-6/12-potential :He))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:mu-alpha1 s1 s2 n1 n2 temp)
				(mu2-1 'MCC s1 s2 p1 p2 temp))))))

(define-test MU2-tc-vs-gkf/HS
  "Comparison between the code in TRANSPORT-COEFFICIENTS and
  GAS-KINETIC-FORMULARY"
  (let ((n1 0.9e22)
	(n2 0.1e22)
	(s1LJ (make-species-lennard-jones-6/12-potential :Ar))
	(temp 300d0)
	(lisp-unit:*epsilon* 1e-6))
    (let ((s2LJ (make-species-lennard-jones-6/12-potential :Ar))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((s1 (make-hard-sphere-potential
		  (mass s1LJ) (sigma s1LJ) (species S1LJ)))
	     (s2 (make-hard-sphere-potential
		  (mass s2LJ) (sigma s2LJ) (species S2LJ)))
	     (coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:mu-alpha1 s1 s2 n1 n2 temp)
				(mu2-1 'MCC s1 s2 p1 p2 temp ))))
    (let ((s2LJ (make-species-lennard-jones-6/12-potential :He))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((s1 (make-hard-sphere-potential
		  (mass s1LJ) (sigma s1LJ) (species S1LJ)))
	     (s2 (make-hard-sphere-potential
		  (mass s2LJ) (sigma s2LJ) (species S2LJ)))
	     (coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:mu-alpha1 s1 s2 n1 n2 temp)
				(mu2-1 'MCC s1 s2 p1 p2 temp))))))

#|
(define-test mu-binary%
  (let ((m-1 1e-26)
	(m-2 1e-26)
	(Temp 300)
	(sigma-1 1e-10)
	(sigma-2 1e-10)
	(n-1 1e22)
	(n-2 1e22)
	(lisp-unit:*epsilon* 1e-4))
    (bind-m*-alpha-beta
      (bind-omega11-alpha-beta
	(bind-omega22-alpha-beta
	  (assert-number-equal
	   (gkf:mu-hs (* +NA+ m-1)
		      temp sigma-1)
	   (mu-binary% n-1 n-2 m-1 m-2 temp omega11-1-1 omega11-1-2 omega11-2-1 omega11-2-2
		       omega22-1-1 omega22-1-2 omega22-2-1 omega22-2-2)
	   "mu"))))))

|#