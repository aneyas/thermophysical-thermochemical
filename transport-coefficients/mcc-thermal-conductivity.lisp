;; Mirko Vukovic
;; Time-stamp: <2012-05-27 21:17:26 mcc-thermal-conductivity.lisp>
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


(export '(k-binary k-alpha1))

(defun k-alpha (P-alpha m-alpha m-beta
		     Phi-alpha Phi-beta
		     nu6-alpha-beta nu6-beta-alpha)
  "Thermal conductivity of species alpha


Sharipov&Kalempa, 2002 (75)"
  (/ (* 2.5 (/ (* P-alpha +k+)
	       m-alpha)
	(+ Phi-beta (* (sqrt (/ m-alpha m-beta))
		       nu6-alpha-beta)))
     (mcc-tc-coeff-denom Phi-alpha Phi-beta
	  nu6-alpha-beta nu6-beta-alpha)))

(define-test k-alpha
  (let ((c (* 5/2 +kb+))
	(lisp-unit:*epsilon* 1e-4))
    (assert-numerical-equal
     (- c)
     (k-alpha 1 1 1 0 0 1 1) "simple-PHIs=0")
    (assert-numerical-equal
     c
     (k-alpha 1 1 1 1 1 0 0) "simple=nus=0")
    (assert-numerical-equal
     (- (* 2 c))
     (k-alpha 2 1 1 0 0 1 1) "P-alpha")
    (assert-numerical-equal
     (/ c -2)
     (k-alpha 1 4 1 0 0 1 1) "m-alpha")
    (assert-numerical-equal
     (/ c -2)
     (k-alpha 1 4 1 0 0 1 1) "m-beta")
    (assert-numerical-equal
     (* 2 c)
     (k-alpha 1 1 1 1/2 2 0 0) "PHI-beta")
    (assert-numerical-equal
     c
     (k-alpha 1 1 1 1 2 0 0) "PHI-alpha")))



(defmethod lambda-ig2-1 ((model (eql 'mcc))
			 species-1 species-2 p-1 p-2 temp
			 &key)
  (let ((m-1 (* (mass species-1) +amu+))
	(m-2 (* (mass species-2) +amu+))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters species-1 species-1))
	(c-1-2 (make-collision-parameters species-1 species-2))
	(c-2-1 (make-collision-parameters species-2 species-1))
	(c-2-2 (make-collision-parameters species-2 species-2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			(m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2 1 2) (1 1 2 2)
	     (nu5 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta
		  omega-22-alpha-beta
		  omega-12-alpha-beta
		  omega-13-alpha-beta)
	     (nu6 m*-alpha-beta m-alpha m-beta n-beta
		  omega-11-alpha-beta
		  omega-22-alpha-beta
		  omega-12-alpha-beta
		  omega-13-alpha-beta))
	  (with-f-alpha-beta*
	      ((1 2) (2 1)
	       (phi nu5-alpha-alpha
		    nu6-alpha-alpha
		    nu5-alpha-beta)
	       (k-alpha p-alpha m-alpha m-beta
			phi-alpha-beta phi-beta-alpha
			nu6-alpha-beta nu6-beta-alpha))
	    (values (+ k-alpha-1-2 k-alpha-2-1)
		    k-alpha-1-2 k-alpha-2-1)))))))


(define-test MU2-tc-vs-gkf/lj
  "Comparison between the code in TRANSPORT-COEFFICIENTS and
  GAS-KINETIC-FORMULARY"
  (let ((n1 0.9e22)
	(n2 0.1e22)
	(s1 (make-species-lennard-jones-6/12-potential :Ar))
	(temp 300d0)
	(lisp-unit:*epsilon* 1e-4))
    (let ((s2 (make-species-lennard-jones-6/12-potential :Ar))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:k-alpha1 s1 s2 n1 n2 temp)
				(lambda-ig2-1 'MCC s1 s2 p1 p2 temp ))))
    (let ((s2 (make-species-lennard-jones-6/12-potential :He))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:k-alpha1 s1 s2 n1 n2 temp)
				(lambda-ig2-1 'MCC s1 s2 p1 p2 temp))))))


(define-test LAMBDA2-tc-vs-gkf/HS
  "Comparison between the code in TRANSPORT-COEFFICIENTS and
  GAS-KINETIC-FORMULARY"
  (let ((n1 0.9e22)
	(n2 0.1e22)
	(s1LJ (make-species-lennard-jones-6/12-potential :Ar))
	(temp 300d0)
	(lisp-unit:*epsilon* 1e-4))
    (let ((s2LJ (make-species-lennard-jones-6/12-potential :Ar))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((s1 (make-hard-sphere-potential
		  (mass s1LJ) (sigma s1LJ) (species S1LJ)))
	     (s2 (make-hard-sphere-potential
		  (mass s2LJ) (sigma s2LJ) (species S2LJ)))
	     (coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:k-alpha1 s1 s2 n1 n2 temp)
				(lambda-ig2-1 'MCC s1 s2 p1 p2 temp ))))
    (let ((s2LJ (make-species-lennard-jones-6/12-potential :He))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let* ((s1 (make-hard-sphere-potential
		  (mass s1LJ) (sigma s1LJ) (species S1LJ)))
	     (s2 (make-hard-sphere-potential
		  (mass s2LJ) (sigma s2LJ) (species S2LJ)))
	     (coll-param (make-collision-parameters s1 s2)))
	(assert-numerical-equal (gkf:k-alpha1 s1 s2 n1 n2 temp)
				(lambda-ig2-1 'MCC s1 s2 p1 p2 temp))))))