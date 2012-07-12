;; Mirko Vukovic
;; Time-stamp: <2012-05-27 22:27:50 mcc-thermal-diffusivity.lisp>
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

(defun at-alpha (Phi-alpha Phi-beta
			m-alpha m-beta
			nu6-alpha-beta nu6-beta-alpha)
  "Sharipov&Kalempa, 2002 (77)

Used in thermal diffusion ratio (73)

I replaced the dependence on nu6-1-2 with nu6-alpha-beta.  That seems
to make more sense"
  (/ (+ Phi-beta (* (sqrt (/ m-alpha m-beta))
		    nu6-alpha-beta))
     (mcc-tc-coeff-denom Phi-alpha Phi-beta
			 nu6-alpha-beta nu6-beta-alpha)))


(define-test at-alpha
  (let ((lisp-unit:*epsilon* 1e-4))
    ;; First batch focuses on PHI-beta in the numerator and
    ;; PHI-alpha/beta in the denumerator by setting all the collision
    ;; frequences to zero.  m-beta =1 to prevent /0 in numerator
    (assert-numerical-equal
     1
     (at-alpha 1 1 0 1 0 0) "simple")
    (assert-numerical-equal
     1
     (at-alpha 1 2 0 1 0 0) "PHI-beta-cancellation/simple")
    (assert-numerical-equal
     1/2
     (at-alpha 2 2 0 1 0 0) "PHI-alpha/simple")
    ;; Second batch adds the mass-fraction terms in the numerator
    (assert-numerical-equal
     2
     (at-alpha 1 1 1 1 1 0) "m-alpha/m-beta * nu6")
    (assert-numerical-equal
     3
     (at-alpha 1 1 4 1 1 0) "m-alpha")
    (assert-numerical-equal
     2
     (at-alpha 1 1 4 4 1 0) "m-beta")
    (assert-numerical-equal
     3
     (at-alpha 1 1 4 4 2 0) "nu6 in numerator")
    ;; Third batch focuses on the nu6 terms in the denominator
    (assert-numerical-equal
     1
     (at-alpha 2 1 0 1 1 1) "nu6 terms")
    (assert-numerical-equal
     -1
     (at-alpha 2 1 0 1 3 1) "nu6-alpha-beta")
    (assert-numerical-equal
     -1
     (at-alpha 2 1 0 1 1 3) "nu6-alpha-beta")))


(defmethod alpha-T-1 ((model (eql 'mcc))
		      (s1 lennard-jones-6/12-potential)
		      (s2 lennard-jones-6/12-potential)
		      p-1 p-2 temp &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass  s2) +amu+))
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
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1) (2)
	     (nu2 m*-alpha-beta m-alpha n-beta
		  omega-11-alpha-beta
		  omega-12-alpha-beta))
	  (with-f-alpha-beta
	      ((1 2 1 2) (1 1 2 2)
	       (nu5 m*-alpha-beta m-alpha m-beta n-beta
		    omega-11-alpha-beta
		    omega-22-alpha-beta
		    omega-12-alpha-beta
		    omega-13-alpha-beta)
	       (nu6 m*-alpha-beta m-alpha m-beta
		    n-beta omega-11-alpha-beta
		    omega-22-alpha-beta
		    omega-12-alpha-beta
		    omega-13-alpha-beta))
	    (with-f-alpha-beta*
		((1 2) (2 1)
		 (phi nu5-alpha-alpha
		      nu6-alpha-alpha
		      nu5-alpha-beta)
		 (at-alpha phi-alpha-beta phi-beta-alpha
			   m-alpha m-beta
			   nu6-alpha-beta nu6-beta-alpha))
	      (let ((n (+ n-1 n-2)))
		(* -1.25 (/ (* n nu2-1-2)
			    n-2)
		   (- at-alpha-1-2 (* (expt (/ m-1 m-2) 2)
				      at-alpha-2-1)))))))))))


(defmethod alpha-T-1 ((model (eql 'mcc))
		      (s1 hard-sphere-potential)
		      (s2 hard-sphere-potential)
		      p-1 p-2 temp &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass  s2) +amu+))
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
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1) (2)
	     (nu2 m*-alpha-beta m-alpha n-beta
		  omega-11-alpha-beta
		  omega-12-alpha-beta))
	  (with-f-alpha-beta
	      ((1 2 1 2) (1 1 2 2)
	       (nu5 m*-alpha-beta m-alpha m-beta n-beta
		    omega-11-alpha-beta
		    omega-22-alpha-beta
		    omega-12-alpha-beta
		    omega-13-alpha-beta)
	       (nu6 m*-alpha-beta m-alpha m-beta
		    n-beta omega-11-alpha-beta
		    omega-22-alpha-beta
		    omega-12-alpha-beta
		    omega-13-alpha-beta))
	    (with-f-alpha-beta*
		((1 2) (2 1)
		 (phi nu5-alpha-alpha
		      nu6-alpha-alpha
		      nu5-alpha-beta)
		 (at-alpha phi-alpha-beta phi-beta-alpha
			   m-alpha m-beta
			   nu6-alpha-beta nu6-beta-alpha))
	      (let ((n (+ n-1 n-2)))
		(* -1.25 (/ (* n nu2-1-2)
			    n-2)
		   (- at-alpha-1-2 (* (expt (/ m-1 m-2) 2)
				      at-alpha-2-1)))))))))))

(define-test alphat-tc-vs-gkf/lj
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
      (assert-numerical-equal (gkf::alpha-t s1 s2 n1 n2 temp)
			      (alpha-t-1 'MCC s1 s2 p1 p2 temp)))
    (let ((s2 (make-species-lennard-jones-6/12-potential :He))
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (assert-numerical-equal (gkf:alpha-t s1 s2 n1 n2 temp)
			      (alpha-t-1 'MCC s1 s2 p1 p2 temp)))))