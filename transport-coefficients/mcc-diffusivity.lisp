;; Mirko Vukovic
;; Time-stamp: <2012-07-18 13:51:45 mcc-diffusivity.lisp>
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

(in-package :transport-coefficients)

(export '(DD-alpha-beta D12b-1 *D12-default-method* CEC MCC))


;;; McCormack formulation, as given by Sharipov & Kalempa, 2002



(defun D (m-alpha m-beta
		     Phi-alpha Phi-beta
		     nu2-alpha-beta nu2-beta-alpha
		     nu6-alpha-beta nu6-beta-alpha)
  "Sharipov&Kalempa, 2002 (76)"
  (/ (- (* nu2-alpha-beta Phi-beta)
	(* nu2-beta-alpha nu6-alpha-beta
	   (sqrt (/ m-beta m-alpha))))
     (mcc-tc-coeff-denom Phi-alpha Phi-beta nu6-alpha-beta nu6-beta-alpha)))



(define-test d
  (let ((lisp-unit:*epsilon* 1e-4))
    (assert-numerical-equal
     1
     (d 1 1 1 1 1 0 0 0) "simple-nu2-beta-alpha/nu6=0")
    (assert-numerical-equal
     2
     (d 1 1 1 1 2 0 0 0) "nu2-alpha-beta")
    (assert-numerical-equal
     2
     (d 1 1 1/2 2 1 0 0 0) "Phi-beta")
    (assert-numerical-equal
     1
     (d 1 1 1 2 1 0 0 0) "Phi-alpha")
    (assert-numerical-equal
     0
     (d 1 1 1 1 1 1 1 0) "Phi-alpha")
    (assert-numerical-equal
     -1
     (d 1 1 1 1 1 2 1 0) "nu2-beta-alpha")
    (assert-numerical-equal
     -3
     (d 1 1 1 1 1 2 2 0) "nu2-beta-alpha")
    (assert-numerical-equal
     -7
     (d 1 4 1 1 1 2 2 0) "m-beta")
    (assert-numerical-equal
     -3
     (d 4 4 1 1 1 2 2 0) "m-alpha")
    (assert-numerical-equal
     3
     (d 4 4 1 1 1 2 2 1) "nu6-beta-alpha")))


(defun Delta-mcc (m-alpha m-beta
		  D-alpha D-beta
		  nu1-alpha-beta nu2-alpha-beta)
  "Delta (second order correction term to the first order D12.

Delta is definded in Ferziger & Kapor 7.3-39,40.

In S&K02, Delta is not explicitly defined, but it's form is seen by
comparing eq. 72 with that of Ferziger & Kapor"
  (* (/ 5.0 8.0)
     (/ nu2-alpha-beta nu1-alpha-beta)
     (+ D-alpha
	(* (/ m-alpha m-beta)
	   D-beta))))

(define-test delta-mcc
  (assert-number-equal 5/8
		       (delta-mcc 1 1 1 0 1 1)
		       "Leading coefficient")
  (assert-number-equal 10/8
		       (delta-mcc 1 1 1 0 1 2)
		       "nu^(2)")
  (assert-number-equal 5/16
		       (delta-mcc 1 1 1 0 2 1)
		       "nu^(1)")
  (assert-number-equal 10/8
		       (delta-mcc 1 1 2 0 1 1)
		       "D1")
  (assert-number-equal 5/8
		       (delta-mcc 1 1 0 1 1 1)
		       "D2-basic")
  (assert-number-equal 10/8
		       (delta-mcc 1 1 0 2 1 1)
		       "D2-scaling")
  (assert-number-equal 10/8
		       (delta-mcc 2 1 0 1 1 1)
		       "m1")
  (assert-number-equal 5/16
		       (delta-mcc 1 2 0 1 1 1)
		       "m2"))

(defun D12-mcc (P-beta n-alpha n-beta
		m-alpha m-beta
		D-alpha D-beta
		nu1-alpha-beta nu2-alpha-beta)
  "Sharipov&Kalempa, 2002 (72), generalized to alpha-beta"
  (let ((n (+ n-alpha n-beta)))
    (/ (/ P-beta
	  (* n m-alpha nu1-alpha-beta))
       (- 1 (Delta-mcc m-alpha m-beta d-alpha d-beta
		       nu1-alpha-beta nu2-alpha-beta)))))

(define-test D12-mcc
  ;; All but the last test exercise only the factors that do not
  ;; involve Delta-mcc.  For that we set nu^2 to zero.  Other
  ;; arguments to Delta-mcc then don't matter.
  (let ((lisp-unit:*epsilon* 1e-6))
    (assert-number-equal 1
			 (D12-mcc 1 0.5 0.5 1 1 1 1 1 0)
			 "base case")
    (assert-number-equal 2
			 (D12-mcc 2 0.5 0.5 1 1 1 1 1 0)
			 "P")
    (assert-number-equal 1/2
			 (D12-mcc 1 1.3 0.7 1 1 1 1 1 0)
			 "n1+n2")
    (assert-number-equal 1/2
			 (D12-mcc 1 0.5 0.5 1 1 1 1 2 0)
			 "nu^(1)")
    (assert-number-equal (/ 1 (- 1 5/8))
			 (D12-mcc 1 0.5 0.5 1 1 1 0 1 1))))
  
  
(defmethod D12b-1
    ((method (eql 'MCC))
     (s1 lennard-jones-6/12-potential)
     (s2 lennard-jones-6/12-potential) p-1 p-2 
     temp &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass s2) +amu+))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2) (m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2) (2 1)
	     (nu1 m*-alpha-beta m-alpha n-beta omega-11-alpha-beta)
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
		 (d m-alpha m-beta phi-alpha-beta
		    phi-beta-alpha nu2-alpha-beta
		    nu2-beta-alpha nu6-alpha-beta
		    nu6-beta-alpha)
		 (d12-mcc p-beta n-alpha n-beta
			  m-alpha m-beta
			  d-alpha-beta d-beta-alpha
			  nu1-alpha-beta nu2-alpha-beta))
	      (values d12-mcc-1-2 d12-mcc-2-1))))))))


(defmethod D12b-1
    ((method (eql 'MCC))
     (s1 hard-sphere-potential)
     (s2 hard-sphere-potential) p-1 p-2 
     temp &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass s2) +amu+))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2) (m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha ((1 2)
		       (p n-alpha temp))
	  (with-f-alpha-beta
	      ((1 2) (2 1)
	       (nu1 m*-alpha-beta m-alpha n-beta omega-11-alpha-beta)
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
		   (d m-alpha m-beta phi-alpha-beta
		      phi-beta-alpha nu2-alpha-beta
		      nu2-beta-alpha nu6-alpha-beta
		      nu6-beta-alpha)
		   (d12-mcc p-beta n-alpha n-beta
		       m-alpha m-beta
		       d-alpha-beta d-beta-alpha
		       nu1-alpha-beta nu2-alpha-beta))
		(values d12-mcc-1-2 d12-mcc-2-1)))))))))


(defmethod delta-1 ((method (eql 'mcc))
		    (s1 lennard-jones-6/12-potential)
		    (s2 lennard-jones-6/12-potential) 
		    p-1 p-2
		    temp
		    &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass s2) +amu+))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2) (m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2) (2 1)
	     (nu1 m*-alpha-beta m-alpha n-beta omega-11-alpha-beta)
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
		 (d m-alpha m-beta phi-alpha-beta
		    phi-beta-alpha nu2-alpha-beta
		    nu2-beta-alpha nu6-alpha-beta
		    nu6-beta-alpha)
		 (delta-mcc m-alpha m-beta
			    d-alpha-beta d-beta-alpha
			    nu1-alpha-beta
			    nu2-alpha-beta))
	      (values delta-mcc-1-2 delta-mcc-2-1))))))))

(defmethod delta-1 ((method (eql 'mcc))
		    (s1 hard-sphere-potential)
		    (s2 hard-sphere-potential) 
		    p-1 p-2
		    temp
		    &key)
  (let ((m-1 (* (mass s1) +amu+))
	(m-2 (* (mass s2) +amu+))
	(n-1 (n p-1 temp))
	(n-2 (n p-2 temp))
	(c-1-1 (make-collision-parameters s1 s1))
	(c-1-2 (make-collision-parameters s1 s2))
	(c-2-1 (make-collision-parameters s2 s1))
	(c-2-2 (make-collision-parameters s2 s2)))
    (with-f-alpha-beta ((1 2 1 2) (1 1 2 2) (m* m-alpha m-beta))
      (with-f-alpha-beta ((1 2 1 2) (1 1 2 2)
			  (omega-11 c-alpha-beta temp)
			  (omega-12 c-alpha-beta temp)
			  (omega-13 c-alpha-beta temp)
			  (omega-22 c-alpha-beta temp))
	(with-f-alpha-beta
	    ((1 2) (2 1)
	     (nu1 m*-alpha-beta m-alpha n-beta omega-11-alpha-beta)
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
		 (d m-alpha m-beta phi-alpha-beta
		    phi-beta-alpha nu2-alpha-beta
		    nu2-beta-alpha nu6-alpha-beta
		    nu6-beta-alpha)
		 (delta-mcc m-alpha m-beta
			    d-alpha-beta d-beta-alpha
			    nu1-alpha-beta
			    nu2-alpha-beta))
	      (values delta-mcc-1-2 delta-mcc-2-1))))))))

(define-test D12-cck/mcc
  "Comparison between first order cck and second order mcc binary
diffusion coefficient.  I expect only very rough agreement."
  (let ((n1 0.9e22)
	(n2 0.1e22)
	(s1 +Ar/LJ6-12+)
	(temp 300d0)
	(lisp-unit:*epsilon* 5e-2))
    (let ((s2 +Ar/LJ6-12+)
	  #+ (or) (lisp-unit:*epsilon* 5e-2)
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let ((coll-params (make-collision-parameters s1 s2)))
	(assert-number-equal (d12-1 'cck coll-params p temp)
			     (d12b-1 'MCC s1 s2 p1 p2 temp))))
    (let ((s2 +H2/LJ6-12+)
	  (p (p (+ n1 n2) temp))
	  (p1 (p n1 temp))
	  (p2 (p n2 temp)))
      (let ((coll-params (make-collision-parameters s1 s2)))
	(assert-number-equal (d12-1 'cck coll-params p temp)
			     (d12b-1 'MCC s1 s2 p1 p2 temp))))))




(define-test D12-MCC-species-swap
  "Is D12 consistent when we swap the species"
  (let ((p-Ar 5d0)
	(p-H2 5d0)
	(temp 300d0)
	(lisp-unit:*epsilon* 1e-12))
    (let ((n-Ar (n p-Ar temp))
	  (n-H2 (n p-H2 temp)))
      (multiple-value-bind (D12/Ar-H2/0 D12/Ar-H2/1)
	  (D12b-1 'MCC +Ar/LJ6-12+ +H2/LJ6-12+ n-Ar n-H2 temp)
	(multiple-value-bind (D12/H2-Ar/0 D12/H2-Ar/1)
	    (D12b-1 'MCC +H2/LJ6-12+ +Ar/LJ6-12+ n-H2 n-Ar temp)
	  (assert-number-equal D12/Ar-H2/0 D12/H2-Ar/1)
	  (assert-number-equal D12/Ar-H2/1 D12/H2-Ar/0))))))

