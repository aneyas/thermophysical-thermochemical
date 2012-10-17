;; Mirko Vukovic
;; Time-stamp: <2012-07-18 13:36:49 cck+mcc-diffusivity-tests.lisp>
;; 
;; Copyright 2012 Mirko Vukovic
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

(in-package #:transport-coefficients)

;;; This file has multiple tests comparing the CCK & MCC diffusion
;;; coefficients.  Since the MCC is accurate to second order, I use
;;; the Delta correction (Ferziger & Capor, 7.3-40) to get the CCK
;;; second order correction.

;;; Lennard-Jones potential tests

(define-symbol-macro CP (* C P))
(define-symbol-macro |(1-C)P| (* (- 1d0 C) P))

(defgeneric D12-Kee-3.100 (coll-pot pressure temp)
  (:documentation "Eq. 3.100 for diffusion from Kee et al

COLL-POT is either a Lennar-Jones-6/12-potential or a Hard-sphere-potential
TEMP -- Kelvin
PRESSURE -- Pascal")
  (:method ((coll-pot lennard-jones-6/12-potential) pressure temp)
    (* 0.018834945 ;; obtained from the let of (D12-1 'cck ...) function
       (/ (sqrt (/ (expt temp 3) (* 1e-0 (mass coll-pot))))
	  (* pressure (expt (sigma coll-pot) 2)
	     (omega*-11 (/ temp (epsilon/k coll-pot)))))))
  (:method ((coll-pot hard-sphere-potential) pressure temp)
    (* 0.018834945
       (/ (sqrt (/ (expt temp 3) (* 1e-0 (mass coll-pot))))
	  (* pressure (expt (sigma coll-pot) 2))))))


(define-test D12-KEE-3.100
  "Test the function against a value that was calculated semi-manually
using Eq. 3.100 of Kee et al"
  ;; we limit the epsilon because the numeric value was calculated
  ;; using the numerical coefficient of 0.0188 while d12-kee-3.100
  ;; uses a more precise value
  (assert-number-equal 0.01854679270047735d0;; 0.5854129255295771d0
		  (d12-kee-3.100 +Ar-Ar/LJ6-12+ 100d0 300d0)))

(define-test D12-1-cck
  "Compare the Ar/Ar D12 diffusion as calculated with the
Chapman-Cowling/Kihara approximation vs eq. 3.100 of Kee et al"
  (let ((temp 300d0)
	(p 100d0))
    (let ((d12 (d12-kee-3.100 +Ar-Ar/LJ6-12+ p temp)))
      (assert-number-equal d12
			   (d12-1 'cck +Ar-Ar/LJ6-12+ p temp)))))


(define-test d12b-1-mcc/cck
  "A rough test.  cck is calculated to first order while mcc to second"
  (let ((temp 300d0)
	(p 100d0)
	(c 1e-9))
    (let ((d12 (d12-1 'cck +Ar-Ar/LJ6-12+ p temp))
	  (lisp-unit:*epsilon* 3e-3))
      (assert-number-equal d12
			   (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+
				   CP |(1-C)P| temp)))))


(define-test delta-cck/mcc
  "Comparison of delta for the Chapman-Enskog and McCormack models for
several values of C.  Incredibly enough, they agree admirably to at
least 10^{-6}"
  (let ((temp 300d0)
	(p 100d0)
	(lisp-unit:*epsilon* 1e-6))
    (let ((c 1e-6))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 1e-1))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 0.3))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 0.5))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 0.8))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 0.9))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))
    (let ((c 0.99999))
      (assert-number-equal (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)
			   (delta-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				    CP |(1-C)P| temp)))))

(define-test D12b-1-mcc/Ar-Ar
  "Test Ar/Ar diffusion calculated with the McCormack collision model
against the second order CCK for several relative Ar concentrations.
We correct the CCK with the Delta coefficient."
  (let ((lisp-unit:*epsilon* 1e-2))
    (let ((temp 300d0)
	  (p 100d0))
      (let ((d12 (d12-1 'cck +Ar-Ar/LJ6-12+ p temp)))
	(let ((C 0.5))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+
				       CP |(1-C)P|  temp)))
	(let ((C 0.1))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.001))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.9))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.999))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +Ar/LJ6-12+ CP |(1-C)P| temp)))))))

(define-test D12b-1-mcc/Ar-H2
  "Test Ar/H2 diffusion calculated with the McCormack collision model
against the second order CCK for several relative Ar concentrations.
We correct the CCK with the Delta coefficient."
  (let ((lisp-unit:*epsilon* 1e-2))
    (let ((temp 300d0)
	  (p 100d0))
      (let ((d12 (d12-1 'cck +Ar-H2/LJ6-12+ p temp)))
	(let ((C 0.5))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+
				       CP |(1-C)P|  temp)))
	(let ((C 0.1))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.001))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.9))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+ CP |(1-C)P| temp)))
	(let ((C 0.999))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +H2/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/LJ6-12+ +H2/LJ6-12+ CP |(1-C)P| temp)))))))




;;; Hard-sphere potential tests


(define-test D12-1-cck/hs
  "Compare the Ar/Ar D12 diffusion as calculated with the
Chapman-Cowling/Kihara approximation vs eq. 3.100 of Kee et al.

We use the hard-sphere potential"
  (let ((lisp-unit:*epsilon* 1e-2))
    (let ((temp 300d0)
	  (p 100d0))
      (let ((d12 (d12-kee-3.100 +Ar-H2/HS+ p temp)))
	(let ((C 0.5))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/LJ6-12+ +Ar/LJ6-12+
						  CP |(1-C)P| temp)))
			       (d12-1 'cck +Ar-H2/HS+ p temp)))))))

(define-test D12b-1-mcc/Ar-H2/HS
  "Test Ar/H2 diffusion calculated with the McCormack collision model
against the first order CCK for several relative Ar/H2 concentrations.

We use the hard-sphere potential.  The MCC diffusion disagrees against
the first order cck as the Ar fraction increases.

We test to about 1% accuracy.  Still most of the tests fail.
"
  (let ((lisp-unit:*epsilon* 1e-2))
    (let ((temp 300d0)
	  (p 100d0))
      (let ((d12 (d12-kee-3.100 +Ar-H2/HS+ p temp)))
	(let ((C 0.001))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/HS+ +H2/HS+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/HS+ +H2/HS+ CP |(1-C)P| temp)
			       "C=0.001"))
	(let ((C 0.1))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/HS+ +H2/HS+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/HS+ +H2/HS+ CP |(1-C)P| temp)
			       "C=0.1"))
	(let ((C 0.5))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/HS+ +H2/HS+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/HS+ +H2/HS+ CP |(1-C)P| temp)
			       "C=0.5"))
	(let ((C 0.9))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/HS+ +H2/HS+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/HS+ +H2/HS+ CP |(1-C)P| temp)
			       "C=0.9"))
	(let ((C 0.999))
	  (assert-number-equal (/ d12
				  (- 1d0 (delta-1 'cck +Ar/HS+ +H2/HS+
						  CP |(1-C)P| temp)))
			       (d12b-1 'mcc +Ar/HS+ +H2/HS+ CP |(1-C)P| temp)
			       "C=0.999"))))))


(defun D12-cck/mcc (s1 s2 p1 p2 temp)
  "Comparison of D12 calculated via CCK and MCC"
  (let ((c (make-collision-parameters s1 s2))
	(p (+ p1 p2)))
    (format t "D12-cck: ~a~%"(d12-1 'cck c temp p))
    (multiple-value-bind (d12-mcc d21-mcc)
	(d12b s1 s2 temp p1 p2)
      (format t "D12-mcc: ~a (~a+~a)~%" (+ d12-mcc d21-mcc)
	      d12-mcc d21-mcc))))