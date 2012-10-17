;; Mirko Vukovic
;; Time-stamp: <2012-07-18 13:52:47 binary-diffusion-plots.lisp>
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

(in-package #:tcu)

(define-plot d12/c-effects
  (let ((c-arr (lseq 1e-3 0.999))
	(temp 300d0)
	(p 100d0))
    (let ((p1 (gmap (lambda (c)
		      (* p c)) c-arr))
	  (p2 (gmap (lambda (c)
		      (* (- 1d0 c) p)) c-arr)))
      (let ((d12-hs-cck (d12-1 'tc:cck +Ar-H2/hs+ p temp))
	    (d12-lj-cck (d12-1 'tc:cck +Ar-H2/LJ6-12+ p temp))
	    (d12-hs-mcc (gpmap (d12b-1 'tc:mcc +Ar/hs+ +H2/hs+ @0p1 @1p2 temp) p1 p2 ))
	    (d12-lj-mcc (gpmap (d12b-1 'tc:mcc +Ar/lj6-12+ +H2/LJ6-12+ @0p1 @1p2 temp)
			       p1 p2 )))
	(set-to ((xlabel "Relative Ar concentration")
		 (ylabel "[m^2/s]")
		 (title "Concentration effects on binary diffusion coefficient in Ar/H2 mixture"))
	  (with-h-bar (d12-hs-cck :tag 1)
	    (with-h-bar (d12-lj-cck :tag 2)
	      (plot-xys c-arr
			(list (list d12-hs-mcc :title "Hard-sphere potential")
			      (list d12-lj-mcc :title "Lennard-Jones potential"))))))))))


(define-plot delta/c-effects
  (let ((c-arr (lseq 1e-3 0.999))
	(temp 300d0)
	(p 100d0))
    (let ((p1 (gmap (lambda (c)
		      (* p c)) c-arr))
	  (p2 (gmap (lambda (c)
		      (* (- 1d0 c) p)) c-arr)))
      (let ((delta-hs-cck (gpmap (delta-1 'tc:cck +Ar/hs+ +H2/hs+ @0p1 @1p2 temp ) p1 p2 ))
	    (delta-lj-cck (gpmap (delta-1 'tc:cck +Ar/lj6-12+ +H2/LJ6-12+ @0p1 @1p2 temp)
				 p1 p2 ))
	    (delta-hs-mcc (gpmap (delta-1 'tc:mcc +Ar/hs+ +H2/hs+ @0p1 @1p2 temp ) p1 p2 ))
	    (delta-lj-mcc (gpmap (delta-1 'tc:mcc +Ar/lj6-12+ +H2/LJ6-12+ @0p1 @1p2 temp )
				 p1 p2 )))
	(set-to ((xlabel "Relative Ar concentration")
		 (ylabel "[-]")
		 (title "Concentration effects on second order correction coefficient in Ar/H2 mixture"))
	  (plot-xys c-arr
		    (list (list delta-hs-cck :title "HS/CCK")
			  (list delta-lj-cck :title "LJ/CCK")
			  (list delta-hs-mcc :with :points :title "HS/MCC")
			  (list delta-lj-mcc :with :points :title "LJ/MCC"))))))))