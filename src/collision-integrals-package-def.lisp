;; Mirko Vukovic
;; Time-stamp: <2011-08-28 16:29:14 collision-integrals-package-def.lisp>
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


(defpackage :collision-integrals
  (:nicknames :omega-xx)
  (:use :cl  :lisp-unit)
  (:import-from :alexandria :with-input-from-file)
  (:import-from :my-utils :polyeval :^2)
  (:import-from :physics-constants :+boltzmann-constant-sp+)
  (:import-from :mv-grid+gsll :tabular-data :make-table
		:init-table-interp :interp-table)
  (:export :omega-11* :omega-22*)
  (:documentation "Formulas for collision integrals used in gas transport coefficient calculations"))

;; export is done in files in each module, either a setup file, or in
;; individual files

(in-package :omega-xx)
(define-symbol-macro +kb+ +boltzmann-constant-sp+)
(defconstant +pi+ (coerce pi 'single-float))



