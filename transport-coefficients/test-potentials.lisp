;; Mirko Vukovic
;; Time-stamp: <2012-07-17 16:44:30 test-potentials.lisp>
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

(export '(+Ar/LJ6-12+ +Ar-Ar/LJ6-12+ +H2/LJ6-12+ +Ar-H2/LJ6-12+
	  +Ar/HS+ +H2/HS+ +Ar-H2/HS+))


(defparameter +Ar/LJ6-12+ (make-species-lennard-jones-6/12-potential :Ar))
(defparameter +Ar-Ar/LJ6-12+ (make-collision-parameters +Ar/LJ6-12+ +Ar/LJ6-12+))
(defparameter +H2/LJ6-12+ (make-species-lennard-jones-6/12-potential :H2))
(defparameter +Ar-H2/LJ6-12+ (make-collision-parameters +Ar/LJ6-12+ +H2/LJ6-12+))

(defparameter +Ar/HS+ (make-hard-sphere-potential (mass +Ar/LJ6-12+) (sigma +Ar/LJ6-12+) :Ar))
(defparameter +H2/HS+ (make-hard-sphere-potential (mass +H2/LJ6-12+) (sigma +H2/LJ6-12+) :H2))
(defparameter +Ar-H2/HS+ (make-collision-parameters +Ar/HS+ +H2/HS+))


