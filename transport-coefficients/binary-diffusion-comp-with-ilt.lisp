(in-package :tcu)

;; We compare binary diffusivity calculation with those of Ivchenko,
;; Loyalka & Tompson, ZAMP, 58-72 (2002) in tables I and II


(defun conversion (c m1 m2)
  (/ (expt (+ (* c m1)
	      (* (- 1 c) m2))
	   2)
     (* m1 m2)))

(defun d12-calc (g1 g2 p c temp)
  (let ((p1 (* c p))
	(p2 (* (- 1 c) p))
	(s1 (species g1))
	(s2 (species g2))
	(m1 (mass g1))
	(m2 (mass g2)))
    (let ((cck-1 (* 1e4 (d12-1 'tc:cck (make-collision-parameters g1 g2)
			p temp)))
	  (delta (delta-1 'tc:cck g1 g2 p1 p2 temp))
	  (mcc (* 1e4 (d12b-1 'tc:mcc g1 g2 p1 p2 temp)))
	  (conversion (conversion c m1 m2 )))
      (format t "~a-~a, ~10tCCK, ~15t1:~20t~5,4f~35t~5,4f~%" 
	      s1 s2 cck-1 (/ cck-1 conversion))
      (format t "~a-~a, ~10tCCK, ~15t2:~20t~5,4f~35t~5,4f~%" 
	      s1 s2 (/ cck-1
		       (- 1 delta))
	      (/ (/ cck-1
		    (- 1 delta))
		 conversion))
      (format t "~a-~a, ~10tMCC, ~15t1:~20t~5,4f~35t~5,4f~%" 
	      s1 s2 mcc (/ mcc conversion)))))

(let ((n2 (make-species-lennard-jones-6/12-potential :N2))
      (ar (make-species-lennard-jones-6/12-potential :Ar))
      (o2 (make-species-lennard-jones-6/12-potential :O2))
      (co2 (make-species-lennard-jones-6/12-potential :CO2))
      (h2 (make-species-lennard-jones-6/12-potential :H2))
      (p 101325.0)
      (temp 293d0))
  (d12-calc n2 h2 p 0.99 temp)
  (d12-calc n2 ar p 0.5 temp)
  (d12-calc n2 o2 p 0.5 temp)
  (d12-calc n2 co2 p 0.5 temp))
