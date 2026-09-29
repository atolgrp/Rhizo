;;;; package.lisp --- Package definition for the Voronoi CA isotropy study.

(defpackage #:voronoi-ca
  (:use #:common-lisp)
  (:nicknames #:vca)
  (:export
   ;; configuration
   #:*width* #:*height* #:*n-seeds* #:*seed*
   #:seed-rng
   ;; simulation
   #:build-tessellation
   #:build-neighbour-table
   #:seed-front
   #:step-front
   ;; measurement
   #:measurement-ids
   #:run-trial
   ;; experiment driver
   #:run-experiment
   #:main
   ;; output
   #:write-tessellation-ppm))
