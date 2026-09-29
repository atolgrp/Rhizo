;;;; voronoi-ca.asd --- System definition.

(asdf:defsystem #:voronoi-ca
  :description "Does a Voronoi tessellation suppress lattice anisotropy in a
cellular-automaton diffusion front? A controlled numerical experiment."
  :author "Apostolos Almpanis"
  :license "MIT"
  :version "1.0.0"
  :serial t
  :components ((:module "src"
                :serial t
                :components ((:file "package")
                             (:file "utils")
                             (:file "grid")
                             (:file "voronoi")
                             (:file "propagation")
                             (:file "ppm")
                             (:file "measurement")
                             (:file "experiment")))))
