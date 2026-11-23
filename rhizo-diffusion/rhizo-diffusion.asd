(asdf:defsystem "rhizo-diffusion"
  :description "Tessellation-agnostic graph diffusion engine (CA-grown Voronoi + front-propagation, and GIS-polygon + weighted Dijkstra diffusion), rewritten from the RhizoCA notebook with explicit per-simulation state so runs are parallel-safe and testable."
  :author "atol"
  :license "MIT"
  :depends-on ("alexandria" "bordeaux-threads" "lparallel" "csv-parser" "cl-json" "cl-heap" "zpng")
  :pathname "src/"
  :serial t
  :components ((:file "package")
               (:file "tessellation")
               (:file "voronoi-ca")
               (:file "nearest-seed")
               (:file "world")
               (:file "diffusion")
               (:file "metrics")
               (:file "io")
               (:file "parallel")
               (:file "geojson")
               (:file "gis-tessellation")
               (:file "weighted-diffusion")
               (:file "render")
               (:file "shape-metrics")
               (:file "validation")
               (:file "reaction-network")
               (:file "flux-stencil")
               (:file "mechanistic-world")
               (:file "mechanistic-step")
               (:file "accuracy"))
  :in-order-to ((asdf:test-op (asdf:test-op "rhizo-diffusion/tests"))))

(asdf:defsystem "rhizo-diffusion/tests"
  :description "FiveAM test suite for rhizo-diffusion."
  :depends-on ("rhizo-diffusion" "fiveam")
  :pathname "tests/"
  :serial t
  :components ((:file "package")
               (:file "tessellation-tests")
               (:file "nearest-seed-tests")
               (:file "world-tests")
               (:file "diffusion-tests")
               (:file "metrics-tests")
               (:file "parallel-tests")
               (:file "gis-tests")
               (:file "weighted-diffusion-tests")
               (:file "render-tests")
               (:file "shape-metrics-tests")
               (:file "validation-tests")
               (:file "reaction-network-tests")
               (:file "mechanistic-tests")
               (:file "flux-stencil-tests")
               (:file "accuracy-tests"))
  :perform (asdf:test-op (op c)
             (declare (ignore op c))
             (uiop:symbol-call :fiveam :run!
                                (uiop:find-symbol* :rhizo-diffusion-suite :rhizo-diffusion.tests))))
