(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(test id->color-is-deterministic-and-varies
  (is (equal (id->color 42) (id->color 42)))
  (is (equal (id->color "A") (id->color "A")))
  (is (not (equal (id->color 1) (id->color 2)))))

(test render-voronoi-grid-png-writes-a-file
  (let* ((grid (make-voronoi-grid 20 20 15 :seed 1))
         (path (merge-pathnames "test-voronoi.png" (uiop:temporary-directory))))
    (render-voronoi-grid-png grid path)
    (is (probe-file path))
    (is (plusp (with-open-file (s path :element-type '(unsigned-byte 8)) (file-length s))))
    (delete-file path)))

(test render-front-propagation-png-writes-a-file
  (let* ((grid (make-voronoi-grid 20 20 15 :seed 2))
         (tess (voronoi-grid-tessellation grid))
         (start (first (cell-ids tess)))
         (world (make-world tess :seed 2))
         (path (merge-pathnames "test-front.png" (uiop:temporary-directory))))
    (seed-molecules-at! world (list start))
    (step-front-propagation! world)
    (render-front-propagation-png grid world path :source-ids (list start))
    (is (probe-file path))
    (delete-file path)))
