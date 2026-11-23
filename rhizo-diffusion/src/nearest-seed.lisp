(in-package :rhizo-diffusion)

;;; Nearest-seed Voronoi grid: the grid used for reaction-diffusion.
;;;
;;; Every pixel belongs to the seed closest to it, by exact Euclidean
;;; distance. This is a pixelated true Voronoi diagram: cell boundaries lie
;;; on the perpendicular bisectors between seeds, to within a pixel. That is
;;; what makes the two-point finite-volume flux consistent - the segment
;;; between two seeds is perpendicular to their shared face (PLAN.md, "Why
;;; the CA engine is the right engine").
;;;
;;; The CA-grown grid (voronoi-ca.lisp) stays as it is, for the isotropy
;;; study. Its alternating 4/8-neighbour growth makes boundaries slightly
;;; octagonal rather than true bisectors, which shows up as about 2% excess
;;; diffusivity (measured 2026-11-23, scripts/fv-accuracy.lisp).
;;;
;;; For the same SEED, MAKE-NEAREST-SEED-GRID and MAKE-VORONOI-GRID draw the
;;; same seed points (both call RANDOM-DISTINCT-POINTS first on a random
;;; state seeded with SEED), so the two grids can be compared cell for cell.
;;;
;;; Lookup is by square bins holding about one seed each, searched in rings
;;; around the pixel's own bin: O(pixels), not O(pixels x seeds).

(defun seed-bin-size (width height n-seeds)
  "Return the side, in pixels, of the square bins used to look up seeds:
chosen so that a WIDTH x HEIGHT grid holds about one of N-SEEDS seeds per
bin on average. Always at least 1."
  (max 1 (ceiling (sqrt (/ (* width height) (max 1 n-seeds))))))

(defun seed-bins (seed-points bin-size width height)
  "Return a fresh 2D array of (ceiling WIDTH BIN-SIZE) x (ceiling HEIGHT
BIN-SIZE) bins. Bin (i j) holds the list of indices, in increasing order, of
the seeds in SEED-POINTS (a vector of (x . y) integer pixels) that lie in
it."
  (let ((bins (make-array (list (ceiling width bin-size) (ceiling height bin-size))
                          :initial-element '())))
    (loop for k from (1- (length seed-points)) downto 0
          do (destructuring-bind (x . y) (aref seed-points k)
               (push k (aref bins (floor x bin-size) (floor y bin-size)))))
    bins))

(defun ring-bins (bx by r n-bx n-by)
  "Return the list of in-range bin coordinates (i . j) at Chebyshev
distance exactly R from bin (BX, BY), on an N-BX x N-BY array of bins."
  (loop for i from (max 0 (- bx r)) to (min (1- n-bx) (+ bx r))
        nconc (loop for j from (max 0 (- by r)) to (min (1- n-by) (+ by r))
                    when (= r (max (abs (- i bx)) (abs (- j by))))
                      collect (cons i j))))

(defun nearest-seed-index (x y seed-points bins bin-size)
  "Return the index of the seed in SEED-POINTS nearest to pixel (X, Y), by
exact integer squared distance; ties go to the lowest index. BINS and
BIN-SIZE come from SEED-BINS. Searches rings of bins of growing radius R
around the pixel's own bin. Any seed not yet searched lies in a ring beyond
R and so is at least R * BIN-SIZE + 1 pixels away, so the search stops as
soon as the best squared distance is at most (R * BIN-SIZE)^2."
  (let* ((bx (floor x bin-size)) (by (floor y bin-size))
         (n-bx (array-dimension bins 0)) (n-by (array-dimension bins 1))
         (max-r (max bx by (- n-bx 1 bx) (- n-by 1 by)))
         (best nil) (best-d2 0))
    (loop for r from 0 to max-r
          do (loop for (i . j) in (ring-bins bx by r n-bx n-by)
                   do (dolist (k (aref bins i j))
                        (destructuring-bind (sx . sy) (aref seed-points k)
                          (let ((d2 (+ (expt (- x sx) 2) (expt (- y sy) 2))))
                            (when (or (null best) (< d2 best-d2) (and (= d2 best-d2) (< k best)))
                              (setf best k best-d2 d2))))))
          until (and best (<= best-d2 (expt (* r bin-size) 2))))
    best))

(defun nearest-seed-ids (width height seed-points)
  "Return a fresh (array fixnum (WIDTH HEIGHT)) in which each pixel holds
the index of its nearest seed in SEED-POINTS (see NEAREST-SEED-INDEX)."
  (let* ((bin-size (seed-bin-size width height (length seed-points)))
         (bins (seed-bins seed-points bin-size width height))
         (ids (make-array (list width height) :element-type 'fixnum)))
    (dotimes (x width ids)
      (dotimes (y height)
        (setf (aref ids x y) (nearest-seed-index x y seed-points bins bin-size))))))

(defun make-nearest-seed-grid (width height n-points &key seed)
  "Return a VORONOI-GRID on a WIDTH x HEIGHT pixel grid with N-POINTS cells,
each pixel assigned to its nearest seed. SEED is the integer seed of this
call's own random state (generated if not given) and is recorded in
VORONOI-GRID-SEED. Each cell's CELL-POSITION is its seed, and FACE-DISTANCE
is seed-to-seed distance.

Face geometry is, for now, the CA grid's: 8-connected adjacency, EDGE-WEIGHT
1.0, staircase pixel-count edge lengths (BUILD-ADJACENCY, BUILD-GEOMETRY).
PLAN.md Phase 1 step 4 replaces the edge lengths."
  (let* ((seed (or seed (random most-positive-fixnum (make-random-state t))))
         (seed-points (coerce (random-distinct-points width height n-points
                                                      (sb-ext:seed-random-state seed))
                              'simple-vector))
         (ids (nearest-seed-ids width height seed-points))
         (positions (positions-table (alexandria:iota (length seed-points))
                                     (lambda (id) (aref seed-points id)))))
    (multiple-value-bind (cell-areas edge-lengths) (build-geometry ids width height)
      (%make-voronoi-grid
       :tessellation (make-fv-tessellation (build-adjacency ids width height)
                                           :cell-areas cell-areas
                                           :edge-lengths edge-lengths
                                           :cell-positions positions)
       :width width :height height :ids ids :seed seed :seed-points seed-points))))
