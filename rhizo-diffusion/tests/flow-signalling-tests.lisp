(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Flow-driven water content and directional signalling (PLAN.md, Phase 5).

(defun um-loam ()
  (van-genuchten-soil (sandy-clay-loam-parameters) :length-scale 1d4))

(defun tracer-flow-world (tess theta0 c0 &key (dt 10d0))
  "Test helper: a world with :WATER and one dissolved tracer at pore-water
concentration C0 in water content THETA0, no reactions."
  (let ((world (make-mechanistic-world tess :species (list (water-species)
                                                           (make-species :name :tracer :diffusivity 490d0))
                                            :reactions '() :dt dt :diffusion-scheme :imex
                                            :initial-state `((:water . ,theta0) (:tracer . ,(* theta0 c0))))))
    (set-water-content! world (constantly theta0) (no-tortuosity))))

(test rain-dilutes-a-tracer-and-water-carries-it-out
  (let* ((soil (um-loam))
         (tess (make-rectangular-tessellation 1 10 :dx 50d0 :dy 50d0))
         (h0 -1d5)
         (theta0 (vg-water-content soil h0))
         (world (tracer-flow-world tess theta0 2d0))
         (problem (make-richards-problem tess :soil soil :cell-index (mechanistic-world-cell-index world)
                                              :boundaries '((:top :flux 1d0) (:bottom :free-drainage))))
         (state (make-richards-state problem h0))
         (balance0 (pool-balance world '(:tracer) '()))
         (end (run-with-flow! world problem state 200 (no-tortuosity)))
         (index (mechanistic-world-cell-index world))
         (concentration (lambda (id) (let ((i (cell-id->index index id)))
                                       (/ (aref (mechanistic-world-state world) i 1)
                                          (aref (mechanistic-world-state world) i 0))))))
    (is (= 2000d0 (richards-state-time end)))
    (is (< (abs (- (pool-balance world '(:tracer) '()) balance0)) (* 1d-10 balance0))
        "tracer in the soil plus tracer carried out is conserved")
    (is (plusp (species-exported world :tracer)) "drainage carries tracer out of the bottom")
    (is (< (funcall concentration 0) (funcall concentration 9)) "clean rain dilutes from the top")
    (is (every (lambda (id) (<= (funcall concentration id) (* 2d0 (+ 1 1d-12)))) (cell-ids tess))
        "no new maximum")
    (is (every (lambda (id) (= (aref (mechanistic-world-state world) (cell-id->index index id) 0)
                               (aref (richards-state-water-contents end) (cell-id->index index id))))
               (cell-ids tess))
        "the world's water is the flow solution's")))

(test still-water-transport-is-plain-diffusion
  ;; hydrostatic heads (h + z constant) and no boundary flow: no water moves
  (let* ((soil (um-loam))
         (tess (make-rectangular-tessellation 3 3 :dx 40d0 :dy 40d0))
         (theta0 (vg-water-content soil -2d5))
         (flowing (tracer-flow-world tess theta0 0d0))
         (still (tracer-flow-world tess theta0 0d0))
         (index (mechanistic-world-cell-index flowing)))
    (dolist (world (list flowing still))
      (setf (aref (mechanistic-world-state world) (cell-id->index index 4) 1) 5d0))
    (let* ((problem (make-richards-problem tess :soil soil :cell-index index))
           (state (make-richards-state problem (lambda (id) (+ -2d5 (cdr (cell-position tess id)))))))
      ;; the hydrostatic heads vary with depth, and so does theta: give the
      ;; comparison world the same water content
      (set-water-content! still (lambda (id) (vg-water-content soil (+ -2d5 (cdr (cell-position tess id)))))
                          (no-tortuosity))
      (run-with-flow! flowing problem state 20 (no-tortuosity))
      (run-mechanistic! still 20)
      (let ((difference (reduce #'max (map 'list (lambda (a b) (abs (- a b)))
                                           (state-column (mechanistic-world-state flowing) 1)
                                           (state-column (mechanistic-world-state still) 1)))))
        (is (< difference 1d-9) "largest difference ~a" difference)))))

(test advective-steady-source-solves-advection-diffusion-decay
  ;; D lap c - q dc/dx - k c = 0 away from the source, checked by finite
  ;; differences; at q = 0 the K0 profile
  (let ((s 3d0) (d 50d0) (k 0.02d0) (q 1.5d0) (h 0.01d0))
    (flet ((c (x y) (advective-point-source-steady-state s d k q x y)))
      (dolist (p '((30d0 . 10d0) (-25d0 . 5d0) (5d0 . -40d0)))
        (destructuring-bind (x . y) p
          (let* ((c0 (c x y))
                 (lap (/ (+ (c (+ x h) y) (c (- x h) y) (c x (+ y h)) (c x (- y h)) (* -4 c0)) (* h h)))
                 (dcdx (/ (- (c (+ x h) y) (c (- x h) y)) (* 2 h)))
                 (residual (- (* d lap) (* q dcdx) (* k c0))))
            (is (< (abs residual) (* 1d-4 (* k c0))) "at ~a: residual ~a vs k c ~a" p residual (* k c0))))))
    (is (< (abs (- (advective-point-source-steady-state s d k 0d0 12d0 5d0)
                   (point-source-steady-state s d k 13d0)))
           1d-14))
    (is (> (advective-point-source-steady-state s d k q 20d0 0d0)
           (advective-point-source-steady-state s d k q 0d0 20d0)
           (advective-point-source-steady-state s d k q -20d0 0d0))
        "downstream > across > upstream")))

(test advective-transient-tends-to-the-steady-state
  (let ((s 3d0) (d 50d0) (k 0.02d0) (q 1.5d0) (capacity 0.4d0))
    (is (< (abs (- (advective-point-source-transient s d k capacity 0d0 12d0 5d0 3600d0)
                   (point-source-transient s d k capacity 13d0 3600d0)))
           1d-15)
        "no flow: the radial transient")
    (let ((steady (advective-point-source-steady-state s d k q 30d0 10d0))
          (late (advective-point-source-transient s d k capacity q 30d0 10d0 1d5)))
      (is (< (abs (- late steady)) (* 1d-3 steady))))))

(test directional-calling-distance-is-longest-downstream
  (let ((s 1d3) (d 50d0) (k 0.002d0) (threshold 1d0))
    (flet ((range (q angle)
             (directional-calling-distance
              (lambda (along across) (advective-point-source-steady-state s d k q along across))
              threshold angle (sqrt (/ d k)))))
      (let ((still (steady-calling-distance s d k threshold)))
        (dolist (angle (list 0d0 (/ pi 2) pi))
          (is (< (abs (- (range 0d0 angle) still)) (* 1d-9 still)))))
      (is (> (range 1d0 0d0) (range 1d0 (/ pi 2)) (range 1d0 pi))))))

(test flow-edges-are-directional
  (let* ((parameters (with-parameter-values (ahl-literature-parameters) '((:water-content . 0.3d0))))
         (positions (random-patch-positions 80 1000d0 1000d0 :seed 4))
         (scenario '(:degrader-density 1d9))
         (still (flow-signalling-edges positions parameters 0d0 scenario :steady))
         (flowing (flow-signalling-edges positions parameters 2d0 scenario :steady)))
    (is (plusp (length still)))
    (is (zerop (one-way-fraction still)) "no flow: every link goes both ways")
    (is (plusp (one-way-fraction flowing)) "flow: some links go one way only")
    ;; flow is along +y: a one-way link points downstream
    (dolist (e flowing)
      (unless (member (cons (cdr e) (car e)) flowing :test #'equal)
        (is (> (cdr (nth (cdr e) positions)) (cdr (nth (car e) positions))))))))

(test flow-graph-connectivity-reports-both-axes
  (let* ((positions '((0d0 . 50d0) (40d0 . 50d0) (95d0 . 50d0) (50d0 . 0d0) (50d0 . 95d0)))
         (edges '((0 . 1) (1 . 2) (2 . 1) (3 . 1) (1 . 4)))
         (c (flow-graph-connectivity positions edges 100d0 100d0 10d0)))
    (is (getf c :spanning-across) "0-1-2 links left and right")
    (is (getf c :spanning-along) "3-1-4 links top and bottom")
    (is (= 5/5 (getf c :largest-weak-fraction)))
    (is (= 2/5 (getf c :largest-strong-fraction)) "only 1 <-> 2 goes both ways")
    (is (= 3/5 (getf c :one-way-fraction)))
    (is (getf c :spanning-along-directed) "a directed path from the top to the bottom: 3 -> 1 -> 4")))
