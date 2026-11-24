(in-package :rhizo-diffusion)

;;; The signalling graph (PLAN.md, Phase 4).
;;;
;;; Nodes are microbial patches (and, later, root segments). There is an
;;; edge from i to j when signal released at i reaches j above j's
;;; detection threshold before it degrades: the graph is directed, since
;;; patches differ in strength and threshold. A reach predicate - a function
;;; of the two node indices and their distance - decides each edge, so the
;;; same graph code serves the analytic calling distance (uniform water) and
;;; fields measured in a simulation.
;;;
;;; Connectivity is reported both ways: weakly connected components (signal
;;; can pass in at least one direction along every link) and strongly
;;; connected ones (every member can relay to every other). A component
;;; spans the domain when it reaches from the left boundary strip to the
;;; right one.
;;;
;;; With equal patches in uniform water the graph is a random geometric
;;; graph of connection radius r_c, the calling distance. Its percolation
;;; point is known: in the plane, discs of radius r_c / 2 percolate at
;;; filling factor eta_c = 1.12808737(6) (Mertens & Moore 2012, Phys Rev E
;;; 86:061109, doi:10.1103/PhysRevE.86.061109), i.e. at mean degree
;;; n pi r_c^2 = 4 eta_c = 4.512.

(defparameter *disc-percolation-filling-factor* 1.12808737d0
  "Critical filling factor eta_c of fully penetrable discs in the plane
\(Mertens & Moore 2012, Phys Rev E 86:061109). A constant of geometry, not
simulation state.")

(defun random-patch-positions (n width height &key seed)
  "Return N points (x . y), double-floats uniform on [0, WIDTH] x [0,
HEIGHT], drawn from a random state seeded with SEED - never the global
one, so the same SEED gives the same patches."
  (let ((state (sb-ext:seed-random-state seed)))
    (loop repeat n
          collect (cons (random (float width 1d0) state) (random (float height 1d0) state)))))

(defun signalling-edges (positions reach)
  "Return the directed edges (i . j), i /= j, among the nodes at POSITIONS
\(a list of points, node i at position i) for which (FUNCALL REACH i j r) is
true, r the distance between them."
  (let ((points (coerce positions 'simple-vector)))
    (loop for i below (length points)
          nconc (loop for j below (length points)
                      when (and (/= i j)
                                (funcall reach i j (point-distance (aref points i) (aref points j))))
                        collect (cons i j)))))

(defun group-components (n representative)
  "Return the components of N nodes as lists of node indices in increasing
order, sorted by their first member; node i belongs to the component
labelled (FUNCALL REPRESENTATIVE i)."
  (let ((groups (make-hash-table)))
    (loop for i from (1- n) downto 0
          do (push i (gethash (funcall representative i) groups)))
    (sort (loop for members being the hash-values of groups collect members) #'< :key #'first)))

(defun weak-components (n edges)
  "Return the weakly connected components of the directed graph on N nodes
with EDGES (pairs (i . j)), by union-find; see GROUP-COMPONENTS for the
order."
  (let ((parent (make-array n)))
    (dotimes (i n) (setf (aref parent i) i))
    (labels ((root (i) (if (= (aref parent i) i)
                           i
                           (setf (aref parent i) (root (aref parent i))))))
      (loop for (i . j) in edges
            do (setf (aref parent (root i)) (root j)))
      (group-components n #'root))))

(defun strong-components (n edges)
  "Return the strongly connected components of the directed graph on N
nodes with EDGES, by Kosaraju's algorithm; see GROUP-COMPONENTS for the
order."
  (let ((out (make-array n :initial-element '()))
        (in (make-array n :initial-element '()))
        (seen (make-array n :initial-element nil))
        (label (make-array n :initial-element nil))
        (order '()))
    (loop for (i . j) in edges
          do (push j (aref out i)) (push i (aref in j)))
    (labels ((visit (i)
               (unless (aref seen i)
                 (setf (aref seen i) t)
                 (mapc #'visit (aref out i))
                 (push i order)))
             (assign (i root)
               (unless (aref label i)
                 (setf (aref label i) root)
                 (dolist (j (aref in i)) (assign j root)))))
      (dotimes (i n) (visit i))
      (dolist (i order) (assign i i))
      (group-components n (lambda (i) (aref label i))))))

(defun largest-component-size (components)
  "Return the number of nodes in the largest of COMPONENTS (0 if none)."
  (reduce #'max components :key #'length :initial-value 0))

(defun spanning-p (component positions x-min x-max margin)
  "Return true if COMPONENT (node indices into POSITIONS) has a node within
MARGIN of the left boundary X-MIN and one within MARGIN of the right
boundary X-MAX."
  (let ((xs (mapcar (lambda (i) (car (nth i positions))) component)))
    (and (some (lambda (x) (<= x (+ x-min margin))) xs)
         (some (lambda (x) (>= x (- x-max margin))) xs))))

(defun critical-connection-radius (n area)
  "Return the connection radius at which a random geometric graph of N
nodes on AREA percolates: the r with mean degree (N / AREA) pi r^2 = 4
eta_c (*DISC-PERCOLATION-FILLING-FACTOR*), for a large domain."
  (sqrt (/ (* 4 *disc-percolation-filling-factor* area) (* pi n))))

;;; Percolation against water content (PLAN.md, Phase 4 gate). With equal
;;; patches in uniform water the connection radius is the calling distance
;;; r_c(theta) (RANGE-VERSUS-WATER), so the graph percolates once r_c
;;; reaches the percolation radius. Measured two ways: on finite random
;;; layouts (PERCOLATION-WATER-CONTENT, one per seed), and predicted for an
;;; infinite domain from the continuum threshold
;;; (PREDICTED-PERCOLATION-WATER-CONTENT).

(defun percolation-at-range (n width height r-connect &key seed)
  "Return the connectivity of N equal patches placed by RANDOM-PATCH-POSITIONS
\(SEED) on a WIDTH x HEIGHT domain, each reaching every other patch closer
than R-CONNECT, as a plist :LARGEST-WEAK-FRACTION and
:LARGEST-STRONG-FRACTION (largest component over N), :SPANNING (a weak
component reaches within R-CONNECT of both the left and right edges) and
:MEAN-DEGREE (edges per node)."
  (let* ((positions (random-patch-positions n width height :seed seed))
         (edges (signalling-edges positions (lambda (i j r) (declare (ignore i j)) (< r r-connect))))
         (weak (weak-components n edges)))
    (list :largest-weak-fraction (/ (largest-component-size weak) n)
          :largest-strong-fraction (/ (largest-component-size (strong-components n edges)) n)
          :spanning (if (some (lambda (c) (spanning-p c positions 0d0 width r-connect)) weak) t nil)
          :mean-degree (/ (length edges) n))))

(defun range-at-water (parameters theta scenario measure)
  "Return the calling distance MEASURE (:STEADY or :TRANSIENT) at water
content THETA (RANGE-VERSUS-WATER)."
  (getf (first (range-versus-water parameters (list theta) :scenario scenario)) measure))

(defun predicted-percolation-water-content (parameters n area &key scenario (measure :steady)
                                                                (theta-dry 0.02d0)
                                                                (theta-wet 0.45d0))
  "Return the water content at which the calling distance MEASURE of equal
patches in SCENARIO reaches the CRITICAL-CONNECTION-RADIUS of N patches on
AREA: the percolation point of an infinite domain at that density. By
bisection between THETA-DRY and THETA-WET, over which the range must rise
through the critical radius."
  (let ((r-critical (critical-connection-radius n area)))
    (bisect-root (lambda (theta) (- (range-at-water parameters theta scenario measure) r-critical))
                 theta-dry theta-wet)))

(defun percolation-water-content (parameters &key n width height seed scenario (measure :steady)
                                               (theta-dry 0.02d0) (theta-wet 0.45d0)
                                               (tolerance 1d-6))
  "Return the water content at which the layout of N patches drawn with
SEED on WIDTH x HEIGHT first spans (PERCOLATION-AT-RANGE with the calling
distance MEASURE as connection radius), to within TOLERANCE, by bisection
between THETA-DRY (must not span) and THETA-WET (must span). Signals an
error otherwise. The range must rise with theta between them."
  (flet ((spans (theta)
           (getf (percolation-at-range n width height
                                       (range-at-water parameters theta scenario measure)
                                       :seed seed)
                 :spanning)))
    (unless (and (not (spans theta-dry)) (spans theta-wet))
      (error "PERCOLATION-WATER-CONTENT: theta ~a and ~a do not bracket spanning (seed ~a)"
             theta-dry theta-wet seed))
    (loop with lo = theta-dry and hi = theta-wet
          while (> (- hi lo) tolerance)
          do (let ((mid (/ (+ lo hi) 2d0)))
               (if (spans mid) (setf hi mid) (setf lo mid)))
          finally (return (/ (+ lo hi) 2d0)))))
