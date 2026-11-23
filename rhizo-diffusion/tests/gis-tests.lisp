(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Synthetic fixture: a 2x3 grid of unevenly-sized rectangles.
;;;
;;;   y=0 +---+-------+---+  y=3
;;;       | A |   B   | C |
;;;   y=2 +---+-------+---+
;;;       | D |   E   | F |
;;;   y=... +---+-------+---+
;;;
;;; (A,D,C,F are 1 wide; B,E are 3 wide; top row 2 tall, bottom row 1 tall)
;;;
;;; Ground truth, worked out by hand from the coordinates in the fixture:
;;;   adjacent (share a full edge):     A-B B-C A-D B-E C-F D-E E-F
;;;   NOT adjacent (corner-touch only): A-E B-D B-F C-E A-C A-F C-D D-F
;;;   centroids: A(0.5,1) B(2.5,1) C(4.5,1) D(0.5,2.5) E(2.5,2.5) F(4.5,2.5)

(defparameter *fixture-path*
  (asdf:system-relative-pathname "rhizo-diffusion" "tests/fixtures/synthetic-regions.geojson"))

(defun region-by-letter (regions letter)
  (find letter regions :key #'region-id :test #'string=))

(test parse-geojson-regions-basic
  (let ((regions (parse-geojson-regions *fixture-path*)))
    (is (= 6 (length regions)))
    (is (set-equal '("A" "B" "C" "D" "E" "F") (mapcar #'region-id regions) :test #'string=))))

(test parse-geojson-regions-centroids
  (let ((regions (parse-geojson-regions *fixture-path*)))
    (flet ((centroid= (letter x y)
             (let ((c (region-centroid (region-by-letter regions letter))))
               (and (< (abs (- (car c) x)) 1d-9) (< (abs (- (cdr c) y)) 1d-9)))))
      (is (centroid= "A" 0.5d0 1.0d0))
      (is (centroid= "B" 2.5d0 1.0d0))
      (is (centroid= "C" 4.5d0 1.0d0))
      (is (centroid= "D" 0.5d0 2.5d0))
      (is (centroid= "E" 2.5d0 2.5d0))
      (is (centroid= "F" 4.5d0 2.5d0)))))

(defun fixture-tessellation-and-ids ()
  (let* ((regions (parse-geojson-regions *fixture-path*))
         (tess (make-gis-tessellation regions)))
    (values tess (loop for l in '("A" "B" "C" "D" "E" "F")
                        collect (cons l (region-id (region-by-letter regions l)))))))

(test gis-adjacency-matches-ground-truth
  (multiple-value-bind (tess ids) (fixture-tessellation-and-ids)
    (flet ((adjacent-p (l1 l2) (member (cdr (assoc l2 ids :test #'string=))
                                       (neighbors-of tess (cdr (assoc l1 ids :test #'string=)))
                                       :test #'equal)))
      (dolist (pair '(("A" "B") ("B" "C") ("A" "D") ("B" "E") ("C" "F") ("D" "E") ("E" "F")))
        (is (adjacent-p (first pair) (second pair))
            "~a and ~a should be adjacent (share a full edge)" (first pair) (second pair)))
      (dolist (pair '(("A" "E") ("B" "D") ("B" "F") ("C" "E") ("A" "C") ("A" "F") ("C" "D") ("D" "F")))
        (is (not (adjacent-p (first pair) (second pair)))
            "~a and ~a should NOT be adjacent (corner-touch at most)" (first pair) (second pair))))))

(test gis-edge-weight-is-centroid-distance
  (multiple-value-bind (tess ids) (fixture-tessellation-and-ids)
    (flet ((id (l) (cdr (assoc l ids :test #'string=))))
      (is (< (abs (- (edge-weight tess (id "A") (id "B")) 2.0d0)) 1d-9))
      (is (< (abs (- (edge-weight tess (id "A") (id "D")) 1.5d0)) 1d-9)))))

(test gis-face-distance-is-centroid-distance
  ;; Centroids A(0.5,1) B(2.5,1) D(0.5,2.5) E(2.5,2.5) - see fixture header.
  (multiple-value-bind (tess ids) (fixture-tessellation-and-ids)
    (flet ((id (l) (cdr (assoc l ids :test #'string=))))
      (is (< (abs (- (face-distance tess (id "A") (id "B")) 2.0d0)) 1d-9))
      (is (< (abs (- (face-distance tess (id "B") (id "A")) 2.0d0)) 1d-9))
      (is (< (abs (- (face-distance tess (id "A") (id "D")) 1.5d0)) 1d-9))
      (is (< (abs (- (face-distance tess (id "D") (id "E")) 2.0d0)) 1d-9)))))

(test tessellation-from-geojson-end-to-end
  (multiple-value-bind (tess regions) (tessellation-from-geojson *fixture-path*)
    (is (= 6 (cell-count tess)))
    (is (= 6 (length regions)))))

;;; Ground truth for the geometry tests below, worked out by hand from the
;;; same fixture coordinates as the adjacency ground truth above:
;;;   areas:       A=2 B=6 C=2 D=1 E=3 F=1  (each region's width * height)
;;;   edge lengths: A-B=2 B-C=2 A-D=1 B-E=3 C-F=1 D-E=1 E-F=1

(test gis-cell-area-matches-ground-truth
  (multiple-value-bind (tess ids) (fixture-tessellation-and-ids)
    (flet ((id (l) (cdr (assoc l ids :test #'string=))))
      (dolist (expected '(("A" . 2.0d0) ("B" . 6.0d0) ("C" . 2.0d0)
                           ("D" . 1.0d0) ("E" . 3.0d0) ("F" . 1.0d0)))
        (is (< (abs (- (cell-area tess (id (car expected))) (cdr expected))) 1d-9)
            "~a: expected area ~a, got ~a" (car expected) (cdr expected)
            (cell-area tess (id (car expected))))))))

(test gis-edge-length-matches-ground-truth
  (multiple-value-bind (tess ids) (fixture-tessellation-and-ids)
    (flet ((id (l) (cdr (assoc l ids :test #'string=))))
      (dolist (expected '(("A" "B" . 2.0d0) ("B" "C" . 2.0d0) ("A" "D" . 1.0d0)
                           ("B" "E" . 3.0d0) ("C" "F" . 1.0d0) ("D" "E" . 1.0d0) ("E" "F" . 1.0d0)))
        (destructuring-bind (l1 l2 . expected-length) expected
          (is (< (abs (- (edge-length tess (id l1) (id l2)) expected-length)) 1d-9)
              "~a-~a: expected length ~a, got ~a" l1 l2 expected-length
              (edge-length tess (id l1) (id l2)))
          (is (< (abs (- (edge-length tess (id l2) (id l1)) expected-length)) 1d-9)
              "~a-~a should be symmetric" l2 l1))))))

(test lookups-work-with-a-freshly-constructed-equal-string-id
  ;; String cell ids come from GeoJSON properties, so a caller will
  ;; naturally look things up with their OWN string literal - a different
  ;; object from whatever cl-json originally allocated while parsing, but
  ;; with the same characters. Every adjacency table/lookup here must be
  ;; keyed by content (:TEST 'EQUAL), never object identity (EQL) - this
  ;; caught a real bug where NEIGHBORS-OF/EDGE-WEIGHT/adjacency tables
  ;; used the default EQL and silently returned nothing for a fresh string.
  (let* ((tess (tessellation-from-geojson *fixture-path*))
         (fresh-a (copy-seq "A")) ; guaranteed a different string object than any "A" cl-json produced
         (fresh-b (copy-seq "B")))
    (is (not (eq fresh-a (region-id (region-by-letter (parse-geojson-regions *fixture-path*) "A")))))
    (is (= 2 (length (neighbors-of tess fresh-a))))
    (is (member fresh-b (neighbors-of tess fresh-a) :test #'string=))
    (is (< (abs (- (edge-weight tess fresh-a fresh-b) 2.0d0)) 1d-9))))
