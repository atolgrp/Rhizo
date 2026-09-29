;;;; experiment.lisp --- The experiment driver.
;;;;
;;;; Two conditions on one code path:
;;;;
;;;;   voronoi  -- n_seeds << width*height, so regions span many lattice cells
;;;;   lattice  -- n_seeds = width*height, one seed per cell, which degenerates
;;;;               exactly to the square lattice
;;;;
;;;; Using the degenerate tessellation as the control means the two conditions
;;;; differ only in seed density; the automaton, the measurement code and the
;;;; data path are identical, so nothing else can explain a difference.

(in-package #:voronoi-ca)

(defparameter *seed* 20210728
  "Base RNG seed. Trial N uses *SEED* + N, so every run is reproducible.")

(defparameter *output-directory* #p"results/"
  "Where CSV and image output is written.")

(defun ensure-output-directory ()
  (ensure-directories-exist *output-directory*))

(defun condition-seed-count (condition width height)
  "Seed count for CONDITION, which is :VORONOI or :LATTICE."
  (ecase condition
    (:lattice (* width height))
    (:voronoi (round (* width height) 10))))

(defun run-experiment (&key (width 200) (height 200) (replicates 5)
                            (conditions '(:lattice :voronoi))
                            (seeds nil)
                            (image t)
                            (verbose t))
  "Run REPLICATES trials of each condition and write one CSV per condition.

SEEDS overrides the per-condition seed count with an explicit alist, e.g.
'((:voronoi . 4000) (:lattice . 40000)). Returns the list of CSV paths."
  (ensure-output-directory)
  (let ((paths '()))
    (dolist (condition conditions (nreverse paths))
      (let* ((n-seeds (or (cdr (assoc condition seeds))
                          (condition-seed-count condition width height)))
             (path (merge-pathnames
                    (format nil "~(~a~)-~dx~d.csv" condition width height)
                    *output-directory*)))
        (when verbose
          (format t "~&~a  ~dx~d, ~d seeds, ~d replicates~%"
                  (string-downcase condition) width height n-seeds replicates)
          (finish-output))
        (with-open-file (stream path :direction :output
                                     :if-exists :supersede
                                     :if-does-not-exist :create)
          (write-csv-row stream (append '("replicate" "seed" "regions")
                                        +probe-names+))
          (dotimes (r replicates)
            (let* ((trial-seed (+ *seed* r))
                   ;; one illustrative snapshot per condition, from trial 0,
                   ;; taken inside the run at the moment the front reaches ring 1
                   (snapshot (when (and image (zerop r))
                               (merge-pathnames
                                (format nil "~(~a~)-~dx~d.ppm" condition width height)
                                *output-directory*))))
              (seed-rng trial-seed)
              (multiple-value-bind (times regions)
                  (run-trial :width width :height height :n-seeds n-seeds
                             :snapshot snapshot
                             :snapshot-scale (max 1 (floor 600 width)))
                (write-csv-row stream
                               (append (list r trial-seed regions) times))
                (finish-output stream)
                (when verbose
                  (format t "  replicate ~d/~d  (~d regions)  ~{~a~^ ~}~%"
                          (1+ r) replicates regions times)
                  (finish-output))))))
        (push path paths)))))

(defun main ()
  "Command-line entry point.

Usage: sbcl --script run.lisp [width] [replicates]"
  (let* ((args (rest #+sbcl sb-ext:*posix-argv* #-sbcl '()))
         (width (if (first args) (parse-integer (first args)) 200))
         (replicates (if (second args) (parse-integer (second args)) 5)))
    (format t "~&Voronoi CA isotropy experiment~%")
    (format t "lattice ~dx~d, ~d replicates per condition, base seed ~d~%~%"
            width width replicates *seed*)
    (let ((start (get-internal-real-time)))
      (run-experiment :width width :height width :replicates replicates)
      (format t "~&~%done in ~,1f s~%"
              (/ (- (get-internal-real-time) start)
                 internal-time-units-per-second)))
    (format t "CSV written to ~a~%" (namestring *output-directory*))
    (format t "Now run:  python3 analysis/analyse.py~%")))
