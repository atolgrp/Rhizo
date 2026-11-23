(in-package :rhizo-diffusion)

;;; A WORLD holds everything one simulation instance needs, explicitly.
;;;
;;; The original notebook kept this in a set of global `defparameter`
;;; specials (*molecules*, *color-dict*, *neighbors-dict*, the shared
;;; global *random-state*, ...). That is exactly what let the two bugs
;;; found while testing it slip in unnoticed (a re-run silently re-used
;;; stale state), and it is exactly what makes running many replicates in
;;; parallel unsafe (two threads mutating the same global list). A WORLD
;;; is passed explicitly to every function that needs it and owns its own
;;; random-state, so two worlds - even running concurrently on different
;;; threads - can never interfere with each other.

(defstruct (world (:constructor %make-world))
  tessellation   ; anything satisfying the tessellation protocol
  molecules      ; hash-table (eql): cell-ids currently on the active front
  frontier       ; hash-table (eql): cell-ids the front has ever reached
  random-state   ; this world's own random-state - never the global one
  seed           ; integer seed used to create random-state, for reproducibility
  (tick 0 :type fixnum))

(defun make-world (tessellation &key seed)
  "Create a fresh, empty WORLD over TESSELLATION. SEED is an integer used
to seed this world's own random-state (generated and recorded if not
supplied) - the run is fully reproducible from WORLD-SEED, and independent
worlds never share random state."
  (let* ((seed (or seed (random most-positive-fixnum (make-random-state t))))
         (random-state (sb-ext:seed-random-state seed)))
    (%make-world :tessellation tessellation
                 :molecules (make-hash-table :test 'equal)
                 :frontier (make-hash-table :test 'equal)
                 :random-state random-state
                 :seed seed
                 :tick 0)))
