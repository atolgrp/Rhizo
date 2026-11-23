(in-package :rhizo-diffusion)

;;; Front-propagation diffusion over any tessellation.
;;;
;;; This is the diffusion mode the original notebook actually used for its
;;; isotropy measurement (update-next-voronoi-generation-front-prop /
;;; prop-time). The old code also had a second, random-walk mode
;;; ("molecules hopping between sources") that was never part of the
;;; validated isotropy result - it is deliberately not ported here to keep
;;; this rewrite scoped to what's proven and needed.
;;;
;;; Note this file never mentions grids, points, or anything CA-specific -
;;; only CELL-IDS/NEIGHBORS-OF from the tessellation protocol. It works
;;; unchanged over a CA-grown VORONOI-GRID today and over any future
;;; tessellation built a different way.

(defun seed-molecules-at! (world cell-ids)
  "Start WORLD's front at CELL-IDS: mark them active and already-visited,
so the front doesn't immediately spread back onto its own starting point."
  (dolist (id cell-ids)
    (setf (gethash id (world-molecules world)) t)
    (setf (gethash id (world-frontier world)) t))
  world)

(defun front-active-p (world cell-id)
  "True if CELL-ID has been reached by WORLD's diffusing front by now."
  (and (gethash cell-id (world-frontier world)) t))

(defun step-front-propagation! (world)
  "Advance WORLD's front by one tick. Every cell on the current front
spreads to each neighbor not yet visited; those neighbors become the new
front. A cell spreads exactly once - FRONTIER (never cleared) is what
stops it from spreading back the way it came."
  (let ((tess (world-tessellation world))
        (next-front (make-hash-table :test 'equal)))
    (loop for id being the hash-keys of (world-molecules world) do
      (dolist (nb (neighbors-of tess id))
        (unless (gethash nb (world-frontier world))
          (setf (gethash nb next-front) t)
          (setf (gethash nb (world-frontier world)) t))))
    (setf (world-molecules world) next-front)
    (incf (world-tick world))
    world))
