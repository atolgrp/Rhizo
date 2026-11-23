(in-package :rhizo-diffusion)

;;; Replicate runner: N independent simulation runs in parallel.
;;;
;;; This is the direct replacement for the old (dotimes (i 10) (prop-check
;;; 301 301 9000)) sweep - except that loop was sequential (and, we found,
;;; buggy: point-density argument silently ignored, setup work redundantly
;;; doubled). Because every replicate builds and owns its own WORLD/
;;; VORONOI-GRID with its own random-state (see world.lisp, voronoi-ca.lisp),
;;; there is no state shared between replicates, so nothing needs a lock -
;;; this is an embarrassingly-parallel batch, and lparallel:pmap is enough.

(defun call-with-kernel (fn &key (workers 4))
  "Call FN with no arguments and return its values, with LPARALLEL:*KERNEL*
bound: to the current kernel if there is one, otherwise to a fresh kernel of
WORKERS threads that is ended when FN returns. Starting threads costs about
a millisecond, so a loop over many ticks should hold one kernel for all of
them (see RUN-MECHANISTIC!) rather than start one per tick."
  (if lparallel:*kernel*
      (funcall fn)
      (let ((lparallel:*kernel* (lparallel:make-kernel workers)))
        (unwind-protect (funcall fn)
          (lparallel:end-kernel :wait t)))))

(defun pmap-cell-range! (n update-fn &key (workers 4) (chunk-size 256))
  "Call (UPDATE-FN i) for every I in 0..N-1 (a dense cell index - see
MAKE-CELL-INDEX in tessellation.lisp), split into contiguous chunks run
across WORKERS threads via LPARALLEL. This is the single-node,
shared-memory half of the mechanistic engine's parallelism story (see
mechanistic-step.lisp REACTION-SUBSTEP!/DIFFUSION-SUBSTEP!) - RUN-REPLICATES
above stays the mechanism for scaling out to a cluster, as a job array of
independent whole-simulation runs.

UPDATE-FN must only ever write to index I itself - reading anything else
(e.g. a neighbor cell's state from the previous tick) is fine, since that's
exactly what DIFFUSION-SUBSTEP! does, but two chunks writing the same index
would race. If LPARALLEL:*KERNEL* is already bound (e.g. this call happens
inside a REPLICATE-FN passed to RUN-REPLICATES, which already has one),
that kernel is reused; otherwise a temporary one sized to WORKERS is
started and torn down, so a single mechanistic run outside a replicate
batch still parallelizes its own per-tick cell updates."
  (call-with-kernel
   (lambda ()
     (lparallel:pmap
      nil
      (lambda (start)
        (loop for i from start below (min n (+ start chunk-size))
              do (funcall update-fn i)))
      (loop for start from 0 below n by chunk-size collect start)))
   :workers workers))

(defun run-replicates (n replicate-fn &key (workers 4) base-seed)
  "Run N replicates of REPLICATE-FN concurrently. REPLICATE-FN takes one
integer seed and should build whatever it needs from it (typically a
fresh VORONOI-GRID and WORLD), run to completion, and return a result.

Returns (VALUES RESULTS BASE-SEED SEEDS): the list of per-replicate
results in order, the top-level seed used to derive all of them, and the
individual seeds - so a whole concurrent batch is reproducible from one
recorded integer."
  (let* ((base-seed (or base-seed (random most-positive-fixnum (make-random-state t))))
         (seed-rs (sb-ext:seed-random-state base-seed))
         (seeds (loop repeat n collect (random most-positive-fixnum seed-rs)))
         (kernel (lparallel:make-kernel workers)))
    (unwind-protect
         (let ((lparallel:*kernel* kernel))
           (values (lparallel:pmap 'list replicate-fn seeds) base-seed seeds))
      (lparallel:end-kernel :wait t))))
