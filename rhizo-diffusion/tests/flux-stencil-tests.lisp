(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

(defun stencil-faces (stencil)
  "List of (i j conductance) for every face stored in STENCIL."
  (let ((offsets (flux-stencil-offsets stencil)))
    (loop for i below (1- (length offsets))
          nconc (loop for k from (aref offsets i) below (aref offsets (1+ i))
                      collect (list i (aref (flux-stencil-neighbors stencil) k)
                                    (aref (flux-stencil-conductances stencil) k))))))

(test flux-stencil-of-two-cells
  ;; area 2, face length 3, distance 4: conductance 0.75 each way
  (let* ((tess (two-cell-tessellation 2.0d0 3.0d0 4.0d0 :hop-weight 1.0d0))
         (stencil (make-flux-stencil tess (make-cell-index tess))))
    (is (equalp #(0 1 2) (flux-stencil-offsets stencil)))
    (is (equalp #(1 0) (flux-stencil-neighbors stencil)))
    (is (equalp #(0.75d0 0.75d0) (flux-stencil-conductances stencil)))
    (is (equalp #(2.0d0 2.0d0) (flux-stencil-areas stencil)))
    (is (= 0.75d0 (stencil-total-conductance stencil 0)))))

(test flux-stencil-is-symmetric-and-matches-the-protocol
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 60 60 25 :seed 3)))
         (index (make-cell-index tess))
         (stencil (make-flux-stencil tess index))
         (faces (stencil-faces stencil)))
    (dolist (face faces)
      (destructuring-bind (i j g) face
        (is (find (list j i g) faces :test #'equal) "face ~a-~a has no mirror" i j)
        (is (= g (face-conductance tess (index->cell-id index i) (index->cell-id index j))))))
    (dotimes (i (cell-index-count index))
      (let ((id (index->cell-id index i)))
        (is (= (cell-area tess id) (aref (flux-stencil-areas stencil) i)))
        (is (< (abs (- (stencil-total-conductance stencil i)
                       (loop for nb in (neighbors-of tess id) sum (face-conductance tess id nb))))
               1d-12))))))

(test flux-stencil-leaves-out-corner-only-contacts
  ;; CA-grown grids list corner-only contacts as neighbours, with zero
  ;; edge length; they carry no flux, so the stencil drops them.
  ;; This grid has 7 corner-only contacts (14 directed), found by search.
  (let* ((tess (voronoi-grid-tessellation (make-voronoi-grid 40 40 40 :seed 4)))
         (stencil (make-flux-stencil tess (make-cell-index tess)))
         (pairs (loop for a in (cell-ids tess)
                      append (loop for b in (neighbors-of tess a) collect (cons a b))))
         (with-face (count-if (lambda (p) (plusp (edge-length tess (car p) (cdr p)))) pairs)))
    (is (< with-face (length pairs)) "this grid should have some corner-only contacts")
    (is (= with-face (length (flux-stencil-neighbors stencil))))
    (is (every #'plusp (flux-stencil-conductances stencil)))))

(test stencil-net-flux-is-conductance-times-difference
  (let* ((tess (two-cell-tessellation 1.0d0 3.0d0 4.0d0))
         (stencil (make-flux-stencil tess (make-cell-index tess)))
         (index (make-cell-index tess))
         (state (make-array '(2 2) :element-type 'double-float :initial-element 0d0)))
    (setf (aref state (cell-id->index index "0") 1) 10d0)
    (is (= -7.5d0 (stencil-net-flux stencil state (cell-id->index index "0") 1)))
    (is (= 7.5d0 (stencil-net-flux stencil state (cell-id->index index "1") 1)))
    (is (= 0d0 (stencil-net-flux stencil state 0 0)) "other columns are untouched")))

(test diffusion-substep-matches-the-protocol-formula
  ;; Oracle: one explicit step computed straight from the tessellation
  ;; protocol, as DIFFUSION-SUBSTEP! did before the stencil existed.
  (let* ((tess (voronoi-grid-tessellation (make-nearest-seed-grid 50 50 30 :seed 8)))
         (d 0.7d0) (dt 0.4d0)
         (species (list (make-species :name :c :diffusivity d)
                        (make-species :name :fixed :diffusivity 0d0)))
         (world (make-mechanistic-world
                 tess :species species :reactions '() :dt dt
                      :initial-state (field-initial-state
                                      tess (list (cons :c (lambda (p) (sin (* 0.1d0 (car p)))))
                                                 (cons :fixed (lambda (p) (cdr p)))))))
         (index (mechanistic-world-cell-index world))
         (before (alexandria:copy-array (mechanistic-world-state world))))
    (flet ((value (id s) (aref before (cell-id->index index id) s)))
      (diffusion-substep! world dt)
      (dolist (id (cell-ids tess))
        (let ((expected (+ (value id 0)
                           (* dt (/ d (cell-area tess id))
                              (loop for nb in (neighbors-of tess id)
                                    sum (* (face-conductance tess id nb)
                                           (- (value nb 0) (value id 0))))))))
          (is (< (abs (- expected (cell-species-value world id :c))) 1d-13))
          (is (= (value id 1) (cell-species-value world id :fixed))))))))

(test run-mechanistic-equals-repeated-steps
  (flet ((fresh ()
           (make-mechanistic-world
            (two-cell-tessellation 1.0d0 1.0d0 1.0d0)
            :species (list (make-species :name :c :diffusivity 0.5d0)) :reactions '() :dt 0.1d0
            :initial-state (lambda (id) (if (string= id "0") '((:c . 10d0)) '((:c . 0d0)))))))
    (let ((a (fresh)) (b (fresh)))
      (run-mechanistic! a 25)
      (dotimes (k 25) (step-mechanistic! b))
      (is (= 25 (mechanistic-world-tick a)))
      (is (equalp (mechanistic-world-state a) (mechanistic-world-state b))))))

(test call-with-kernel-reuses-or-makes-a-kernel
  (is (null lparallel:*kernel*) "precondition: no kernel bound in the test runner")
  (is (equal '(1 2) (multiple-value-list (call-with-kernel (lambda () (values 1 2))))))
  (is (call-with-kernel (lambda () lparallel:*kernel*)) "a temporary kernel is bound inside")
  (is (null lparallel:*kernel*) "and gone afterwards")
  (let ((outer (lparallel:make-kernel 2)))
    (unwind-protect
         (let ((lparallel:*kernel* outer))
           (is (eq outer (call-with-kernel (lambda () lparallel:*kernel*)))))
      (let ((lparallel:*kernel* outer)) (lparallel:end-kernel :wait t)))))
