(in-package :rhizo-diffusion)

;;; Tagged signal: each patch's AHL as its own species (PLAN.md, Phase 4).
;;;
;;; An edge of the signalling graph asks how much of patch i's signal
;;; reaches patch j - which a single AHL field cannot say. But transport is
;;; linear in the sources: diffusion, hydrolysis and sorption act on each
;;; molecule alike. So the AHL released by patch k is tracked as its own
;;; species TAG-k (with its sorbed pool), produced only where the immobile
;;; marker NODE-k holds patch k's producers. All tags share the water
;;; content, diffusivity, decay and sorption, so their sum is the total AHL
;;; field and each tag is exactly patch k's contribution.
;;;
;;; Two couplings act on the total, and stay exact: quorum feedback
;;; (production alpha + beta H(total pore-water AHL), when FEEDBACK) and
;;; Michaelis-Menten biodegradation, whose saturation is shared (each tag is
;;; degraded at kcat rho_d c_k / (Km + c_total)). Shared ledger pools
;;; TAGS-PRODUCED and TAGS-LOST keep mass balance pass/fail.

(defun tag-name (k)
  "Return the species name of patch K's tagged AHL, e.g. :TAG-3."
  (intern (format nil "TAG-~d" k) :keyword))

(defun node-name (k)
  "Return the name of patch K's immobile producer marker, e.g. :NODE-3."
  (intern (format nil "NODE-~d" k) :keyword))

(defun tagged-species (n parameters)
  "Return the species of a tagged world with N patches: TAG-k (free-water
diffusivity PARAMETERS' :DIFFUSIVITY) and its sorbed pool, NODE-k, the
ledger pools TAGS-PRODUCED and TAGS-LOST, uniform DEGRADER biomass and the
WATER-SPECIES."
  (let ((d0 (float (parameter-value parameters :diffusivity) 1d0)))
    (append (loop for k below n collect (make-species :name (tag-name k) :diffusivity d0))
            (loop for k below n collect (make-species :name (signal-pool (tag-name k) :sorbed)))
            (loop for k below n collect (make-species :name (node-name k)))
            (mapcar (lambda (name) (make-species :name name)) '(:tags-produced :tags-lost :degrader))
            (list (water-species)))))

(defun total-tag-concentration (state species-index tags)
  "Return the total pore-water concentration of the species TAGS in STATE."
  (/ (loop for tag in tags sum (state-value state species-index tag))
     (state-value state species-index :water)))

(defun tagged-reactions (n parameters &key feedback)
  "Return the reactions of a tagged world with N patches and PARAMETERS:
for each patch k, production of TAG-k at NODE-k's producer density times
alpha + beta H(total pore-water AHL) when FEEDBACK, else alpha alone;
hydrolysis of the dissolved tag; biodegradation at kcat rho_d c_k / (Km +
c_total); and linear sorption toward rho_b K_d c_k."
  (let ((tags (loop for k below n collect (tag-name k))))
    (flet ((p (name) (parameter-value parameters name)))
      (loop for k below n
            for tag = (tag-name k)
            for node = (node-name k)
            nconc
            (list
             (make-reaction
              :name (signal-pool tag :production)
              :stoichiometry `((,tag . 1d0) (:tags-produced . 1d0))
              :rate-fn (let ((node node))
                         (lambda (state index)
                           (* (state-value state index node)
                              (+ (p :basal-production)
                                 (if feedback
                                     (* (p :induced-production)
                                        (hill-activation (total-tag-concentration state index tags)
                                                         (p :induction-threshold)
                                                         (p :hill-coefficient)))
                                     0d0))))))
             (make-reaction
              :name (signal-pool tag :decay)
              :stoichiometry `((,tag . -1d0) (:tags-lost . 1d0))
              :rate-fn (first-order-decay-rate (p :hydrolysis-rate) tag))
             (make-reaction
              :name (signal-pool tag :biodegradation)
              :stoichiometry `((,tag . -1d0) (:tags-lost . 1d0))
              :rate-fn (let ((tag tag))
                         (lambda (state index)
                           (let ((c-total (total-tag-concentration state index tags)))
                             (/ (* (p :biodegradation-kcat) (state-value state index :degrader)
                                   (max 0d0 (concentration-of state index (pore-water tag))))
                                (+ (p :biodegradation-km) (max 0d0 c-total)))))))
             (make-reaction
              :name (signal-pool tag :sorption)
              :stoichiometry `((,tag . -1d0) (,(signal-pool tag :sorbed) . 1d0))
              :rate-fn (linear-reversible-rate (p :desorption-rate) (pore-water tag)
                                               (signal-pool tag :sorbed)
                                               (sorption-coefficient parameters))))))))

(defun make-tagged-ahl-world (tessellation node-cells &key (parameters (ahl-literature-parameters))
                                                        dt (producer-density 5d14)
                                                        (degrader-density 0d0) feedback
                                                        water-content
                                                        (relative-diffusivity :millington-quirk)
                                                        (diffusion-scheme :imex))
  "Return a fresh tagged MECHANISTIC-WORLD over TESSELLATION, patch k
holding PRODUCER-DENSITY producers per litre of soil in cell (NTH k
NODE-CELLS), DEGRADER-DENSITY degraders everywhere, with PARAMETERS,
FEEDBACK (see TAGGED-REACTIONS), and water content and relation as for
MAKE-AHL-WORLD."
  (let ((n (length node-cells)))
    (set-water-content!
     (make-mechanistic-world
      tessellation :species (tagged-species n parameters)
                   :reactions (tagged-reactions n parameters :feedback feedback)
                   :dt dt :diffusion-scheme diffusion-scheme
                   :initial-state (lambda (id)
                                    (acons :degrader degrader-density
                                           (loop for cell in node-cells
                                                 for k from 0
                                                 when (equal id cell)
                                                   collect (cons (node-name k) producer-density)))))
     (water-content-field water-content parameters)
     (relative-diffusivity-model relative-diffusivity parameters))))

(defun tagged-edges (world node-cells threshold)
  "Return the directed edges (i . j), i /= j, of the signalling graph in the
tagged WORLD: patch i's tag reaches patch j's cell (NTH j NODE-CELLS) at a
pore-water concentration of at least THRESHOLD."
  (let* ((state (mechanistic-world-state world))
         (index (mechanistic-world-cell-index world))
         (species-index (mechanistic-world-species-index world))
         (water (gethash :water species-index))
         (rows (mapcar (lambda (cell) (cell-id->index index cell)) node-cells)))
    (loop for i below (length node-cells)
          for tag = (gethash (tag-name i) species-index)
          nconc (loop for row in rows
                      for j from 0
                      when (and (/= i j) (>= (/ (aref state row tag) (aref state row water)) threshold))
                        collect (cons i j)))))

(defun tagged-pool-balance (world)
  "Return all tags and their sorbed pools plus TAGS-LOST minus TAGS-PRODUCED
in WORLD: invariant in a conserving run."
  (let ((n (count-if (lambda (sp) (search "NODE-" (symbol-name (species-name sp))))
                     (mechanistic-world-species-list world))))
    (pool-balance world
                  (append (loop for k below n collect (tag-name k))
                          (loop for k below n collect (signal-pool (tag-name k) :sorbed))
                          '(:tags-lost))
                  '(:tags-produced))))
