(in-package :rhizo-diffusion)

;;; Soil water and solute diffusion (PLAN.md, Phase 3).
;;;
;;; Each cell carries a volumetric water content theta. A dissolved solute
;;; is stored in the water (storage capacity theta, mechanistic-world.lisp)
;;; and diffuses through it with an effective diffusivity
;;;
;;;     D_s = D_0 xi(theta),
;;;
;;; D_0 its diffusivity in free water and xi the relative diffusivity (what
;;; Chou et al. 2012 call the tortuosity factor; it includes theta itself),
;;; so the diffusive flux is -D_s grad C with C the pore-water
;;; concentration. Several relations for xi are in use and they disagree,
;;; so the relation is a plug-in: a function of theta, built by one of the
;;; constructors below (RELATIVE-DIFFUSIVITY-MODEL picks one by name).
;;;
;;; Chou, Wu, Zeng & Chang (2012, Water Resour Res 48:W10539,
;;; doi:10.1029/2011WR011653) measured bromide diffusion in a sand, a sandy
;;; clay loam and a clay. Millington-Quirk with its standard constants
;;; overestimated all three and had the highest error; a quadratic in theta
;;; (Mullins & Sommer 1986) had the largest Akaike weight for the sand and
;;; the sandy clay loam. (PLAN.md cites this study as "Chou, Hu & Dane".)

(defun millington-quirk (porosity)
  "Return the Millington-Quirk relative diffusivity xi(theta) =
theta^(10/3) / POROSITY^2 (Millington & Quirk 1961, Trans Faraday Soc
57:1200; the solute form evaluated by Chou et al. 2012, their eq. 5 with
gamma = 1, n = 10/3, m = 2). The baseline relation of PLAN.md, Phase 3."
  (let ((phi2 (expt (float porosity 1d0) 2)))
    (lambda (theta) (/ (expt (float theta 1d0) 10/3) phi2))))

(defun papendick-campbell (coefficient exponent)
  "Return the power-law relative diffusivity xi(theta) = COEFFICIENT
theta^EXPONENT (Papendick & Campbell 1980: 2.8 theta^3; Chou et al. 2012,
eq. 4)."
  (lambda (theta) (* coefficient (expt (float theta 1d0) exponent))))

(defun mullins-sommer-threshold (delta mu)
  "Return the water content -MU / DELTA at which the Mullins-Sommer
relative diffusivity falls to zero."
  (- (/ mu delta)))

(defun mullins-sommer (delta mu)
  "Return the quadratic relative diffusivity xi(theta) = DELTA theta^2 + MU
theta (Mullins & Sommer 1986; Chou et al. 2012, eq. 7), taken as zero at
and below its threshold (MULLINS-SOMMER-THRESHOLD), where the water films
no longer connect and diffusion stops."
  (let ((threshold (mullins-sommer-threshold delta mu)))
    (lambda (theta)
      (if (<= theta threshold)
          0d0
          (+ (* delta theta theta) (* mu theta))))))

(defun no-tortuosity ()
  "Return the relative diffusivity xi(theta) = theta: D_s = theta D_0, so
the pore-water diffusivity is that in free water - a liquid medium such as
agar, and the free-water assumption of PLAN.md Phase 2. Not a soil
relation."
  (lambda (theta) (float theta 1d0)))

(defun soil-parameters ()
  "Return a fresh list of the soil MODEL-PARAMETERs: :BULK-DENSITY,
:POROSITY and :WATER-CONTENT, used by sorption and the storage capacity,
and the constants of the relative-diffusivity relations."
  (list
   (make-model-parameter
    :bulk-density 1.3d0 "kg/L"
    :note "Typical mineral-soil value, not from a measurement; set to the bulk density of the soil modelled.")
   (make-model-parameter
    :porosity (- 1d0 (/ 1.3d0 2.65d0)) "1"
    :note "1 - rho_b / rho_s from the flagged bulk density and an assumed particle density of 2.65 g/cm^3; set to the soil modelled.")
   (make-model-parameter
    :water-content 0.3d0 "1"
    :source "Scenario choice: volumetric water content (Phase 3 makes it a field)")
   (make-model-parameter
    :mullins-sommer-delta 2.1679d0 "1"
    :source "Chou et al. 2012, Water Resour Res 48:W10539, doi:10.1029/2011WR011653, Table 4: Mullins-Sommer quadratic fitted to their sandy clay loam"
    :note "With MU, threshold water content 0.2505. Sand: delta 0.2038, mu -0.0244 (threshold 0.120).")
   (make-model-parameter
    :mullins-sommer-mu -0.5431d0 "1"
    :source "Chou et al. 2012, Table 4: Mullins-Sommer quadratic fitted to their sandy clay loam")
   (make-model-parameter
    :papendick-campbell-coefficient 2.8d0 "1"
    :source "Papendick & Campbell 1980, as given by Chou et al. 2012 (eq. 4, c = 2.8)")
   (make-model-parameter
    :papendick-campbell-exponent 3d0 "1"
    :source "Papendick & Campbell 1980, as given by Chou et al. 2012 (eq. 4, exponent 3.0)")))

(defun relative-diffusivity-model (name parameters)
  "Return the relative diffusivity xi(theta) called NAME - :MILLINGTON-QUIRK,
:PAPENDICK-CAMPBELL, :MULLINS-SOMMER or :NO-TORTUOSITY - with its constants from PARAMETERS
\\(SOIL-PARAMETERS). Signals an error for another NAME."
  (flet ((p (key) (parameter-value parameters key)))
    (ecase name
      (:no-tortuosity (no-tortuosity))
      (:millington-quirk (millington-quirk (p :porosity)))
      (:papendick-campbell (papendick-campbell (p :papendick-campbell-coefficient)
                                               (p :papendick-campbell-exponent)))
      (:mullins-sommer (mullins-sommer (p :mullins-sommer-delta) (p :mullins-sommer-mu))))))

;;; Water content in a world. A world that tracks water has an immobile
;;; species :WATER holding each cell's theta, which rate laws read through
;;; (:PER name :WATER) references (reaction-network.lisp). Every pool is an
;;; amount per soil volume; every species with a positive free-water
;;; diffusivity counts as dissolved: its storage capacity is theta and its
;;; diffusivity D_0 xi(theta).

(defun water-species ()
  "Return the immobile species :WATER, each cell's volumetric water
content."
  (make-species :name :water :diffusivity 0d0))

(defun set-water-content! (world theta-field relation &key (blocked-p (constantly nil)))
  "Set the water content of every cell of WORLD to (FUNCALL THETA-FIELD
cell-id), which must be positive: the :WATER species, and for each
dissolved species (free-water diffusivity D_0 > 0) its capacity theta and
diffusivity D_0 (FUNCALL RELATION theta), or zero in cells satisfying
BLOCKED-P (a no-flux wall, e.g. a root). Amounts are not touched: when
water leaves, solute stays behind and its pore-water concentration rises.
Changes only WORLD. Returns WORLD."
  (let* ((index (mechanistic-world-cell-index world))
         (water (or (gethash :water (mechanistic-world-species-index world))
                    (error "SET-WATER-CONTENT!: the world has no :WATER species")))
         (state (mechanistic-world-state world)))
    (dolist (sp (mechanistic-world-species-list world))
      (let ((d0 (species-diffusivity sp)))
        (when (plusp d0)
          (set-species-capacity! world (species-name sp) theta-field)
          (set-species-diffusivity! world (species-name sp)
                                    (lambda (id)
                                      (if (funcall blocked-p id)
                                          0d0
                                          (* d0 (funcall relation (funcall theta-field id)))))))))
    (dotimes (i (cell-index-count index) world)
      (let ((theta (float (funcall theta-field (index->cell-id index i)) 1d0)))
        (unless (plusp theta)
          (error "SET-WATER-CONTENT!: water content ~a in cell ~s is not positive"
                 theta (index->cell-id index i)))
        (setf (aref state i water) theta)))))
