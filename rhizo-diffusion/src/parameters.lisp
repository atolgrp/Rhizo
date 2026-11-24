(in-package :rhizo-diffusion)

;;; Model parameters with their sources.
;;;
;;; CLAUDE.md, rule 3: physical and biological parameters carry their
;;; literature reference, and a parameter without one is flagged as such,
;;; not quietly defaulted. A MODEL-PARAMETER holds a value, its unit, and a
;;; SOURCE string; a parameter whose SOURCE is NIL is "flagged" - it still
;;; has a value so a run can go ahead, but FLAGGED-PARAMETERS lists it and
;;; every report prints it (PARAMETER-TABLE-ROWS).
;;;
;;; A source string says what kind of value it is: a literature value
;;; (author, year, where in the paper), a value derived from literature
;;; values (with the derivation in NOTE), or a scenario choice (a property of
;;; the simulated experiment, such as a cell density, not of the chemistry).
;;;
;;; Parameter sets are plain lists, built fresh by functions such as
;;; AHL-LITERATURE-PARAMETERS (signalling.lisp), never global state.

(defstruct (model-parameter (:constructor make-model-parameter (name value unit &key source note)))
  name    ; keyword
  value   ; a number, in UNIT
  unit    ; a string, e.g. "um^2/s"
  source  ; a string, or NIL for a flagged (unsourced) parameter
  note)   ; a string or NIL: derivation, caveats, analogue compound, ...

(defun find-parameter (parameters name)
  "Return the MODEL-PARAMETER named NAME in the list PARAMETERS. Signals
an error if there is none."
  (or (find name parameters :key #'model-parameter-name)
      (error "FIND-PARAMETER: no parameter ~s" name)))

(defun parameter-value (parameters name)
  "Return the value of the parameter named NAME in PARAMETERS. Signals an
error if there is none."
  (model-parameter-value (find-parameter parameters name)))

(defun flagged-parameters (parameters)
  "Return the names, in order, of the parameters in PARAMETERS that have no
source."
  (mapcar #'model-parameter-name (remove-if #'model-parameter-source parameters)))

(defun with-parameter-values (parameters overrides)
  "Return a fresh list of PARAMETERS with the values in OVERRIDES, an alist
\(name . value), replaced. Each replaced parameter keeps its unit and source,
and its note records the value the source gave. PARAMETERS is not
modified. Signals an error for an unknown name."
  (dolist (entry overrides)
    (find-parameter parameters (car entry)))
  (mapcar (lambda (p)
            (let ((entry (assoc (model-parameter-name p) overrides)))
              (if entry
                  (make-model-parameter
                   (model-parameter-name p) (cdr entry) (model-parameter-unit p)
                   :source (model-parameter-source p)
                   :note (format nil "overridden; the source gives ~a~@[ (~a)~]"
                                 (model-parameter-value p) (model-parameter-note p)))
                  p)))
          parameters))

(defun parameter-table-rows (parameters)
  "Return one list (name value unit source note) per parameter, for
reports. A flagged parameter's source reads \"NO SOURCE (flagged)\"."
  (mapcar (lambda (p)
            (list (model-parameter-name p) (model-parameter-value p) (model-parameter-unit p)
                  (or (model-parameter-source p) "NO SOURCE (flagged)")
                  (model-parameter-note p)))
          parameters))
