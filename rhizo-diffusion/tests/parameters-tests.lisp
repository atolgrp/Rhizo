(in-package :rhizo-diffusion.tests)
(in-suite rhizo-diffusion-suite)

;;; Model parameters carry their literature source (CLAUDE.md rule 3).

(defun example-parameters ()
  (list (make-model-parameter :d 490d0 "um^2/s" :source "Hense et al. 2012")
        (make-model-parameter :kd 33.6d0 "L/kg" :source "Sheng et al. 2018" :note "C10-HSL analogue")
        (make-model-parameter :rho-b 1.3d0 "kg/L" :note "typical mineral soil")))

(test parameter-value-looks-up-by-name
  (let ((parameters (example-parameters)))
    (is (= 490d0 (parameter-value parameters :d)))
    (is (= 1.3d0 (parameter-value parameters :rho-b)))
    (is (string= "Sheng et al. 2018" (model-parameter-source (find-parameter parameters :kd))))
    (signals error (parameter-value parameters :no-such-parameter))))

(test flagged-parameters-are-those-without-a-source
  (is (equal '(:rho-b) (flagged-parameters (example-parameters))))
  (is (null (flagged-parameters (butlast (example-parameters))))))

(test with-parameter-values-returns-a-fresh-list
  (let* ((parameters (example-parameters))
         (changed (with-parameter-values parameters '((:kd . 84.3d0)))))
    (is (= 84.3d0 (parameter-value changed :kd)))
    (is (= 33.6d0 (parameter-value parameters :kd)) "the original list is not modified")
    (is (= 490d0 (parameter-value changed :d)))
    (is (search "33.6" (model-parameter-note (find-parameter changed :kd)))
        "an overridden parameter records the value its source gave")
    (is (equal (flagged-parameters parameters) (flagged-parameters changed)))
    (signals error (with-parameter-values parameters '((:no-such-parameter . 1d0))))))

(test parameter-table-rows-mark-missing-sources
  (let ((rows (parameter-table-rows (example-parameters))))
    (is (= 3 (length rows)))
    (is (equal '(:d 490d0 "um^2/s" "Hense et al. 2012" nil) (first rows)))
    (is (search "NO SOURCE" (fourth (third rows))))))
