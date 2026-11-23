(in-package :rhizo-diffusion)

(defun write-csv (path header rows)
  "Write HEADER (a list of column names) then ROWS (a list of lists) to
PATH as CSV. Always :SUPERSEDE, never :APPEND - a fresh header is written
every time. (The old code got bitten by exactly the opposite choice: a
header-creation call was commented out in one code path while a later
function still opened that file in append-only mode, so a run from a
clean directory errored with no file to append to.)"
  (with-open-file (stream path :direction :output :if-exists :supersede :if-does-not-exist :create)
    (csv-parser:write-csv-line stream header)
    (dolist (row rows)
      (csv-parser:write-csv-line stream row))))
