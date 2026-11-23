From: Apostolos Almpanis
To: neuss at femlisp.org
Subject: [Femlisp 2.2.1] :OUTPUT NIL crashes SOLVE with an unrelated type error

Hello,

I have been evaluating Femlisp as the discretisation layer for a porous-media
reactive-transport library in Common Lisp, and hit a small bug while learning
the API. Patch and reasoning below; everything was verified on SBCL 2.2.9 with
femlisp-20250622-git (VERSION: Femlisp-2.2.1) from Quicklisp dist 2026-01-01.

--------------------------------------------------------------------------
SYMPTOM
--------------------------------------------------------------------------

Passing :OUTPUT NIL to a blackboard -- the obvious way to ask for no output --
signals a type error from deep inside the solver, naming neither OUTPUT nor the
blackboard:

    (solve (blackboard :problem (cdr-model-problem 2)
                       :success-if '(> :nr-levels 2)
                       :output nil))

    ;=> The value NIL is not of type REAL

:OUTPUT 0 works and is the correct way to silence output, so there is an easy
workaround once you know it. The difficulty is finding it: the error points
nowhere near the cause, and NIL is what most people will try first.

--------------------------------------------------------------------------
CAUSE
--------------------------------------------------------------------------

The backtrace bottoms out in OUTPUT-P:

    0: ((:METHOD FL.ITERATION::OUTPUT-P (<ITERATION> T)) ...)

In src/strategy/gps.lisp, SOLVE on a blackboard binds *OUTPUT-DEPTH* from the
blackboard's :OUTPUT entry, using :KEEP as the not-supplied sentinel:

    (let* ((bb-output (getbb blackboard :output :keep))
           (*output-depth* (if (eq bb-output :keep)
                               *output-depth*
                               bb-output)))
      (solve strategy blackboard))

So :OUTPUT NIL binds *OUTPUT-DEPTH* to NIL. A nested iteration whose own OUTPUT
slot is unbound then reaches the fallback branch of OUTPUT-P, in
src/iteration/iterate.lisp:

    (or (member *output-depth* '(t :all :infinity))
        (<= *iteration-depth* *output-depth*))

With *OUTPUT-DEPTH* = NIL the MEMBER fails and NIL is passed to <= as its second
argument, which is the type error.

Two details suggest NIL is meant to be a legal value rather than user error:

  * The slot is documented as "A boolean indicating if output is to be done."
    NIL is the other boolean.

  * The :KEEP sentinel exists precisely so that "not supplied" can be told apart
    from a supplied value -- which is only necessary if NIL is expected to be
    supplied. Had NIL been the not-supplied marker, :KEEP would be redundant.

The same crash is reachable from several application entry points that bind
*OUTPUT-DEPTH* directly from a user-supplied argument, so it is not specific to
the blackboard path:

    src/applications/cdr/hom-cdr.lisp:160      (lret* ((*output-depth* output)
    src/applications/cdr/bl-cdr.lisp:80        (lret* ((*output-depth* output)
    src/applications/elasticity/elahom.lisp:132     (*output-depth* output)
    src/applications/navier-stokes/hom-ns.lisp:96   (*output-depth* output)
    src/applications/articles/heisig-neuss-2017.lisp:143 (*output-depth* output)

Each of these will signal the same error if called with OUTPUT = NIL. That is
the main argument for fixing OUTPUT-P rather than gps.lisp: one change covers
every path.

Incidentally, iterate.lisp already carries a disabled probe on this slot,

    #+(or) (when (null slot) (error "output-slot is NIL"))

which suggests a NIL in this area has caused trouble before.

--------------------------------------------------------------------------
PROPOSED FIX
--------------------------------------------------------------------------

Guard the comparison so that anything which is neither a real nor one of the
unlimited sentinels means "no output", which is the intent for NIL:

--- a/src/iteration/iterate.lisp
+++ b/src/iteration/iterate.lisp
@@ -60,7 +60,7 @@
 (defgeneric output-p (iteration blackboard)
   (:documentation "Determines if status output is done.")
   (:method ((iter <iteration>) blackboard)
-    "Checks if @var{*iteration-depth*} is smaller or equal
-@var{*output-depth*}."
+    "Checks if @var{*iteration-depth*} is smaller or equal
+@var{*output-depth*}.  A @var{*output-depth*} of NIL means no output."
     (declare (ignore blackboard))
     (if (slot-boundp iter 'output)
 	(let ((slot (slot-value iter 'output)))
 	  #+(or) (when (null slot) (error "output-slot is NIL"))
 	  slot)
 	(or (member *output-depth* '(t :all :infinity))
-	    (<= *iteration-depth* *output-depth*)))))
+	    (and (realp *output-depth*)
+		 (<= *iteration-depth* *output-depth*))))))

Three lines, no new dependency, and it keeps the existing shape of the code.

--------------------------------------------------------------------------
VERIFICATION
--------------------------------------------------------------------------

I checked that the change is behaviour-preserving rather than assuming it, by
evaluating both versions of the branch across every documented value of
*OUTPUT-DEPTH* crossed with *ITERATION-DEPTH* in {0,1,2}:

  *output-depth*   original      patched     verdict
  ----------------------------------------------------
  0, 1, 2          T / NIL       identical   same  (9 cases)
  T                T             T           same  (3 cases)
  :ALL             T             T           same  (3 cases)
  :INFINITY        T             T           same  (3 cases)
  NIL              TYPE-ERROR    NIL         FIXED (3 cases)

18 cases unchanged, 3 fixed, 0 unintended differences. End to end, with the
method redefined in a live image:

  before patch   :output nil -> TYPE-ERROR: The value NIL is not of type REAL
  after patch    :output nil -> OK (<ANSATZ-SPACE-VECTOR>)

and the values that already worked still do: :output 0, :output 1, :output T,
and :output absent all solve as before.

The script that produces all of the above is attached (output-nil-bug.lisp); it
is self-contained and runs from a bare SBCL with Quicklisp.

I also applied the attached patch file to a clean tree with `patch -p1`, cleared
the FASL cache, and let Femlisp recompile from the patched source. All five
cases then solve: :output NIL, 0, 1, T, and absent.

--------------------------------------------------------------------------
STILL PRESENT UPSTREAM, AS FAR AS I CAN SEE
--------------------------------------------------------------------------

The GitHub mirror at gitGNU/gnu_femlisp still carries the unguarded form on
master (src/iteration/iterate.lisp line 71, with no occurrence of REALP or
NUMBERP in that file), so this does not appear to be already fixed. That is a
mirror rather than your working tree, though, so if you have since changed this
please disregard.

--------------------------------------------------------------------------
ALTERNATIVE, IF YOU PREFER
--------------------------------------------------------------------------

If NIL should instead be rejected as invalid input, the equivalent improvement
would be to signal it where it enters, in gps.lisp, with a message that names
the argument -- e.g. "output must be a non-negative integer or one of T, :ALL,
:INFINITY". That would still be a large improvement over the current error, but
it contradicts the slot's "boolean" documentation and would leave the five
application entry points above untouched, so I would favour the patch.

Thanks for Femlisp, and for keeping it maintained -- it is the only thing of
its kind in Common Lisp, and having it work out of the box from Quicklisp made
the evaluation far easier than I expected.

Best regards,
Apostolos Almpanis
