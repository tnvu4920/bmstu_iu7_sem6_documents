;;;; This file contains structures and functions for the maintenance of
;;;; basic information about defined types. Different object systems
;;;; can be supported simultaneously.

;;;; This software is part of the SBCL system. See the README file for
;;;; more information.
;;;;
;;;; This software is derived from the CMU CL system, which was
;;;; written at Carnegie Mellon University and released into the
;;;; public domain. The software is in the public domain and is
;;;; provided with absolutely no warranty. See the COPYING and CREDITS
;;;; files for more information.

(in-package "SB-KERNEL")

(!begin-collecting-cold-init-forms)
;;; Has the type system been properly initialized? (I.e. is it OK to
;;; use it?)
(define-load-time-global *type-system-initialized* nil)
(!cold-init-forms (setq *type-system-initialized* nil))

(in-package "SB-PCL")

(declaim (#+sb-xc-host special #-sb-xc-host global *the-class-t*
                  *the-class-slot-object*
                  *the-class-structure-object*
                  *the-class-standard-object*
                  *the-class-function*
                  *the-class-funcallable-standard-object*
                  *the-class-system-class*
                  *the-class-slot-class*
                  *the-class-condition-class*
                  *the-class-structure-class*
                  *the-class-standard-class*
                  *the-class-funcallable-standard-class*
                  *the-class-forward-referenced-class*
                  *the-class-method*
                  *the-class-standard-method*
                  *the-class-standard-reader-method*
                  *the-class-standard-writer-method*
                  *the-class-standard-boundp-method*
                  *the-class-global-reader-method*
                  *the-class-global-writer-method*
                  *the-class-global-boundp-method*
                  *the-class-standard-generic-function*
                  *the-class-standard-direct-slot-definition*
                  *the-class-standard-effective-slot-definition*

                  *the-eslotd-standard-class-slots*
                  *the-eslotd-funcallable-standard-class-slots*))

(in-package "SB-KERNEL")


;;;; the CLASSOID structure

;;; The CLASSOID structure is a supertype of all classoid types.
;;; Its definition occurs in 'early-classoid.lisp'
#+sb-xc-host
(defmethod make-load-form ((self classoid) &optional env)
  (declare (ignore env))
  `(find-classoid ',(classoid-name self)))


;;;; basic LAYOUT stuff

;;; a vector of conses, initialized by genesis
;;;
;;; In each cons, the car is the symbol naming the layout, and the
;;; cdr is the layout itself.
;;; If #+metaspace then the cdr is actually of type WRAPPER,
;;; and if #-metaspace then the wrapper is a LAYOUT.
(defvar *!initial-wrappers*)

;;; a table mapping class names to layouts for classes we have
;;; referenced but not yet loaded. This is initialized from an alist
;;; created by genesis describing the layouts that genesis created at
;;; cold-load time.
(define-load-time-global *forward-referenced-wrappers*
    ;; FIXME: why is the test EQUAL and not EQ? Aren't the keys all symbols?
    (make-hash-table :test 'equal))
#-sb-xc-host
(!cold-init-forms
  ;; *forward-referenced-wrappers* is protected by *WORLD-LOCK*
  ;; so it does not need a :synchronized option.
 (setq *forward-referenced-wrappers* (make-hash-table :test 'equal))
 (dovector (x *!initial-wrappers*)
   (let ((expected (hash-layout-name (car x)))
         (actual (wrapper-clos-hash (cdr x))))
     (unless (= actual expected) (bug "XC layout hash calculation failed")))
   (setf (gethash (car x) *forward-referenced-wrappers*) (cdr x))))

;;; FIXME: This lock is only seized in the classoid/layout/class
;;; system, and is now a misnomer.
#-sb-xc-host
(define-load-time-global **world-lock** nil)
#-sb-xc-host
(!cold-init-forms
 (setq **world-lock** (sb-thread:make-mutex :name "World Lock")))

(defmacro with-world-lock (() &body body)
  #+sb-xc-host `(progn ,@body)
  #-sb-xc-host `(sb-thread:with-recursive-lock (**world-lock**) ,@body))

;;; The LAYOUT structure itself is defined in 'early-classoid.lisp'

#+sb-xc-host
(progn
(defun make-layout (hash classoid &rest keys)
  (macrolet ((make (&rest extra)
               `(apply #'host-make-wrapper
                       (cdr (assq (classoid-name classoid) *popular-structure-types*))
                       hash classoid ,@extra :allow-other-keys t keys)))
    #-metaspace (make)
    #+metaspace (let* ((layout (%make-layout))
                       (wrapper (make :friend layout)))
                  (setf (layout-friend layout) wrapper)
                  wrapper)))
;; The target reconstructs wrappers using FOP-LAYOUT but the host uses MAKE-LOAD-FORM.
(defmethod cl:make-load-form ((wrapper wrapper) &optional env)
  (declare (ignore env))
  (labels ((externalize (wrapper &aux (classoid (wrapper-classoid wrapper))
                                      (name (classoid-name classoid)))
             (when (or (wrapper-invalid wrapper)
                       (not name)
                       (typep classoid 'undefined-classoid))
               (sb-c:compiler-error "can't dump ~S" wrapper))
             `(xc-load-wrapper ',name
                               ,(wrapper-depthoid wrapper)
                               (vector ,@(map 'list #'externalize (wrapper-inherits wrapper)))
                               ,(wrapper-length wrapper)
                               ,(wrapper-bitmap wrapper))))
    (externalize wrapper)))
(defun xc-load-wrapper (name depthoid inherits length bitmap)
  (let ((classoid (find-classoid name)))
    (aver (and classoid (not (undefined-classoid-p classoid))))
    (let ((wrapper (classoid-wrapper classoid)))
      (unless (and (= (wrapper-depthoid wrapper) depthoid)
                   (= (length (wrapper-inherits wrapper)) (length inherits))
                   (every #'eq (wrapper-inherits wrapper) inherits)
                   (= (wrapper-length wrapper) length)
                   (= (wrapper-bitmap wrapper) bitmap))
        (error "XC can't reload layout for ~S with ~S vs ~A"
               name (list depthoid inherits length bitmap) wrapper))
      wrapper)))
) ; end PROGN

(defmethod print-object ((wrapper wrapper) stream)
  (print-unreadable-object (wrapper stream :type t :identity t)
    (format stream
            "~@[(ID=~d) ~]for ~S~@[, INVALID=~S~]"
            (layout-id wrapper)
            (wrapper-proper-name wrapper)
            (wrapper-invalid wrapper))))

(eval-when (#-sb-xc :compile-toplevel :load-toplevel :execute)
  (defun wrapper-proper-name (wrapper)
    (classoid-proper-name (wrapper-classoid wrapper))))

;;; Return the layout currently installed in the classoid named NAME.
;;; If there is none, then make a layout referring for an undefined classoid.
;;; NB: for #+metaspace this returns a WRAPPER, not a LAYOUT.
(declaim (ftype (sfunction (symbol) wrapper) find-layout))
(defun find-layout (name)
  (binding* ((classoid (find-classoid name nil) :exit-if-null) ; threadsafe
             (wrapper (classoid-wrapper classoid) :exit-if-null))
    (return-from find-layout wrapper))
  (let ((table *forward-referenced-wrappers*))
    (with-world-lock ()
      (let ((classoid (find-classoid name nil)))
        (or (and classoid (classoid-wrapper classoid))
            (values (ensure-gethash name table
                                    (make-layout
                                     (hash-layout-name name)
                                     (or classoid
                                         (make-undefined-classoid name))))))))))

;;; If LAYOUT's slot values differ from the specified slot values in
;;; any interesting way, then give a warning and return T.
(declaim (ftype (function (simple-string
                           wrapper
                           simple-string
                           index
                           simple-vector
                           layout-depthoid
                           layout-bitmap))
                warn-if-altered-layout))
(defun warn-if-altered-layout (old-context old-layout context
                               length inherits depthoid bitmap)
  (let ((name (wrapper-proper-name old-layout))
        (old-inherits (wrapper-inherits old-layout)))
    (or (when (mismatch old-inherits inherits :key #'wrapper-proper-name)
          (warn "change in superclasses of class ~S:~%  ~
                       ~A superclasses: ~S~%  ~
                       ~A superclasses: ~S"
                      name
                      old-context
                      (map 'list #'wrapper-proper-name old-inherits)
                      context
                      (map 'list #'wrapper-proper-name inherits))
          t)
        (let ((diff (mismatch old-inherits inherits)))
          (when diff
            (warn "in class ~S:~%  ~
                    ~@(~A~) definition of superclass ~S is incompatible with~%  ~
                    ~A definition."
                   name
                   old-context
                   (wrapper-proper-name (svref old-inherits diff))
                   context)
            t))
        (let ((old-length (wrapper-length old-layout)))
          (unless (= old-length length)
            (warn "change in instance length of class ~S:~%  ~
                   ~A length: ~W~%  ~
                   ~A length: ~W"
                  name
                  old-context old-length
                  context length)
            t))
        (let ((old-bitmap (wrapper-bitmap old-layout)))
          (unless (= old-bitmap bitmap)
            (warn "change in placement of raw slots of class ~S ~
between the ~A definition and the ~A definition"
                  name old-context context)
            t))
        (unless (= (wrapper-depthoid old-layout) depthoid)
          (warn "change in the inheritance structure of class ~S~%  ~
                 between the ~A definition and the ~A definition"
                name old-context context)
          t))))

(defun load-layout (name depthoid inherits length bitmap flags)
  (let* ((layout
          (or (binding* ((classoid (find-classoid name nil) :exit-if-null))
                (classoid-wrapper classoid))
              (let ((table *forward-referenced-wrappers*))
                (with-world-lock ()
                 (let ((classoid (find-classoid name nil)))
                   (or (and classoid (classoid-wrapper classoid))
                       (ensure-gethash
                        name table
                        (make-layout
                         (hash-layout-name name)
                         (or classoid (make-undefined-classoid name))
                         :depthoid depthoid :inherits inherits
                         :length length :bitmap bitmap :flags flags))))))))
         (classoid
          (or (find-classoid name nil) (wrapper-classoid layout))))
    (if (or (eq (wrapper-invalid layout) :uninitialized)
            (not *type-system-initialized*))
        (setf (wrapper-classoid layout) classoid)
        ;; There was an old layout already initialized with old
        ;; information, and we'll now check that old information
        ;; which was known with certainty is consistent with current
        ;; information which is known with certainty.
        (when (warn-if-altered-layout "current" layout "compile time"
                                    length inherits depthoid bitmap)
          (error "The loaded code expects an incompatible layout for class ~S."
                 (wrapper-proper-name layout))))
    layout))

(defun classoid-lock (classoid)
  #+sb-xc-host (declare (ignore classoid))
  #-sb-xc-host
  (or (classoid-%lock classoid)
      (let* ((lock (sb-thread:make-mutex :name "classoid lock"))
             (oldval (cas (classoid-%lock classoid) nil lock)))
        (if (eq oldval nil) lock oldval))))

(defun add-subclassoid (super sub wrapper)
  (with-system-mutex ((classoid-lock super))
    (let ((table (classoid-subclasses super))
          (count 0))
      (cond ((hash-table-p table)
             (setf (gethash sub table) wrapper))
            ((dolist (cell table)
               (when (eq (car cell) sub)
                 (return (setf (cdr cell) wrapper)))
               (incf (truly-the fixnum count))))
            ((<= count 7)
             (setf (classoid-subclasses super) (acons sub wrapper table))
             ;; Is this barrier really necessary? mutex release is a release barrier.
             ;; I was probably struggling with crashes in the 'classoid-typep.impure'
             ;; test, and was just trying anything and everything.
             (sb-thread:barrier (:write))
             wrapper)
            (t
             ;; Upgrade to a hash-table
             (let ((new #+sb-xc-host (make-hash-table :test 'eq)
                        #-sb-xc-host (make-hash-table :hash-function #'type-hash-value
                                                      :test 'eq)))
               (loop for (key . val) in table do (setf (gethash key new) val))
               (setf (gethash sub new) wrapper)
               (setf (classoid-subclasses super) new)
               (sb-thread:barrier (:write))
               wrapper))))))

;;; Mnemonic device: the argument order is as GETHASH (1st = key, 2nd = table).
;;; But the 2nd arg is the superclassoid, *not* its subclassoid table,
;;; because the mutex is stored in the classoid, not the table.
(defun get-subclassoid (sub super)
  (sb-thread:barrier (:read))
  (let ((table (classoid-subclasses super)))
  (when table
    (cond ((listp table) (cdr (assq sub table)))
          (t
           ;; our hash-table are not safe to read with a single writer
           (with-system-mutex ((classoid-lock super))
             (values (gethash sub table))))))))

;;; Mnemonic device: it's like REMHASH (1st = key, 2nd = table)
(defun remove-subclassoid (sub super)
  (sb-thread:barrier (:read))
  (when (classoid-subclasses super)
    (with-system-mutex ((classoid-lock super))
      (let ((table (classoid-subclasses super)))
        (cond ((listp table)
               (setf (classoid-subclasses super)
                     (delete sub table :key #'car :test #'eq)))
              (t
               ;; There's no reason to demote a table to a list ever.
               (remhash sub table))))))
  nil)

(defmacro do-subclassoids (((classoid-var wrapper-var) super) &body body)
  (let ((f (make-symbol "FUNCTION")))
    `(dx-flet ((,f (,classoid-var ,wrapper-var) ,@body))
       (call-with-subclassoids #',f (the classoid ,super)))))

(defun call-with-subclassoids (function super &aux (table (classoid-subclasses super)))
  ;; Uses of DO-SUBCLASSOIDS don't need to acquire the classoid lock on SUPER.
  ;; Even if there are readers or writers, hash-table iteration is safe.
  ;; This was not always so - iteration could overrun the k/v array because it always
  ;; re-fetched the scan limit, which could see a higher limit than corresponded
  ;; to the k/v vector that it had gotten initially.
  ;; If you're doing concurrent modification of the class heterarchy, there are no
  ;; real guarantees. We dont' always hold a lock at a wider scope than the table lock,
  ;; but sometimes we do, such as in REGISTER-LAYOUT.
  (if (listp table)
      (loop for (key . value) in table do (funcall function key value))
      (maphash (lambda (key value) (funcall function key value))
               table))
  nil)

;;; Record LAYOUT as the layout for its class, adding it as a subtype
;;; of all superclasses. This is the operation that "installs" a
;;; layout for a class in the type system, clobbering any old layout.
;;; However, this does not modify the class namespace; that is a
;;; separate operation (think anonymous classes.)
;;; -- If INVALIDATE, then all the layouts for any old definition
;;;    and subclasses are invalidated, and the SUBCLASSES slot is cleared.
;;; -- If DESTRUCT-LAYOUT, then this is some old layout, and is to be
;;;    destructively modified to hold the same type information.
(macrolet ((set-bitmap-from-layout (to-layout from-layout)
             `(let ((to-index
                     (+ (type-dd-length sb-vm:layout)
                        (calculate-extra-id-words (layout-depthoid ,to-layout))))
                    (from-index
                     (+ (type-dd-length sb-vm:layout)
                        (calculate-extra-id-words (layout-depthoid ,from-layout)))))
                (dotimes (i (bitmap-nwords ,from-layout))
                  (%raw-instance-set/word ,to-layout (+ to-index i)
                        (%raw-instance-ref/word ,from-layout (+ from-index i)))))))
(defun register-layout (wrapper &key (invalidate t) modify)
  (declare (type wrapper wrapper) (type (or wrapper null) modify))
  (with-world-lock ()
    (let* ((classoid (wrapper-classoid wrapper))
           (classoid-wrapper (classoid-wrapper classoid)))

      ;; Attempting to register ourselves with a temporary undefined
      ;; class placeholder is almost certainly a programmer error. (I
      ;; should know, I did it.) -- WHN 19990927
      (aver (not (undefined-classoid-p classoid)))

      ;; This assertion dates from classic CMU CL. The rationale is
      ;; probably that calling REGISTER-LAYOUT more than once for the
      ;; same LAYOUT is almost certainly a programmer error.
      (aver (not (eq classoid-wrapper wrapper)))

      ;; Figure out what classes are affected by the change, and issue
      ;; appropriate warnings and invalidations.
      (when classoid-wrapper
        (%modify-classoid classoid)
        (do-subclassoids ((subclass subclass-wrapper) classoid) ; under WORLD-LOCK
            (%modify-classoid subclass)
            (when invalidate
              (%invalidate-layout subclass-wrapper)))
        (when invalidate
          (%invalidate-layout classoid-wrapper)
          (setf (classoid-subclasses classoid) nil)))

      (if modify
          #+sb-xc-host (error "Why mutate a layout in XC host?")
          #-sb-xc-host
          ;; Destructively modifying a layout is not threadsafe at all.
          ;; Use at your own risk (interactive use only).
          (let ((inherits (wrapper-inherits wrapper))
                (depthoid (wrapper-depthoid wrapper)) ; "new" depthoid
                (extra-id-words ; "old" extra words
                 (calculate-extra-id-words (wrapper-depthoid modify)))
                (layout (wrapper-friend wrapper))
                (id ; read my ID before screwing with the depthoid
                 (layout-id modify)))
            (aver (logtest +structure-layout-flag+ (layout-flags layout)))
            (aver (= (length inherits) depthoid))
            ;; DEPTHOID implies the number of words of "extra" IDs preceding the bitmap.
            ;; Layout alteration is forbidden if it would affect the number of such words.
            ;; So MUTABLE-LAYOUT-P should have checked that this is OK, but assert it
            ;; again to be certain. Heap corruption is the greater evil versus a minor
            ;; inconvenience of not offering the RECKLESSLY-CONTINUE restart.
            (aver (= (calculate-extra-id-words depthoid) extra-id-words))
            #-64-bit (setf (wrapper-depthoid modify) (wrapper-depthoid wrapper)
                           (wrapper-length modify) (wrapper-length wrapper))
            (setf (layout-flags (wrapper-friend modify)) (layout-flags layout)
                  (wrapper-info modify) (wrapper-info wrapper))
            ;; Zero out the inherited ID values one word at a time.
            ;; This makes self-ID transiently disappear, but what else can we do?
            ;; It's may be in the wrong slot anyway, depending on whether depthoid changed.
            ;; The calculation of the min word count of 3 or 6 is done as
            ;;   (/ (- (1+ layout-id-vector-fixed-capacity) 2) number-of-ids-per-word)
            ;; which is surely more confusing than spelling it as 3 or 6.
            (dotimes (i (+ extra-id-words #+64-bit 3 #-64-bit 6))
              (%raw-instance-set/word (wrapper-friend modify)
                                      (+ (get-dsd-index sb-vm:layout id-word0) i)
                                      0))
            (set-layout-inherits modify inherits t id)
            (let ((dst (wrapper-friend modify))
                  (src (wrapper-friend wrapper)))
              (set-bitmap-from-layout dst src))
            (setf (wrapper-invalid modify) nil
                  (classoid-wrapper classoid) modify))
          (setf (wrapper-invalid wrapper) nil
                (classoid-wrapper classoid) wrapper))

      (dovector (super-wrapper (wrapper-inherits wrapper))
        (let ((super (wrapper-classoid super-wrapper)))
          (when (and (eq (classoid-state super) :sealed)
                     (not (get-subclassoid classoid super)))
            (warn "unsealing sealed class ~S in order to subclass it"
                  (classoid-name super))
            (setf (classoid-state super) :read-only))
          (add-subclassoid super classoid (or modify wrapper))))))

  (values)))

;;; Arrange the inherited layouts to appear at their expected depth,
;;; ensuring that hierarchical type tests succeed. Layouts with
;;; DEPTHOID >= 0 (i.e. hierarchical classes) are placed first,
;;; at exactly that index in the INHERITS vector. Then, non-hierarchical
;;; layouts are placed in remaining elements. Then, any still-empty
;;; elements are filled with their successors, ensuring that each
;;; element contains a valid layout.
;;;
;;; *** FIXME *** the preceding comment seems dubious, and I'm not sure whether
;;; to fix the code or the comment or both. The code works as-is, but is too hairy.
;;; I fail to see how "still-empty" elements can exist after filling in mandatory
;;; elements. It seems to anticipate being able to create a type whose INHERITS vector
;;; length exceeds depthoid, or, say, a type at depthoid 5 which inherits STREAM but
;;; might lack an entry at depth index 1 for example. As to why I think the comment
;;; is false: FILE-STREAM and STRING-STREAM each have depthoid 4, but their INHERITS
;;; vector has length 2. So they don't store elements that would be at index 2 and 3.
;;; Length less than depthoid is opposite of what the fill-in logic supports.
;;; How, in practice, could a user achieve such weird states as need this logic?
;;; If impossible, then simplify it.

;;; This reordering may destroy CPL ordering, so the inherits should
;;; not be read as being in CPL order.
(defun order-layout-inherits (layouts)
  (declare (simple-vector layouts))
  (let ((length (length layouts))
        (max-depth -1))
    (dotimes (i length)
      (let ((depth (wrapper-depthoid (svref layouts i))))
        (when (> depth max-depth)
          (setf max-depth depth))))
    (let* ((new-length (max (1+ max-depth) length))
           ;; KLUDGE: 0 here is the "uninitialized" element.  We need
           ;; to specify it explicitly for portability purposes, as
           ;; elements can be read before being set [ see below, "(EQL
           ;; OLD-LAYOUT 0)" ].  -- CSR, 2002-04-20
           (inherits (make-array new-length :initial-element 0)))
      (dotimes (i length)
        (let* ((layout (svref layouts i))
               (depth (wrapper-depthoid layout)))
          (unless (eql depth -1)
            (let ((old-layout (svref inherits depth)))
              (unless (or (eql old-layout 0) (eq old-layout layout))
                (error "layout depth conflict: ~S~%" layouts)))
            (setf (svref inherits depth) layout))))
      (do ((i 0 (1+ i))
           (j 0))
          ((>= i length))
        (declare (type index i j))
        (let* ((layout (svref layouts i))
               (depth (wrapper-depthoid layout)))
          (when (eql depth -1)
            (loop (when (eql (svref inherits j) 0)
                    (return))
                  (incf j))
            (setf (svref inherits j) layout))))
      (do ((i (1- new-length) (1- i)))
          ((< i 0))
        (declare (type fixnum i))
        (when (eql (svref inherits i) 0)
          (setf (svref inherits i) (svref inherits (1+ i)))))
      inherits)))

;;;; class precedence lists

;;; Topologically sort the list of objects to meet a set of ordering
;;; constraints given by pairs (A . B) constraining A to precede B.
;;; When there are multiple objects to choose, the tie-breaker
;;; function is called with both the list of object to choose from and
;;; the reverse ordering built so far.
(defun topological-sort (objects constraints tie-breaker)
  (declare (list objects constraints)
           (function tie-breaker))
  (let ((obj-info (make-hash-table :size (length objects)))
        (free-objs nil)
        (result nil))
    (loop for (obj1 . obj2) in constraints do
       (incf (first (ensure-gethash obj2 obj-info (list 0))))
       (push obj2 (rest (ensure-gethash obj1 obj-info (list 0)))))
    (dolist (obj objects)
      (let ((info (gethash obj obj-info)))
        (when (or (not info) (zerop (first info)))
          (push obj free-objs))))
    (loop
     (flet ((next-result (obj)
              (push obj result)
              (dolist (successor (rest (gethash obj obj-info)))
                (let* ((successor-info (gethash successor obj-info))
                       (count (1- (first successor-info))))
                  (setf (first successor-info) count)
                  (when (zerop count)
                    (push successor free-objs))))))
       (cond ((endp free-objs)
              (dohash ((obj info) obj-info)
                (unless (zerop (first info))
                  (error "Topological sort failed due to constraint on ~S."
                         obj)))
              (return (nreverse result)))
             ((endp (rest free-objs))
              (next-result (pop free-objs)))
             (t
              (let ((obj (funcall tie-breaker free-objs result)))
                (setf free-objs (remove obj free-objs))
                (next-result obj))))))))


;;; standard class precedence list computation
(defun std-compute-class-precedence-list (class)
  (let ((classes nil)
        (constraints nil))
    (labels ((note-class (class)
               (unless (member class classes)
                 (push class classes)
                 (let ((superclasses (classoid-direct-superclasses class)))
                   (do ((prev class)
                        (rest superclasses (rest rest)))
                       ((endp rest))
                     (let ((next (first rest)))
                       (push (cons prev next) constraints)
                       (setf prev next)))
                   (dolist (class superclasses)
                     (note-class class)))))
             (std-cpl-tie-breaker (free-classes rev-cpl)
               (dolist (class rev-cpl (first free-classes))
                 (let* ((superclasses (classoid-direct-superclasses class))
                        (intersection (intersection free-classes
                                                    superclasses)))
                   (when intersection
                     (return (first intersection)))))))
      (note-class class)
      (topological-sort classes constraints #'std-cpl-tie-breaker))))


;;; Return the layout for an object. This is the basic operation for
;;; finding out the "type" of an object, and is used for generic
;;; function dispatch. The standard doesn't seem to say as much as it
;;; should about what this returns for built-in objects. For example,
;;; it seems that we must return NULL rather than LIST when X is NIL
;;; so that GF's can specialize on NULL.
;;; x86-64 has a vop that implements this without even needing to place
;;; the vector of layouts in the constant pool of the containing code.
#-(or sb-xc-host (and compact-instance-header x86-64))
(progn
(declaim (inline wrapper-of))
(defun wrapper-of (x)
  (declare (optimize (speed 3) (safety 0)))
  (cond ((%instancep x) (%instance-wrapper x))
        ((funcallable-instance-p x) (%fun-wrapper x))
        ;; Compiler can dump literal layouts, which handily sidesteps
        ;; the question of when cold-init runs L-T-V forms.
        ((null x) #.(find-layout 'null))
        (t
         ;; Note that WIDETAG-OF is slightly suboptimal here and could be
         ;; improved - we've already ruled out some of the lowtags.
         (layout-friend
          (svref (load-time-value **primitive-object-layouts** t)
                 (widetag-of x)))))))

#-sb-xc-host
(progn
(declaim (inline classoid-of))
(defun classoid-of (object)
  "Return the class of the supplied object, which may be any Lisp object, not
   just a CLOS STANDARD-OBJECT."
  (wrapper-classoid (wrapper-of object))))


;;;; classoid namespace

(eval-when (#-sb-xc :compile-toplevel :load-toplevel :execute)
  (defun (setf find-classoid) (new-value name)
    (declare (type (or null classoid) new-value))
    (aver new-value)
    (with-world-lock ()
        (let ((cell (find-classoid-cell name :create t)))
          (ecase (info :type :kind name)
            ((nil))
            (:forthcoming-defclass-type
             ;; FIXME: Currently, nothing needs to be done in this case.
             ;; Later, when PCL is integrated tighter into SBCL, this
             ;; might need more work.
             nil)
            (:instance
             (aver cell)
             (let ((old-value (classoid-cell-classoid cell)))
               (aver old-value)
               ;; KLUDGE: The reason these clauses aren't directly
               ;; parallel is that we need to use the internal
               ;; CLASSOID structure ourselves, because we don't
               ;; have CLASSes to work with until PCL is built. In
               ;; the host, CLASSes have an approximately
               ;; one-to-one correspondence with the target
               ;; CLASSOIDs (as well as with the target CLASSes,
               ;; modulo potential differences with respect to
               ;; conditions).
               #+sb-xc-host
               (let ((old (cl:class-of old-value))
                     (new (cl:class-of new-value)))
                 (unless (eq old new)
                   (bug "Trying to change the metaclass of ~S from ~S to ~S in the ~
                            cross-compiler."
                        name (cl:class-name old) (cl:class-name new))))
               #-sb-xc-host
               (let ((old (classoid-of old-value))
                     (new (classoid-of new-value)))
                 (unless (eq old new)
                   (warn "Changing meta-class of ~S from ~S to ~S."
                         name (classoid-name old) (classoid-name new))))))
            (:primitive
             (error "Cannot redefine standard type ~
                     ~/sb-impl:print-type-specifier/." name))
            (:defined
             (warn "redefining DEFTYPE type to be a class: ~
                    ~/sb-ext:print-symbol-with-prefix/" name)
             (clear-info :type :expander name)
             (clear-info :type :source-location name)))

          (remhash name *forward-referenced-wrappers*)
          (%note-type-defined name)
          ;; FIXME: I'm unconvinced of the need to handle either of these.
          ;; Package locks preclude the latter, and in the former case,
          ;; once you've made some random thing into a :PRIMITIVE kind of type,
          ;; you've painted yourself into a corner - those types
          ;; elicit vociferous complaints if you try to redefine them.
          ;;
          ;; we need to handle things like
          ;;   (setf (find-class 'foo) (find-class 'integer))
          ;; and
          ;;   (setf (find-class 'integer) (find-class 'integer))
          (cond ((built-in-classoid-p new-value)
                 ;; But I can't figure out how to get assertions to pass
                 ;; without violation what would otherwise be invariants
                 ;; of the internal representation of types. This sucks.
                 (setf (info :type :kind name)
                       (or (info :type :kind name) :defined)))
                (t
                 (setf (info :type :kind name) :instance)))
          (setf (classoid-cell-classoid cell) new-value)
          (unless (eq (info :type :compiler-layout name)
                      (classoid-wrapper new-value))
            (setf (info :type :compiler-layout name)
                  (classoid-wrapper new-value)))))
    new-value)

  (defun %clear-classoid (name cell)
    (ecase (info :type :kind name)
      ((nil))
      (:defined)
      (:primitive
       (error "Attempt to remove :PRIMITIVE type: ~
              ~/sb-impl:print-type-specifier/" name))
      ((:forthcoming-defclass-type :instance)
       (when cell
         ;; Note: We cannot remove the classoid cell from the table,
         ;; since compiled code may refer directly to the cell, and
         ;; getting a different cell for a classoid with the same name
         ;; just would not do.

         ;; Remove the proper name of the classoid, if this was it.
         (let* ((classoid (classoid-cell-classoid cell))
                (proper-name (classoid-name classoid)))
           (when (eq proper-name name)
             (setf (classoid-name classoid) nil)))

         ;; Clear the cell.
         (setf (classoid-cell-classoid cell) nil
               (classoid-cell-pcl-class cell) nil))
       (clear-info :type :kind name)
       (clear-info :type :documentation name)
       (clear-info :type :compiler-layout name)
       (values-specifier-type-cache-clear)))))

(defun find-classoid-cell (name &key create)
  (let ((real-name (uncross name)))
    (cond ((info :type :classoid-cell real-name))
          (create
           (get-info-value-initializing :type :classoid-cell real-name
                                        (make-classoid-cell real-name))))))

;;; Return the classoid with the specified NAME. If ERRORP is false,
;;; then NIL is returned when no such class exists.
(defun find-classoid (name &optional (errorp t))
  (declare (type symbol name))
  (let ((cell (find-classoid-cell name)))
    (cond ((and cell (classoid-cell-classoid cell)))
          (errorp
           (error 'simple-type-error
                  :datum nil
                  :expected-type 'class
                  :format-control "Class not yet defined: ~S"
                  :format-arguments (list name))))))

;;; Called when we are about to define NAME as a class meeting some
;;; predicate (such as a meta-class type test.) The first result is
;;; always of the desired class. The second result is any existing
;;; LAYOUT for this name.
;;;
;;; Again, this should be compiler-only, but easier to make this
;;; thread-safe.
(defun insured-find-classoid (name predicate constructor)
  (declare (type function predicate)
           (type (or function symbol) constructor))
  (let ((table *forward-referenced-wrappers*))
    (with-system-mutex ((hash-table-lock table))
      (let* ((old (find-classoid name nil))
             (res (if (and old (funcall predicate old))
                      old
                      (funcall constructor :name name)))
             (found (or (gethash name table)
                        (when old (classoid-wrapper old)))))
        (when found
          (setf (wrapper-classoid found) res))
        (values res found)))))

;;; If the classoid has a proper name, return the name, otherwise return
;;; the classoid.
(defun classoid-proper-name (classoid)
  (declare (type classoid classoid))
  (let ((name (classoid-name classoid)))
    (if (and name (eq (find-classoid name nil) classoid))
        name
        classoid)))

;;;; CLASS type operations

;; CLASSOID-ENUMERABLE-P is referenced during compile by !DEFINE-TYPE-CLASS.
;; But don't redefine it when building the target since we've already
;; got a perfectly good definition loaded for the host.
(eval-when (#-sb-xc :compile-toplevel :load-toplevel :execute)
  ;; Actually this definition makes very little sense because
  ;;     (TYPE-ENUMERABLE (FIND-CLASSOID 'CHARACTER)) => T
  ;; but (TYPE-ENUMERABLE (SPECIFIER-TYPE 'CHARACTER)) => NIL.
  ;; You should never see the CLASSOID used as a type though,
  ;; at least not from parsing and set operations.
  ;; On a related note, (TYPE-ENUMERABLE (FIND-CLASSOID 'NULL))
  ;; should probably be T, but you'll never see that type either.
  ;; Perhaps a better definition of this function would be
  ;;   (if (classoid-translation x) (bug "enumerable-p classoid?") nil)
  (defun classoid-enumerable-p (x) (eq (classoid-name x) 'character)))
(define-type-class classoid :enumerable #'classoid-enumerable-p
                    :might-contain-other-types nil)

(defun classoid-inherits-from (sub super-or-name)
  (declare (type classoid sub)
           (type (or symbol classoid) super-or-name))
  (let ((super (if (symbolp super-or-name)
                   (find-classoid super-or-name)
                   super-or-name)))
    (find (classoid-wrapper super)
          (wrapper-inherits (classoid-wrapper sub)))))

;;; We might be passed classoids with invalid layouts; in any pairwise
;;; class comparison, we must ensure that both are valid before
;;; proceeding.
(defun %ensure-classoid-valid (classoid layout error-context)
  (declare (ignorable error-context)) ; not used on host
  (aver (eq classoid (wrapper-classoid layout)))
  (or (not (wrapper-invalid layout))
      ;; Avoid accidentally reaching code that can't work.
      #+sb-xc-host (bug "(TYPEP x 'STANDARD-CLASSOID) can't be tested")
      #-sb-xc-host
      (if (typep classoid 'standard-classoid)
          (let ((class (classoid-pcl-class classoid)))
            (cond
              ((sb-mop:class-finalized-p class)
               (sb-pcl::%force-cache-flushes class)
               t)
              ((sb-pcl::class-has-a-forward-referenced-superclass-p class)
               (when error-context
                 (bug "~@<Invalid class ~S with forward-referenced superclass ~
                       ~S in ~A.~%~:@>"
                      class
                      (sb-pcl::class-has-a-forward-referenced-superclass-p class)
                      error-context))
               nil)
              (t
               (sb-mop:finalize-inheritance class)
               t)))
          (bug "~@<Don't know how to ensure validity of ~S (not a STANDARD-CLASSOID) ~
                for ~A.~%~:@>"
               classoid (or error-context 'subtypep)))))

(defun %ensure-both-classoids-valid (class1 class2 &optional errorp)
  (do ((layout1 (classoid-wrapper class1) (classoid-wrapper class1))
       (layout2 (classoid-wrapper class2) (classoid-wrapper class2))
       (i 0 (+ i 1)))
      ((and (not (wrapper-invalid layout1)) (not (wrapper-invalid layout2)))
       t)
    (aver (< i 2))
    (unless (and (%ensure-classoid-valid class1 layout1 errorp)
                 (%ensure-classoid-valid class2 layout2 errorp))
      (return-from %ensure-both-classoids-valid nil))))

;;; Simple methods for TYPE= and SUBTYPEP should never be called when
;;; the two classes are equal, since there are EQ checks in those
;;; operations.
(define-type-method (classoid :simple-=) (type1 type2)
  (aver (not (eq type1 type2)))
  (values nil t))

(define-type-method (classoid :simple-subtypep) (class1 class2)
  (aver (not (eq class1 class2)))
  (with-world-lock () ; FIXME: why such coarse lock granularity here?
    (if (%ensure-both-classoids-valid class1 class2)
        (let ()
          (if (get-subclassoid class1 class2)
              (values t t)
              (if (and (typep class1 'standard-classoid)
                       (typep class2 'standard-classoid)
                       (or (sb-pcl::class-has-a-forward-referenced-superclass-p
                            (classoid-pcl-class class1))
                           (sb-pcl::class-has-a-forward-referenced-superclass-p
                            (classoid-pcl-class class2))))
                  ;; If there's a forward-referenced class involved we don't know for sure.
                  ;; (There are cases which we /could/ figure out, but that doesn't seem
                  ;; to be required or important, really.)
                  (values nil nil)
                  (values nil t))))
        (values nil nil))))

;;; When finding the intersection of a sealed class and some other
;;; class (not hierarchically related) the intersection is the union
;;; of the currently shared subclasses.
(defun sealed-class-intersection2 (sealed other)
  (declare (type classoid sealed other))
  (let ((s-sub (classoid-subclasses sealed))
        (o-sub (classoid-subclasses other)))
    (if (and s-sub o-sub)
        ;; FIXME: should we put more locking here?
        ;; [contrast with define-type-method (classoid :simple-subtypep)]
        (collect ((res *empty-type* type-union))
          (do-subclassoids ((subclass wrapper) sealed)
            (declare (ignore wrapper))
            (when (get-subclassoid subclass other)
              (res (specifier-type subclass))))
          (res))
        *empty-type*)))

(define-type-method (classoid :simple-intersection2) (class1 class2)
  (declare (type classoid class1 class2))
  (with-world-lock ()
    (%ensure-both-classoids-valid class1 class2 "type intersection")
    (cond ((eq class1 class2)
           class1)
          ;; If one is a subclass of the other, then that is the
          ;; intersection.
          ((get-subclassoid class1 class2) class1)
          ((get-subclassoid class2 class1) class2)
          ;; Otherwise, we can't in general be sure that the
          ;; intersection is empty, since a subclass of both might be
          ;; defined. But we can eliminate it for some special cases.
          ((or (structure-classoid-p class1)
               (structure-classoid-p class2))
           ;; No subclass of both can be defined.
           *empty-type*)
          ((eq (classoid-state class1) :sealed)
           ;; checking whether a subclass of both can be defined:
           (sealed-class-intersection2 class1 class2))
          ((eq (classoid-state class2) :sealed)
           ;; checking whether a subclass of both can be defined:
           (sealed-class-intersection2 class2 class1))
          ;; If exactly one of CLASS{1,2} is a CONDITION-CLASSOID,
          ;; there can be no intersection: sub-/superclass relations
          ;; between CONDITION-CLASSOIDs and other CLASSOIDs are not
          ;; possible and a CONDITION-CLASSOIDs cannot be changed into
          ;; different CLASSOIDs.
          ((let ((c1 (condition-classoid-p class1))
                 (c2 (condition-classoid-p class2)))
             (or (and c1 (not c2)) (and (not c1) c2)))
           *empty-type*)
          (t
           ;; uncertain, since a subclass of both might be defined
           nil))))

;;; KLUDGE: we need this to deal with the special-case INSTANCE and
;;; FUNCALLABLE-INSTANCE types (which used to be CLASSOIDs until CSR
;;; discovered that this was incompatible with the MOP class
;;; hierarchy).  See NAMED :COMPLEX-SUBTYPEP-ARG2
(declaim (type cons **non-instance-classoid-types**))
(defglobal **non-instance-classoid-types**
  '(symbol system-area-pointer weak-pointer code-component
    #-(or x86 x86-64 arm64) lra
    fdefn random-class))

(defun classoid-non-instance-p (classoid)
  (declare (type classoid classoid))
  (member classoid **non-instance-classoid-types**
          :key #'find-classoid))

;;; KLUDGE: we need this because of the need to represent
;;; intersections of two classes, even when empty at a given time, as
;;; uncanonicalized intersections because of the possibility of later
;;; defining a subclass of both classes.  The necessity for changing
;;; the default return value from SUBTYPEP to NIL, T if no alternate
;;; method is present comes about because, unlike the other places we
;;; use INVOKE-COMPLEX-SUBTYPEP-ARG1-METHOD, in HAIRY methods and the
;;; like, classes are in their own hierarchy with no possibility of
;;; mixtures with other type classes.
(define-type-method (classoid :complex-subtypep-arg2) (type1 class2)
  (if (and (intersection-type-p type1)
           (> (count-if #'classoid-p (intersection-type-types type1)) 1))
      (values nil nil)
      (invoke-complex-subtypep-arg1-method type1 class2 nil t)))

(define-type-method (classoid :negate) (type) (make-negation-type type))

(define-type-method (classoid :unparse) (type)
  (classoid-proper-name type))

;;;; built-in classes

;;; The BUILT-IN-CLASSES list is a data structure which configures the
;;; creation of all the built-in classes. It contains all the info
;;; that we need to maintain the mapping between classes, compile-time
;;; types and run-time type codes. These options are defined:
;;;
;;; :TRANSLATION (default none)
;;;     When this class is "parsed" as a type specifier, it is
;;;     translated into the specified internal type representation,
;;;     rather than being left as a class. This is used for types
;;;     which we want to canonicalize to some other kind of type
;;;     object because in general we want to be able to include more
;;;     information than just the class (e.g. for numeric types.)
;;;
;;; :STATE (default :SEALED)
;;;     The value of CLASS-STATE which we want on completion,
;;;     indicating whether subclasses can be created at run-time.
;;;
;;; :HIERARCHICAL-P (default T unless any of the inherits are non-hierarchical)
;;;     True if we can assign this class a unique inheritance depth.
;;;
;;; :CODES (default none)
;;;     Run-time type codes which should be translated back to this
;;;     class by CLASS-OF. Unspecified for abstract classes.
;;;
;;; :INHERITS (default this class and T)
;;;     The class-precedence list for this class, with this class and
;;;     T implicit.
;;;
;;; :DIRECT-SUPERCLASSES (default to head of CPL)
;;;     List of the direct superclasses of this class.
;;;
;;; NB: not to be confused with SB-PCL::*BUILT-IN-CLASSES*
;;; (note the difference in spelling, to help keep things unconfusing)
#+sb-xc-host
(defvar *builtin-classoids*
   `((t :state :read-only :translation t)
     (character :codes (,sb-vm:character-widetag)
                :translation (character-set)
                :prototype-form (code-char 42))
     (symbol :codes (,sb-vm:symbol-widetag)
             :predicate symbolp
             :prototype-form '*)

     (system-area-pointer :codes (,sb-vm:sap-widetag)
                          :predicate system-area-pointer-p
                          :prototype-form (int-sap 0))
     (weak-pointer :codes (,sb-vm:weak-pointer-widetag)
                   :predicate weak-pointer-p
                   :prototype-form (make-weak-pointer 0))
     (code-component :codes (,sb-vm:code-header-widetag)
                     :predicate code-component-p
                     :prototype-form (fun-code-header #'identity))
     #-(or x86 x86-64 arm64) (lra :codes (,sb-vm:return-pc-widetag)
                                  :predicate lra-p
                                  ;; Make the PROTOTYPE slot unbound.
                                  :prototype-form sb-pcl:+slot-unbound+)
     (fdefn :codes (,sb-vm:fdefn-widetag)
            :predicate fdefn-p
            :prototype-form (find-or-create-fdefn 'sb-mop:class-prototype))
     (random-class ; used for unknown type codes
            ;; Make the PROTOTYPE slot unbound.
            :prototype-form sb-pcl:+slot-unbound+)
     (function
      :codes (,sb-vm:closure-widetag ,sb-vm:simple-fun-widetag)
      :predicate functionp
      :state :read-only
      :prototype-form #'identity)

     (number :translation number :prototype-form 0)
     (complex
      :translation complex
      :inherits (number)
      :codes (,sb-vm:complex-widetag)
      :prototype-form ,(complex 0 1))
     (complex-single-float
      :translation (complex single-float)
      :inherits (complex number)
      :codes (,sb-vm:complex-single-float-widetag)
      :prototype-form ,(complex $0f0 $0f0))
     (complex-double-float
      :translation (complex double-float)
      :inherits (complex number)
      :codes (,sb-vm:complex-double-float-widetag)
      :prototype-form ,(complex $0d0 $0d0))
     #+long-float
     (complex-long-float
      :translation (complex long-float)
      :inherits (complex number)
      :codes (,sb-vm:complex-long-float-widetag)
      :prototype-form ,(complex $0L0 $0L0))
     #+sb-simd-pack
     (simd-pack
      :translation simd-pack
      :codes (,sb-vm:simd-pack-widetag)
      :prototype-form (%make-simd-pack-ub64 42 42))
     #+sb-simd-pack-256
     (simd-pack-256
      :translation simd-pack-256
      :codes (,sb-vm:simd-pack-256-widetag)
      :prototype-form
      ;; KLUDGE: doesn't work without AVX2 support from the CPU
      ;; (%make-simd-pack-256-ub64 42 42 42 42)
      sb-pcl:+slot-unbound+)
     (real :translation real :inherits (number) :prototype-form 0)
     (float :translation float :inherits (real number) :prototype-form $0f0)
     (single-float
      :translation single-float
      :inherits (float real number)
      :codes (,sb-vm:single-float-widetag)
      :prototype-form $0f0)
     (double-float
      :translation double-float
      :inherits (float real number)
      :codes (,sb-vm:double-float-widetag)
      :prototype-form $0d0)
     #+long-float
     (long-float
      :translation long-float
      :inherits (float real number)
      :codes (,sb-vm:long-float-widetag)
      :prototype-form $0L0)
     (rational
      :translation rational :inherits (real number) :prototype-form 0)
     (ratio
      :translation (and rational (not integer))
      :inherits (rational real number)
      :codes (,sb-vm:ratio-widetag)
      :prototype-form 1/42)
     (integer
      :translation integer :inherits (rational real number) :prototype-form 0)
     (fixnum
      :translation (integer ,most-negative-fixnum ,most-positive-fixnum)
      :inherits (integer rational real number)
      :codes ,(mapcar #'symbol-value sb-vm::fixnum-lowtags)
      :prototype-form 42)
     (bignum
      :translation (and integer (not fixnum))
      :inherits (integer rational real number)
      :codes (,sb-vm:bignum-widetag)
      :prototype-form ,(1+ most-positive-fixnum))

     (array :translation array :codes (,sb-vm:complex-array-widetag)
            :hierarchical-p nil
            :prototype-form (make-array nil :adjustable t))
     (simple-array
      :translation simple-array :codes (,sb-vm:simple-array-widetag)
      :inherits (array)
      :prototype-form (make-array nil))
     (sequence
      :translation (or cons (member nil) vector extended-sequence)
      :state :read-only
      :depth 1)
     (vector
      :translation vector :codes (,sb-vm:complex-vector-widetag)
      :direct-superclasses (array sequence)
      :inherits (array sequence)
      :prototype-form (make-array 0 :adjustable t))
     (simple-vector
      :translation simple-vector :codes (,sb-vm:simple-vector-widetag)
      :direct-superclasses (vector simple-array)
      :inherits (vector simple-array array sequence)
      :prototype-form (make-array 0))
     (bit-vector
      :translation bit-vector :codes (,sb-vm:complex-bit-vector-widetag)
      :inherits (vector array sequence)
      :prototype-form (make-array 0 :element-type 'bit :fill-pointer t))
     (simple-bit-vector
      :translation simple-bit-vector :codes (,sb-vm:simple-bit-vector-widetag)
      :direct-superclasses (bit-vector simple-array)
      :inherits (bit-vector vector simple-array
                 array sequence)
      :prototype-form #*)
     (string
      :translation string
      :direct-superclasses (vector)
      :inherits (vector array sequence)
      :prototype-form "")
     (simple-string
      :translation simple-string
      :direct-superclasses (string simple-array)
      :inherits (string vector simple-array array sequence)
      :prototype-form "")
     (vector-nil
      :translation (vector nil)
      :inherits (vector array sequence)
      :prototype-form (make-array 0 :element-type 'nil :fill-pointer t))
     ;; This name is imperfect. It should be SIMPLE-RANK1-ARRAY-NIL
     ;; to clearly convey that the dimensions are '(*) and not '*.
     (simple-array-nil
      :translation (simple-array nil (*))
      :codes (,sb-vm:simple-array-nil-widetag)
      :direct-superclasses (vector-nil)
      :inherits (vector-nil vector simple-array array sequence)
      :prototype-form (make-array 0 :element-type 'nil))
     (base-string
      :translation base-string
      :codes (,sb-vm:complex-base-string-widetag)
      :direct-superclasses (string)
      :inherits (string vector array sequence)
      :prototype-form (make-array 0 :element-type 'base-char :fill-pointer t))
     (simple-base-string
      :translation simple-base-string
      :codes (,sb-vm:simple-base-string-widetag)
      :direct-superclasses (base-string simple-string)
      :inherits (base-string simple-string string vector simple-array
                 array sequence)
      :prototype-form (make-array 0 :element-type 'base-char))
     #+sb-unicode
     (character-string
      :translation (vector character)
      :codes (,sb-vm:complex-character-string-widetag)
      :direct-superclasses (string)
      :inherits (string vector array sequence)
      :prototype-form (make-array 0 :element-type 'character :fill-pointer t))
     #+sb-unicode
     (simple-character-string
      :translation (simple-array character (*))
      :codes (,sb-vm:simple-character-string-widetag)
      :direct-superclasses (character-string simple-string)
      :inherits (character-string simple-string string vector simple-array
                 array sequence)
      :prototype-form (make-array 0 :element-type 'character))
     (list
      :translation (or cons (member nil))
      :inherits (sequence)
      :prototype-form 'nil)
     (cons
      :codes (,sb-vm:list-pointer-lowtag)
      :translation cons
      :inherits (list sequence)
      :prototype-form (cons nil nil))
     (null
      :translation (member nil)
      :inherits (symbol list sequence)
      :direct-superclasses (symbol list)
      :prototype-form 'nil)

     (sb-pcl::slot-object
      :translation (or structure-object standard-object condition)
      :predicate slot-object-p
      :hierarchical-p nil
      :state :read-only
      :prototype-form (make-defstruct-description t 'arbitrary))

     ;; KLUDGE: the length must match the subsequent defstruct.
     (pathname :depth 1
               :predicate pathnamep
               :length ,(+ 7 sb-vm:instance-data-start)
               :prototype-form (make-trivial-default-pathname))
     (logical-pathname :depth 2
                       :predicate logical-pathname-p
                       :length ,(+ 7 sb-vm:instance-data-start)
                       :prototype-form (make-trivial-default-logical-pathname)
                       :inherits (pathname))

     ;; These last few are strange. STREAM has only T as an ancestor,
     ;; so you'd think it would be at depth 1. FILE- and STRING-STREAM
     ;; each have STREAM and T as ancestors, so you'd think they'd be at depth
     ;; 1 greater than STREAM, instead of 2 greater. But changing any of
     ;; these to the "obvious" value makes various type checks go wrong.
     ;;
     ;; Essentially the hardwiring corresponds to the indices of the
     ;; respective types in the inherits vector for FD-STREAM.
     ;;  * (wrapper-inherits (find-layout 'fd-stream))
     ;;  #(#<LAYOUT for T {50300003}>
     ;;    #<LAYOUT for STRUCTURE-OBJECT {50300103}>
     ;;    #<LAYOUT for STREAM {50301003}>
     ;;    #<LAYOUT for ANSI-STREAM {50301183}>
     ;;    #<LAYOUT for FILE-STREAM {50303303}>)

     (stream
      :predicate streamp
      :state :read-only
      :depth 2)
     (file-stream
      :predicate file-stream-p
      :state :read-only
      :depth 4
      :inherits (stream))
     (string-stream
      :predicate string-stream-p
      :state :read-only
      :depth 4
      :inherits (stream))
     ,@(loop for x across sb-vm:*specialized-array-element-type-properties*
             unless (member (sb-vm:saetp-specifier x) '(t character base-char nil bit))
             collect
             ;; I'm not sure if it's an accident that there are distinct SB-KERNEL
             ;; versus SB-VM symbols for the specialized arrays. The former are types
             ;; in the language, and the latter are primitive object types,
             ;; but istm they should be designated by the same symbols.
             `(,(intern (string (sb-vm:saetp-primitive-type-name x)) *package*)
               :translation (simple-array ,(sb-vm:saetp-specifier x) (*))
               :codes (,(sb-vm:saetp-typecode x))
               :direct-superclasses (vector simple-array)
               :inherits (vector simple-array array sequence)
               :prototype-form
               (logically-readonlyize
                (make-array 0 :element-type ',(sb-vm:saetp-specifier x)))))))

(eval-when (#-sb-xc-host :compile-toplevel)
  (defun compute-builtin-classoids ()
    (mapcar (lambda (x)
              (let* ((name (car x))
                     (classoid (find-classoid name))
                     (translation (built-in-classoid-translation classoid))
                     (predicate
                      (if (member name '(t random-class))
                          'error
                          (or (getf (cdr x) :predicate)
                              (sb-c::backend-type-predicate translation)))))
                (assert predicate)
                ;; destructuring-bind will see the first :translation
                ;; keyword; we don't need to delete the other one.
                (list* name :predicate predicate :translation translation (cdr x))))
            *builtin-classoids*)))

;;; The read interceptor has to be disabled to avoid infinite recursion on CTYPEs
(eval-when (:compile-toplevel) (setq sb-cold::*choke-on-host-irrationals* nil))
#-sb-xc-host
(define-load-time-global *builtin-classoids* nil)
#-sb-xc-host
(!cold-init-forms
 (setq *builtin-classoids* '#.(compute-builtin-classoids)))
(eval-when (:compile-toplevel) (setq sb-cold::*choke-on-host-irrationals* t))

;;; See also src/code/type-init.lisp where we finish setting up the
;;; translations for built-in types.
(!cold-init-forms
 (dolist (x *builtin-classoids*)
   #-sb-xc-host (/show0 "at head of loop over *BUILTIN-CLASSOIDS*")
   (destructuring-bind
       (name &key
               (translation nil trans-p)
               predicate
               inherits
               codes
               state
               depth
               (length 0)
               prototype-form
               (hierarchical-p t) ; might be modified below
               (direct-superclasses (if inherits
                                        (list (car inherits))
                                        '(t))))
       x
     (declare (ignorable codes state translation trans-p predicate))
     ;; instance metatypes and T don't need a prototype, everything else does
     (unless (or prototype-form depth (eq name 't))
       (error "Missing prototype in ~S" x))
     (let* ((pred-fn (if (fboundp predicate) (symbol-function predicate) #'error))
            (inherits-list (if (eq name t)
                               ()
                               (cons t (reverse inherits))))
            (classoid
              (acond #-sb-xc-host ; genesis dumps some classoid literals
                     ((find-classoid name nil)
                      (%instance-set it (get-dsd-index built-in-classoid predicate)
                                     pred-fn)
                      ;; Unseal it so that REGISTER-LAYOUT doesn't warn
                      (setf (classoid-state it) nil)
                      it)
                     (t
                      (setf (classoid-cell-classoid
                             (find-classoid-cell name :create t))
                            (!make-built-in-classoid
                             :%bits (pack-ctype-bits classoid name)
                             :name name
                             :translation #+sb-xc-host (if trans-p :initializing nil)
                                          #-sb-xc-host translation
                             :allow-other-keys t :predicate pred-fn
                             :direct-superclasses
                             (if (eq name t)
                                 nil
                                 (mapcar #'find-classoid
                                         direct-superclasses))))))))
       (setf (info :type :kind name) :primitive)
       #+sb-xc-host
       (unless trans-p
         (setf (info :type :builtin name) classoid))
       #-sb-xc-host (setf (info :type :builtin name) (or translation classoid))
       (let* ((inherits-vector
                (map 'simple-vector
                     (lambda (x)
                       (let ((super-layout
                               (classoid-wrapper (find-classoid x))))
                         (when (minusp (wrapper-depthoid super-layout))
                           (setf hierarchical-p nil))
                         super-layout))
                     inherits-list))
              (depthoid (if hierarchical-p
                            (or depth (length inherits-vector))
                            -1)))
         (register-layout (load-layout name
                                       depthoid
                                       inherits-vector
                                       length
                                       +layout-all-tagged+
                                       0) ; flags
                          :invalidate nil)))))
 (/show0 "done with loop over *!BUILTIN-CLASSOIDS*"))

;;; Now that we have set up the class hierarchy, seal the sealed
;;; classes. This must be done after the subclasses have been set up.
(!cold-init-forms
  (dolist (x *builtin-classoids*)
    (destructuring-bind (name &key (state :sealed) &allow-other-keys) x
      (setf (classoid-state (find-classoid name)) state))))

;;;; class definition/redefinition

;;; This is to be called whenever we are altering a class.
#+sb-xc-host
(defun %modify-classoid (classoid) (bug "MODIFY-CLASSOID ~S" classoid))
#-sb-xc-host
(defun %modify-classoid (classoid)
  (clear-type-caches)
  (awhen (classoid-state classoid)
    ;; FIXME: This should probably be CERROR.
    (warn "making ~(~A~) class ~S writable" it (classoid-name classoid))
    (setf (classoid-state classoid) nil)))

;;; Mark LAYOUT as invalid. This is called only on CONDITION and STRUCTURE
;;; subtypes when redefining incompatibly. PCL objects use invalidate-wrapper.
;;; Remove class from all superclasses
;;; too (might not be registered, so might not be in subclasses of the
;;; nominal superclasses.)  We set the layout-clos-hash slots to 0 to
;;; invalidate the wrappers for specialized dispatch functions, which
;;; use those slots as indexes into tables.
(defun %invalidate-layout (wrapper)
  (declare (type wrapper wrapper))
  #+sb-xc-host (error "Can't invalidate layout ~S" wrapper)
  #-sb-xc-host
  (progn
    (setf (wrapper-invalid wrapper) t)
    ;; Ensure that the INVALID slot conveying ancillary data describing the
    ;; invalidity reason is published before causing the invalid layout trap.
    (sb-thread:barrier (:write))
    #+metaspace (setf (layout-clos-hash (wrapper-friend wrapper)) 0)
    (setf (wrapper-clos-hash wrapper) 0)
    (let ((inherits (wrapper-inherits wrapper))
          (classoid (wrapper-classoid wrapper)))
      (%modify-classoid classoid)
      (dovector (super inherits)
        (remove-subclassoid classoid (wrapper-classoid super)))))
  (values))

;;;; cold loading initializations

;;; FIXME: It would be good to arrange for this to be called when the
;;; cross-compiler is being built, not just when the target Lisp is
;;; being cold loaded. Perhaps this could be moved to its own file
;;; late in the build-order.lisp-expr sequence, and be put in
;;; !COLD-INIT-FORMS there?
(defun !class-finalize ()
  (dohash ((name wrapper) *forward-referenced-wrappers*)
    (let ((class (find-classoid name nil)))
      (cond ((not class)
             (error "How is there no classoid for ~S ?" name))
            ((eq (classoid-wrapper class) wrapper)
             (remhash name *forward-referenced-wrappers*))
            (t
             (error "Something strange with forward layout for ~S:~%  ~S"
                    name wrapper))))))

(!defun-from-collected-cold-init-forms !classes-cold-init)
