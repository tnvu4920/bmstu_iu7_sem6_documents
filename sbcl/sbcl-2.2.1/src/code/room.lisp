;;;; heap-grovelling memory usage stuff

;;;; This software is part of the SBCL system. See the README file for
;;;; more information.
;;;;
;;;; This software is derived from the CMU CL system, which was
;;;; written at Carnegie Mellon University and released into the
;;;; public domain. The software is in the public domain and is
;;;; provided with absolutely no warranty. See the COPYING and CREDITS
;;;; files for more information.

(in-package "SB-VM")


;;;; type format database

(defstruct (room-info (:constructor make-room-info (mask name kind))
                      (:copier nil))
    ;; the mask applied to HeaderValue to compute object size
    (mask 0 :type (and fixnum unsigned-byte) :read-only t)
    ;; the name of this type
    (name nil :type symbol :read-only t)
    ;; kind of type (how to reconstitute an object)
    (kind nil
          :type (member :other :closure :instance :list :code :fdefn)
          :read-only t))
(declaim (freeze-type room-info))

(defun room-info-type-name (info)
    (if (specialized-array-element-type-properties-p info)
        (saetp-primitive-type-name info)
        (room-info-name info)))

(defconstant tiny-boxed-size-mask #xFF)
(defun compute-room-infos ()
  (let ((infos (make-array 256 :initial-element nil))
        (default-size-mask (mask-field (byte 23 0) -1)))
    (dolist (obj *primitive-objects*)
      (let ((widetag (primitive-object-widetag obj))
            (lowtag (primitive-object-lowtag obj))
            (name (primitive-object-name obj)))
        (when (and (eq lowtag 'other-pointer-lowtag)
                   (not (member widetag '(t nil))))
          (setf (svref infos (symbol-value widetag))
                (case name
                 (fdefn  (make-room-info 0 name :fdefn))
                 (symbol (make-room-info tiny-boxed-size-mask name :other))
                 (t      (make-room-info default-size-mask name :other)))))))

    (let ((info (make-room-info tiny-boxed-size-mask 'array-header :other)))
      (dolist (code (list #+sb-unicode complex-character-string-widetag
                          complex-base-string-widetag simple-array-widetag
                          complex-bit-vector-widetag complex-vector-widetag
                          complex-array-widetag))
        (setf (svref infos code) info)))

    (setf (svref infos bignum-widetag)
          ;; Lose 1 more bit than n-widetag-bits because fullcgc robs 1 bit,
          ;; not that this is expected to work concurrently with gc.
          (make-room-info (ash most-positive-word (- (1+ n-widetag-bits)))
                          'bignum :other))
    (setf (svref infos filler-widetag)
          (make-room-info (ash most-positive-word (- (1+ n-widetag-bits)))
                          'filler :other))

    (setf (svref infos closure-widetag)
          (make-room-info short-header-max-words 'closure :closure))

    (dotimes (i (length *specialized-array-element-type-properties*))
      (let ((saetp (aref *specialized-array-element-type-properties* i)))
        (when (saetp-specifier saetp) ;; SIMPLE-ARRAY-NIL is a special case.
          (setf (svref infos (saetp-typecode saetp)) saetp))))

    ;; This one is here for completeness only- LENGTH does not imply size.
    (setf (svref infos simple-array-nil-widetag)
          (make-room-info 0 'simple-array-nil :other))

    (setf (svref infos code-header-widetag)
          (make-room-info 0 'code :code))

    (setf (svref infos instance-widetag)
          (make-room-info 0 'instance :instance))

    (setf (svref infos funcallable-instance-widetag)
          (make-room-info short-header-max-words 'funcallable-instance :closure))

    (let ((cons-info (make-room-info 0 'cons :list)))
      ;; A cons consists of two words, both of which may be either a
      ;; pointer or immediate data.  According to the runtime this means
      ;; either a fixnum, a character, an unbound-marker, a single-float
      ;; on a 64-bit system, or a pointer.
      (dotimes (i (ash 1 (- n-widetag-bits n-fixnum-tag-bits)))
        (setf (svref infos (ash i n-fixnum-tag-bits)) cons-info))

      (dotimes (i (ash 1 (- n-widetag-bits n-lowtag-bits)))
        (setf (svref infos (logior (ash i n-lowtag-bits) instance-pointer-lowtag))
              cons-info)
        (setf (svref infos (logior (ash i n-lowtag-bits) list-pointer-lowtag))
              cons-info)
        (setf (svref infos (logior (ash i n-lowtag-bits) fun-pointer-lowtag))
              cons-info)
        (setf (svref infos (logior (ash i n-lowtag-bits) other-pointer-lowtag))
              cons-info))

      (setf (svref infos character-widetag) cons-info)

      (setf (svref infos unbound-marker-widetag) cons-info)

      ;; Single-floats are immediate data on 64-bit systems.
      #+64-bit (setf (svref infos single-float-widetag) cons-info))

    infos))

(define-load-time-global *room-info* (compute-room-infos))
(declaim (type (simple-vector 256) *room-info*))

(defconstant-eqx +heap-spaces+
  '((:dynamic   "Dynamic space"   dynamic-usage)
    #+immobile-space
    (:immobile  "Immobile space"  sb-kernel::immobile-space-usage)
    (:read-only "Read-only space" sb-kernel::read-only-space-usage)
    (:static    "Static space"    sb-kernel::static-space-usage))
  #'equal)

(defconstant-eqx +stack-spaces+
  '((:control-stack "Control stack" sb-kernel::control-stack-usage)
    (:binding-stack "Binding stack" sb-kernel::binding-stack-usage))
  #'equal)

(defconstant-eqx +all-spaces+ (append +heap-spaces+ +stack-spaces+) #'equal)

(defconstant-eqx +heap-space-keywords+ (mapcar #'first +heap-spaces+) #'equal)
(deftype spaces () `(member . ,+heap-space-keywords+))


;;;; MAP-ALLOCATED-OBJECTS

;;; Return the lower limit and current free-pointer of SPACE as fixnums
;;; whose raw bits (at the register level) represent a pointer.
;;; This makes it "off" by a factor of (EXPT 2 N-FIXNUM-TAG-BITS) - and/or
;;; possibly negative - if you look at the value in Lisp,
;;; but avoids potentially needing a bignum on 32-bit machines.
;;; 64-bit machines have no problem since most current generation CPUs
;;; use an address width that is narrower than 64 bits.
;;; This function is private because of the wacky representation.
(defun %space-bounds (space)
  (macrolet ((bounds (a b) `(values (%make-lisp-obj ,a) (%make-lisp-obj ,b))))
    (ecase space
      (:static
       (bounds static-space-start
               (sap-int *static-space-free-pointer*)))
      (:read-only
       (bounds read-only-space-start
               (sap-int *read-only-space-free-pointer*)))
      #+immobile-space
      (:fixed
       (bounds fixedobj-space-start
               (sap-int *fixedobj-space-free-pointer*)))
      #+immobile-space
      (:variable
       (bounds varyobj-space-start
               (sap-int *varyobj-space-free-pointer*)))
      (:dynamic
       (bounds (current-dynamic-space-start)
               (sap-int (dynamic-space-free-pointer)))))))

;;; Return the total number of bytes used in SPACE.
(defun space-bytes (space)
  (if (eq space :immobile)
      (+ (space-bytes :immobile-fixed)
         (space-bytes :immobile-variable))
      (multiple-value-bind (start end) (%space-bounds space)
        (ash (- end start) n-fixnum-tag-bits))))

(defun instance-length (instance) ; excluding header, not aligned to even
  ;; Add 1 if expressed length PLUS header (total number of words) would be
  ;; an even number, and the hash state bits indicate hashed-and-moved.
  (+ (%instance-length instance)
     ;; Compute 1 or 0 depending whether the instance was physically extended
     ;; by one word for the stable hash value. Extension occurs when and only when
     ;; the hash state is hashed-and-moved, and the apparent total number of words
     ;; inclusive of header (and exclusive of extension) is even. ANDing the least
     ;; significant bit of the payload size with HASH-SLOT-PRESENT arrives at the
     ;; desired boolean value. If apparent size is odd in hashed-and-moved state,
     ;; the physical size undergoes no change.
     (let ((header-word (instance-header-word instance)))
       (logand (ash header-word (- instance-length-shift))
               (ash header-word (- hash-slot-present-flag))
               1))))

;;; Macros not needed after this file (and avoids a redefinition warning this way)
(eval-when (:compile-toplevel)
(defmacro widetag@baseptr (sap)
  #+big-endian `(sap-ref-8 ,sap ,(1- n-word-bytes))
  #+little-endian `(sap-ref-8 ,sap 0))

(defmacro lispobj@baseptr (sap widetag)
  `(%make-lisp-obj
    (logior (sap-int ,sap)
            (logand (deref (extern-alien "widetag_lowtag" (array char 256)) ,widetag)
                    lowtag-mask)))))

;;; This uses the funny fixnum representation of ADDRESS. I'd like to change this
;;; to take a SAP but god forbid people are using it?
;;; DO NOT USE THIS! It is soon to be removed
(defun reconstitute-object (address)
  (let ((sap (descriptor-sap address)))
    (lispobj@baseptr sap (widetag@baseptr sap))))

;;; Iterate over all the objects in the contiguous block of memory
;;; with the low address at START and the high address just before
;;; END, calling FUN with the object, the object's type code, and the
;;; object's total size in bytes, including any header and padding.
;;; START and END are untagged, aligned memory addresses interpreted
;;; as FIXNUMs (unlike SAPs or tagged addresses, these will not cons).
(defun map-objects-in-range (fun start end &optional (strict-bound t))
  (declare (type function fun))
  (declare (dynamic-extent fun))
  (let ((start (descriptor-sap start))
        (end (descriptor-sap end)))
    (loop
     (if (sap>= start end) (return))
     (binding* ((widetag (widetag@baseptr start))
                (obj (lispobj@baseptr start widetag))
                ((typecode size)
                 ;; PRIMITIVE-OBJECT-SIZE works on conses, but they're exceptions already
                 ;; because of absence of a widetag, so may as well not call the sizer.
                 (if (listp obj)
                     (values list-pointer-lowtag (* 2 n-word-bytes))
                     (values widetag (primitive-object-size obj)))))
       ;; SIZE is surely a fixnum. Non-fixnum would imply at least
       ;; a 512MB object if 32-bit words, and is inconceivable if 64-bit.
       ;; But check to be sure.
       (aver (not (logtest (the fixnum size) lowtag-mask)))
       (unless (= typecode filler-widetag)
         (funcall fun obj typecode size))
       (setq start (sap+ start size))))
    (when strict-bound
     ;; If START is not eq to END, then we have blown past our endpoint.
      #+sb-devel
     (unless (sap= start end)
       ;; don't make things go more wrong than they already are.
       (alien-funcall (extern-alien "printf" (function void system-area-pointer))
                      (vector-sap #.(format nil "map-objects-in-range failure~%")))
       (alien-funcall (extern-alien "ldb_monitor" (function void))))
     #-sb-devel
     (aver (sap= start end)))))

;;; Access to the GENCGC page table for better precision in
;;; MAP-ALLOCATED-OBJECTS
#+gencgc
(progn
  (define-alien-type (struct page)
      (struct page
              ;; To cut down the size of the page table, the scan_start_offset
              ;; - a/k/a "start" - is measured in 4-byte integers regardless
              ;; of word size. This is fine for 32-bit address space,
              ;; but if 64-bit then we have to scale the value. Additionally
              ;; there is a fallback for when even the scaled value is too big.
              (start #+64-bit (unsigned 32) #-64-bit signed)
              ;; On platforms with small enough GC pages, this field
              ;; will be a short. On platforms with larger ones, it'll
              ;; be an int. It should probably never be an int.
              (words-used (unsigned
                           #.(if (typep gencgc-page-words '(unsigned-byte 16))
                                 16
                                 32)))
              (flags (unsigned 8))
              (gen (signed 8))))
  #+immobile-space
  (progn
    (define-alien-type (struct immobile-page)
        ;; ... and yet another place for Lisp to become out-of-sync with C.
        (struct immobile-page
                (flags (unsigned 8))
                (obj-spacing (unsigned 8))
                (obj-size (unsigned 8))
                (generations (unsigned 8))
                (free-index (unsigned 32))
                (page-link (unsigned 16))
                (prior-free-index (unsigned 16))))
    (define-alien-variable "fixedobj_pages" (* (struct immobile-page))))
  (declaim (inline find-page-index))
  (define-alien-routine ("ext_find_page_index" find-page-index)
    long (index unsigned))
  (define-alien-variable "next_free_page" sb-kernel::page-index-t)
  (define-alien-variable "page_table" (* (struct page))))

#+immobile-space
(progn
(deftype immobile-subspaces () '(member :fixed :variable))
(declaim (ftype (sfunction (function &rest immobile-subspaces) null)
                map-immobile-objects))
(defun map-immobile-objects (function &rest subspaces) ; Perform no filtering
  (declare (dynamic-extent function))
  (do-rest-arg ((subspace) subspaces)
    (multiple-value-bind (start end) (%space-bounds subspace)
      (map-objects-in-range function start end)))))

#|
MAP-ALLOCATED-OBJECTS is fundamentally unsafe to use if the user-supplied
function allocates anything. Consider what can happens when NEXT-FREE-PAGE
points to a partially filled page, and one more object is created extending
an allocation region that began on the formerly "last" page:

   0x10027cfff0: 0x00000000000000d9     <-- this was Lisp's view of
   0x10027cfff8: 0x0000000000000006         the last page (page 1273)
   ---- page boundary ----
   0x10027d0000: 0x0000001000005ecf     <-- next_free_page moves here (page 1274)
   0x10027d0008: 0x00000000000000ba
   0x10027d0010: 0x0000000000000040
   0x10027d0018: 0x0000000000000000

Lisp did not think that the page starting at 0x10027d0000 was allocated,
so it believes the stopping point is page 1273.  When we read the bytes-used
on that page, we see a totally full page, but do not consider adjoining any
additional pages into the contiguous block.
However the object, a vector, that started on page 1273 ends on page 1274,
causing MAP-OBJECTS-IN-RANGE to assert that it overran 0x10027d0000.

We could try a few things to mitigate this:
* Try to "chase" the value of next-free-page.  This is literally impossible -
  it's a moving target, and it's extremely likely to exhaust memory doing so,
  especially if the supplied lambda is an interpreted function.
  (Each object scanned causes consing of more bytes, and we never
  "catch up" to the moving next-free-page)

* If the page that we're looking at is full but the FINALLY clause is hit,
  don't stop looking for more pages in that one case. Instead keep looking
  for the end of the contiguous block, but stop as soon as any potential
  stopping point is found; don't chase next-free-page.  This is tricky
  as well and just about as infeasible.

* Pass a flag to MAP-OBJECTS-IN-RANGE specifying that it's OK to
  surpass the expected bound - silently accept our fate.
  This is what we do since it's simple, and seems to work.
|#

;;; Iterate over all the objects allocated in each of the SPACES, calling FUN
;;; with the object, the object's type code, and the object's total size in
;;; bytes, including any header and padding. As a special case, if exactly one
;;; space named :ALL is requested, then map over the known spaces.
(defun map-allocated-objects (fun &rest spaces)
  (declare (type function fun)
           ;; KLUDGE: rest-arg and self calls do not play nice and it'll get consed
           (optimize (sb-c::recognize-self-calls 0)))
  (declare (dynamic-extent fun))
  (when (and (= (length spaces) 1) (eq (first spaces) :all))
    (return-from map-allocated-objects
     (map-allocated-objects fun
                            :read-only :static
                            #+immobile-space :immobile
                            :dynamic)))
  ;; You can't specify :ALL and also a list of spaces. Check that up front.
  (do-rest-arg ((space) spaces) (the spaces space))
  (flet ((do-1-space (space)
           (ecase space
            (:static
             ;; Static space starts with NIL, which requires special
             ;; handling, as the header and alignment are slightly off.
             (multiple-value-bind (start end) (%space-bounds space)
               (declare (ignore start))
               (funcall fun nil symbol-widetag (* sizeof-nil-in-words n-word-bytes))
               (let ((start (+ (logandc2 nil-value lowtag-mask)
                               (ash (- sizeof-nil-in-words 2) word-shift))))
                 (map-objects-in-range fun
                                       (ash start (- n-fixnum-tag-bits))
                                       end))))

            ((:read-only #-gencgc :dynamic)
             ;; Read-only space (and dynamic space on cheneygc) is a block
             ;; of contiguous allocations.
             (multiple-value-bind (start end) (%space-bounds space)
                                  (map-objects-in-range fun start end)))
            #+immobile-space
            (:immobile
             ;; Filter out filler objects. These either look like cons cells
             ;; in fixedobj subspace, or code without enough header words
             ;; in varyobj subspace. (cf 'filler_obj_p' in gc-internal.h)
             (dx-flet ((filter (obj type size)
                         (unless (= type list-pointer-lowtag)
                           (funcall fun obj type size))))
               (map-immobile-objects #'filter :fixed))
             (dx-flet ((filter (obj type size)
                         (unless (and (code-component-p obj)
                                      (code-obj-is-filler-p obj))
                           (funcall fun obj type size))))
               (map-immobile-objects #'filter :variable))))))
    (do-rest-arg ((space) spaces)
      (if (eq space :dynamic)
          (without-gcing #+cheneygc (do-1-space space)
                         #+gencgc (walk-dynamic-space fun #b1111111 0 0))
          (do-1-space space)))))

;;; Using the mask bits you can make many different match conditions resulting
;;; from a product of {boxed,unboxed,code,any} x {large,non-large,both}
;;; e.g. mask = #b10011" constraint = "#b10010"
;;; matches "large & (unboxed | code)"
;;;
;;; I think, when iterating over only code, that if we grab the code_allocator_lock
;;; and free_pages_lock, that this can be made reliable (both crash-free and
;;; guaranteed to visit all chosen objects) despite other threads running.
;;; As things are it is only "maybe" reliable, regardless of the parameters.
#+gencgc
(defun walk-dynamic-space (fun generation-mask
                               page-type-mask page-type-constraint)
  (declare (function fun)
           (type (unsigned-byte 7) generation-mask)
           (type (unsigned-byte 5) page-type-mask page-type-constraint))
       ;; Dynamic space on gencgc requires walking the GC page tables
       ;; in order to determine what regions contain objects.

       ;; We explicitly presume that any pages in an allocation region
       ;; that are in-use have a BYTES-USED of GENCGC-PAGE-BYTES
       ;; (indicating a full page) or an otherwise-valid BYTES-USED.
       ;; We also presume that the pages of an open allocation region
       ;; after the first page, and any pages that are unallocated,
       ;; have a BYTES-USED of zero.  GENCGC seems to guarantee this.

       ;; Our procedure is to scan forward through the page table,
       ;; maintaining an "end pointer" until we reach a page where
       ;; BYTES-USED is not GENCGC-PAGE-BYTES or we reach
       ;; NEXT-FREE-PAGE.  We then MAP-OBJECTS-IN-RANGE if the range
       ;; is not empty, and proceed to the next page (unless we've hit
       ;; NEXT-FREE-PAGE).

       ;; FIXME: WITHOUT-GCING prevents a GC flip, but doesn't prevent
       ;; closing allocation regions and opening new ones.  This may
       ;; prove to be an issue with concurrent systems, or with
       ;; spectacularly poor timing for closing an allocation region
       ;; in a single-threaded system.
  (close-thread-alloc-region)
  (do ((initial-next-free-page next-free-page)
       (base (int-sap (current-dynamic-space-start)))
       (start-page 0)
       (end-page 0)
       (end-page-bytes-used 0))
      ((> start-page initial-next-free-page))
         ;; The type constraint on page indices is probably too generous,
         ;; but it does its job of producing efficient code.
    (declare (type (integer 0 (#.(/ (ash 1 n-machine-word-bits) gencgc-page-bytes)))
                   start-page end-page))
    (setq end-page start-page)
    (loop (setq end-page-bytes-used
                (ash (slot (deref page-table end-page) 'words-used) word-shift))
          ;; See 'page_ends_contiguous_block_p' in gencgc.c
          (when (or (< end-page-bytes-used gencgc-page-bytes)
                    (= (slot (deref page-table (1+ end-page)) 'start) 0))
            (return))
          (incf end-page))
    (let ((start (sap+ base (truly-the signed-word
                                       (logand (* start-page gencgc-page-bytes)
                                               most-positive-word))))
          (end (sap+ base (truly-the signed-word
                                     (logand (+ (* end-page gencgc-page-bytes)
                                                end-page-bytes-used)
                                             most-positive-word)))))
      (when (sap> end start)
        ;; The bits in the 5-bit flag field have fixed positions,
        ;; but the position of the field itself depends on endianness.
        (let ((flags (ldb #+little-endian (byte 5 0) #+big-endian (byte 5 3)
                          (slot (deref page-table start-page) 'flags))))
          ;; The GEN slot is declared as (SIGNED 8) which does not satisfy the
          ;; type restriction on the first argument to LOGBITP.
          ;; Masking it to 3 bits fixes that, and allows using the other 5 bits
          ;; for something potentially.
          (when (and (logbitp (logand (slot (deref page-table start-page) 'gen) 7)
                              generation-mask)
                     (= (logand flags page-type-mask) page-type-constraint))
            (map-objects-in-range fun
                                  (%make-lisp-obj (sap-int start))
                                  (%make-lisp-obj (sap-int end))
                                  (< start-page initial-next-free-page))))))
    (setq start-page (1+ end-page))))

;; Users are often surprised to learn that a just-consed object can't
;; necessarily be seen by MAP-ALLOCATED-OBJECTS, so close the region
;; to update the page table.
;; Since we're in WITHOUT-GCING, there can be no interrupts.
;; Moreover it's probably not safe in the least to walk any thread's
;; allocation region, unless the observer and observed aren't consing.
(defun close-thread-alloc-region ()
  #+gencgc (alien-funcall (extern-alien "close_thread_region" (function void)))
  nil)

;;;; MEMORY-USAGE

#+immobile-space
(progn
(declaim (ftype (function (immobile-subspaces) (values t t t &optional))
                immobile-fragmentation-information))
(defun immobile-fragmentation-information (subspace)
  (binding* (((start free-pointer) (%space-bounds subspace))
             (used-bytes (ash (- free-pointer start) n-fixnum-tag-bits))
             (holes '())
             (hole-bytes 0))
    (map-immobile-objects
     (lambda (obj type size)
       (let ((address (logandc2 (get-lisp-obj-address obj) lowtag-mask)))
         (when (case subspace
                 (:fixed (= type list-pointer-lowtag))
                 (:variable (hole-p address)))
           (push (cons address size) holes)
           (incf hole-bytes size))))
     subspace)
    (values holes hole-bytes used-bytes)))

(defun show-fragmentation (&key (subspaces '(:fixed :variable))
                                (stream *standard-output*))
  (dolist (subspace subspaces)
    (format stream "~(~A~) subspace fragmentation:~%" subspace)
    (multiple-value-bind (holes hole-bytes total-space-used)
        (immobile-fragmentation-information subspace)
      (loop for (start . size) in holes
            do (format stream "~2@T~X..~X ~8:D~%" start (+ start size) size))
      (format stream "~2@T~18@<~:D hole~:P~> ~8:D (~,2,2F% of ~:D ~
                      bytes used)~%"
              (length holes) hole-bytes
              (/ hole-bytes total-space-used) total-space-used))))

(defun sb-kernel::immobile-space-usage ()
  (binding* (((nil fixed-hole-bytes fixed-used-bytes)
              (immobile-fragmentation-information :fixed))
             ((nil variable-hole-bytes variable-used-bytes)
              (immobile-fragmentation-information :variable))
             (total-used-bytes (+ fixed-used-bytes variable-used-bytes))
             (total-hole-bytes (+ fixed-hole-bytes variable-hole-bytes)))
    (values total-used-bytes total-hole-bytes)))
) ; end PROGN

;;; Return a list of 3-lists (bytes object type-name) for the objects
;;; allocated in Space.
(defun type-breakdown (space)
  (declare (muffle-conditions compiler-note))
  (let ((sizes (make-array 256 :initial-element 0 :element-type 'word))
        (counts (make-array 256 :initial-element 0 :element-type 'word)))
    (map-allocated-objects
     (lambda (obj type size)
       (declare (word size) (optimize (speed 3)) (ignore obj))
       (incf (aref sizes type) size)
       (incf (aref counts type)))
     space)

    (let ((totals (make-hash-table :test 'eq)))
      (dotimes (i 256)
        (let ((total-count (aref counts i)))
          (unless (zerop total-count)
            (let* ((total-size (aref sizes i))
                   (name (room-info-type-name (aref *room-info* i)))
                   (found (ensure-gethash name totals (list 0 0 name))))
              (incf (first found) total-size)
              (incf (second found) total-count)))))

      (collect ((totals-list))
        (maphash (lambda (k v)
                   (declare (ignore k))
                   (totals-list v))
                 totals)
        (sort (totals-list) #'> :key #'first)))))

;;; Handle the summary printing for MEMORY-USAGE. Totals is a list of lists
;;; (space-name . totals-for-space), where totals-for-space is the list
;;; returned by TYPE-BREAKDOWN.
(defun print-summary (spaces totals)
  (let ((summary (make-hash-table :test 'eq))
        (space-count (length spaces)))
    (dolist (space-total totals)
      (dolist (total (cdr space-total))
        (push (cons (car space-total) total)
              (gethash (third total) summary))))

    (collect ((summary-totals))
      (maphash (lambda (k v)
                 (declare (ignore k))
                 (let ((sum 0))
                   (declare (unsigned-byte sum))
                   (dolist (space-total v)
                     (incf sum (first (cdr space-total))))
                   (summary-totals (cons sum v))))
               summary)

      (format t "~2&Summary of space~P: ~(~{~A ~}~)~%" space-count spaces)
      (let ((summary-total-bytes 0)
            (summary-total-objects 0))
        (declare (unsigned-byte summary-total-bytes summary-total-objects))
        (dolist (space-totals
                 (mapcar #'cdr (sort (summary-totals) #'> :key #'car)))
          (let ((total-objects 0)
                (total-bytes 0)
                name)
            (declare (unsigned-byte total-objects total-bytes))
            (collect ((spaces))
              (dolist (space-total space-totals)
                (let ((total (cdr space-total)))
                  (setq name (third total))
                  (incf total-bytes (first total))
                  (incf total-objects (second total))
                  (spaces (cons (car space-total) (first total)))))
              (format t "~%~A:~%    ~:D bytes, ~:D object~:P"
                      name total-bytes total-objects)
              (unless (= 1 space-count)
                (dolist (space (spaces))
                  (format t ", ~D% ~(~A~)"
                          (round (* (cdr space) 100) total-bytes) (car space))))
              (format t ".~%")
              (incf summary-total-bytes total-bytes)
              (incf summary-total-objects total-objects))))
        (format t "~%Summary total:~%    ~:D bytes, ~:D objects.~%"
                summary-total-bytes summary-total-objects)))))

(declaim (ftype (sfunction (index &key (:comma-interval (and (integer 1) index))) index)
                decimal-with-grouped-digits-width))
(defun decimal-with-grouped-digits-width (value &key (comma-interval 3))
  (let ((digits (length (write-to-string value :base 10))))
    (+ digits (floor (1- digits) comma-interval))))

;;; Report object usage for a single space.
(defun report-space-total (space-info cutoff)
  (declare (list space-info) (type (or single-float null) cutoff))
  (destructuring-bind (space . types) space-info
    (format t "~2&Breakdown for ~(~A~) space:~%" space)
    (let* ((total-bytes (reduce #'+ (mapcar #'first types)))
           (bytes-width (decimal-with-grouped-digits-width total-bytes))
           (total-objects (reduce #'+ (mapcar #'second types)))
           (objects-width (decimal-with-grouped-digits-width total-objects))
           (cutoff-point (if cutoff
                             (truncate (* (float total-bytes) cutoff))
                             0))
           (reported-bytes 0)
           (reported-objects 0))
      (declare (unsigned-byte total-objects total-bytes cutoff-point
                              reported-objects reported-bytes))
      (flet ((type-usage (bytes objects name &optional note)
               (format t "  ~V:D bytes for ~V:D ~(~A~) object~2:*~P~*~
                          ~:[~; ~:*(~A)~]~%"
                       bytes-width bytes objects-width objects name note)))
        (loop for (bytes objects name) in types do
             (when (<= bytes cutoff-point)
               (type-usage (- total-bytes reported-bytes)
                           (- total-objects reported-objects)
                           "other")
               (return))
             (incf reported-bytes bytes)
             (incf reported-objects objects)
             (type-usage bytes objects name))
        (terpri)
        (type-usage total-bytes total-objects space "space total")))))

;;; Print information about the heap memory in use. PRINT-SPACES is a
;;; list of the spaces to print detailed information for.
;;; COUNT-SPACES is a list of the spaces to scan. For either one, T
;;; means all spaces (i.e. :STATIC, :DYNAMIC and :READ-ONLY.) If
;;; PRINT-SUMMARY is true, then summary information will be printed.
;;; The defaults print only summary information for dynamic space. If
;;; true, CUTOFF is a fraction of the usage in a report below which
;;; types will be combined as OTHER.
(defun memory-usage (&key print-spaces (count-spaces '(:dynamic #+immobile-space :immobile))
                          (print-summary t) cutoff)
  (declare (type (or single-float null) cutoff))
  (let* ((spaces (if (eq count-spaces t) +heap-space-keywords+ count-spaces))
         (totals (mapcar (lambda (space)
                           (cons space (type-breakdown space)))
                         spaces)))

    (dolist (space-total totals)
      (when (or (eq print-spaces t)
                (member (car space-total) print-spaces))
        (report-space-total space-total cutoff)))

    (when print-summary (print-summary spaces totals)))

  (values))

;;; Print a breakdown by instance type of all the instances allocated
;;; in SPACE. If TOP-N is true, print only information for the
;;; TOP-N types with largest usage.
(defun instance-usage (space &key (top-n 15))
  (declare (type spaces space) (type (or fixnum null) top-n))
  (format t "~2&~@[Top ~W ~]~(~A~) instance types:~%" top-n space)
  (let ((totals (make-hash-table :test 'eq))
        (total-objects 0)
        (total-bytes 0))
    (declare (unsigned-byte total-objects total-bytes))
    (map-allocated-objects
     (lambda (obj type size)
       (declare (optimize (speed 3)))
       (when (or (eql type instance-widetag)
                 (eql type funcallable-instance-widetag))
         (incf total-objects)
         (let* ((layout (if (eql type funcallable-instance-widetag)
                            (%fun-wrapper obj)
                            (%instance-wrapper obj)))
                (classoid (wrapper-classoid layout))
                (found (ensure-gethash classoid totals (cons 0 0)))
                (size size))
           (declare (fixnum size))
           (incf total-bytes size)
           (incf (the fixnum (car found)))
           (incf (the fixnum (cdr found)) size))))
     space)
    (let* ((sorted (sort (%hash-table-alist totals) #'> :key #'cddr))
           (interesting (if top-n
                            (subseq sorted 0 (min (length sorted) top-n))
                            sorted))
           (bytes-width (decimal-with-grouped-digits-width total-bytes))
           (objects-width (decimal-with-grouped-digits-width total-objects))
           (totals-label (format nil "~:(~A~) instance total" space))
           (types-width (reduce #'max interesting
                                :key (lambda (info)
                                       (let ((type (first info)))
                                         (length
                                          (typecase type
                                            (string
                                             type)
                                            (classoid
                                             (with-output-to-string (stream)
                                               (sb-ext:print-symbol-with-prefix
                                                stream (classoid-name type))))))))
                                :initial-value (length totals-label)))
           (printed-bytes 0)
           (printed-objects 0))
      (declare (unsigned-byte printed-bytes printed-objects))
      (flet ((type-usage (type objects bytes)
               (etypecase type
                 (string
                  (format t "  ~V@<~A~> ~V:D bytes, ~V:D object~:P.~%"
                          (1+ types-width) type bytes-width bytes
                          objects-width objects))
                 (classoid
                  (format t "  ~V@<~/sb-ext:print-symbol-with-prefix/~> ~
                             ~V:D bytes, ~V:D object~:P.~%"
                          (1+ types-width) (classoid-name type) bytes-width bytes
                          objects-width objects)))))
        (loop for (type . (objects . bytes)) in interesting
              do (incf printed-bytes bytes)
                 (incf printed-objects objects)
                 (type-usage type objects bytes))
        (terpri)
        (let ((residual-objects (- total-objects printed-objects))
              (residual-bytes (- total-bytes printed-bytes)))
          (unless (zerop residual-objects)
            (type-usage "Other types" residual-objects residual-bytes)))
        (type-usage totals-label total-objects total-bytes))))
  (values))

;;;; PRINT-ALLOCATED-OBJECTS

;;; This function is sheer madness.  You're better off using
;;; LIST-ALLOCATED-OBJECTS and then iterating over that, to avoid
;;; seeing all the junk created while doing this thing.
(defun print-allocated-objects (space &key (percent 0) (pages 5)
                                      type larger smaller count
                                      (stream *standard-output*))
  (declare (type (integer 0 99) percent) (type index pages)
           (type stream stream) (type spaces space)
           (type (or index null) type larger smaller count))
  (multiple-value-bind (start end) (%space-bounds space)
    (let* ((space-start (ash start n-fixnum-tag-bits))
           (space-end (ash end n-fixnum-tag-bits))
           (space-size (- space-end space-start))
           (pagesize sb-c:+backend-page-bytes+)
           (start (+ space-start (round (* space-size percent) 100)))
           (printed-conses (make-hash-table :test 'eq))
           (pages-so-far 0)
           (count-so-far 0)
           (last-page 0))
      (declare (type word last-page start)
               (fixnum pages-so-far count-so-far pagesize))
      (labels ((note-conses (x)
                 (unless (or (atom x) (gethash x printed-conses))
                   (setf (gethash x printed-conses) t)
                   (note-conses (car x))
                   (note-conses (cdr x)))))
        (map-allocated-objects
         (lambda (obj obj-type size)
           (let ((addr (get-lisp-obj-address obj)))
             (when (>= addr start)
               (when (if count
                         (> count-so-far count)
                         (> pages-so-far pages))
                 (return-from print-allocated-objects (values)))

               (unless count
                 (let ((this-page (* (the (values word t)
                                       (truncate addr pagesize))
                                     pagesize)))
                   (declare (type word this-page))
                   (when (/= this-page last-page)
                     (when (< pages-so-far pages)
                       ;; FIXME: What is this? (ERROR "Argh..")? or
                       ;; a warning? or code that can be removed
                       ;; once the system is stable? or what?
                       (format stream "~2&**** Page ~W, address ~X:~%"
                               pages-so-far addr))
                     (setq last-page this-page)
                     (incf pages-so-far))))

               (when (and (or (not type) (eql obj-type type))
                          (or (not smaller) (<= size smaller))
                          (or (not larger) (>= size larger)))
                 (incf count-so-far)
                 (case type
                   (#.code-header-widetag
                    (let ((dinfo (%code-debug-info obj)))
                      (format stream "~&Code object: ~S~%"
                              (if dinfo
                                  (sb-c::compiled-debug-info-name dinfo)
                                  "No debug info."))))
                   (#.symbol-widetag
                    (format stream "~&~S~%" obj))
                   (#.list-pointer-lowtag
                    (unless (gethash obj printed-conses)
                      (note-conses obj)
                      (let ((*print-circle* t)
                            (*print-level* 5)
                            (*print-length* 10))
                        (format stream "~&~S~%" obj))))
                   (t
                    (fresh-line stream)
                    (let ((str (write-to-string obj :level 5 :length 10
                                                :pretty nil)))
                      (unless (eql type instance-widetag)
                        (format stream "~S: " (type-of obj)))
                      (format stream "~A~%"
                              (subseq str 0 (min (length str) 60))))))))))
         space))))
  (values))

;;;; LIST-ALLOCATED-OBJECTS, LIST-REFERENCING-OBJECTS

(defun list-allocated-objects (space &key type larger smaller count
                                     test)
  (declare (type (or (eql :all) spaces) space)
           (type (or (unsigned-byte 8) null) type)
           (type (or index null) larger smaller count)
           (type (or function null) test))
  (declare (dynamic-extent test))
  (when (eql count 0)
    (return-from list-allocated-objects nil))
  ;; This function was pretty much random as to what subset of the heap it
  ;; visited- it might see half the heap, 1/10th of the heap, who knows, because
  ;; it stopped based on hitting a sentinel cons cell made just prior to the loop.
  ;; That stopping condition was totally wrong because allocation does not occur
  ;; linearly.  Taking 2 passes (first count, then store) stands a chance of
  ;; getting a reasonable point-in-time view as long as other threads are not consing
  ;; like crazy. If the user-supplied TEST function conses at all, then the result is
  ;; still very arbitrary - including possible duplication of objects if we visit
  ;; something and then see it again after GC transports it higher. The only way to
  ;; allow consing in the predicate would be to use dedicated "arenas" for new
  ;; allocations, that being a concept which we do not now - and may never - support.
  (sb-int:dx-flet ((wantp (obj widetag size)
                     (and (or (not type) (eql widetag type))
                          (or (not smaller) (<= size smaller))
                          (or (not larger) (>= size larger))
                          (or (not test) (funcall test obj)))))
    ;; Unless COUNT is smallish, always start by counting. Don't just trust the user
    ;; because s/he might specify :COUNT huge-num which is acceptable provided that
    ;; huge-num is an INDEX which could either exhaust the heap, or at least be
    ;; wasteful if but a tiny handful of objects would actually satisfy WANTP.
    (let* ((output (make-array
                    (if (typep count '(integer 0 100000))
                        count
                        (let ((n 0))
                          (map-allocated-objects
                           (lambda (obj widetag size)
                             (when (wantp obj widetag size) (incf n)))
                           space)
                          n))))
           (index 0))
      (block done
       (map-allocated-objects
        (lambda (obj widetag size)
          (when (wantp obj widetag size)
            (setf (aref output index) obj)
            (when (= (incf index) (length output))
              (return-from done))))
        space))
      (let ((list
             (cond ((= index (length output)) ; easy case
                    (coerce output 'list))
                   (t ; didn't fill the array
                    (collect ((res))
                      (dotimes (i index (res))
                        (res (svref output i))))))))
        (fill output 0) ; assist GC a bit
        list))))

;;; Calls FUNCTION with all objects that have (possibly conservative)
;;; references to them on current stack.
;;; This is for use by SB-INTROSPECT. (Other consumers, at your own risk)
;;; Note that we do not call MAKE-LISP-OBJ in the errorp=nil mode, as it
;;; potentially uses FORMAT and MAKE-UNPRINTABLE-OBJECT with each invocation.
;;; And see the cautionary remarks above that function regarding its dangerous
;;; nature (more so on precise GC).  On conservative GC we should be OK here
;;; because we know that there's a stack reference.
(defun map-stack-references (function)
  (declare (type function function))
  (declare (dynamic-extent function))
  (macrolet ((iter (step limit test)
               `(do ((sp (current-sp) (sap+ sp (,step n-word-bytes)))
                     (limit (sb-di::descriptor-sap ,limit))
                     (seen nil))
                    ((,test sp limit))
                  (let ((word (sap-ref-word sp 0)))
                    ;; Explicitly skip non-pointer words. The callable that
                    ;; SB-INTROSPECT provides ignores immediate values anyway.
                    (when (and (is-lisp-pointer word)
                               (not (zerop (sb-di::valid-lisp-pointer-p (int-sap word)))))
                      (let ((obj (%make-lisp-obj word)))
                        (unless (memq obj seen)
                          (push obj seen)
                          (funcall function obj))))))))
    #+stack-grows-downward-not-upward (iter + *control-stack-end* sap>)
    #-stack-grows-downward-not-upward (iter - *control-stack-start* sap<)))

(declaim (inline symbol-extra-slot-p))
(defun symbol-extra-slot-p (x)
  (> (logand (get-header-data x) tiny-boxed-size-mask)
     (1- symbol-size)))

;;; Invoke FUNCTOID (a macro or function) on OBJ and any values in MORE.
;;; Note that neither OBJ nor items in MORE undergo ONCE-ONLY treatment.
;;; The fact that FUNCTOID can be a macro allows treatment of its first argument
;;; as a generalized place in the manner of SETF, allowing read/write access.
;;; CLAUSES are used to modify the output of this macro. See example uses
;;; for more detail.
;;; HIGH EXPERIMENTAL: PROCEED AT YOUR OWN RISK.
(defmacro do-referenced-object ((obj functoid &rest more) &rest alterations
                                &aux (n-matched-alterations 0))
  (labels ((make-case (type &rest actions)
             (apply #'make-case* type
                    (mapcar (lambda (action) `(,functoid ,action ,@more))
                            actions)))
           (make-case* (type &rest actions)
             (let* ((found (assoc type alterations :test 'equal))
                    (alteration (or (cdr found) '(:extend))))
               (when found
                 (incf n-matched-alterations))
               (ecase (car alteration)
                 (:override (list `(,type ,@(cdr alteration))))
                 (:extend   (list `(,type ,@actions ,@(cdr alteration))))
                 (:delete))))) ; no clause
    (prog1
      `(typecase ,obj
         ;; Until the compiler can learn how to efficiently compile jump tables
         ;; by widetag, test in descending order of popularity.
         ;; These two are in fact generally the most frequently occurring type.
         ,.(make-case 'cons `(car ,obj) `(cdr ,obj))
         ,.(make-case* 'instance
            `(progn
               (,functoid (%instance-layout ,obj) ,@more)
               (do-instance-tagged-slot (.i. ,obj nil)
                 (,functoid (%instance-ref ,obj .i.) ,@more))))
         (function
          (typecase ,obj
            ,.(make-case* 'closure
               `(,functoid (%closure-fun ,obj) ,@more)
               `(do-closure-values (.o. ,obj :include-extra-values t)
                  ;; FIXME: doesn't allow setf, but of course there is
                  ;; no closure-index-set anyway, so .O. might be unused
                  ;; if functoid is a macro that does nothing.
                  (,functoid .o. ,@more)))
            ,.(make-case* 'funcallable-instance
               `(progn
                  ;; As for INSTANCE, allow the functoid to see the access form
                  (,functoid (%fun-layout ,obj) ,@more)
                  (,functoid (%funcallable-instance-fun ,obj) ,@more)
                  ;; Unfortunately for FUNCALLABLE-INSTANCEs, the relation
                  ;; between layout bitmap indices and indices as given to
                  ;; FUNCALLABLE-INSTANCE-INFO is not so obvious, and it's
                  ;; both tricky and unnecessary to generalize iteration.
                  ;; So just hardcode the few cases that exist.
                  #+compact-instance-header
                  (ecase (wrapper-bitmap (%fun-wrapper ,obj))
                    (-1 ; external trampoline, all slots are tagged
                     ;; In this case, the trampoline word is scanned, with no ill effect.
                     (loop for .i. from 0
                           to (- (get-closure-length ,obj) funcallable-instance-info-offset)
                           do (,functoid (%funcallable-instance-info ,obj .i.) ,@more)))
                    (#b0110 ; internal trampoline, 2 raw slots, 1 tagged slot
                     ;;   ^ ---- trampoline
                     ;;  ^------ implementation function
                     ;; ^------- (FUNCALLABLE-INSTANCE-INFO 0)
                     ;; and the rest of the words are raw.
                     (,functoid (%funcallable-instance-info ,obj 0) ,@more)))
                  #-compact-instance-header
                  (progn
                    (aver (eql (wrapper-bitmap (%fun-wrapper ,obj)) -4))
                    ;;            v ----trampoline
                    ;; = #b1...1100
                    ;;           ^----- layout
                    ;;          ^------ implementation function
                    ;;         ^------- (FUNCALLABLE-INSTANCE-INFO 0)
                    (loop for .i. from 0
                          to (- (get-closure-length ,obj) funcallable-instance-info-offset)
                          do (,functoid (%funcallable-instance-info ,obj .i.) ,@more)))))
            .,(make-case 'function))) ; in case there was code provided for it
         (t
          ;; TODO: the generated code is pretty horrible. OTHER-POINTER-LOWTAG
          ;; is known at this point, but tested N times.
          (typecase ,obj
            ,.(make-case* 'simple-vector
               `(dotimes (.i. (length ,obj))
                  (,functoid (data-vector-ref ,obj .i.) ,@more)))
            ;; Fancy arrays aren't highly popular, but this case must precede ARRAY
            ;; because ARRAY weeds out all other arrays, namely the ones that
            ;; hold no pointers: simple-string, simple-bit-vector, etc.
            ,.(make-case '(satisfies array-header-p)
               `(%array-data ,obj)
               `(%array-displaced-p ,obj)
               `(%array-displaced-from ,obj))
            ,.(make-case 'array)
            ,.(make-case* 'symbol
               `(,functoid (%primitive sb-c:fast-symbol-global-value ,obj) ,@more)
               `(,functoid (symbol-%info ,obj) ,@more)
               `(,functoid (symbol-name ,obj) ,@more)
               `(,functoid (symbol-package ,obj) ,@more)
               `(when (symbol-extra-slot-p ,obj)
                  (,functoid (symbol-extra ,obj) ,@more)))
            ,.(make-case 'fdefn
               `(fdefn-name ,obj)
               `(fdefn-fun ,obj)
               ;; While it looks like we could easily allow a pointer to a movable object
               ;; in the fdefn-raw-addr slot, it is not exactly trivial- at a bare minimum,
               ;; translating the raw-addr to a lispobj might have to be pseudoatomic,
               ;; since we don't know what object to pin when reconstructing it.
               ;; For simple-funs in dynamic space, it doesn't have to be pseudoatomic
               ;; because a reference to the interior of code pins the code.
               ;; Closure trampolines would be fine as well. That leaves funcallable instances
               ;; as the pain point. Those could go on pages of code as well, but see the
               ;; comment in conservative_root_p() in gencgc as to why that alone
               ;; would be inadequate- we require a properly tagged descriptor
               ;; to enliven any object other than code.
               #+immobile-code
               `(%make-lisp-obj
                 (alien-funcall (extern-alien "fdefn_callee_lispobj" (function unsigned unsigned))
                                (logandc2 (get-lisp-obj-address ,obj) lowtag-mask))))
            ,.(make-case* 'code-component
               `(loop for .i. from 2 below (code-header-words ,obj)
                      do (,functoid (code-header-ref ,obj .i.) ,@more)))
            ,.(make-case '(or float (complex float) bignum
                           #+sb-simd-pack simd-pack
                           #+sb-simd-pack-256 simd-pack-256
                           system-area-pointer)) ; nothing to do
            ,.(make-case 'weak-pointer `(weak-pointer-value ,obj))
            ,.(make-case 'ratio `(%numerator ,obj) `(%denominator ,obj))
            ;; Visitor won't be invoked on (COMPLEX float)
            ,.(make-case '(complex rational) `(%realpart ,obj) `(%imagpart ,obj))
            ;; Caller can do anything in the fallback case.
            ,.(make-case 't))))
      (when (> (length alterations) n-matched-alterations)
        (error "DO-REFERENCED-OBJECT usage error")))))

;;; Return T if and only if THIS references THAT.
;;; code-components are considered to reference their embedded
;;; simple-funs for this purpose; if THIS is a simple-fun, it is ignored.
(defun references-p (this that)
  (declare (optimize (sb-c::aref-trapping 0)))
  (macrolet ((test (x) `(when (eq ,x that) (go win))))
    (tagbody
       (do-referenced-object (this test)
         (code-component
          :extend
          (dotimes (i (code-n-entries this))
            (let ((f (%code-entry-point this i)))
              (when (eq f that)
                (go win)))))
         (t
          :extend
          (case (widetag-of this)
            (#.value-cell-widetag
             (test (value-cell-ref this)))
            (#.filler-widetag)
            (t
             (bug "Unknown object type #x~x addr=~x"
                  (widetag-of this)
                  (get-lisp-obj-address this))))))
       (return-from references-p nil)
     win
       (return-from references-p t))))

;;; This interface allows one either to be agnostic of the referencing space,
;;; or specify exactly one space, but not specify a list of spaces.
;;; An upward-compatible change would be to assume a list, and call ENSURE-LIST.
(defun map-referencing-objects (fun space object)
  (declare (type (or (eql :all) spaces) space))
  (declare (dynamic-extent fun))
  (let (list)
    (map-allocated-objects
     (lambda (referer widetag size)
       (declare (ignore widetag size))
       ;; Don't count a self-reference as a reference
       (when (and (neq referer object)
                  (references-p referer object))
         (push referer list)))
     space)
    (mapc (%coerce-callable-to-fun fun) list)))

(defun list-referencing-objects (space object)
  (collect ((res))
    (map-referencing-objects
     (lambda (obj) (res obj)) space object)
    (res)))

;;;; ROOM

(defun room-minimal-info ()
  (multiple-value-bind (names name-width
                        used-bytes used-bytes-width
                        overhead-bytes)
      (loop for (nil name function) in +all-spaces+
            for (space-used-bytes space-overhead-bytes)
               = (multiple-value-list (funcall function))
            collect name into names
            collect space-used-bytes into used-bytes
            collect space-overhead-bytes into overhead-bytes
            maximizing (length name) into name-maximum
            maximizing space-used-bytes into used-bytes-maximum
            finally (return (values
                             names name-maximum
                             used-bytes (decimal-with-grouped-digits-width
                                         used-bytes-maximum)
                             overhead-bytes)))
    (loop for name in names
          for space-used-bytes in used-bytes
          for space-overhead-bytes in overhead-bytes
          do (format t "~V@<~A usage is:~> ~V:D bytes~@[ (~:D bytes ~
                        overhead)~].~%"
                     (+ name-width 10) name used-bytes-width space-used-bytes
                     space-overhead-bytes)))
  #+sb-thread
  (format t "Control and binding stack usage is for the current thread ~
             only.~%")
  (format t "Garbage collection is currently ~:[enabled~;DISABLED~].~%"
          *gc-inhibit*))

(defun room-intermediate-info ()
  (room-minimal-info)
  (memory-usage :count-spaces '(:dynamic #+immobile-space :immobile)
                :print-spaces t
                :cutoff 0.05f0
                :print-summary nil))

(defun room-maximal-info ()
  (let ((spaces '(:dynamic #+immobile-space :immobile :static)))
    (room-minimal-info)
    (memory-usage :count-spaces spaces)
    (dolist (space spaces)
      (instance-usage space :top-n 10))))

(defun room (&optional (verbosity :default))
  "Print to *STANDARD-OUTPUT* information about the state of internal
  storage and its management. The optional argument controls the
  verbosity of output. If it is T, ROOM prints out a maximal amount of
  information. If it is NIL, ROOM prints out a minimal amount of
  information. If it is :DEFAULT or it is not supplied, ROOM prints out
  an intermediate amount of information."
  (fresh-line)
  (ecase verbosity
    ((t)
     (room-maximal-info))
    ((nil)
     (room-minimal-info))
    (:default
     (room-intermediate-info)))
  (values))

#+nil ; for debugging
(defun show-dynamic-space-code (&optional (stream *standard-output*)
                                &aux (n-code-bytes 0)
                                     (total-pages next-free-page)
                                     (pages
                                      (make-array total-pages :element-type 'bit
                                                  :initial-element 0)))
  (flet ((dump-page (page-num)
           (format stream "~&Page ~D~%" page-num)
           (let ((where (+ dynamic-space-start (* page-num gencgc-page-bytes)))
                 (seen-filler nil))
             (loop
               (let* ((obj (let ((sap (int-sap where)))
                             (lispobj@baseptr sap (widetag@baseptr sap))))
                      (size (primitive-object-size obj)))
                 (when (code-component-p obj)
                   (incf n-code-bytes size))
                 (when (if (and (consp obj) (eq (car obj) 0) (eq (cdr obj) 0))
                           (if seen-filler
                               (progn (write-char #\. stream) nil)
                               (setq seen-filler t))
                           (progn (setq seen-filler nil) t))
                   (let ((*print-pretty* nil))
                     (format stream "~&  ~X ~4X ~S " where size obj)))
                 (incf where size)
                 (loop for index from page-num to (find-page-index (1- where))
                       do (setf (sbit pages index) 1)))
               (let ((next-page (find-page-index where)))
                 (cond ((= (logand where (1- gencgc-page-bytes)) 0)
                        (format stream "~&-- END OF PAGE --~%")
                        (return next-page))
                       ((eq next-page page-num))
                       (t
                        (setq page-num next-page seen-filler nil))))))))
    (let ((i 0))
      (loop while (< i total-pages)
            do (let ((type (ldb (byte 2 0) (slot (deref page-table i) 'flags))))
                 (if (= type 3)
                     (setq i (dump-page i))
                     (incf i)))))
    (let* ((n-pages (count 1 pages))
           (tot (* n-pages gencgc-page-bytes))
           (waste (- tot n-code-bytes)))
      (format t "~&Used-bytes=~D Pages=~D Waste=~D (~F%)~%"
              n-code-bytes n-pages waste
              (* 100 (/ waste tot))))))

#+nil ; for debugging
(defun show-immobile-spaces (which)
  (flet ((show (obj type size)
           (declare (ignore type size))
           (let ((*print-pretty* nil) (*print-length* 3))
             (format t "~x: ~s~%" (get-lisp-obj-address obj) obj))))
    (when (or (eq which :fixed) (eq which :both))
      (format t "Fixedobj space~%==============~%")
      (map-objects-in-range #'show
        (%make-lisp-obj fixedobj-space-start)
        (%make-lisp-obj (sap-int *fixedobj-space-free-pointer*))))
    (when (or (eq which :variable) (eq which :both))
      (format t "Varyobj space~%=============~%")
      (map-objects-in-range #'show
        (%make-lisp-obj varyobj-space-start)
        (%make-lisp-obj (sap-int *varyobj-space-free-pointer*))))))

#+gencgc
(defun generation-of (object)
  (with-pinned-objects (object)
    (let* ((addr (get-lisp-obj-address object))
           (page (find-page-index addr)))
      (cond ((>= page 0) (slot (deref page-table page) 'gen))
            #+immobile-space
            ((immobile-space-addr-p addr)
             ;; SIMPLE-FUNs don't contain a generation byte
             (when (simple-fun-p object)
               (setq addr (get-lisp-obj-address (fun-code-header object))))
             (let ((sap (int-sap (logandc2 addr lowtag-mask))))
               (logand (if (fdefn-p object) (sap-ref-8 sap 1) (sap-ref-8 sap 3))
                       #xF)))))))

;;; Show objects in a much simpler way than print-allocated-objects.
;;; Probably don't use this for generation 0 of dynamic space. Other spaces are ok.
;;; (And this is removed from the image; it's meant for developers)
#+gencgc
(defun show-generation-objs (gen space)
  (let ((*print-pretty* nil))
    (map-allocated-objects
     (lambda (obj type size)
       (declare (ignore type))
       (when (= (generation-of obj) gen)
         (format t "~x ~3x ~s~%" (get-lisp-obj-address obj) size obj)))
     space)))

;;; Unfortunately this is a near total copy of the test in gc.impure.lisp
(defun !ensure-genesis-code/data-separation ()
  #+gencgc
  (let* ((n-bits
          (progn
            (close-thread-alloc-region)
            (+ next-free-page 50)))
         (code-bits (make-array n-bits :element-type 'bit :initial-element 0))
         (data-bits (make-array n-bits :element-type 'bit :initial-element 0))
         (total-code-size 0))
    (map-allocated-objects
     (lambda (obj type size)
       (declare ((and fixnum (integer 1)) size))
       ;; M-A-O disables GC, therefore GET-LISP-OBJ-ADDRESS is safe
       (let ((obj-addr (get-lisp-obj-address obj))
             (array (cond ((= type code-header-widetag)
                           (incf total-code-size size)
                           code-bits)
                          (t
                           data-bits)))
             (other-array (cond ((= type code-header-widetag)
                                 data-bits)
                                (t
                                 code-bits))))
         ;; This is not the most efficient way to update the bit arrays,
         ;; but the simplest and clearest for sure. (The loop could avoided
         ;; if the current page is the same as the previously seen page)
         (loop for index from (find-page-index obj-addr)
               to (find-page-index (truly-the word
                                              (+ (logandc2 obj-addr lowtag-mask)
                                                 (1- size))))
               do (cond ((= (sbit other-array index) 1)
                         (format t "~&broken on page index ~d base ~x~%"
                                 index
                                 (+ dynamic-space-start (* index gencgc-page-bytes)))
                         (alien-funcall (extern-alien "ldb_monitor" (function void))))
                        (t
                         (setf (sbit array index) 1))))))
     :dynamic)))
;;; Because pseudo-static objects can not move nor be freed,
;;; this is a valid test that genesis separated code and data.
(!ensure-genesis-code/data-separation)

(defun hexdump (thing &optional (n-words nil wordsp)
                                ;; pass NIL explicitly if T crashes on you
                                (decode t))
  (multiple-value-bind (obj addr count)
      (if (typep thing 'word) ; ambiguous in the edge case, but assume it's
          ;; an address (though you might be trying to dump a bignum's data)
          (values nil thing (if wordsp n-words 1))
          (values
           thing
           (logandc2 (get-lisp-obj-address thing) lowtag-mask)
           (if wordsp
               n-words
               (if (and (typep thing 'code-component) (plusp (code-n-entries thing)))
                   ;; Display up through the first fun header
                   (+ (code-header-words thing)
                      (ash (%code-fun-offset thing 0) (- word-shift))
                      simple-fun-insts-offset)
                   ;; at most 16 words
                   (min 16 (ash (primitive-object-size thing) (- word-shift)))))))
    (with-pinned-objects (obj)
      (dotimes (i count)
        (let ((word (sap-ref-word (int-sap addr) (ash i word-shift))))
          (multiple-value-bind (lispobj ok fmt)
              (cond ((and (typep thing 'code-component)
                          (< 1 i (code-header-words thing)))
                     (values (code-header-ref thing i) t))
                    #+compact-symbol
                    ((and (typep thing '(and symbol (not null)))
                          (= i symbol-name-slot))
                     (values (list (sb-impl::symbol-package-id thing)
                                   (symbol-name thing))
                             t
                             "{~{~A,~S~}}"))
                    (decode
                     (make-lisp-obj word nil)))
            (let ((*print-lines* 1)
                  (*print-pretty* t))
              (format t "~x: ~v,'0x~:[~; = ~@?~]~%"
                      (+ addr (ash i word-shift))
                      (* 2 n-word-bytes)
                      word ok (or fmt "~S") lispobj))))))))
#+sb-thread
(defun show-tls-map ()
  (let ((list
         (sort (list-allocated-objects
                :all
                :type symbol-widetag
                :test (lambda (x) (plusp (sb-kernel:symbol-tls-index x))))
               #'<
               :key #'sb-kernel:symbol-tls-index))
        (prev 0))
    (dolist (x list)
      (let ((n  (ash (sb-kernel:symbol-tls-index x) (- word-shift))))
        (when (and (> n primitive-thread-object-length)
                   (> n (1+ prev)))
          (format t "(unused)~%"))
        (format t "~5d = ~s~%" n x)
        (setq prev n)))))

#+gencgc
(flet ((print-it (obj type size)
         (declare (ignore type size))
         (let ((*print-level* 2) (*print-lines* 1))
           (format t "~x ~s~%" (get-lisp-obj-address obj) obj))))
(defun print-all-code ()
  (walk-dynamic-space #'print-it #x7f #b01111 #b00111))
(defun print-large-code ()
  (walk-dynamic-space #'print-it #x7f #b11111 #b10111))
(defun print-large-unboxed ()
  (walk-dynamic-space #'print-it #x7f #b11111 #b10010))
;;; Use this only if you know what you're doing. It can fail because a page
;;; that needs to continue onto the next page will cause the "overrun" check
;;; to fail.
(defun print-page-contents (page)
  (let* ((start
          (+ (current-dynamic-space-start) (* gencgc-page-bytes page)))
         (end
          (+ start gencgc-page-bytes)))
    (map-objects-in-range #'print-it (%make-lisp-obj start) (%make-lisp-obj end)))))

(defun map-code-objects (fun)
  (dx-flet ((filter (obj type size)
              (declare (ignore size))
              (when (= type code-header-widetag)
                (funcall fun obj)))
            (nofilter (obj type size)
              (declare (ignore type size))
              (funcall fun obj)))
    #+cheneygc (map-allocated-objects #'filter :all)
    #+gencgc
    (without-gcing
      #+immobile-code
      (map-objects-in-range #'nofilter
                            (ash varyobj-space-start (- n-fixnum-tag-bits))
                            (%make-lisp-obj (sap-int *varyobj-space-free-pointer*)))
      (alien-funcall (extern-alien "close_code_region" (function void)))
      (walk-dynamic-space #'nofilter
                          #b1111111 ; all generations
                          #b111 #b111)))) ; type mask and constraint

(export 'code-from-serialno)
(defun code-from-serialno (serial)
  (dx-flet ((visit (obj)
              (when (= (%code-serialno obj) serial)
                (return-from code-from-serialno obj))))
    (map-code-objects #'visit)))

(defun show-all-layouts ()
  (let ((l (list-allocated-objects :all :test #'sb-kernel::wrapper-p))
        zero trailing-raw trailing-tagged vanilla)
    (dolist (x l)
      (let ((m (wrapper-bitmap x)))
        (cond ((eql m +layout-all-tagged+) (push x vanilla))
              ((eql m 0) (push x zero))
              ((minusp m) (push x trailing-tagged))
              (t (push x trailing-raw)))))
    (flet ((legend (newline str list)
             (when newline (terpri))
             (let ((s (format nil str (length list))))
               (format t "~A~%~A~%" s (make-string (length s) :initial-element #\-)))))
      (when zero
        (legend nil "Zero bitmap (~d):" zero)
        (dolist (x zero) (format t "~a~%" (wrapper-classoid-name x))))
      (when trailing-raw
        (legend t "Trailing raw (~d):" trailing-raw)
        (dolist (x trailing-raw)
          (let ((m (wrapper-bitmap x)))
            (format t "~30a 0...~v,'0b~%"
                    (wrapper-classoid-name x)
                    (acond ((wrapper-info x) (1+ (dd-length it))) (t 32))
                    m))))
      (when trailing-tagged
        (legend t "Trailing tagged (~d):" trailing-tagged)
        (dolist (x trailing-tagged)
          (let ((m (wrapper-bitmap x)))
            (format t "~30a 1...~b~%"
                    (wrapper-classoid-name x)
                    (acond ((wrapper-info x) (ldb (byte (dd-length it) 0) m))
                           (t (ldb (byte 32 0) m)))))))
      (legend t "Default: (~d) [not shown]" vanilla))))

#+ubsan
(defun find-poisoned-vectors (&aux result)
  (dolist (v (list-allocated-objects :all :type simple-vector-widetag)
             result)
    (when (dotimes (i (length v))
            (declare (optimize (sb-c::aref-trapping 0)))
            (let ((val (svref v i)))
              (when (= (get-lisp-obj-address val) no-tls-value-marker-widetag)
                (return t))))
      (push (make-weak-pointer v) result)
      (let* ((origin (vector-extra-data v))
             (code (sb-di::code-header-from-pc (ash origin -3)))
             (*print-array* nil))
        (format t "g~d ~a ~a~%" (sb-kernel:generation-of v) v code)))))

(in-package "SB-C")
;;; As soon as practical in warm build it makes sense to add
;;; cold-allocation-patch-points into the weak hash-table.
;;; FIXME: I suspect that this wants to be just a weak vector
;;; (all code objects that have any allocation profiling compiled in),
;;; and not a hash-table, and that the list of fixups in the component
;;; can be attached to the debug info (in the manner of debug funs).
;;; When this was first implemented, weak-vectors weren't a thing. Maybe?
(defvar *!cold-allocation-patch-point*)
(loop for (code . points) in *!cold-allocation-patch-point*
      do (setf (gethash code *allocation-patch-points*) points))
