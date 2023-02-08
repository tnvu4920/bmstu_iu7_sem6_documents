#+gencgc (setf (extern-alien "gc_allocate_dirty" char) 1)

(setq *compile-print* nil)
sb-ext::(declaim (muffle-conditions compiler-note))
(require :asdf)

(in-package :asdf)

(defun keywordize (x)
  (intern (string-upcase x) :keyword))

sb-ext::(declaim (muffle-conditions sb-kernel:redefinition-warning))
(defun wrapping-source-registry ()
  '(:source-registry (:tree #p"SYS:CONTRIB;") :ignore-inherited-configuration))
sb-ext::(declaim (unmuffle-conditions sb-kernel:redefinition-warning))

(defun setup-asdf-contrib ()
  ;;(setf *resolve-symlinks* nil)
  (let* ((sbcl-top (merge-pathnames (getenv-pathname "SBCL_TOP" :ensure-directory t)))
         (src-contrib (subpathname sbcl-top "contrib/"))
         (asdf-cache (subpathname sbcl-top "obj/asdf-cache/"))
         (source-registry '(:source-registry :ignore-inherited-configuration))
         (output-translations `(:output-translations (,(namestring src-contrib)
                                                      ,(namestring asdf-cache))
                                :ignore-inherited-configuration))
         (src.pat (wilden src-contrib))
         (src.dir.pat (merge-pathnames* *wild-inferiors* src-contrib))
         (out.pat (wilden asdf-cache)))
    (ensure-directories-exist asdf-cache)
    (setf (logical-pathname-translations "SYS")
          `(("CONTRIB;**;*.*.*" ,src.pat))) ;; this makes recursive tree search work.
    (initialize-source-registry source-registry)
    (initialize-output-translations output-translations)
    (setf (logical-pathname-translations "SYS")
          (labels ((typepat (type base)
                     `(,(format nil "CONTRIB;**;*.~:@(~A~).*" type)
                       ,(make-pathname :type (string-downcase type) :defaults base)))
                   (outpat (type) (typepat type out.pat))
                   (srcpat (type) (typepat type src.pat))
                   (outpats (&rest types) (mapcar #'outpat types))
                   (srcpats (&rest types) (mapcar #'srcpat types)))
            `(,@(srcpats :lisp :asd)
              ,@(outpats :fasl :sbcl-warnings :build-report
                               :out :exe :lisp-temp :o :c :test-report :html)
              ("CONTRIB;**;" ,src.dir.pat)
              #|("CONTRIB;**;*.*.*" ,src.pat)|#)))
    (setf *central-registry* nil)))

(defun build-asdf-contrib (system)
  (setq *features*
        (append '(:sb-building-contrib) sb-impl:+internal-features+ *features*))
  (setup-asdf-contrib)
  (let* ((name (string-downcase system))
         (sbcl-top (merge-pathnames (getenv-pathname "SBCL_TOP" :ensure-directory t)))
         (out-contrib (subpathname sbcl-top "obj/sbcl-home/contrib/"))
         (cache-module (subpathname sbcl-top (format nil "obj/asdf-cache/~a/" name)))
         (system (find-system name))
         (system.fasl (output-file 'compile-bundle-op system))
         (module.fasl (subpathname out-contrib (strcat name ".fasl")))
         (module-setup.lisp (subpathname cache-module "module-setup.lisp"))
           (module-setup.fasl (subpathname cache-module "module-setup.fasl"))
         (dependencies (mapcar 'keywordize (component-sideway-dependencies system)))
           (input-fasls (list module-setup.fasl system.fasl)))
    (ensure-directories-exist out-contrib)
    (ensure-directories-exist cache-module)
    (with-open-file (o module-setup.lisp
                       :direction :output :if-exists :rename-and-delete)
      (format o "(provide :~A)~%~{(require ~(~S~))~%~}" name dependencies))
    (compile-file module-setup.lisp :output-file module-setup.fasl)
    (operate 'compile-bundle-op system)
    (let ((s (find-symbol "DUMP/RESTORE-INTERESTING-TYPES" "SB-C")))
      (when s (funcall s 'write)))
    (concatenate-files input-fasls module.fasl)))

(defun test-asdf-contrib (system)
  (setq *features*
        (append '(:sb-testing-contrib) sb-impl:+internal-features+ *features*))
  (setup-asdf-contrib)
  (asdf:test-system system))
