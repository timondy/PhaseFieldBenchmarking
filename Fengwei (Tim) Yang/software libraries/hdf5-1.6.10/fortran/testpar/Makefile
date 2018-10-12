## HDF5 Fortran Parallel Library Test Makefile(.in)
##

##
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.
##
top_srcdir=../..
top_builddir=../..
srcdir=.
##------------------------------------------------------------ -*- makefile -*-
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.


## The following section of this makefile comes from the
## `./config/commence' file which was generated with config.status
## from `./config/commence.in'.
##-----------------------------------------------------------------------------

## Things that Make needs
.SUFFIXES:
.SUFFIXES: .c .f90 .o .lo


## Directories to search
VPATH=$(srcdir) $(top_builddir)/src $(top_srcdir)/src

## Programs
SHELL=/bin/sh
CC=mpicc
F9X=mpif90
F9XMODEXT=mod
F9XMODFLAG=-I
FSEARCH_DIRS=
## H5_CPPFLAGS, H5_CFLAGS, and H5_FFLAGS are used for extra flags used
## to build HDF5 (flags that don't need to be included in h5cc,
## for instance).
CFLAGS= -fPIC -ansi -O -fomit-frame-pointer -finline-functions 
CPPFLAGS=-D_BSD_SOURCE -D_POSIX_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE   
FFLAGS=-I.   -I. 
LDFLAGS=
ROOT=/home/f/fy/fy44/Programs/hdf5-1.6.10/fortran
LIBS=-lm 
AR=ar
RANLIB=ranlib
PERL=perl
RM=rm -f
CP=cp
INSTALL=/usr/bin/install -c
INSTALL_PROGRAM=${INSTALL}
INSTALL_DATA=${INSTALL} -m 644
PARALLEL=yes
RUNSERIAL=
RUNPARALLEL=/usr/bin/mpirun -np $${NPROCS:=2}
RUNTEST=$(RUNSERIAL)
TIME=time
TRACE=:

## Installation points
prefix=/home/f/fy/fy44/Programs/hdf5-1.6.10
exec_prefix=${prefix}
bindir=$(DESTDIR)${exec_prefix}/bin
libdir=$(DESTDIR)${exec_prefix}/lib
includedir=$(DESTDIR)${prefix}/include
docdir=$(DESTDIR)${prefix}/doc
PUB_LIB=$(LIB)

## Shared libraries
LT_STATIC_EXEC=
DYNAMIC_DIRS=
LT=$(top_builddir)/libtool
LT_CCOMPILE=$(LT) --mode=compile $(CC)
LT_LINK_CLIB=$(LT) --mode=link $(CC) $(DYNAMIC_DIRS) $(LT_STATIC_EXEC) -rpath $(libdir)
LT_LINK_CEXE=$(LT) --mode=link $(CC) $(DYNAMIC_DIRS) -rpath $(bindir)
LT_RUN=$(LT) --mode=execute
LT_INSTALL_PROG=$(LT) --mode=install $(INSTALL_PROGRAM)
LT_INSTALL_LIB=$(LT) --mode=install $(INSTALL_DATA)
LT_UNINSTALL=$(LT) --mode=uninstall $(RM)

# Fortran compilation/linking stuff
LT_FCOMPILE=$(LT) --mode=compile $(F9X)
LT_LINK_FLIB=$(LT) --mode=link $(F9X) -static $(DYNAMIC_DIRS)
LT_LINK_FEXE=$(LT) --mode=link $(F9X) $(LT_STATIC_EXEC) $(DYNAMIC_DIRS)
FORTMODFILES=h5a h5d h5e h5f h5g h5i h5p h5r h5s h5t h5z h5fortran_flags h5fortran_types h5global h5lib hdf5

## Optional variables. We must declare them here because Irix pmake
## complains if it sees a reference to a variable which has never been
## defined. The main makefile is free to redefine these to something else.
DOCDIR=$(docdir)
EXAMPLEDIR=$(docdir)/hdf5/examples/fortran
LIB=
LIB_SRC=
LIB_OBJ=
PUB_HDR=
PUB_PROGS=
PROGS=
TEST_PROGS=
TEST_PROGS_PARA=
TEST_FLAGS=
TEST_SCRIPTS=
TEST_SCRIPTS_PARA=
AUX_LIB=
EXAMPLE_PROGS=
SUBDIRS=
LIBHDF5=

## The default is to build the library and/or programs.  We must build
## them sequentially.
all:
	$(MAKE) lib
	$(MAKE) progs
	$(MAKE) tests

## The following rules insure that the Makefile is up-to-date by rerunning
## various autoconf components (although not all versions of make assume
## that the makefile is implicitly a target).  We use time stamp files to
## keep track of the most recent update of H5config.h.in and H5config.h
## because autoheader and config.status don't update the modification time
## if the contents don't change.
## Invoke it by "gmake reconfigure".
##
## Graphically, the dependencies are:
##
##		         configure.in
##                           |  |
##         +-----------------+  +------------+
##         |                                 |
##      stamp1                           configure
##   (H5config.h.in)                         |
##       |  |                                |
##       |  +---------------+    +-----------+
##       |                  |    |
##       |               config.status
##       |                     |
##       |  +------------------+
##       |  |
##      stamp2
##    (H5config.h)                   Makefile.in et al
##         |                                 |
##         +------------------+  +-----------+
##                            |  |
##                          Makefile
##
## A side effect of updating stamp1 is to generate H5config.h.in and a
## side effect of updating stamp2 is to generate H5config.h.  When using
## a version of make that doesn't treat the makefile as the initial target
## the user may want to occassionally type `make Makefile' in any source
## directory.
## The `Makefile' target has been renamed to `reconfigure' so that the
## autoconf and make depend processes do not start up automatically.
## One must do `make reconfigure' explicitedly to start the process.
## (When srcdir is used and if more than one machines are running,
## this automatic Makefile/autoconf can get things unstable.)
##
STAMP1=$(top_builddir)/config/stamp1
STAMP2=$(top_builddir)/config/stamp2

MAKEFILE_PARTS=$(srcdir)/Makefile.in		\
	       $(top_srcdir)/config/commence.in	\
	       $(top_srcdir)/config/conclude.in	\
	       $(top_srcdir)/config/depend1.in	\
	       $(top_srcdir)/config/depend2.in	\
	       $(top_srcdir)/config/depend3.in	\
	       $(top_srcdir)/config/depend4.in	\
	       $(top_srcdir)/config/dependN.in

$(STAMP1): $(top_srcdir)/configure.in
	touch $(STAMP1)
	-cd $(top_srcdir); autoheader

$(STAMP2): $(STAMP1) $(top_builddir)/config.status
	touch $(STAMP2)
	-cd $(top_builddir); \
	   CONFIG_FILES= CONFIG_HEADERS=src/H5config.h ./config.status

$(top_srcdir)/configure: $(top_srcdir)/configure.in
	-cd $(top_srcdir); autoconf

$(top_builddir)/config.status: $(top_srcdir)/configure $(STAMP1)
	-cd $(top_builddir); ./config.status --recheck

# rerun the autoconf process if any configure components have changed. 
reconfigure: $(MAKEFILE_PARTS) $(STAMP2)
	-cd $(top_builddir); CONFIG_HEADERS= ./config.status

##-----------------------------------------------------------------------------
## The following section of this makefile comes from the middle of
## `Makefile.in' from this directory.  It was generated by running
## `config.status'.
##-----------------------------------------------------------------------------

hdf5_builddir=$(top_builddir)/src

## Add the include directory to the C preprocessor flags the the hdf5 library
## to the library list.
LT_LINK_LIB=$(LT) --mode=link $(F9X) -static -rpath $(libdir)
LIBFORTRAN=../src/libhdf5_fortran.la
LIBHDF5=$(hdf5_builddir)/libhdf5.la
FTLIB = ../test/libh5test_fortran.la
CTLIB = ../../test/.libs/libh5test.a
FSEARCH_DIRS =-I../src 

## These are our main targets
TEST_PROGS_PARA=ptesthdf5_fortran

## Temporary files
MOSTLYCLEAN=ParaEg[123].h5f
DISTCLEAN=go

## Test source files
TEST_PFORTRAN_SRC=ptesthdf5_fortran.f90 thdf5.f90 thyperslab_wr.f90
TEST_PFORTRAN_OBJ=$(TEST_PFORTRAN_SRC:.f90=.lo)
TEST_SRC=$(TEST_PFORTRAN_SRC)
TEST_OBJ=$(TEST_SRC:.f90=.lo)
TEST_HDR=

## How to build the tests... They all depend on the hdf5 library
$(TEST_PROGS_PARA): $(LIBHDF5)

$(TEST_OBJ): $(TEST_HDR)

## Put in dependencies to the THDF5 module so that things will be built
## in the correct order.
thyperslab_wr.lo:	    $(srcdir)/thyperslab_wr.f90 thdf5.lo
ptesthdf5_fortran.lo:	    $(srcdir)/ptesthdf5_fortran.f90 thdf5.lo

ptesthdf5_fortran: thdf5.lo thyperslab_wr.lo ptesthdf5_fortran.lo
	@$(LT_LINK_FEXE) $(FFLAGS) -o $@ ptesthdf5_fortran.lo thdf5.lo thyperslab_wr.lo  $(FTLIB) $(CTLIB) $(LIBFORTRAN) $(LIBHDF5) $(LIBS)

##------------------------------------------------------------ -*- makefile -*-
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.


## The following section of this makefile comes from the
## `./config/conclude' file which was generated with config.status
## from `./config/conclude.in'.
##-----------------------------------------------------------------------------

## This is the target for the library described in the main body of the
## makefile.
##
lib: $(LIB)
$(LIB) __no_library__: $(LIB_OBJ)
	@$(LT_LINK_FLIB) -o $@ $(FFLAGS) $(LIB_OBJ) $(LDFLAGS) $(LIBS)

progs: $(LIB) $(PROGS)

## Build a tags file in this directory.
TAGS: $(LIB_SRC)
	$(RM) $@
	-etags $(LIB_SRC)

## Runs each test in order, passing $(TEST_FLAGS) to the program.
tests: $(TEST_PROGS) $(TEST_PROGS_PARA) $(LIB)
check test _test: tests
	@echo "===Fortran tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-s
	$(MAKE) _check-p
	@echo "===Fortran tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

## serial tests
check-s: tests
	@echo "===Serial Fortran tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-s
	@echo "===Serial Fortran tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

## parallel tests
check-p: tests
	@echo "===Parallel Fortran tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-p
	@echo "===Parallel Fortran tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

## setup for non-repeating tests.
TEST_PROG_CHKEXE=$(TEST_PROGS:=.chkexe_)
TEST_PROG_PARA_CHKEXE=$(TEST_PROGS_PARA:=.chkexe_)
TEST_SCRIPT_CHKSH=$(TEST_SCRIPTS:=.chkexe_)
TEST_SCRIPT_PARA_CHKSH=$(TEST_SCRIPTS_PARA:=.chkexe_)
# .chkexe files are used to mark tests that have run successfully.
# .chklog files are output from those tests.
# *.clog are from the MPE option.
CHECK_CLEANFILES=*.chkexe *.chklog *.clog

_check-s: $(TEST_PROG_CHKEXE) $(TEST_SCRIPT_CHKSH)

# The dummy.chkexe here prevents the target from being
# empty if there are no tests in the current directory.
# $${log} is the log file.
# $${tname} is the name of test.
$(TEST_PROG_CHKEXE) $(TEST_PROG_PARA_CHKEXE) dummy.chkexe_:
	@if test "X$@" != "X.chkexe_" && test "X$@" != "Xdummy.chkexe_"; then \
	   tname=$(@:.chkexe_=)$(EXEEXT);\
	   log=$(@:.chkexe_=.chklog); \
	   echo "============================"; \
	   if $(top_srcdir)/bin/newer $(@:.chkexe_=.chkexe) $${tname}; then \
	      echo "No need to test $${tname} again."; \
	   else \
	      echo "============================" > $${log}; \
	      if test "X$(HDF_FORTRAN)" = "Xyes"; then \
	         echo "Fortran API: Testing $(HDF5_DRIVER) $${tname} $(TEST_FLAGS)"; \
	         echo "Fortran API: $(HDF5_DRIVER) $${tname} $(TEST_FLAGS) Test Log" >> $${log}; \
	      elif test "X$(HDF_CXX)" = "Xyes"; then \
	         echo "C++ API: Testing $(HDF5_DRIVER) $${tname} $(TEST_FLAGS)"; \
	         echo "C++ API: $(HDF5_DRIVER) $${tname} $(TEST_FLAGS) Test Log" >> $${log};\
	      else \
	         echo "Testing $(HDF5_DRIVER) $${tname} $(TEST_FLAGS)"; \
	         echo "$(HDF5_DRIVER) $${tname} $(TEST_FLAGS) Test Log" >> $${log}; \
	      fi; \
	      echo "============================" >> $${log}; \
	      srcdir="$(srcdir)" \
	         $(TIME) $(RUNTEST) ./$${tname} $(TEST_FLAGS) >> $${log} 2>&1 \
	         && touch $(@:.chkexe_=.chkexe) || \
	         (test $$HDF5_Make_Ignore && echo "*** Error ignored") || \
	         (cat $${log} && false) || exit 1; \
	      echo "" >> $${log}; \
	      echo "Finished testing $${tname} $(TEST_FLAGS)" >> $${log}; \
	      echo "============================" >> $${log}; \
	      echo "Finished testing $${tname} $(TEST_FLAGS)"; \
	      cat $${log}; \
	   fi; \
	fi
 
# The dummysh.chkexe here prevents the target from being
# empty if there are no tests in the current directory.
# $${log} is the log file.
# $${tname} is the name of test.
$(TEST_SCRIPT_CHKSH) $(TEST_SCRIPT_PARA_CHKSH) dummysh.chkexe_:
	@if test "X$@" != "X.chkexe_" && test "X$@" != "Xdummysh.chkexe_"; then \
	   cmd=$(@:.chkexe_=);\
	   tname=`basename $$cmd`;\
	   chkname=`basename $(@:.chkexe_=.chkexe)`;\
	   log=`basename $(@:.chkexe_=.chklog)`; \
	   echo "============================"; \
	   if $(top_srcdir)/bin/newer $${chkname} $$cmd $(SCRIPT_DEPEND); then \
	      echo "No need to test $${tname} again."; \
	   else \
	      echo "============================" > $${log}; \
	      if test "X$(HDF_FORTRAN)" = "Xyes"; then \
	         echo "Fortran API: Testing $${tname} $(TEST_FLAGS)"; \
	         echo "Fortran API: $${tname} $(TEST_FLAGS) Test Log" >> $${log}; \
	      elif test "X$(HDF_CXX)" = "Xyes"; then \
	         echo "C++ API: Testing $${tname} $(TEST_FLAGS)"; \
	         echo "C++ API: $${tname} $(TEST_FLAGS) Test Log" >> $${log}; \
	      else \
	         echo "Testing $${tname} $(TEST_FLAGS)"; \
	         echo "$${tname} $(TEST_FLAGS) Test Log" >> $${log}; \
	      fi; \
	      echo "============================" >> $${log}; \
	      RUNSERIAL="$(RUNSERIAL)" RUNPARALLEL="$(RUNPARALLEL)"           \
	      srcdir="$(srcdir)" \
	         $(TIME) $(SHELL) $$cmd $(TEST_FLAGS) >> $${log} 2>&1 \
	         && touch $${chkname} || \
	         (test $$HDF5_Make_Ignore && echo "*** Error ignored") || \
	         (cat $${log} && false) || exit 1; \
	      echo "" >> $${log}; \
	      echo "Finished testing $${tname} $(TEST_FLAGS)" >> $${log}; \
	      echo "============================" >> $${log}; \
	      echo "Finished testing $${tname} $(TEST_FLAGS)"; \
	      cat $${log}; \
	   fi; \
	   echo "============================"; \
	fi

# Actual execution of check-p.
_check-p: $(LIB) $(PROGS) $(TESTS)
	@if test -n "$(TEST_PROGS_PARA)$(TEST_SCRIPTS_PARA)"; then      \
	   echo "===Parallel tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="; \
	fi
	@if test -n "$(TEST_PROGS_PARA)"; then                                \
	    echo "**** Hint ****";                                            \
	    echo "Parallel test files reside in the current directory"        \
	         "by default.";                                               \
	    echo "Set HDF5_PARAPREFIX to use another directory. E.g.,";       \
	    echo "    HDF5_PARAPREFIX=/PFS/user/me";                          \
	    echo "    export HDF5_PARAPREFIX";                                \
	    echo "    make check";                                            \
	    echo "**** end of Hint ****";                                     \
	fi
	@for test in $(TEST_PROGS_PARA) dummy; do                             \
	   if test $$test != dummy; then                                      \
	      $(MAKE) $(AM_MAKEFLAGS) $$test.chkexe_ \
	      RUNTEST="$(RUNPARALLEL)" || exit 1; \
	   fi;                                                                \
	done
	@for test in $(TEST_SCRIPTS_PARA) dummy; do                           \
	  if test $$test != dummy; then                                      \
	    $(MAKE) $(AM_MAKEFLAGS) $$test.chkexe_ || exit 1; \
	  fi;                                                                 \
	done
	@if test -n "$(TEST_PROGS_PARA)$(TEST_SCRIPTS_PARA)"; then      \
	    echo "===Parallel tests in `echo ${PWD} | sed -e s:.*/::` ended `date`===";\
	fi

## Make installation directories directories if they don't exist.
$(libdir):
	$(top_srcdir)/bin/mkdirs $@

$(includedir):
	$(top_srcdir)/bin/mkdirs $@

$(bindir):
	$(top_srcdir)/bin/mkdirs $@

$(DOCDIR):
	$(top_srcdir)/bin/mkdirs $@

$(EXAMPLEDIR):
	$(top_srcdir)/bin/mkdirs $@

## Install the library, the public header files, and public programs.
## Install h5fc as h5pfc for parallel mode.
install: $(PUB_LIB) $(PUB_HDR) $(PUB_PROGS) $(libdir) $(includedir) $(bindir)
	@for f in X $(PUB_LIB); do					      \
	   if test $$f != X; then					      \
	      (cd .libs && rm -f $(PUB_LIB)i && ln -s ../$(PUB_LIB) $(PUB_LIB)i); \
	      ($(LT_INSTALL_LIB) $$f $(libdir)/. || exit 1);		      \
	   fi;								      \
	done
	@if test -n "$(F9XMODEXT)" && test "X$(F9XMODEXT)" != "Xo"; then      \
	   if test -f "hdf5.$(F9XMODEXT)" || test -f "HDF5.$(F9XMODEXT)"; then \
	      ((cp *.$(F9XMODEXT) $(libdir)/. && chmod 644 $(libdir)/*.$(F9XMODEXT)) || exit 1); \
	   fi;                                                                \
	fi
	@for f in X $(PUB_HDR); do					      \
	   if test $$f != X; then					      \
	      if test -f $$f; then					      \
		 (set -x; $(INSTALL_DATA) $$f $(includedir)/. || exit 1);     \
	      else							      \
	         (set -x; $(INSTALL_DATA) $(srcdir)/$$f $(includedir)/. ||    \
		    exit 1);						      \
	      fi;							      \
	   fi;								      \
	done
	@for f in X $(PUB_PROGS); do					      \
	   if test $$f != X; then					      \
	     if test "X$(PARALLEL)" != "Xno" -a $$f = "h5fc"; then            \
	      ($(LT_INSTALL_PROG) $$f $(bindir)/h5pfc || exit 1);	      \
	     else							      \
	      ($(LT_INSTALL_PROG) $$f $(bindir)/. || exit 1);		      \
	     fi;                                                              \
	   fi;							              \
	done

install-examples: $(EXAMPLE_PROGS) $(EXAMPLEDIR)
	@for f in X $(EXAMPLE_PROGS); do				      \
	  if test $$f != X; then					      \
	    (set -x; $(INSTALL_DATA) $(srcdir)/$$f $(EXAMPLEDIR)/. || exit 1);\
	  fi;								      \
	done

uninstall-examples:
	@if test -n "$(EXAMPLE_PROGS)"; then 				      \
	  set -x; cd $(EXAMPLEDIR) && $(RM) $(EXAMPLE_PROGS);		      \
	fi

## Removes those things that `make install' (would have) installed.
uninstall:
	@if test -n "$(PUB_HDR)"; then  				      \
	  set -x; cd $(includedir) && $(RM) $(PUB_HDR);			      \
	fi
	@if test -n "$(F9XMODEXT)" -a "X$(F9XMODEXT)" != "Xo"; then	      \
	  if test -f "$(libdir)/hdf5.$(F9XMODEXT)" -o -f "$(libdir)/HDF5.$(F9XMODEXT)"; then \
	    for f in $(FORTMODFILES); do                                    \
	      set -x; $(RM) $(libdir)/$$f.$(F9XMODEXT);                     \
	    done                                                            \
	  fi;								      \
	fi
	@for f in X $(PUB_PROGS); do					      \
	  if test $$f != X; then					      \
	    $(LT_UNINSTALL) $(bindir)/$$f;				      \
	  fi;								      \
	done

# General rule for recursive cleaning targets.  Like the rule above,
# but doesn't require building BUILT_SOURCES.
check-clean:
	@$(MAKE) $(AM_MAKEFLAGS) build-$@ || exit 1;
	@for d in X $(SUBDIRS); do                                             \
	    if test $$d != X && test $$d != .; then                           \
	        (set -x; cd $$d && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;     \
	    fi;                                                               \
	done

# clean up files generated by tests so they can be re-run.
build-check-clean:
	$(RM) -rf $(CHECK_CLEANFILES)

## Removes temporary files without removing the final target files.  That is,
## remove things like object files but not libraries or executables.
##
mostlyclean: build-check-clean
	@if test -n "$(LIB_OBJ)"; then                                        \
	    $(RM) $(LIB_OBJ) $(LIB_OBJ:.lo=.o);                               \
	fi
	@if test -n "$(TEST_OBJ)"; then                                       \
	    $(RM) $(TEST_OBJ) $(TEST_OBJ:.lo=.o);                             \
	fi
	@if test -n "$(PROG_OBJ)" || test -n "$(MOSTLYCLEAN)"; then           \
	    $(RM) $(PROG_OBJ) $(PROG_OBJ:.lo=.o) $(MOSTLYCLEAN);              \
	fi

## Like `mostlyclean' except it also removes the final targets: things like
## libraries and executables.  This target doesn't remove any file that
## is part of the HDF5 distribution.
##
clean: mostlyclean
	@if test -n "$(LIB)" -o -n "$(TEST_PROGS)" -o -n "$(TEST_PROGS_PARA)" -o -n "$(PROGS)" -o -n "$(CLEAN)"; then \
	    $(RM) $(LIB) $(TEST_PROGS) $(TEST_PROGS_PARA) $(PROGS) $(CLEAN);                     \
	fi
	-$(RM) *.M *.a *.mod
	-$(RM) -r .libs

## Like `clean' except it also removes files that were created by running
## configure.  If you've unpacked the source and built HDF5 without creating
## any other files, then `make distclean' will leave only the files that were
## in the distribution.
##
distclean: clean
	-$(RM) .depend TAGS *~ core *.core *.bak *.old *.new $(DISTCLEAN)
	@if test -f $(srcdir)/Makefile.in; then				      \
	   (set -x; $(RM) Makefile);					      \
	fi

## Like `distclean' except it deletes all files that can be regenerated from
## the makefile, including those generated from autoheader and autoconf.
##
maintainer-clean: distclean
	-$(RM) *~ core core.* *.core *.bak *.contrib gmon.out

## Implicit rules
.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

.c.lo:
	@$(LT_CCOMPILE) $(CFLAGS) $(CPPFLAGS) -c $<

.f90.o:
	$(F9X) $(FFLAGS) $(FSEARCH_DIRS) -c $<

.f90.lo:
	@$(LT_FCOMPILE) $(FFLAGS) $(FSEARCH_DIRS) -c $<

##-----------------------------------------------------------------------------
## The following section of this makefile contains dependencies between the
## source files and the header files.
##-----------------------------------------------------------------------------
.PHONY: dep depend
dep depend: $(srcdir)/Dependencies

##                                                             -*- makefile -*-
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the files COPYING and Copyright.html.  COPYING can be found at the root
# of the source code distribution tree; Copyright.html can be found at the
# root level of an installed copy of the electronic HDF5 document set and
# is linked from the top-level documents page.  It can also be found at
# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
# access to either file, you may request a copy from help@hdfgroup.org.



## We keep a list of dependencies in `.depend' for each of the source
## files on which it depends.  When one of the source files is modified
## we remove its record from .depend and regenerate its dependencies,
## tacking them onto the end of .depend.  By including the .depend file
## into the makefile, we're telling make that the makefile depends on
## the dependency list in the .depend file.
##
## This is as fast as the `.d' method described in the GNU make manual
## for automatic dependencies, but has the added advantage that all
## dependencies are stored in one place.  The advantage over the
## `makedepend' program is that only those files that are out of date
## have dependency information rebuilt, and the Makefile is not
## modified.
##
## This is also where tracing information is updated.  The $(TRACE)
## program is run on each source file to make sure that the H5TRACE()
## macros are up to date.  If they are then the file is not modified,
## otherwise the file is changed and a backup is saved by appending a
## tilde to the file name.
##
$(srcdir)/Dependencies: .depend
	@if test "$(srcdir)" != "."; then					\
	  (									\
	    echo '##'; 								\
	    echo '## This file is machine generated on GNU systems.';		\
	    echo '## Only temporary changes may be made here.';			\
	    echo '##'; 								\
	    echo '## Copyright by the Board of Trustees of the University of Illinois.'; \
	    echo '## All rights reserved.'; 					\
	    echo '##'; 								\
	    echo '## This file is part of HDF5.  The full HDF5 copyright notice, including'; \
	    echo '## terms governing use, modification, and redistribution, is contained in'; \
	    echo '## the files COPYING and Copyright.html.  COPYING can be found at the root'; \
	    echo '## of the source code distribution tree; Copyright.html can be found at the'; \
	    echo '## root level of an installed copy of the electronic HDF5 document set and'; \
	    echo '## is linked from the top-level documents page.  It can also be found at'; \
	    echo '## http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have'; \
	    echo '## access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu.'; \
	    echo;								\
	  ) >$@;								\
	  $(PERL) -p $(top_srcdir)/bin/distdep .depend >>$@;			\
	else									\
	  echo 'Dependencies cannot be built when $$srcdir == $$builddir';	\
	fi

.depend: $(LIB_SRC) $(TEST_SRC) $(PROG_SRC)
	@touch .depend
	@for dep in $? dummy; do                                              \
      if test $$dep != "dummy" -a -n "$(PERL)"; then                      \
        case "$$dep" in                                                   \
          *.c)                                                            \
            echo Building dependencies for $$dep;                         \
            obj=`basename $$dep .c`.lo;                                   \
            sed '\%^'"$$obj"':%,\%[^\\]$$%d' <$@ >$@- && mv $@- $@;       \
            $(TRACE) $$dep;                                               \
            $(CC) -MM -MG $(CPPFLAGS) $$dep 2>/dev/null >>$@;             \
            $(PERL) -w $(top_srcdir)/bin/dependencies --srcdir=$(srcdir) --top_srcdir=$(top_srcdir) --top_builddir=$(top_builddir) $@; \
            ;;                                                            \
        esac;                                                             \
      fi;                                                                 \
    done

-include .depend

