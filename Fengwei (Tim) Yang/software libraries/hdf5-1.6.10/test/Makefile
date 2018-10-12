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
## HDF5 Library Test Makefile(.in)
##
top_srcdir=..
top_builddir=..
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
##
##
## The following section of this makefile comes from the
## `./config/commence' file which was generated with config.status
## from `./config/commence.in'.
##-----------------------------------------------------------------------------

## Things that Make needs
.SUFFIXES:
.SUFFIXES: .c .o .lo


## Directories to search
VPATH=$(srcdir) $(top_builddir)/src $(top_srcdir)/src

## Programs
SHELL=/bin/sh
CC=mpicc
## H5_CPPFLAGS and H5_CFLAGS are used for extra flags used to build HDF5
## (flags that don't need to be included in h5cc, for instance).
CFLAGS= -ansi -O -fomit-frame-pointer -finline-functions 
CPPFLAGS=-UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE  
H5_CFLAGS= -ansi -O -fomit-frame-pointer -finline-functions
H5_CPPFLAGS=-UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE 
LDFLAGS=
ROOT=/home/f/fy/fy44/Programs/hdf5-1.6.10
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
RUNPARALLEL=/usr/bin/mpirun -np $${NPROCS:=3}
RUNTEST=$(RUNSERIAL)
TIME=time
TRACE=:

## Installation points
## Insert DESTDIR to individual installation directory because they can be
## defined independent of prefix.
prefix=/home/f/fy/fy44/hdf5-1.6.10
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
LT_COMPILE=$(LT) --mode=compile $(CC)
LT_LINK_LIB=$(LT) --mode=link $(CC) -rpath $(libdir) $(DYNAMIC_DIRS)
LT_LINK_EXE=$(LT) --mode=link $(CC) $(LT_STATIC_EXEC) -dlopen self -rpath $(libdir) $(DYNAMIC_DIRS)
LT_RUN=$(LT) --mode=execute
LT_INSTALL_PROG=$(LT) --mode=install $(INSTALL_PROGRAM)
LT_INSTALL_LIB=$(LT) --mode=install $(INSTALL_DATA)
LT_UNINSTALL=$(LT) --mode=uninstall $(RM)

## Optional variables. We must declare them here because Irix pmake
## complains if it sees a reference to a variable which has never been
## defined. The main makefile is free to redefine these to something else.
DOCDIR=$(docdir)
SUBDOCDIR=
EXAMPLEDIR=$(docdir)/hdf5/examples/c
LIB=
LIB_SRC=
LIB_OBJ=
PUB_HDR=
PUB_PROGS=
PUB_DOCS=
PUB_SUBDOCS=
PROGS=
TEST_PROGS=
PRIV_PROGS=
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
	@for sub in   dummy; do					\
	  if test $$sub != dummy; then					\
	    -cd $$sub; $(MAKE) reconfigure;				\
	  fi;								\
	done;

## make help target
help:
	$(top_srcdir)/bin/makehelp

##-----------------------------------------------------------------------------
## The following section of this makefile comes from the middle of
## `Makefile.in' from this directory.  It was generated by running
## `config.status'.
##-----------------------------------------------------------------------------

## Add include directory to the C preprocessor flags and the h5test and hdf5
## libraries to the library list.
CPPFLAGS=-I. -I$(srcdir) -I../src -I$(top_srcdir)/src -UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE 

## These are our main targets. They should be listed in the order to be
## executed, generally most specific tests to least specific tests.
TEST_PROGS=testhdf5 lheap ohdr stab gheap cache hyperslab istore bittests    \
    dtypes dsets cmpd_dset extend external links unlink big mtime fillval    \
    mount flush1 flush2 enum set_extent     \
    ttsafe getname file_handle     \
    ntypes dangle reserved err_vers

# Test scripts.
TEST_SCRIPTS = testlibinfo.sh testcheck_version.sh

TIMINGS=testmeta
PROGS=$(TEST_PROGS) tcheck_version

## The libh5test.a library provides common support code for the tests. We link
## this library statically because some systems can only link executables to
## a single shared library and libhdf5 is much bigger than libh5test.
LIBHDF5=../src/libhdf5.la
LIB=libh5test.la
LIB_SRC=h5test.c testframe.c
LIB_OBJ=$(LIB_SRC:.c=.lo)
PUB_LIB=

## Temporary files.  These files are the ones created by setting the
## HDF5_NOCLEANUP environment variable and running `make test' without
## specifying a file prefix or low-level driver.  Changing the file
## prefix or low-level driver with environment variables will influence 
## the temporary file name in ways that the makefile is not aware of.
MOSTLYCLEAN=cmpd_dset.h5 compact_dataset.h5 dataset.h5 extend.h5 istore.h5   \
            tfile[1-3].h5 th5s1.h5 th5s3.h5 lheap.h5 ohdr.h5 stab[1-2].h5 \
            extern_[1-3].h5 extern_[1-4][ab].raw gheap[0-4].h5      \
            links.h5 links[1-3].h5 big.data big[0-9][0-9][0-9][0-9][0-9].h5  \
            sec2.h5 stdio.h5 dtypes[1-7].h5 tattr.h5 coord.h5         \
            tselect.h5 mtime.h5 unlink.h5 fillval_[0-9].h5 fillval.raw       \
            mount_[0-9].h5 testmeta.h5 ttime.h5 trefer[1-3].h5 tvltypes.h5   \
            tvlstr.h5 tvlstr2.h5 flush.h5 enum1.h5 titerate.h5 ttsafe.h5     \
            tarray1.h5 tgenprop.h5 tmisc1.h5 tmisc2a.h5 tmisc2b.h5 tmisc3.h5 \
            tmisc4a.h5 tmisc4b.h5 tmisc5.h5 tmisc6.h5 tmisc7.h5 tmisc8.h5    \
            tmisc9.h5 tmisc10.h5 tmisc11.h5 tmisc12.h5 tmisc13a.h5           \
            tmisc13b.h5 tmisc14.h5 set_extent_read.h5 set_extent_create.h5   \
            getname.h5 getname1.h5 getname2.h5 getname3.h5 sec2_file.h5      \
            family_file000[0-3][0-9].h5 new_family_v18_000[0-3][0-9].h5      \
            multi_file-[rs].h5 core_file new_move_[ab].h5 ntypes.h5          \
            dangle.h5 test_filters.h5 tstint[1-2].h5 unlink_chunked.h5       \
            err_vers.h5

CLEAN=$(TIMINGS)
DISTCLEAN=testerror.sh testlibinfo.sh testcheck_version.sh

## Source and object files for programs...  The TEST_SRC list contains all the
## source files and is used for things like dependencies, archiving, etc.  The
## other source lists are for the individual tests, the files of which may
## overlap with other tests.

TEST_SRC=big.c bittests.c cache.c cmpd_dset.c dsets.c dtypes.c extend.c      \
    external.c fillval.c flush1.c flush2.c gheap.c h5test.c hyperslab.c      \
    istore.c lheap.c links.c mount.c mtime.c ohdr.c stab.c tarray.c tattr.c  \
    tcheck_version.c tconfig.c tcoords.c testhdf5.c testmeta.c tfile.c       \
    tgenprop.c th5s.c theap.c titerate.c tmeta.c tmisc.c trefer.c trefstr.c  \
    tselect.c tskiplist.c ttime.c ttst.c tvltypes.c tvlstr.c                 \
    unlink.c enum.c ttsafe.c ttsafe_dcreate.c ttsafe_error.c                 \
    ttsafe_cancel.c ttsafe_acreate.c set_extent.c getname.c file_handle.c    \
    ntypes.c dangle.c reserved.c err_vers.c

TEST_OBJ=$(TEST_SRC:.c=.lo)

VFD_LIST=sec2 stdio core split multi family

## Private header files (not to be installed)...
PRIVATE_HDR=testhdf5.h ttsafe.h

## Additional targets
.PHONY: timings _timings
timings _timings: $(TIMINGS)
	@for timing in $(TIMINGS) dummy; do				      \
	   if test $$timing != dummy; then				      \
	      echo "Running $$timing $(TEST_FLAGS)"; 			      \
	      $(RUNTEST) ./$$timing $(TEST_FLAGS) || exit 1; 		      \
	   fi;								      \
	done;

# The flush1 test must run before the flush2 test
flush2.chkexe_: flush1.chkexe_

## Programs have to be built before they can be tested!
##
check test _test: $(PROGS)

## How to build the tests...  They all depend on the test and hdf5 libraries.
$(PROGS): $(LIB) $(LIBHDF5)

TESTHDF5_OBJ=testhdf5.lo tarray.lo tattr.lo tconfig.lo tfile.lo               \
    tgenprop.lo th5s.lo tcoords.lo theap.lo titerate.lo tmeta.lo tmisc.lo     \
    trefer.lo trefstr.lo tselect.lo tskiplist.lo ttime.lo ttst.lo             \
    tvltypes.lo tvlstr.lo

TTS_OBJ=ttsafe.lo ttsafe_dcreate.lo ttsafe_error.lo             \
    ttsafe_cancel.lo ttsafe_acreate.lo

testhdf5: $(TESTHDF5_OBJ)
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ $(TESTHDF5_OBJ) $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

lheap: lheap.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ lheap.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

cache: cache.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ cache.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

ohdr: ohdr.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ ohdr.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

stab: stab.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ stab.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

gheap: gheap.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ gheap.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

dsets: dsets.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ dsets.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

bittests: bittests.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ bittests.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

dtypes: dtypes.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ dtypes.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

hyperslab: hyperslab.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ hyperslab.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

istore: istore.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ istore.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

cmpd_dset: cmpd_dset.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ cmpd_dset.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

extend: extend.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ extend.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

external: external.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ external.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

big: big.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ big.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

links: links.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ links.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

mtime: mtime.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ mtime.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

unlink: unlink.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ unlink.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

fillval: fillval.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ fillval.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

mount: mount.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ mount.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

flush1: flush1.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ flush1.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

flush2: flush2.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ flush2.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

enum: enum.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ enum.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

tcheck_version: tcheck_version.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ tcheck_version.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

ttsafe: $(TTS_OBJ)
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ $(TTS_OBJ) $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

set_extent: set_extent.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ set_extent.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

testmeta: testmeta.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ testmeta.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

getname: getname.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ getname.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

file_handle: file_handle.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ file_handle.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

ntypes: ntypes.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ ntypes.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

dangle: dangle.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ dangle.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

reserved: reserved.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ reserved.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

err_vers: err_vers.lo
	@$(LT_LINK_EXE) $(CFLAGS) -o $@ err_vers.lo $(LIB) $(LIBHDF5) $(LDFLAGS) $(LIBS)

## Nothing is installed from test directory when 'make install' is invoked, but
## this will force all test targets to be made, as other directories may need
## the test library built in order to 'make install' (i.e., perform directory).
install: all

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
##
##
## The following section of this makefile comes from the
## `./config/conclude' file which was generated with config.status
## from `./config/conclude.in'.
##-----------------------------------------------------------------------------

## This is the target for the library described in the main body of the
## makefile.
##
lib: $(LIB)
$(LIB) __no_library__: $(LIB_OBJ)
	@$(LT_LINK_LIB) -o $@ $(CFLAGS) $(LIB_OBJ) $(LDFLAGS) $(LIBS)
	@for d in X $(SUBDIRS); do                                      \
	    if test $$d != X; then                                      \
	        (set -x; cd $$d && $(MAKE) lib) || exit 1;              \
	    fi;                                                         \
	done

progs: $(LIB) $(PROGS)
	@for d in X $(SUBDIRS); do                                      \
	    if test $$d != X; then                                      \
	        (set -x; cd $$d && $(MAKE) progs) || exit 1;            \
	    fi;                                                         \
	done

## Build a tags file in this directory.
##
TAGS: $(LIB_SRC)
	$(RM) $@
	-etags $(LIB_SRC)

## Runs each test in order, passing $(TEST_FLAGS) to the program.
## Since tests are done in a shell loop, "make -i" does apply inside it.
## Set HDF5_Make_Ignore to a non-blank string to ignore errors inside the loop.
## The timestamps give a rough idea how much time the tests use.
##
tests: $(TEST_PROGS) $(TEST_PROGS_PARA) $(TEST_SCRIPTS) $(LIB)
check test _test: tests
	@echo "===Tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-s
	$(MAKE) _check-p
	@echo "===Tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

## serial tests
check-s: tests
	@echo "===Serial tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-s
	@echo "===Serial tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

## parallel tests
check-p: tests
	@echo "===Parallel tests in `echo ${PWD} | sed -e s:.*/::` begin `date`==="
	$(MAKE) _check-p
	@echo "===Parallel tests in `echo ${PWD} | sed -e s:.*/::` ended `date`==="

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


## Run test with different Virtual File Driver
check-vfd:
	@for vfd in $(VFD_LIST) dummy; do				      \
	    if test $$vfd != dummy; then                                      \
		echo "============================";                          \
		echo "Testing Virtual File Driver $$vfd";		      \
		echo "============================";                          \
                $(MAKE) $(AM_MAKEFLAGS) check-clean || exit 1;                \
		HDF5_DRIVER=$$vfd $(MAKE) check;			      \
	    fi;								      \
	done

## Make installation directories directories if they don't exist.
##
$(libdir):
	-$(top_srcdir)/bin/mkdirs $@

$(includedir):
	-$(top_srcdir)/bin/mkdirs $@

$(bindir):
	-$(top_srcdir)/bin/mkdirs $@

$(DOCDIR):
	-$(top_srcdir)/bin/mkdirs $@

$(EXAMPLEDIR):
	-$(top_srcdir)/bin/mkdirs $@

## Install the library, the public header files, and public programs.
##
install: $(AUX_LIB) $(PUB_LIB) $(PUB_HDR) $(PUB_PROGS) $(libdir) $(includedir) $(bindir)
	@for f in X $(PUB_LIB); do					      \
	   if test $$f != X; then					      \
	      ($(LT_INSTALL_LIB) $$f $(libdir)/. || exit 1);		      \
	   fi;								      \
	done
	@if test -f libhdf5.settings; then				      \
	   (set -x; $(INSTALL_DATA) libhdf5.settings $(libdir)/. || exit 1);  \
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
	     if test "X$(PARALLEL)" != "Xno"; then                            \
	       if test $$f = "h5cc"; then                     		      \
	         ($(LT_INSTALL_PROG) $$f $(bindir)/h5pcc || exit 1);          \
	       elif test $$f = "h5fc"; then                                   \
	         ($(LT_INSTALL_PROG) $$f $(bindir)/h5pfc || exit 1);          \
	       else							      \
	         ($(LT_INSTALL_PROG) $$f $(bindir)/. || exit 1);	      \
	       fi;							      \
	     else                                                             \
	       ($(LT_INSTALL_PROG) $$f $(bindir)/. || exit 1);		      \
	     fi;							      \
	   fi;							              \
	done
	@for d in X $(SUBDIRS); do                                            \
	    if test $$d != X; then                                            \
	       (set -x; cd $$d && $(MAKE) $@) || exit 1;                      \
	    fi;                                                               \
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

## Install the documents.
##
install-doc: $(PUB_DOCS) $(DOCDIR)
	@if test -n "$(DOCDIR)"; then                                         \
	    if test -d "$(DOCDIR)"; then                                      \
	       :;                                                             \
	    else                                                              \
	       mkdir $(DOCDIR) && chmod 755 $(DOCDIR);                        \
	    fi;                                                               \
	fi
	@if test -n "$(SUBDOCDIR)"; then                                      \
	    if test -d "$(SUBDOCDIR)" ; then                                  \
	       :;                                                             \
	    else                                                              \
	       mkdir $(SUBDOCDIR) && chmod 755 $(SUBDOCDIR);                  \
	    fi;                                                               \
	fi
	@for d in X $(SUBDIRS); do                                            \
	   if test $$d != X; then                                             \
	      (set -x; cd $$d && $(MAKE) $@) || exit 1;                       \
	   fi;                                                                \
	done
	@for f in X $(PUB_DOCS); do                                           \
	   if test $$f != X; then                                             \
	      if test -f $$f; then                                            \
		 (set -x; $(INSTALL_DATA) $$f $(DOCDIR)/. || exit 1);         \
	      else                                                            \
	         (set -x; $(INSTALL_DATA) $(srcdir)/$$f $(DOCDIR)/. ||        \
		    exit 1);                                                  \
	      fi;                                                             \
	   fi;                                                                \
	done
	@for f in X $(PUB_SUBDOCS); do                                        \
	   if test $$f != X; then                                             \
	      if test -f $$f; then                                            \
		 (set -x; $(INSTALL_DATA) $$f $(SUBDOCDIR)/. || exit 1);      \
	      else                                                            \
	         (set -x; $(INSTALL_DATA) $(srcdir)/$$f $(SUBDOCDIR)/. ||     \
		    exit 1);                                                  \
	      fi;                                                             \
	   fi;                                                                \
	done

## Removes those things that `make install' (would have) installed.
##
uninstall:
	@for f in libhdf5.settings $(LIB); do				      \
          $(LT_UNINSTALL) $(libdir)/$$f;				      \
	done
	@if test -n "$(PUB_HDR)"; then                                        \
	  set -x; cd $(includedir) && $(RM) $(PUB_HDR);			      \
	fi
	@for f in X $(PUB_PROGS); do					      \
	    if test $$f != X; then					      \
	        $(LT_UNINSTALL) $(bindir)/$$f;				      \
	    fi;								      \
	done

## Removes those things that `make install-doc' (would have) installed.
##
uninstall-doc:
	@for d in X $(SUBDIRS); do					      \
	  if test $$d != X; then					      \
	    (set -x; cd $$d && $(MAKE) $@) || exit 1;                         \
	  fi;								      \
	done
	@if test -n "$(PUB_DOCS)"; then 				      \
	  set -x; cd $(DOCDIR) && $(RM) $(PUB_DOCS);			      \
	fi
	@if test -n "$(PUB_SUBDOCS)"; then				      \
	  set -x; cd $(DOCDIR) && $(RM) $(PUB_SUBDOCS);			      \
	fi

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
	@for d in X $(SUBDIRS); do                                            \
	    if test $$d != X; then                                            \
	        (set -x; cd $$d; $(MAKE) $@) || exit 1;                       \
	    fi;                                                               \
	done
	@if test -n "$(LIB_OBJ)"; then                                        \
	    $(RM) $(LIB_OBJ) $(LIB_OBJ:.lo=.o);                              \
	fi
	@if test -n "$(TEST_OBJ)"; then                                       \
	    $(RM) $(TEST_OBJ) $(TEST_OBJ:.lo=.o);                            \
	fi
	@if test -n "$(PROG_OBJ)" || test -n "$(MOSTLYCLEAN)"; then           \
	    $(RM) $(PROG_OBJ) $(PROG_OBJ:.lo=.o) $(MOSTLYCLEAN);             \
	fi

## Like `mostlyclean' except it also removes the final targets: things like
## libraries and executables.  This target doesn't remove any file that
## is part of the HDF5 distribution.
##
clean: mostlyclean
	@for d in X $(SUBDIRS); do                                            \
	    if test $$d != X; then                                            \
	        (set -x; cd $$d; $(MAKE) $@) || exit 1;                       \
	    fi;                                                               \
	done
	@if test -n "$(LIB)" -o -n "$(TEST_PROGS)" -o -n "$(TEST_PROGS_PARA)" -o -n "$(PROGS)" -o -n "$(CLEAN)"; then \
	    $(RM) $(LIB) $(TEST_PROGS) $(TEST_PROGS_PARA) $(PROGS) $(CLEAN);                     \
	fi;
	-$(RM) -r .libs

## Like `clean' except it also removes files that were created by running
## configure.  If you've unpacked the source and built HDF5 without creating
## any other files, then `make distclean' will leave only the files that were
## in the distribution.
##
distclean: clean
	-$(RM) .depend TAGS *~ core *.core *.bak *.old *.new $(DISTCLEAN)
	@for d in X $(SUBDIRS); do					      \
	  if test $$d != X; then					      \
	    (set -x; cd $$d; $(MAKE) $@) || exit 1;			      \
	  fi;								      \
	done
	@if test -f $(srcdir)/Makefile.in; then				      \
	   (set -x; $(RM) Makefile);					      \
	fi

## Like `distclean' except it deletes all files that can be regenerated from
## the makefile, including those generated from autoheader and autoconf.
##
maintainer-clean: distclean
	@for d in X $(SUBDIRS); do                                            \
	    if test $$d != X; then                                            \
	        (set -x; cd $$d; $(MAKE) $@) || exit 1;                       \
	    fi;                                                               \
	done
	-$(RM) *~ core core.* *.core *.bak *.contrib gmon.out

## Implicit rules
##
.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

.c.lo:
	@$(LT_COMPILE) $(CFLAGS) $(CPPFLAGS) -c $<

##-----------------------------------------------------------------------------
## The following section of this makefile contains dependencies between the
## source files and the header files.
##-----------------------------------------------------------------------------
.PHONY: dep depend
dep depend: $(srcdir)/Dependencies

##                                                             -*- makefile -*-
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
	    echo '#'; 								\
	    echo '# This file is machine generated on GNU systems.';		\
	    echo '# Only temporary changes may be made here.';			\
	    echo '#'; 								\
            echo '# Copyright by The HDF Group.';				\
            echo '# Copyright by the Board of Trustees of the University of Illinois.'; \
            echo '# All rights reserved.';					\
            echo '#';								\
            echo '# This file is part of HDF5.  The full HDF5 copyright notice, including'; \
            echo '# terms governing use, modification, and redistribution, is contained in'; \
            echo '# the files COPYING and Copyright.html.  COPYING can be found at the root'; \
            echo '# of the source code distribution tree; Copyright.html can be found at the'; \
            echo '# root level of an installed copy of the electronic HDF5 document set and'; \
            echo '# is linked from the top-level documents page.  It can also be found at'; \
            echo '# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have'; \
            echo '# access to either file, you may request a copy from help@hdfgroup.org.'; \
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

