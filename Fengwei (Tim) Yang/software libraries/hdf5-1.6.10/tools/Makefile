## Tools HDF5 Makefile(.in)
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

CONFIG = ordered

## All subdirectories except for the ``lib'' subdirectory.
##
SUBDIRS=h5diff h5ls misc h5import h5repack h5jam h5dump

## Add include directory to the C preprocessor flags, add -lh5tools and
## -lhdf5 to the list of libraries.
##
CPPFLAGS=-I. -I$(srcdir) -I$(top_builddir)/src -I$(top_srcdir)/src  \
         -I$(top_srcdir)/tools/lib -UH5_DEBUG_API -DNDEBUG -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_POSIX_SOURCE -D_BSD_SOURCE 

## Test programs and scripts.
##
TEST_PROGS=
TEST_SCRIPTS=

## These are our main targets: library and tools.
##
LIBTOOLS=$(top_builddir)/tools/lib/libh5tools.la
LIBHDF5=$(top_builddir)/src/libhdf5.la

PROGS=$(PUB_PROGS) $(TEST_PROGS)

## Source and object files for the tests
##
TEST_SRC=
TEST_OBJ=$(TEST_SRC:.c=.lo)

lib ::
	(cd lib && $(MAKE) $@) || exit 1;

progs: $(LIBTOOLS) $(LIBHDF5)
	@:; for d in X $(SUBDIRS); do                                    \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@) || exit 1;                                 \
	  fi;                                                                 \
	done

dep depend:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@) || exit 1;                                 \
	  fi;                                                                 \
	done

check check-s check-p check-clean tests test _test: $(PROGS)
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@) || exit 1;                                 \
	  fi;                                                                 \
	done

uninstall:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@) || exit 1;                                 \
	  fi;                                                                 \
	done

install:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@) || exit 1;                                 \
	  fi;                                                                 \
	done

install-doc:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@);                                           \
	  fi;                                                                 \
	done

.PHONY: all lib progs test _test install uninstall clean           \
	mostlyclean distclean maintainer-clean

clean mostlyclean:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@);                                           \
	  fi;                                                                 \
	done

distclean:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@);                                           \
	  fi;                                                                 \
	done
	-$(RM) Makefile

maintainer-clean:
	@:; for d in X lib $(SUBDIRS); do                                \
	  if test $$d != X; then                                              \
	    (cd $$d && $(MAKE) $@);                                           \
	  fi;                                                                 \
	done

