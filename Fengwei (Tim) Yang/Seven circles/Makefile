#--------------------------------------------------------------------------
#
# This is an updated version of the main Makefile for the Phasefield 
# solvers using Paramesh, developed at the University of Leeds.
#
# Assumed environment variables that the user must have set:
#        HDF5_DIR     : The directory where a v1.6 HDF5 is installed
#        PARAMESH_DIR : The directory where paramesh is installed
# 
# Rewritten : Chris Goodyer        March 2012.
#             added debug option   May 2012
#--------------------------------------------------------------------------
export cur-dir := $(shell pwd)

# Set the location of the paramesh top directory
export paramesh_dir = $(PARAMESH_DIR)
export phase_dir = $(cur-dir)
export hdf5_dir=$(HDF5_DIR)

# Define the fortran compiler
export FC = mpif90
export CC = mpicc 

export MPICH_F90 = ifort
export MPICH_CC = icc

#-----------------------------------------------

# Set the desired compilation flags

#export ADD_LIB = -L$(hdf5_dir)/lib -L../libs -lhdf5 -lhdf5_hl -lz -lgfortran -lhdf5_fortran
#export ADD_LIB = -L$(hdf5_dir)/lib -L../libs -lhdf5 -lhdf5_hl -lz -lhdf5_fortran
###export ADD_LIB = -L$(hdf5_dir)/lib -L$(paramesh_dir)/libs -L$(phase_dir)/libs -lhdf5 -lhdf5_hl -lz

export ADD_LIB = -L$(hdf5_dir)/lib -L$(paramesh_dir)/libs -L$(phase_dir)/libs -lhdf5 -lhdf5_hl

#export FFLAGS +=-Wno-unused-variable

# ifort intel 8.0
###export FFLAGS = -C -g -r8 -I$(paramesh_dir)/headers -warn all
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers -i8
###export CFLAGS = -C -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8
# Optimised
###export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers -i8
###export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8

# UnOptimised Arc1
#export FFLAGS = -C -g -r8 -I$(paramesh_dir)/headers -traceback
#export CFLAGS = -C -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8

# Optimised Arc1
###export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers -xSSE4.2 
###export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8 -xSSE4.2
#export ADD_LIB += -limf

# Optimised CSLIN
###export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers -fpp -fpe1
export FFLAGS = -O3 -I$(paramesh_dir)/headers
export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include  
#export ADD_LIB += -limf
#export LDFLAGS += -lsvml

#FFLAGS += -gen-dep

.PHONY:	default
default:	MG-SRC PFBM-SRC

debug: export FFLAGS += -g -warn all -warn nounused -fpe1
debug: export CFLAGS += -g -w3
debug: default 

PFdebug: export FFLAGS += -g -warn all -warn nounused -fpe1
PFdebug: export CFLAGS += -g -w3
PFdebug: MG-SRC PF-SRC 

OMP: export FFLAGS += -openmp
OMP: default


FULLDIRS = MG-SRC PFBM-SRC
#FOO = ${VAR:-default}
#--------------------------------------------------------------------------

#multi

.PHONY: all
all: src 
#multi

.PHONY: clean
clean:
	$(RM) -r *~ libs/libmodules.a 
	for dir in $(FULLDIRS); do \
	  $(MAKE) -C $(phase_dir)/$$dir -f Makefile.intel clean; \
	done

# An example target to match an application directory name other than Tests
# in which the users application files are located.
.PHONY: src
src: 
	$(MAKE) -C $(phase_dir)/$@ -f Makefile.intel

# Compile the generic components in the MG-SRC directory into the appropriate library
.PHONY: MG-SRC
MG-SRC: 
	$(MAKE) -C $(phase_dir)/$@ -f Makefile.intel

# Compile the cell migration code
.PHONY: PFBM-SRC
PFBM-SRC:
	$(MAKE) -C $(phase_dir)/$@ -f Makefile.intel
	
.PHONY: DOCS
DOCS:
	$(MAKE) -C $(phase_dir)/$@ -f Makefile.latex

.PHONY: DOCSclean
DOCSclean:
	$(MAKE) -C $(phase_dir)/DOCS -f Makefile.latex clean

#.PHONY: libs
#libs:
#	mkdir $@

#.PHONY: paralib
#paralib: $(paramesh_dir)/libs/libparamesh.a $(paramesh_dir)/libs/libmpi_paramesh.a
#	echo '--- New libraries -- removing executable ----'
#	rm -f pf.intel multi.intel
	
	
