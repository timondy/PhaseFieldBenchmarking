#--------------------------------------------------------------------------
#
# This is an updated version of the main Makefile for the Phasefield 
# solvers using Paramesh, developed at the University of Leeds.
# It is used to compile the Extractor program that uses the generated HDF5
# files and outputs a selection of other formats.
#
# Assumed environment variables that the user must have set:
#        HDF5_DIR     : The directory where a v1.6 HDF5 is installed
#        PARAMESH_DIR : The directory where paramesh is installed
# 
# Rewritten : Chris Goodyer        March 2012.
#
#--------------------------------------------------------------------------
export cur-dir := $(shell pwd)

# Set the location of the paramesh top directory
export paramesh_dir = $(PARAMESH_DIR)
export phase_dir = $(cur-dir)/../
export extractor_dir = $(cur-dir)
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
export ADD_LIB = -L$(hdf5_dir)/lib -L$(paramesh_dir)/libs -L$(phase_dir)/libs -lhdf5 -lhdf5_hl -lz

#export FFLAGS +=-Wno-unused-variable

# ifort intel 8.0
export FFLAGS = -C -g -r8 -I$(paramesh_dir)/headers -warn all
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers -i8
export CFLAGS = -C -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8
# Optimised
export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers -i8
export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8

# UnOptimised Arc1
#export FFLAGS = -C -g -r8 -I$(paramesh_dir)/headers -traceback
#export CFLAGS = -C -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8

# Optimised Arc1
export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers -xSSE4.2 
export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include -r8 -i8 -xSSE4.2
#export ADD_LIB += -limf

# Optimised CSLIN
export FFLAGS = -O3 -r8 -I$(paramesh_dir)/headers -fpp
export CFLAGS = -Wl,-lm -O3 -I$(paramesh_dir)/headers -I$(hdf5_dir)/include  
#export ADD_LIB += -limf

# Run as a pseudo-2D case?
#export FFLAGS += -DPSEUDO2D

#export FFLAGS += -warn all -g -pg
#export CFLAGS += -warn all -g -pg

#--------------------------------------------------------------------------

.PHONY:	default
default:	paralib src

.PHONY: all
all: src

.PHONY: clean
clean:
	for dir in src; do \
	  $(MAKE) -C $(extractor_dir)/$$dir -f Makefile.intel clean; \
	done

# An example target to match an application directory name other than Tests
# in which the users application files are located.
.PHONY: src
src: 
	$(MAKE) -C $(extractor_dir)/$@ -f Makefile.intel

.PHONY: paralib
paralib: $(paramesh_dir)/libs/libparamesh.a $(paramesh_dir)/libs/libmpi_paramesh.a
	echo '--- New libraries -- removing executable ----'
	rm -f Extractor
