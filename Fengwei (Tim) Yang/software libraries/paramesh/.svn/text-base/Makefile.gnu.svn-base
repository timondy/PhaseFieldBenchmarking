
#--------------------------------------------------------------------------
#
# This Gmake file will compile the PARAMESH library and create a
# set of library files to which you can link. To use it, make sure
# it is in the PARAMESH root directory.
# It works by running gmake on the Makefile.gnu files which appear
# in the headers, source and mpi_source sub-directories.
# To simply create these PARAMESH library files, type
#     gmake -f Makefile.gnu
# when in the PARAMESH root directory. The library files will
# then be found in a newly created sub-directory called libs.
#
# If you type
#     gmake -f Makefile.gnu Tests
# it will also compile and link the test programs in the Tests
# sub-directory. There is a file called Makefile.gnu inside Tests
# which is used.
#
# To compile and link application files in a sub-directory called
# User_applic you could type
#     gmake -f Makefile.gnu User_applic
# provided you copy Makefile.gnu from Tests to User_applic, and modify
# it appropriately.
#
#
# Written : Ernest Mamikonyan        April 2002.
#
#--------------------------------------------------------------------------
export cur-dir := $(shell pwd)

# Set the location of the paramesh top directory
export paramesh_dir = $(cur-dir)
#export hdf5_dir=/usr/not-backed-up/jrg/local/hdf5-1.6.8/hdf5
#export hdf5_dir=/usr/

#export hdf5_dir=$(HOME)/hdf5
#export hdf5_dir=$(HOME)/SYSHDF5
export hdf5_dir=/usr/local/hdf5

# Define the fortran compiler

# CEG works on feng-gps2
export FC = mpif90 -f90=gfortran44
export CC = mpicc 

# CEG Intel for feng-gps2
export FC = mpif90 -f90=gfortran44
export CC = mpicc 

#CEG For Everest
#export FC = mpif90 
#export CC = mpicc
#export hdf5_dir=$(HOME)/PHASEFIELD/hdf5


#-----------------------------------------------

# Set the desired compilation flags
# SGIs
#export FFLAGS = -cpp -O3 -r8 -I$(paramesh_dir)/headers 

# NAG95 linux - you must use -float-store. 
#export FFLAGS = -O2 -r8 -dusty -w -I$(paramesh_dir)/headers
#export FFLAGS = -O2 -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -g -C -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O2 -float-store -nan -r8  -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -float-store -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -float-store -r8 -C=all -gline -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -I$(paramesh_dir)/headers
#export FFLAGS = -O0 -fdefault-real-8 -fdefault-double-8 -Wall -I$(paramesh_dir)/headers
#export CFLAGS = -O4 -g -Wall -I$(paramesh_dir)/headers

# CEG debug works on feng-gps2
export FFLAGS = -g -fdefault-real-8 -fdefault-double-8 -Wall -I$(paramesh_dir)/headers -I$(hdf5_dir)/include
export CFLAGS = -Wl,-lm -O0 -g -Wall -I$(paramesh_dir)/headers -I$(hdf5_dir)/include
# CEG optimised for feng-gps2
export FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -I$(paramesh_dir)/headers -I$(hdf5_dir)/include
export CFLAGS = -Wl,-lm -O3 -Wall -I$(paramesh_dir)/headers -I$(hdf5_dir)/include

#export ADD_LIB += $(HOME)/PhaseField/hdf5/libs/libhdf5.a


export FFLAGS +=-Wno-unused-variable 


#export ADD_LIB = -L/home/macneice/Autopack/1.3.2/lib -lautopack-myrinet-O



# ifort intel 8.0
#export FFLAGS = -g -check all -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers

# CEG TRying for Everest
#export FFLAGS = -g -check all -r8 -I$(paramesh_dir)/headers
#export FFLAGS = -O4 -r8 -I$(paramesh_dir)/headers

# Halem
#export FFLAGS = -cpp -g -r8 -I$(paramesh_dir)/headers
#export ADD_LIB = -lmpi
#-----------------------------------------------

# Additional libraries to link to. You do not need
# to add the shmem library. This is automatically added
# if you define SHMEM=1 below.
# export ADD_LIB = /usr/lib32/libmpi.so

# FOR ifc, if not using a composite command like mpif90
#export FC = ifc
#export FFLAGS = -O3 -ip -i4 -r8 -I$(paramesh_dir)/headers -I/usr/local/mpich-intel/include
#export ADD_LIB = -L/usr/local/mpich-intel/lib -lmpichf90 -lmpich -lPEPCF90

#-----------------------------------------------

# some compilers can generate make rules to stdout from the source files
# if you have such a compiler, provide the flags, otherwise comment it out
#export MY_CPP := gcc -E -MM -MG  # for the GNU C Preprocessor

#-----------------------------------------------

# SHMEM or MPI ?
# uncomment to use SHMEM
#export SHMEM = 1

#--------------------------------------------------------------------------


.PHONY: all
ifdef SHMEM
all: libs headers source
else
all: libs headers mpi_source source
endif

.PHONY: headers
headers:
	$(MAKE) -C $(paramesh_dir)/$@ -f Makefile.gnu
	cp -f $(paramesh_dir)/headers/libmodules.a $(paramesh_dir)/libs

.PHONY: mpi_source
mpi_source: headers
	$(MAKE) -C $(paramesh_dir)/$@ -f Makefile.gnu
	cp -f $(paramesh_dir)/mpi_source/libmpi_paramesh.a $(paramesh_dir)/libs

.PHONY: source
source: headers
	$(MAKE) -C $(paramesh_dir)/$@ -f Makefile.gnu
	cp -f $(paramesh_dir)/source/libparamesh.a $(paramesh_dir)/libs

.PHONY: clean
clean:
	$(RM) -r *~ libs/libmodules.a 
	for dir in headers mpi_source source Tests ; do \
	  $(MAKE) -C $(paramesh_dir)/$$dir -f Makefile.gnu clean; \
	done

.PHONY: Tests
Tests: all
	$(MAKE) -C $(paramesh_dir)/$@ -f Makefile.gnu

# An example target to match an application directory name other than Tests
# in which the users application files are located.
.PHONY: User_applic
User_applic: all
	$(MAKE) -C $(paramesh_dir)/$@ -f Makefile.gnu

libs:
	mkdir $@

