#LD_SPARSKIT = \
#  ../Unsupported/SPARSKIT2/UNSUPP/BLAS1/blas1.o \
#  ../Unsupported/distdot.o \
#  -L../Unsupported/SPARSKIT2 \
#  -lskit 

LD_SPARSKIT = 


LDFLAGS += -L../libs
ifndef SHMEM
LDFLAGS += -lmpi_paramesh -lmodules
#LDFLAGS += -lmpi_paramesh
endif

# normal LDFLAGS
ifndef T3E
LDFLAGS += -lparamesh -lmodules
#LDFLAGS += -lparamesh
endif

LDFLAGS += $(LD_SPARSKIT)

ifndef SHMEM
LDFLAGS += -lmpi_paramesh -lmodules -lparamesh
#LDFLAGS += -lmpi_paramesh
endif

ifdef T3E
# T3E LDFLAGS
LDFLAGS += -lparamesh -lmodules
endif

#LDFLAGS += -lmodules

ifdef SHMEM
ifdef SGI
# SGI LDFLAGS
LDFLAGS += -lsma
endif
ifdef COMPAQ
# Compaq LDFLAGS
LDFLAGS += -lshmem
endif
endif

LDFLAGS += $(ADD_LIB)

# List all application source files
sources := \
 amr_1blk_bcset.F90 \
 dipole_field.F90 \
 check_data.F90 \
 zero_guardcells.F90 \
 output_tecplot_refmap.F90 \
 convert_rthetaphi_to_xyz.F90 \
 ftheta.F90 \
 report.F90 \
 gtest_neigh_data1.F90 \
 mesh_test.F90

# divb_test_supp.F90

%.o:%.F90
	$(FC) -c $(FFLAGS) $<

objects := $(sources:.F90=.o)

# tests are assumed to be name test_*.F90
# so this strips test_ and .F90 from all files that match test_*.F90
# tests := $(patsubst %.F90,%,$(patsubst test_%,%,$(wildcard test_*.F90)))
tests := \
multigrid \
singular_coord_sph_gcell \
singular_coord_sph_prol \
singular_coord_pol \
derefine_1blk_1 \
derefine_1blk_2 \
guardcell_1blk \
prolong_1blk \
prolong_divb \
prolong_divb2 \
c_to_f_1blk \
c_to_f_1blk_2 \
c_to_f_1blk_3 \
c_to_f_1blk_4 \
c_to_f_1blk_divb \
1blk_guardcell_icoord \
1blk_guardcell_nlayers \
1blk_guardcell_big \
1blk_guardcell_med \
prolong_multilevel_1blk \
checkpoint \
checkpoint1 \
checkpoint_hdf5 \
checkpoint_hdf5_2 \
checkpoint_hdf5_scale \
restrict_1blk \
flux_conserve_1blk \
flux_conserve_1blk_2 \
edges_1blk \
multi_level_1 \
multi_level_2 \
force_consist \
gcell_on \
bcset_example 


# compiles all tests as defined above
.PHONY: all
all: $(wildcard test_*.F90) $(wildcard test_*.c) $(objects)
	for test in $(tests); do $(MAKE) -f Makefile.gnu $$test; done

# compiles one particular test
# this is the target that does the work for the one above
%: test_%.F90 $(objects)
	$(FC) $(FFLAGS) -o test_$@ $^ $(LDFLAGS)

%: test_%.c 
	$(CC) $(CFLAGS) -c test_c_interface.c
	$(FC) $(FFLAGS) -o test_c_interface test_c_interface.o amr_1blk_bcset.o dipole_field.o check_data.o zero_guardcells.o ftheta.o $(LDFLAGS)

$(objects): test_defs.fh

.PHONY: clean
clean:
	$(RM) *.o *.i *~
	for test in $(tests); do $(RM) $(addprefix test_,$$test); done
