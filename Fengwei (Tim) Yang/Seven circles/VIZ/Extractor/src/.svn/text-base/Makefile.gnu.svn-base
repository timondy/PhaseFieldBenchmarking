LDFLAGS += -L$(hdf5_dir)/lib -L../libs -lhdf5 -lhdf5_hl -lz -lgfortran -lparamesh -lmpi_paramesh -lmodules
ifndef SHMEM
#LDFLAGS += -lmpi_paramesh
endif

# normal LDFLAGS
#LDFLAGS += -lparamesh -lmodules -lgfortran -lhdf5 -lhdf5_hl
#LDFLAGS += -lparamesh -lmodules
# T3E LDFLAGS
#LDFLAGS += -lparamesh

ifdef SHMEM
# SGI LDFLAGS
#LDFLAGS += -lsma
# Compaq LDFLAGS
#LDFLAGS += -lshmem
endif

#FFLAGS+=-g

LDFLAGS += $(ADD_LIB)
LDFLAGS += -lparamesh

# List all application source files
sources := 			\
amr_multigrid_user_modules.f90 	\
amr_1blk_bcset.f90		\
amr_multigrid_stat_03.f90 	\
amr_multigrid_control_03.f90  	\
amr_multigrid_user_edits_03.f90 \
ceg_mesh_partitioning.f90	\
non_linear_mg.f90 
csources := profilingLinker.c
#amr_test_refinement.f90 amr_advance_soln.f90 amr_timestep.f90 
%.o:%.f90
	$(FC) -c $(FFLAGS) $<
%.o:%.c
	$(CC) -c $(FFLAGS) $<
objects := $(sources:.f90=.o)
exobjects := $(csources:.c=.o)


# Identify the main program
main := non_linear_mg.f90
mainobject := $(.f90=.o)
exobject := $(.c=.o)

# Set the executable name
CMD := ../Extractor




# compiles the program.
$(CMD): $(mainobject) $(objects) $(exobjects)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

#$(CMD): $(mainobject) $(objects) 
#	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)


.PHONY: clean
clean:
	$(RM) $(CMD) *.o *.i *.mod *~
