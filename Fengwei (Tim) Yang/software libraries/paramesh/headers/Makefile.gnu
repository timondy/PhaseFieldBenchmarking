PHLIB_RPATH := $(paramesh_dir)/headers

sources := \
paramesh_dimensions.F90 \
paramesh_interfaces.F90 \
physicaldata.F90 \
prolong_arrays.F90 \
timings.F90 \
tree.F90 \
amr_mg_common.F90 \
workspace.F90 \
io.F90 \
constants.F90 \
paramesh_comm_data.F90

ifndef SHMEM
sources += mpi_morton.F90 paramesh_mpi_interfaces.F90
endif

.NOTPARALLEL %.lo:%.F90
	libtool  --mode=compile --tag=FC $(FC) -c $(FFLAGS) $<

objects := $(sources:.F90=.lo)

libmodules.a: $(objects)
	libtool  --mode=link --tag=FC $(FC) $(FFLAGS) -o libmodules.la $^ $(LDFLAGS)  -rpath $(PHLIB_RPATH)
	if [ -f $(PHLIB_RPATH)/.libs/libmodules.a ]; then ln -sf $(PHLIB_RPATH)/.libs/libmodules.a $(paramesh_dir)/libs/ ; fi
	if [ -f $(PHLIB_RPATH)/.libs/libmodules.so ]; then ln -sf $(PHLIB_RPATH)/.libs/libmodules.so $(paramesh_dir)/libs/ ; fi

#	$(AR) $(ARFLAGS) $@ $^

ifdef MY_CPP
GNUmakefile.include: $(sources)
	find . -name \*.F90 | xargs $(MY_CPP) > $@
include GNUmakefile.include
else
$(objects): $(wildcard *.fh)
endif

.PHONY: clean
clean:
	$(RM) libmodules.a *.o *.lo *.mod *.d *.la *~ GNUmakefile.include .libs
