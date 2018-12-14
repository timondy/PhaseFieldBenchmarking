# Uintah PhaseField Component

Uintah is distributed under the MIT License at
* http://uintah-build.sci.utah.edu/trac/wiki/License

More information about Uintah are available at
* http://uintah.utah.edu/

This repository is not updated as it is meant to provide a snapshot of the
software used for computing benchmark results.

The PhaseField component implementing these benchmark simulations is implemented
upon revision 58893 from the official SVN repository at
* https://gforge.sci.utah.edu/svn/uintah/trunk


## REQUIREMENTS

* make
  tested with gnu (3.82, 4.1)

* F77 and C++11 compilers:
  tested with gnu (4.8.5, 6.3.0) and intel (18.0.2) compilers

* MPI
  tested with openmpi (1.10.7, 2.0.2) and intelmpi (2018.2.199)

* libxml2
  tested (2.9.1, 2.9.4)

* zlib
  tested (1.2.7, 1.2.8)

* blas/lapack
  tested with netlib blas/lapack (3.4.2, 3.7.0) and intel mkl (2018.2)

* no solver is required, neither hypre nor petsc are required.


## UINTAH CONFIGURATION

Minimal configuration scripts for building Uintah with the PhaseField component
are available in ./config-examples for the following systems:

* debian: Debian stretch (stable) x86_64 GNU/Linux
* centos: CentOS 7 x86_64 GNU/Linux
* arc3:   Advanced Research Computing Node 3 (http://www.arc.leeds.ac.uk)


## BUILD INSTRUCTIONS

Only the sus target needs to be compiled:

1.  Create build directory
    ```
    $ mkdir build
    ```

2.  Copy (and edit) best-suited configuration script to build directory
    ```
    $ cp config-examples/configure.#####-###-##### build/configure.sh
    $ <edit> build/configure.sh
    ```

3.  Change directory and configure
    ```
    $ cd build
    $ ./configure.sh
    ```

4.  Compile
    ```
    $ make [-j#] sus
    ```

* no install method is provided with Uintah


## HOW TO RUN BENCHMARK TESTS

Input files for the PhaseField Benchmark applications are located in
[src/StandAlone/inputs/PhaseField/](./src/StandAlone/inputs/PhaseField/)benchmark##

1.  Create output directory and change directory
    ```
    $ mkdir out
    $ cd out
    ```

2.  Run sus with mpi
    ```
    $ mpirun -np <num_proc> <build_path>/sus <path_to_input_file>
    ```

A new folder with extension .uda will be crated in the currend working directory
Output can be visualized with Hypre
(https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)

*   Interrupted simulations may be restarted with
    ```
    $ mpirun -np <num_proc> <build_path>/sus -restart <path_to_output_file>
    ```

## HOW SCHEDULE BENCHMARK TESTS

Examples of sge job scripts are located in [src/scripts/PhaseField/benchmark/sge/](./src/scripts/PhaseField/benchmark/sge/)

1.  Jobs can be submitted from the output directory with
    ```
    $ qsub <path_to_sge_file>
    ```


## HOW TO POSTPROCESS BENCHMARK TESTS

Bash scripts are located in [src/scripts/PhaseField/benchmark/bash](./src/scripts/PhaseField/benchmark/bash/)
and Matlab scripts are in [src/scripts/PhaseField/benchmark/matlab/](./src/scripts/PhaseField/benchmark/matlab/)

1.  Reorganize output file.
    From the output directory run to create a tree of directories at ./dat and
    copy relevant files from origina uda outputs.
    ```
    $ <bash_script_path>/organize_output
    ```

2.  Join together output files.
    Deepest directories within ./dat will now contain one or more files ending
    with progressive numbers corresponding to the original uda suffixes
    identifying output from the first run and the eventual following restarts
    For example the content of one of these could be:

    - energy_000.dat : contains the energy profile from the initial run
    - energy_001.dat : contains the energy profile from the first restart
    - energy_002.dat : contains the energy profile from the second restart
    - u0_000.dat :     contains the value of the solution at the center from the initial run
    - u0_001.dat :     contains the value of the solution at the center from the first restart
    - u0_002.dat :     contains the value of the solution at the center from the second restart

    a.  For every restart file check the first timestep
    ```
    $ head energy_001.dat energy_002.dat u0_001.dat u0_002.dat
    ```

    b.  Delete all lines after that timestep in the previous file
    ```
    $ <edit> energy_000.dat energy_001.dat u0_000.dat u0_001.dat
    ```

    c.  Join all output together
    ```
    $ cat energy_000.dat energy_001.dat energy_002.dat > energy.dat
    $ cat u0_000.dat u0_001.dat u0_002.dat > energy.dat
    ```

3.  Run Matlab scripts
    This will generate a csv file with computed benchmark values for each
    simulation found in the path

    a.  Start Matlab and add scripts path and cd into output directory
    ```
    $ addpath <matlab_scripts_path>
    $ cd <output_path>
    ```

    b.  Call the appropriate postprocess function
    ```
    $ postprocess##(<path>)
    ```


## PHASE FIELD COMPONENT DOCUMENTATION

The [documentation](https://timondy.github.io/PhaseFieldBenchmarking/Uintah/doc/Components/PhaseField/html/index.html)
([fork](https://jonmatteochurch.github.io/PhaseFieldBenchmarking/Uintah/doc/Components/PhaseField/html/index.html))
of the PhaseField component can be generated with Doxygen

```
$ cd doc/Components/PhaseField/
$ doxygen
```
