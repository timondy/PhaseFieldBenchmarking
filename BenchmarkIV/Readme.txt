This directory contains the data for the benchmark IV from the paper, "High Accuracy Benchmark Problems for Allen-Cahn and Cahn-Hilliard Dynamics". 

Specifically, this contains the MATLAB data files run8.mat and run9.mat containing the arrays E and t, with E(t) the computed energy for the benchmark. These were computed with the spectral spatial discretization, adaptive time stepping method with grid N=512, and local error tolerance sigma=5e-6 and 2.5e-6 for run8 and run9, respectively. 

Also included is a MATLAB program L1diff.m that evaluates the benchmark error between the two runs. Substitute data from your own computation for the run8.mat to evaluate the benchmark for your computation. 

Please direct comments and questions to:
Brian Wetton, wetton@math.ubc.ca 