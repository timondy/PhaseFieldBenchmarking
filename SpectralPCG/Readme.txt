This directory contains the MATLAB code to generate the benchmark results in the paper, "High Accuracy Benchmark Problems for Allen-Cahn and Cahn-Hilliard Dynamics". 

Specifically, this code generates the results for the second benchmark in that paper (2D Cahn Hilliard seven circles initial conditions) with the spectral spatial discretization, DIRK2 variable time stepping approach.

Results that match the reported benchmarks in that paper for this method can be obtained with the parameters:

N = 192, sigma = 1e-4 (epsilon=0.1) 
N = 256, sigma = 1e-5 (epsilon=0.05) 
N = 384, sigma = 1e-5 (epsilon=0.0.025) 

These parameters can be found at the beginning of the main program, ch2d.m, currently set to the first case above. The code was developed for MATLAB version 2108a. 

Please direct comments and questions to:
Brian Wetton, wetton@math.ubc.ca 