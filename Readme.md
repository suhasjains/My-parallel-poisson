#Poisson Equation Solver

Serial Solvers:
* Jacobi solver 
* Gauss-Seidel with successive over relaxation

Parallel Solver:
* Jacobi solver

Code style:
* Written in C with structures.

Usage:
* To compile serial solver: g++ poisson.c -o run, 
* To compile parallel solver on local machine machine/cluster: mpic++ parallel_poisson.c -o run
* To compile parallel solver on Cray supercomputer: CC -O3 parallel_poisson.c -o run
* To run on local machine/cluster: qsub pp.sh
* To run on Cray supercomputer: qsub submit.sh 
* Type plot_output within an octave terminal to see the surface plot of the output or run the python script plot_output.py

Change log:
* Parallelization on distributed memory system using MPI.
* Domain decomposition is done.
* Communication of buffer cells is working now.
* Global residual added now.
* Speedup tested upto 16 processes on local machine.
* Speedup tested upto 720 processes on Cray XC40.

Things to do:
* Parllelize Conjugate gradient solver.
* Non-blocking communication is to be tested.

