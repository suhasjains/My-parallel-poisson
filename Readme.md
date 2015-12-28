#Poisson Equation Solver

Parallel Solver:
* Jacobi solver
* Red-Black Gauss-Seidel with successive over relaxation

Code style:
* Written using OOPS in C++.

Usage:
* To compile parallel solver on local machine machine/cluster: make
* changed - (To compile parallel solver on Cray supercomputer: CC -O3 parallel_poisson.c -o run)
* To run on cluster: qsub pp.sh
* changed - (To run on Cray supercomputer: qsub submit.sh) 
* Type plot_output within an octave terminal to see the surface plot of the output or run the python script plot_output.py

Change log:
* Parallelization on distributed memory system using MPI.
* Domain decomposition is done.
* Communication of buffer cells is working now.
* Global residual added now.
* Speedup tested upto 16 processes on local machine.
* Speedup tested upto 720 processes on Cray XC40.
* Red-Black Gauss-Seidel solver added.

Things to do:
* Parllelize Conjugate gradient solver.
* Non-blocking communication is to be tested.
* update report.pptx
* Mesh importing package to be added

