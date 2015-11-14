#Poisson Equation Solver

Serial Solvers:
* Jacobi solver 
* Gauss-Seidel with successive over relaxation

Parallel Solver:
* Jacobi solver


Written in C with structures.

* Structures are used to store the field variables and for easy specification of boundary conditions.


Usage:
* To compile use: g++ poisson.c, mpic++ parallel_poisson.c
* Type plot_output within an octave terminal to see the surface plot of the output.
* To run use: mpirun -n 
* To run on a cluster/supercomputer use: qsub pp.sh

Things done:
* Parallelization on distributed memory system using MPI.
* Domain decomposition is done.
* Communication of buffer cells is working now.
* Global residual added now.
* Speedup tested upto 16 processes.

Things to do:
* Parllelize Conjugate gradient solver.
* Scalability test with large cores.
* Non-blocking communication is to be tested.

