#Poisson Equation Solver

Solvers:
* Jacobi solver 
* Gauss-Seidel with successive over relaxation



Written in C with structures.

* Structures are used to store the field variables and for easy specification of boundary conditions.


Usage:
* To compile use: g++ poisson.c
* Type plot within an octave terminal to see the surface plot of the output.


Under progress:
* Parallelization on distributed memory system using MPI.
* Domain decomposition is done.


Finally scalability of the parallelization is to be tested.

