#!/bin/bash
#$ -cwd
#$ -N poisson
#$ -pe mpi 4

time mpirun -n 4 run >log
