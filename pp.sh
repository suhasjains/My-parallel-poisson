#!/bin/bash
#$ -cwd
#$ -N poisson
#$ -pe mpi 16

mpirun -n 16 a.out
