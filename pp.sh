#!/bin/bash
#$ -cwd
#$ -N poisson
#$ -pe mpi 4

mpirun -n 4 run > log
