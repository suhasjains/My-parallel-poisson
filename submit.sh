#!/bin/sh
#This job should be redirected to small queue
#PBS -N Poisson144
#PBS -l select=3:ncpus=24       
#PBS -l walltime=00:10:00               
#PBS -l place=scatter
#PBS -e error.log               
#PBS -o output.log              
#PBS -S /bin/sh -V              
#PBS -S /bin/sh@sdb -V

./opt/modules/default/init/sh
cd /mnt/lustre/mec3/mecgrp/Suhas/poisson72/
aprun -j1 -n 72 -N 24 run

