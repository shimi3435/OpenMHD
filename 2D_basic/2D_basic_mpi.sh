#!/bin/bash
#------- qsub option -----------
#PBS -P NIFS25KISC015
#PBS -N 2D_basic_mpi_PS
#PBS -q B_S
#PBS -l select=4:ncpus=24:mpiprocs=1:ompthreads=24
#PBS -l walltime=1:00
#PBS -j oe

#------- Program execution -----------
module load openmpi/5.0.7/rocm6.3.3

cd ${PBS_O_WORKDIR}

mpirun --display-map --map-by numa -x UCX_MAX_RNDV_RAILS=4 ./ap.out
