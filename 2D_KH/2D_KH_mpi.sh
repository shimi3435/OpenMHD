#!/bin/bash
#------- qsub option -----------
#PBS -P NIFS25KISC015
#PBS -N 2D_KH_PS
#PBS -q B_S
#PBS -l select=4:ncpus=24:mpiprocs=1:ompthreads=24
#PBS -l walltime=5:00
#PBS -j oe

#------- Program execution -----------
module load openmpi/5.0.7/rocm6.3.3

export OMP_NUM_THREADS=24

cd ${PBS_O_WORKDIR}

mpirun --display-map --bind-to NUMA -x UCX_MAX_RNDV_RAILS=4 ./ap.out
