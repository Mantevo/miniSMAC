#!/bin/bash

#SBATCH --nodes=16                    # Number of nodes - all cores per node are allocated to the job
#SBATCH --time=10:10:00               # Wall clock time (HH:MM:SS) - once the job exceeds this time, the job will be terminated (default is 5 minutes)
#SBATCH --account=FYXXXXXXX           # WC ID
#SBATCH --job-name=smac2d_mpi_16_omp_8      # Name of job

nodes=$SLURM_JOB_NUM_NODES           # Number of nodes - the number of nodes you have requested (for a list of SLURM environment variables see "man sbatch")
mpi_ranks_per_node=1                # Number MPI processes to run on each node (a.k.a. PPN)
                                     # TLCC2 has 16 cores/node, 8 cores/socket

export OMP_NUM_THREADS=8             # threads per MPI rank
mpiexec --bind-to-socket --bysocket --npernode $mpi_ranks_per_node --n $(($mpi_ranks_per_node*$nodes)) /gscratch1/dwbarne/ms2d_2d_only_partitioned_grids/src/smac2d > temp_mpi_16_omp_8.out 
# and the following one line works for threads as well:
#mpiexec --bind-to-socket --bysockett --npernode $mpi_ranks_per_node --n $(($mpi_ranks_per_node*$nodes)) env OMP_NUM_THREADS=8 /gscratch1/dwbarne/ms2d_2d_only_partitioned_grids/src/smac2d > temp1x8_5kx5k.out 

# Note 1: This will start ($mpi_ranks_per_node * $cores) total MPI processes using $cores per node.  
#           If you want some other number of processes, add "-np N" after the mpiexec, where N is the total you want.
#           Example:  mpiexec -np 24  ......(for a 2 node job, this will load 16 processes on the first node and 8 processes on the second node)
#           If you want a specific number of process to run on each node, (thus increasing the effective memory per core), use the --npernode option.
#           Example: mpiexec -np 24 --npernode 12  ......(for a 2 node job, this will load 12 processes on each node)
# Note 2: The option "--bind-to-core" binds each MPI process to a core.
#           It locks processes to specific cores and locks the memory they allocate
#           to the nearest memory controller, improving efficiency and repeatability.

# To submit your job, do:
# sbatch <script_name>
#
#The slurm output file will by default, be written into the directory you submitted your job from  (slurm-JOBID.out)
