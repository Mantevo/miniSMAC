#!/bin/bash

#SBATCH --nodes=1                    # Number of nodes - all cores per node are allocated to the job
#SBATCH --time=60:10:00               # Wall clock time (HH:MM:SS) - once the job exceeds this time, the job will be terminated (default is 5 minutes)
#SBATCH --account=FY140087           # WC ID
#SBATCH --job-name=smac2d_1           # Name of job

nodes=$SLURM_JOB_NUM_NODES           # Number of nodes - the number of nodes you have requested (for a list of SLURM environment variables see "man sbatch")
cores=1                              # Number MPI processes to run on each node (a.k.a. PPN)
                                     # TLCC2 has 16 cores per node

mpiexec --bind-to-core --npernode $cores --n $(($cores*$nodes)) /gscratch1/dwbarne/ms2d_2d_only_partitioned_grids/src/smac2d > temp1_7kx7k.out 

# Note 1: This will start ($nodes * $cores) total MPI processes using $cores per node.  
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
