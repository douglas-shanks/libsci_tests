#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=test_dgemm
#SBATCH --time=0:15:0
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
echo "Starting job $SLURM_JOB_ID at `date`"
echo "nodes"
echo $SLURM_JOB_NODELIST
echo

module list
export OMP_PROC_BIND=close
export OMP_PLACES=cores

rm pdsygvx
ftn sample_pdsyev_call.f -o pdsygvx
EXE=pdsygvx

export OMP_NUM_THREADS=1
srun -c${OMP_NUM_THREADS} ${EXE}

echo "Finished job at `date`"
