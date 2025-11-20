#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=test_dgemm
#SBATCH --time=0:15:0
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
echo "Starting job $SLURM_JOB_ID at `date`"
echo "nodes"
echo $SLURM_JOB_NODELIST
echo

module list
export OMP_PROC_BIND=close
export OMP_PLACES=cores

rm dgemm
ftn dgemm_test.f90 -o dgemm
EXE=dgemm

for i in 1 2 4 8 16 32 84 128
do
    echo $i
    export OMP_NUM_THREADS=$i
    srun -n1 -c${OMP_NUM_THREADS} ${EXE}
done

echo "Finished job at `date`"
