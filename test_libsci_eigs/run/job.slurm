#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=pdsyevd
#SBATCH --time=00:5:0
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=y02-coe
#SBATCH --partition=standard
#SBATCH --qos=standard

# Turn of the stop message from the Cray compiler
export NO_STOP_MESSAGE=1


export OMP_NUM_THREADS=1
module load PrgEnv-cray

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
module list

EXE=pdsyevd_bench
#ldd $EXE


module swap craype-network-ofi craype-network-ucx 
module swap cray-mpich cray-mpich-ucx
module load libfabric

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

echo " ================ UCX RUN ================ "
time srun --hint=nomultithread  --distribution=block:block ./${EXE}

module restore 
module list
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

export MPICH_OFI_STARTUP_CONNECT=0
echo " ================ OFI RUN ================ "

echo "MPICH_OFI_STARTUP_CONNECT = " $MPICH_OFI_STARTUP_CONNECT
time srun --hint=nomultithread  --distribution=block:block ./${EXE}


#export MPICH_OFI_STARTUP_CONNECT=1
#echo " ================ OFI RUN ================ "

#echo "MPICH_OFI_STARTUP_CONNECT = " $MPICH_OFI_STARTUP_CONNECT
#time srun --hint=nomultithread  --distribution=block:block ./${EXE}
