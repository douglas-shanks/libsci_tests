#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=pdsyevd
#SBATCH --time=0:40:0
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
module load PrgEnv-gnu

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
module list

ulimit -c unlimited
ulimit -s unlimited

EXE=pdsyevd_bench
#ldd $EXE

#export PE_SCI_EXT_LIBPATH="-L/work/y02/y02/dshanks/tools/pdsyev_bench/elpa_gcc/2023.05.001/openmp/lib/"
#export PE_SCI_EXT_LIBNAME="-Wl,--undefined=elpa_get_communicators\ -lelpa_openmp"
export LIBSCI_USE_ELPA=1

module swap craype-network-ofi craype-network-ucx 
module swap cray-mpich cray-mpich-ucx
module load libfabric

echo " ================ UCX RUN ================ "
time srun --hint=nomultithread  --distribution=block:block ./${EXE}

module restore 
module load PrgEnv-gnu
module list

export MPICH_OFI_STARTUP_CONNECT=0
echo " ================ OFI RUN ================ "

echo "MPICH_OFI_STARTUP_CONNECT = " $MPICH_OFI_STARTUP_CONNECT
time srun --hint=nomultithread  --distribution=block:block ./${EXE}


#export MPICH_OFI_STARTUP_CONNECT=1
#echo " ================ OFI RUN ================ "

#echo "MPICH_OFI_STARTUP_CONNECT = " $MPICH_OFI_STARTUP_CONNECT
#time srun --hint=nomultithread  --distribution=block:block ./${EXE}
