#!/bin/bash

module load PrgEnv-gnu
#cray
#module load cpe/22.04

export PE_SCI_EXT_LIBPATH="-L/work/y02/y02/dshanks/tools/pdsyev_bench/elpa_2016/2016.05.004/openmp/lib/"
export PE_SCI_EXT_LIBNAME="-Wl,--undefined=elpa_get_communicators\ -lelpa_openmp"

echo "Build benchmark_pdsyevd"
ftn -g -fopenmp -fallow-argument-mismatch -o pdsyevd_bench benchmark_pdsyevd.f90
echo "Build benchmark_pzheevd"
ftn -g -fopenmp -fallow-argument-mismatch -o pzheevd_bench benchmark_pzheevd.f90

mv *_bench ../run
