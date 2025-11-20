#!/bin/bash

module load PrgEnv-cray
#module load cpe/22.04

echo "Build benchmark_pdsyevd"
ftn -h omp -o pdsyevd_bench benchmark_pdsyevd.f90
echo "Build benchmark_pzheevd"
ftn -h omp -o pzheevd_bench benchmark_pzheevd.f90

mv *_bench ../run
