#!/bin/bash
#
module restore
#module load PrgEnv-cray
module load PrgEnv-gnu
#module load cray-fftw
module load cray-python
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

#export FCFLAGS="-dC -homp -eZ -hipa1"
export FCFLAGS="-fallow-argument-mismatch"

ELPA_VERSION=2016.05.004 #2023.05.001
ELPA_LABEL=elpa
ELPA_NAME=${ELPA_LABEL}-${ELPA_VERSION}
ELPA_ROOT=${PWD}/${ELPA_LABEL}

rm -rf ${ELPA_ROOT}
mkdir -p ${ELPA_ROOT}
cd ${ELPA_ROOT}

wget -q https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${ELPA_VERSION}/${ELPA_NAME}.tar.gz
tar zxf ${ELPA_NAME}.tar.gz
rm ${ELPA_NAME}.tar.gz

cd ${ELPA_NAME}
#mkdir build-serial
#cd build-serial

#export LIBS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"

#CC=cc CXX=CC FC=ftn LDFLAGS=-dynamic ../configure       \
#  --enable-openmp=no --enable-shared=no \
#  --disable-avx512 --disable-detect-mpi-launcher \
#  --without-threading-support-check-during-build \
#  --prefix=${ELPA_ROOT}/${ELPA_VERSION}/serial

#make -j 8
#make -j 8 install
#make -j 8 clean

cd ${ELPA_ROOT}/${ELPA_NAME}
mkdir build-openmp
cd build-openmp

#export LIBS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"

CC=cc CXX=CC FC=ftn LDFLAGS=-dynamic ../configure \
  --enable-openmp=yes --enable-shared=no --enable-allow-thread-limiting \
  --disable-avx512 --disable-detect-mpi-launcher \
  --without-threading-support-check-during-build \
  --prefix=${ELPA_ROOT}/${ELPA_VERSION}/openmp

make -j 8
make -j 8 install
make -j 8 clean
