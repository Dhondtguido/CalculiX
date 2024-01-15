#!/bin/bash

if ! [[ -d build ]]; then
    mkdir build
fi
cd build

INSTALLPATH="/home/guido/PaStiX/parsec_i8"

umask 022

# fixes
sed -i '/-1 == cpu/i return cpu;' parsec/bindthread.c

cmake \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${INSTALLPATH} \
    -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-10.2 \
    -DCUDA_DIR=/usr/local/cuda-10.2 \
    -DCUDA_USE_STATIC_CUDA_RUNTIME=ON \
    -DCMAKE_CUDA_HOST_COMPILER=gcc \
    -DPARSEC_GPU_WITH_CUDA=ON \
    -DHWLOC_DIR=/home/guido/PaStiX/hwloc_i8 \
    ..

make -j8

rm -rf ${INSTALLPATH}
make install
