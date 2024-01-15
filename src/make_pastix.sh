#!/bin/bash
if ! [[ -d build ]]; then
    mkdir build
fi
cd build


cmake   \
        -DBLAS_DIR=/home/guido/OpenBLAS_i8 \
	-DHWLOC_DIR=/home/guido/PaStiX/hwloc_i8 \
	-DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-10.2 \
	-DCMAKE_INSTALL_PREFIX=/home/guido/PaStiX/pastix_i8 \
	-DCMAKE_BUILD_TYPE=Release \
	-DPASTIX_WITH_PARSEC=ON \
	-DPARSEC_DIR=/home/guido/PaStiX/parsec_i8 \
	-DSCOTCH_DIR=/home/guido/PaStiX/scotch_i8 \
	-DPASTIX_WITH_CUDA=ON \
	-DCUDA_DIR=/usr/local/cuda-10.2 \
	-DPASTIX_ORDERING_SCOTCH=ON \
	-DCMAKE_C_COMPILER=gcc \
	-DCMAKE_CXX_COMPILER=g++ \
	-DCMAKE_Fortran_COMPILER=gfortran \
	-DCMAKE_C_FLAGS="-fopenmp" \
        ..

make -j8
make install
