#!/bin/bash

umask 022
./configure --prefix=/home/guido/PaStiX/hwloc_i8 CC=gcc CXX=g++
make -j8
make install
