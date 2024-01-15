#!/bin/bash

umask 022
cd src
ln -sf Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
sed -i '/CFLAGS/ s/$/ -DINTSIZE64/' Makefile.inc

make -j8 scotch
make -j8 esmumps
make prefix=/home/guido/PaStiX/scotch_i8 install
