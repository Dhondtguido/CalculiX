#!/bin/sh

# run a sensitivity analysis for the beam

../src/CalculiX opt1

# rename the file controlling the mesh modification

mv opt1.equ opt2.inp

# calculating the mesh modification

../src/CalculiX opt2

# adding the mesh modification to the coordinates

cgx -b opt2.fbl

# renaming NNall into Nall

perl -pi -e 's/NNall/Nall/g' Nall.msh

# running the sensitivity analysis for the once-optimized-beam

rm -f opt3.inc
mv Nall.msh opt3.inc
../src/CalculiX opt3
