#!/bin/sh
for x in *.f *.c
do
    if [ "$x" = "CalculiX.c" ]; then
	echo $x "is kept"
	continue
    fi
    if [ "$x" = "gauss.f" ]; then
	echo $x "is kept"
	continue
    fi
    if [ "$x" = "xlocal.f" ]; then
	echo $x "is kept"
	continue
    fi
    
    if grep -q $x Makefile.inc
    then
	echo $x "is kept"
    else
	echo $x "is deleted"
	rm -f $x
    fi
done
