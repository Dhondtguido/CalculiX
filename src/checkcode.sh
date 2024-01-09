#!/bin/sh
for x in *.f *.c
	 
do
    if grep -qi proprietary $x
    then
	echo "proprietary occurs in" $x
	exit
    fi
    if grep -qi pratt $x
    then
	echo "pratt occurs in" $x
	exit
    fi
    if grep -qi whitney $x
    then
	echo "whitney occurs in" $x
	exit
    fi
    if grep -qi mtu $x
    then
	echo "mtu occurs in" $x
	exit
    fi
    if grep -qi larsson $x
    then
	echo "larsson occurs in" $x
	exit
    fi
done
