#!/bin/sh

DIR=$PWD

#!- configure for 
./configure --prefix=${PWD} \
 --program-prefix=BRAMS --enable-jules --with-fpcomp=mpif90 --with-fcomp=gfortran \
 --with-cpcomp=mpicc --with-ccomp=gcc --with-chem=RELACS_TUV --with-aer=SIMPLE

