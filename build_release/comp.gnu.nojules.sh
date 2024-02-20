#!/bin/sh

DIR=$PWD

#!- configure for 
./configure --prefix=${PWD} \
 --program-prefix=BRAMS --with-fpcomp=mpif90 --with-fcomp=gfortran \
 --with-cpcomp=mpicc --with-ccomp=gcc --with-chem=RELACS_TUV --with-aer=SIMPLE \
 --with-hdf5lib=$(spack location -i hdf5)

