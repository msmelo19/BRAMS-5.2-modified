#!/bin/sh

DIR=$PWD

#!- configure for 
./configure --prefix=${PWD} \
 --program-prefix=BRAMS --with-fpcomp=mpif90 --with-fcomp=gfortran \
 --with-cpcomp=mpicc --with-ccomp=gcc --with-chem=RELACS_TUV --with-aer=SIMPLE \
 --with-hdf5lib=/home/u111227/spack/v0.18.1/opt/spack/linux-ubuntu20.04-skylake_avx512/gcc-11.2.0/hdf5-1.12.2-7kuvbuflp22jmpcao7nzyyxr6ma2c72g
