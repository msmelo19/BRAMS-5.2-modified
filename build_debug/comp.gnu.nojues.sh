#!/bin/sh

DIR=$PWD

#!- configure for 
./configure --prefix=${PWD} \
 --program-prefix=BRAMS --with-fpcomp=mpif90 --with-fcomp=gfortran \
 --with-cpcomp=mpicc --with-ccomp=gcc --with-chem=RELACS_TUV --with-aer=SIMPLE \
 --with-hdf5lib=/home/rpsouto/usr/local/spack/git/v0.18.1/opt/spack/linux-ubuntu20.04-skylake/gcc-8.4.0/hdf5-1.12.2-wbqymagzc7m2cx7r6hf7tlzefsadzspz
