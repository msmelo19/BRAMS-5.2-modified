#!/bin/sh
cd ./ 
#!- configure for Tupa 

./configure --program-prefix=BRAMS --enable-jules --with-fpcomp=ftn --with-fcomp=pgf90 --with-cpcomp=ftn --with-ccomp=pgcc --with-chem=RELACS_TUV --with-aer=SIMPLE
make
cp  brams-5.2 /Users/saulofreitas/work/testcase_brams5.1/BRAMS-5.2
exit

#!- configure for pgi 
#./configure --program-prefix=BRAMS --enable-jules --with-fpcomp=/opt/hdf5-pgi12.6/bin/h5pfc --with-fcomp=/opt/pgi/osx86-64/12.6/bin/pgf90 \
#--with-cpcomp=/opt/hdf5-pgi12.6/bin/h5pcc --with-ccomp=/opt/pgi/osx86-64/12.6/bin/pgcc --with-chem=RELACS_TUV --with-aer=SIMPLE
make

cp  brams-5.1 /Users/saulofreitas/work/testcase_brams5.1/BRAMS-5.1.1-relacs-tuv
exit
