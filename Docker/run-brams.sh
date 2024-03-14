#!/bin/bash

MPI_NUM_PROCESS=2
RAMSIN=RAMSIN
BRAMS_EXEC=brams-5.2

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/view/lib/

cp /BRAMS-5.2-modified/build_release/${BRAMS_EXEC} /brams-data/.
cd /brams-data
mkdir tmp
export TMPDIR=./tmp

mpirun --allow-run-as-root -np ${MPI_NUM_PROCESS} ${BRAMS_EXEC} -f ./${RAMSIN}
