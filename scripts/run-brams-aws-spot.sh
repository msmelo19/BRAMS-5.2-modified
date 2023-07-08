#!/bin/bash
#SBATCH --job-name=brams-5.2
#SBATCH --nodes=2              # Number of nodes
#SBATCH --ntasks=2             # Number of MPI ranks
#SBATCH --ntasks-per-node=1    # Number of MPI ranks per node

set -e
ulimit -s unlimited

BRAMS_BIN=${1}
RAMSIN_INITIAL_PATH=${2}

source ~/.spack/v0.18.1/env_spack.sh
spack load hdf5@1.12.2

export TMPDIR=./tmp
#export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi2.so

if [ -z ${SLURM_RESTART_COUNT} ] || [ ${SLURM_RESTART_COUNT} -eq 0 ]; then
  echo "INITIAL"
  mpirun -n ${SLURM_NTASKS} ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH} 2>&1 | tee brams-initial-${SLURM_NTASKS}-${SLURM_JOB_ID}.out &
else
  echo "HISTORY"
  ./create-ramsin-history-mode.sh ${RAMSIN_INITIAL_PATH}
  mpirun -n ${SLURM_NTASKS} ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH}-history 2>&1 | tee brams-history-${SLURM_NTASKS}-${SLURM_JOB_ID}.out
fi
