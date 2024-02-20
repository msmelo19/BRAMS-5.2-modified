#!/bin/bash
#SBATCH --job-name=brams-5.2
#SBATCH --nodes=2              # Number of nodes
#SBATCH --ntasks=96             # Number of MPI ranks
#SBATCH --ntasks-per-node=48    # Number of MPI ranks per node

set -e
ulimit -s unlimited

BRAMS_BIN=${1}
RAMSIN_INITIAL_PATH=${2}
LOG_DIR=result/history/2

source ~/.spack/v0.18.1/env_spack.sh
spack load hdf5@1.12.2

export TMPDIR=./tmp
#export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi2.so

if [ -z ${SLURM_RESTART_COUNT} ] || [ ${SLURM_RESTART_COUNT} -eq 0 ]; then
  echo "INITIAL"
  date > ${LOG_DIR}/date-ini.txt
  time mpirun --map-by ppr:48:node ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH} 2>&1 | tee ${LOG_DIR}/brams-initial-${SLURM_NTASKS}-${SLURM_JOB_ID}.out
else
  echo "HISTORY"
  ./create-ramsin-history-mode.sh ${RAMSIN_INITIAL_PATH}
  date > ${LOG_DIR}/date-his-${SLURM_RESTART_COUNT}.txt
  time mpirun --map-by ppr:48:node ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH}-history 2>&1 | tee brams-history-${SLURM_NTASKS}-${SLURM_JOB_ID}-${SLURM_RESTART_COUNT}.out
fi

date > ${LOG_DIR}/date-fim
