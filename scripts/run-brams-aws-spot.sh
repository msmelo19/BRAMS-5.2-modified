#!/bin/bash
#SBATCH --job-name=brams-5.2
#SBATCH --nodes=4               # Number of nodes
#SBATCH --ntasks=192             # Number of MPI ranks
#SBATCH --ntasks-per-node=48    # Number of MPI ranks per node

#######################################################################################################
#
# USAGE: sbatch -p <nome_da_fila> run-brams-aws-spot.sh <executável do brams> <caminho do RAMSIN em modo INITIAL>
#
#######################################################################################################

set -e
ulimit -s unlimited

BRAMS_BIN=${1}
RAMSIN_INITIAL_PATH=${2}
RESULT_DIR=./result/ondemand/${SLURM_NNODES}-nos

source ~/.spack/v0.18.1/env_spack.sh
spack load hdf5@1.12.2

export TMPDIR=./tmp
mkdir -p ${RESULT_DIR}

for i in $(seq 1 3); do
	if [ -z ${SLURM_RESTART_COUNT} ] || [ ${SLURM_RESTART_COUNT} -eq 0 ]; then
	  echo "INITIAL"
	  mpirun --map-by ppr:${SLURM_NTASKS_PER_NODE}:node ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH} 2>&1 | tee ${RESULT_DIR}/brams-initial-${SLURM_NTASKS}-${SLURM_JOB_ID}-${i}.out
	else
	  echo "HISTORY"
	  ./create-ramsin-history-mode.sh ${RAMSIN_INITIAL_PATH}
	  mpirun --map-by ppr:${SLURM_NTASKS_PER_NODE}:node ${BRAMS_BIN} -f ${RAMSIN_INITIAL_PATH}-history 2>&1 | tee brams-history-${SLURM_NTASKS}-${SLURM_JOB_ID}-${i}-${SLURM_RESTART_COUNT}.out
	fi
done

mv slurm-${SLURM_JOB_ID}.out ${RESULT_DIR}
