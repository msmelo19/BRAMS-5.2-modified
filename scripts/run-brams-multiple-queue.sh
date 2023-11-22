#!/bin/bash
#SBATCH --job-name=queue-controller
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=1              # Number of MPI ranks

###############################################################
#
# USAGE: sbatch -p <nome_da_fila> run-brams-multiple-queue.sh <primeira opção de fila> <segunda opção de fila> <script de submissão> <executável do brams> <caminho do RAMSIN em modo INITIAL>
#
###############################################################

PRIMARY_QUEUE=${1}
SECONDARY_QUEUE=${2}
SUBMISSION_SCRIPT=${3}
BRAMS_BIN=${4}
RAMSIN_INITIAL_PATH=${5}

JOB_ID=$(sbatch -p ${PRIMARY_QUEUE} ${SUBMISSION_SCRIPT} ${BRAMS_BIN} ${RAMSIN_INITIAL_PATH} | awk '{print $NF}')

python3 ./queue-changer.py ${JOB_ID} ${SECONDARY_QUEUE} 