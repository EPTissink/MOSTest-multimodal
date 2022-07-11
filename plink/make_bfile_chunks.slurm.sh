#!/bin/bash

#SBATCH --job-name=make_bfile_chunks
#SBATCH --account=p33_norment
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=1
#SBATCH --array=0-906

source /cluster/bin/jobsetup
module purge

source ${HOME}/.profile
source ${HOME}/cluster/bin/activate

module list
which plink2
which python

MAX_PADDING=3
CHUNK_ID=$(printf %0${MAX_PADDING}d ${SLURM_ARRAY_TASK_ID})
DIR=/cluster/projects/p33/users/alexeas/elleke/data
SNP_CHUNKS_DIR=${DIR}/snp_chunks
BFILE_CHUNKS_DIR=${DIR}/bfile_chunks
SNP_CHUNK_FILE=${SNP_CHUNKS_DIR}/UKB34886_T1_QCed_230222.chunk${CHUNK_ID}.txt

BFILE=${DIR}/UKB34886_T1_QCed_230222
BFILE_CHUNK=${BFILE_CHUNKS_DIR}/UKB34886_T1_QCed_230222.chunk${CHUNK_ID}
BFILE_CHUNK_PERM=${BFILE_CHUNKS_DIR}/UKB34886_T1_QCed_230222.chunk${CHUNK_ID}.perm

plink2 --bfile ${BFILE} --extract ${SNP_CHUNK_FILE} --out ${BFILE_CHUNK} --make-bed --threads 1 --memory 8192

PY=/cluster/projects/p33/users/alexeas/elleke/src/permute_bed.py
python ${PY} --bfile ${BFILE_CHUNK} --out ${BFILE_CHUNK_PERM}

