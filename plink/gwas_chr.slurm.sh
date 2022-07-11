#!/bin/bash
#SBATCH --job-name=uinv_gwas
#SBATCH --account=p33_norment
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=1-22

set -o errexit
source ../settings.sh

PHENO=$1 # a file in PHENO_DIR
COVAR=$2 # a file in COVAR_DIR
RES=$3 # prefix of a file to create in OUT_DIR

ROOT_DIR=/cluster/projects/p33/users/elizabept/multimodal/discovery
PHENO_DIR=${ROOT_DIR}/phenos/results/QCed_norm
COVAR_DIR=${ROOT_DIR}/phenos/data
#OUT_DIR=${OUT}
PLINK2=${SOFTWARE}/plink2
THREADS=32
MEMORY=$((8192 * ${THREADS}))

BFILE=${ROOT_DIR}/genos/results/UKB34886_T1_QCed_230222
PHENO=${PHENO_DIR}/${PHENO}
COVAR=${COVAR_DIR}/${COVAR}
RES=${OUT}/${RES}_${SLURM_ARRAY_TASK_ID}

${PLINK2} --bfile ${BFILE} \
	--chr ${SLURM_ARRAY_TASK_ID} \
	--glm hide-covar cols='chrom,pos,ref,nobs,beta,se,tz,p' \
	--vif 500 --covar ${COVAR} \
	--covar-variance-standardize \
	--pheno ${PHENO} \
	--threads ${THREADS} \
	--memory ${MEMORY} \
	--out ${RES}

