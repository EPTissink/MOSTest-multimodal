#!/bin/bash
#SBATCH --job-name=munge_sumstats
#SBATCH --account=p33_norment
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=0-400

set -o errexit
source ../settings.sh

source ${SOFTWARE}/ldsr/bin/activate

which python

GWAS_DIR=/cluster/projects/p33/users/elizabept/multimodal/discovery/plink/results/orig
LDSR_SRC=${DATA}/ldsc
LDSR_OUT=${OUT}/munge
PHENOS=($(ls -1 ${GWAS_DIR}/univ_QCed_orig.*_PCs*.glm.linear))

FILENAME=${PHENOS[${SLURM_ARRAY_TASK_ID}]}
PHENO=$(cut -f2 -d'.' <<< $FILENAME)
STAT="T_STAT,0"

python ${LDSR_SRC}/munge_sumstats.py --sumstats $FILENAME --snp ID --signed-sumstats ${STAT} --N-col OBS_CT --a1 A1 --a2 REF --ignore BETA --merge-alleles ${LDSR_SRC}/w_hm3.snplist --out ${LDSR_OUT}/${PHENO}
