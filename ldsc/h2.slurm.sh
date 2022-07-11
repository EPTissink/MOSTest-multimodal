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

LDSR_SRC=${DATA}/ldsc
LDSR_OUT=${OUT}/h2
LDSR_IN=${OUT}/munge
PHENOS=($(ls -1 ${LDSR_IN}/*_PCs*.sumstats.gz))

FILENAME=${PHENOS[${SLURM_ARRAY_TASK_ID}]}
PHENO=$(cut -f1 -d'.' <<< $(basename $FILENAME))

python ${LDSR_SRC}/ldsc.py --h2 $FILENAME --ref-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --w-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --out ${LDSR_OUT}/${PHENO}.h2
