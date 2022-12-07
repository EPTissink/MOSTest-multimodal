#!/bin/bash
#SBATCH --job-name=rg_sumstats
#SBATCH --account=p33_norment
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=0-1

set -o errexit
source ../settings.sh
source ${SOFTWARE}/ldsr/bin/activate
which python

LDSR_SRC=${DATA}/ldsc
LDSR_OUT=${OUT}/rg
LDSR_IN=${OUT}/munge
readarray -t LIST < ${DATA}/pheno_path.txt
FILENAME1=${LIST[${SLURM_ARRAY_TASK_ID}]}
PHENO=$(cut -f1 -d'.' <<< $(basename $FILENAME1))
FILENAME2=$( echo ${LIST[@]:${SLURM_ARRAY_TASK_ID}+1} | sed 's/ /,/g' )

python ${LDSR_SRC}/ldsc.py --rg $FILENAME1,$FILENAME2 --ref-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --w-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --out ${LDSR_OUT}/${PHENO}.rg
