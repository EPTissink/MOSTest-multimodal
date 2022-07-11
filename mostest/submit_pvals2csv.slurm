#!/bin/bash
#SBATCH --job-name=pval2csv
#SBATCH --account=p33_norment
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=0-50:2

set -o errexit
source ../settings.sh
source ${SOFTWARE}/tsd/bin/activate
BFILE=/cluster/projects/p33/users/elizabept/multimodal/discovery/genos/results/UKB34886_T1_QCed_230222.bim

for MODALITY in T1 diffusion rsfmri multimodal; do
	
	python pvals2csv.py --mat ${OUT}/${MODALITY}/${MODALITY}_mostest_eig${SLURM_ARRAY_TASK_ID} --bim ${BFILE} --n ${DATA}/sample_size.${MODALITY}.txt --out ${OUT}/${MODALITY}/${MODALITY}_eig${SLURM_ARRAY_TASK_ID}

done
