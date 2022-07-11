#!/bin/bash
#SBATCH --job-name=clump
#SBATCH --account=p33_norment
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --array=0-50:2

set -o errexit
source ../settings.sh
source ${SOFTWARE}/tsd/bin/activate
module load plink/1.90b6.2

BFILES=/cluster/projects/p33/users/elizabept/multimodal/discovery/genos/results/LD_ref

for MODALITY in diffusion T1 rsfmri multimodal; do
	
		python3 ${SOFTWARE}/python_convert-master/sumstats.py clump \
			--clump-field PVAL \
			--force  \
			--plink plink \
			--exclude-ranges 6:25119106-33854733 \
			--sumstats ${OUT}/${MODALITY}/${MODALITY}_eig${SLURM_ARRAY_TASK_ID}.most.orig.csv.gz \
			--bfile-chr ${BFILES}/chr@ \
			--clump-p1 5e-8 \
			--out ${OUT}/${MODALITY}/${MODALITY}_mostest_eig${SLURM_ARRAY_TASK_ID}.most.orig.clump
done
