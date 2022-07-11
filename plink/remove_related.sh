#!/bin/bash
#SBATCH --job-name=subsetBfiles
#SBATCH --account=p33
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8GB

source ~/.bashrc
module purge
module load plink2/2.00a2LM

set -o errexit

plink2 	--bfile ${ukb}/genetics/UKB42k_160320/UKB42k_QCed_160320 \
	--king-cutoff 0.05 \
	--keep ${multimodal}/phenos/results/UKB42k_T1_qnorm_ecdf_230222.txt \
	--maf 0.005 \
	--make-bed \
	--out ${multimodal}/phenos/results/UKB34886_T1_QCed_230222 \
	--threads 32
