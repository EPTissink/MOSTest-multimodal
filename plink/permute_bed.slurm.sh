#!/bin/bash

#SBATCH --job-name=permute_bfile
#SBATCH --account=p33_norment
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=4

source /cluster/bin/jobsetup
module purge

module load  Python/3.7.4-GCCcore-8.3.0
source /cluster/projects/p33/users/elizabept/software/py3/bin/activate

module list

python3 permute_bed.py \
	--bfile /cluster/projects/p33/users/elizabept/multimodal/discovery/genos/results/UKB34886_T1_QCed_230222 \
	--out /cluster/projects/p33/users/elizabept/multimodal/discovery/genos/results/UKB34886_T1_QCed_080322_perm
