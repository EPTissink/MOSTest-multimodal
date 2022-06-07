#!/bin/bash

which plink2

ROOT=/cluster/projects/p33/users/alexeas/elleke/pleioprs/geno

BFILE_IN=${ROOT}/all
BFILE_OUT=${ROOT}/all_qc

plink2 --bfile ${BFILE_IN} --maf 0.01 --geno 0.05 --make-bed --out ${BFILE_OUT}

