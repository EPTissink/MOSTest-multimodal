#!/bin/bash

which plink2
which plink

ROOT=/cluster/projects/p33/users/alexeas/elleke/pleioprs/geno
GENO_DIR=/cluster/projects/p33/groups/biostat/genetics/imputed/cur/stage1
BATCHES="norment_i_a0807b norment_ii_fe2200 norment_iii_0c942c"
MERGE_LIST=${ROOT}/merge_list.txt
rm -f ${MERGE_LIST}

for BATCH in ${BATCHES}; do
    BFILE=${GENO_DIR}/${BATCH}/bed/all
    OUT=${ROOT}/${BATCH}
    plink2 --bfile ${BFILE} --keep ${ROOT}/ids2keep.txt --make-bed --out ${OUT}
    echo ${OUT} >> ${MERGE_LIST}
done

BFILE_MERGED=${ROOT}/all
plink --merge-list ${MERGE_LIST} --make-bed --out ${BFILE_MERGED}

