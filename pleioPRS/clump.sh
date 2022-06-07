#!/bin/bash

which plink

ROOT=/cluster/projects/p33/users/alexeas/elleke/pleioprs
BFILE=${ROOT}/geno/all_qc
SUMSTATS_DIR=${ROOT}/sumstats

for TRAIT in ADHD; do
    SUMSTATS=${SUMSTATS_DIR}/${TRAIT}.sumstats.prs.id.txt
    for CLUMP_FIELD in PVAL FDR; do
        echo ${TRAIT} ${CLUMP_FIELD}
        OUT=${ROOT}/sumstats/clump/${TRAIT}.${CLUMP_FIELD}
        plink --bfile ${BFILE} --clump ${SUMSTATS} --clump-field ${CLUMP_FIELD} \
            --clump-kb 250 --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-snp-field ID \
            --memory 8192 --out ${OUT} --silent --threads 8
    done
done
