#!/bin/bash

ROOT=/cluster/projects/p33/users/alexeas/elleke/pleioprs
CLUMP_DIR=${ROOT}/sumstats/clump

for TRAIT in ADHD; do
    for CLUMP_FIELD in PVAL FDR; do
        CLUMPED_FILE=${CLUMP_DIR}/${TRAIT}.${CLUMP_FIELD}.clumped
        for N in 10 100 1000 10000 50000 100000 150000; do
            OUT=${CLUMPED_FILE}.${N}.snps
            tail -n+2 ${CLUMPED_FILE} | awk '{print($3)}' | head -n ${N} > ${OUT}
        done
    done
done
