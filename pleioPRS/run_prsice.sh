#!/bin/bash

which PRSice_linux

ROOT=/cluster/projects/p33/users/alexeas/elleke/pleioprs
SUMSTATS_DIR=${ROOT}/sumstats
CLUMP_DIR=${SUMSTATS_DIR}/clump
PRS_DIR=${ROOT}/prs
TARGET=${ROOT}/geno/all_qc
PHENO=${ROOT}/pheno/pheno.txt
COVAR=${ROOT}/pheno/covar.txt
THREADS=8

for TRAIT in ASD; do
    SUMSTATS=${SUMSTATS_DIR}/${TRAIT}.sumstats.prs.id.txt
    for CLUMP_FIELD in PVAL FDR; do
        CLUMPED_FILE=${CLUMP_DIR}/${TRAIT}.${CLUMP_FIELD}.clumped
        for N in 10 100 1000 10000 50000 100000 150000; do
            echo ${TRAIT} ${CLUMP_FIELD} ${N}
            SNPS2USE=${CLUMPED_FILE}.${N}.snps
            OUT=${PRS_DIR}/${TRAIT}.${CLUMP_FIELD}.${N}
            if [ ${TRAIT} == "ASD" ]; then
                PHENO_COL="AUT"
            else
                PHENO_COL=${TRAIT}
            fi
            PRSice_linux --base ${SUMSTATS} --binary-target T --target ${TARGET} \
                --pheno ${PHENO} --pheno-col ${PHENO_COL} --cov ${COVAR} --snp ID \
                --stat BETA --a1 A1 --a2 A2 --bp BP --chr CHR --pvalue ${CLUMP_FIELD} \
                --no-clump --print-snp --bar-levels 1 --no-full --fastscore \
                --extract ${SNPS2USE} --thread ${THREADS} --out ${OUT}
        done
    done
done
