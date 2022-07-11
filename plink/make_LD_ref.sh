#!/bin/bash
source ../settings.sh


PLINK2=${SOFTWARE}/plink2
THREADS=8
MEMORY=$((8192 * ${THREADS}))

BFILE=${OUT}/UKB34886_T1_QCed_230222

for i in {1..22}; do

	${PLINK2} --bfile ${BFILE} \
		--chr ${i} \
		--keep ${OUT}/LD_ref.txt \
		--threads ${THREADS} \
		--memory ${MEMORY} \
		--make-bed \
		--out ${OUT}/LD_ref/chr${i}
done
