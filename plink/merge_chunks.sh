#!/bin/bash

source ${HOME}/.profile
which plink

FIRST_CHUNK_I=0
LAST_CHUNK_I=906
DIR=/cluster/projects/p33/users/alexeas/elleke/data
BFILE_CHUNKS_DIR=${DIR}/bfile_chunks
MAX_PADDING=3
MERGE_LIST=${DIR}/bfile_chunks_list.txt
BFILE_MERGED=${DIR}/UKB34886_T1_QCed_230222.perm

rm -rf ${MERGE_LIST}
for i in $(seq ${FIRST_CHUNK_I} ${LAST_CHUNK_I}); do
    CHUNK_ID=$(printf %0${MAX_PADDING}d ${i})
    echo "${BFILE_CHUNKS_DIR}/UKB34886_T1_QCed_230222.chunk${CHUNK_ID}.perm" >> ${MERGE_LIST}    
done

plink --merge-list ${MERGE_LIST} --make-bed --memory 40000 --out ${BFILE_MERGED}

