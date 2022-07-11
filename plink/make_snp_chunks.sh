#!/bin/bash

BIM=/cluster/projects/p33/users/elizabept/multimodal/discovery/genos/results/UKB34886_T1_QCed_230222.bim
N_SNP_IN_CHUNK=10000
OUT_PREFIX=/cluster/projects/p33/users/alexeas/elleke/data/snp_chunks/UKB34886_T1_QCed_230222.chunk

split -l ${N_SNP_IN_CHUNK} --numeric-suffixes --suffix-length 3 --additional-suffix ".txt" <(cut -f2 ${BIM}) ${OUT_PREFIX}

