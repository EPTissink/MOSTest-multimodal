#! /usr/bin/env Rscript

rm(list=ls())
library(GenomicRanges)
library(data.table)
library(eulerr)

options(bitmapType="cairo")
Path <- "/path/"
setwd(Path)

########################################### 
# point towards loci txt 
###########################################

t1 <- fread(paste0(Path, "results/T1/T1_mostest_eig16.most.orig.clump.loci.csv"))
rsfmri <- fread(paste0(Path, "results/rsfmri/rsfmri_mostest_eig0.most.orig.clump.loci.csv"))
diffusion <- fread(paste0(Path, "results/diffusion/diffusion_mostest_eig10.most.orig.clump.loci.csv"))
multi <- fread(paste0(Path, "results/multimodal/multimodal_mostest_eig26.most.orig.clump.loci.csv"))

## Remove singleton loci
t1 <- t1[t1$MaxBP-t1$MinBP>0,]
rsfmri <- rsfmri[rsfmri$MaxBP-rsfmri$MinBP>0,]
diffusion <- diffusion[diffusion$MaxBP-diffusion$MinBP>0,]
multi <- multi[multi$MaxBP-multi$MinBP>0,]

## Make GR objects
t1.gr=GRanges(seqnames=t1$CHR, ranges=IRanges(start=t1$MinBP,end=t1$MaxBP),names=t1$locusnum)
rsfmri.gr=GRanges(seqnames=rsfmri$CHR, ranges=IRanges(start=rsfmri$MinBP,end=rsfmri$MaxBP),names=rsfmri$locusnum)
diffusion.gr=GRanges(seqnames=diffusion$CHR, ranges=IRanges(start=diffusion$MinBP,end=diffusion$MaxBP),names=diffusion$locusnum)
multi.gr=GRanges(seqnames=multi$CHR, ranges=IRanges(start=multi$MinBP,end=multi$MaxBP),names=multi$locusnum)


############################
## Modality specific loci
############################

## T1
t1_diffusion = subsetByOverlaps(t1.gr,diffusion.gr)
t1_rsfmri = subsetByOverlaps(t1.gr,rsfmri.gr)
t1_multi = subsetByOverlaps(t1.gr,multi.gr)

t1_overlap = reduce(c(t1_multi,t1_rsfmri,t1_diffusion))
t1_unique = as(subsetByOverlaps(t1.gr,t1_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(t1_unique) = c("CHR","MinBP","MaxBP","locusnum")

t1_overlap = reduce(c(t1_rsfmri,t1_diffusion))
t1_specific = as(subsetByOverlaps(t1.gr,t1_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(t1_specific) = c("CHR","MinBP","MaxBP","locusnum")

## Diffusion
diffusion_t1 = subsetByOverlaps(diffusion.gr,t1.gr)
diffusion_rsfmri = subsetByOverlaps(diffusion.gr,rsfmri.gr)
diffusion_multi = subsetByOverlaps(diffusion.gr,multi.gr)

diffusion_overlap = reduce(c(diffusion_multi,diffusion_rsfmri,diffusion_t1))
diffusion_unique = as(subsetByOverlaps(diffusion.gr,diffusion_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(diffusion_unique) = c("CHR","MinBP","MaxBP","locusnum")

diffusion_overlap = reduce(c(diffusion_rsfmri,diffusion_t1))
diffusion_specific = as(subsetByOverlaps(diffusion.gr,diffusion_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(diffusion_specific) = c("CHR","MinBP","MaxBP","locusnum")

## rsfMRI
rsfmri_diffusion = subsetByOverlaps(rsfmri.gr,diffusion.gr)
rsfmri_t1 = subsetByOverlaps(rsfmri.gr,t1.gr)
rsfmri_multi = subsetByOverlaps(rsfmri.gr,multi.gr)

rsfmri_overlap = reduce(c(rsfmri_multi,rsfmri_t1,rsfmri_diffusion))
rsfmri_unique = as(subsetByOverlaps(rsfmri.gr,rsfmri_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(rsfmri_unique) = c("CHR","MinBP","MaxBP","locusnum")

rsfmri_overlap = reduce(c(rsfmri_t1,rsfmri_diffusion))
rsfmri_specific = as(subsetByOverlaps(rsfmri.gr,rsfmri_overlap,invert=T), "data.frame")[,-c(4,5)]
colnames(rsfmri_specific) = c("CHR","MinBP","MaxBP","locusnum")


############################
## Cross-modality loci
############################

## Unique for multimodal
multi_rsfmri = subsetByOverlaps(multi.gr,rsfmri.gr)
multi_t1 = subsetByOverlaps(multi.gr,t1.gr)
multi_diffusion = subsetByOverlaps(multi.gr,diffusion.gr)

union = reduce(c(multi_rsfmri,multi_t1,multi_diffusion))
multi_unique = as(subsetByOverlaps(multi.gr,union,invert=T), "data.frame")[,-c(4,5)]
colnames(multi_unique) = c("CHR","MinBP","MaxBP","locusnum")

## In >2 modalities
multi_diffusion_t1 = subsetByOverlaps(multi_diffusion,multi_t1)
multi_diffusion_rsfmri = subsetByOverlaps(multi_diffusion,multi_rsfmri)
multi_rsfmri_t1 = subsetByOverlaps(multi_rsfmri,multi_t1)

crossmodal = unique(rbind(values(multi_rsfmri_t1),values(multi_diffusion_t1),values(multi_diffusion_rsfmri)))
crossmodal = multi[multi$locusnum %in% crossmodal$names,c(1,2,5,6)]

## Merge
multimodal = multi[multi$locusnum %in% c(crossmodal$locusnum,multi_unique$locusnum),c(1,2,5,6)]

############################
## Save loci
############################
write.table(t1_specific,file=paste0(Path,"results/T1/T1_specific_loci.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(rsfmri_specific,file=paste0(Path,"results/rsfmri/rsfmri_specific_loci.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(diffusion_specific,file=paste0(Path,"results/diffusion/diffusion_specific_loci.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(multimodal,file=paste0(Path,"results/multimodal/crossmodality_loci.txt"),quote = F,row.names = F,col.names = T,sep = "\t")

##########################
## Venn Diagram
##########################
# A = multimodal
# B = t1
# C = rsfmri
# D = diffusion
pd = euler(c(A = dim(multi_unique)[1],
             B = dim(t1_unique)[1],
             C = dim(rsfmri_unique)[1],
             D = dim(diffusion_unique)[1],
             "A&B" = ((dim(t1_specific)[1])-(dim(t1_unique)[1])),
             "A&C" = (dim(rsfmri_specific)[1]-dim(rsfmri_unique)[1]),
             "A&D" = (dim(diffusion_specific)[1]-dim(diffusion_unique)[1]),
             "B&C" = length(subsetByOverlaps(rsfmri_t1,multi.gr,invert=T)),
             "B&D" = length(subsetByOverlaps(t1_diffusion,multi.gr,invert=T)),
             "C&D" = length(subsetByOverlaps(diffusion_rsfmri,multi.gr,invert=T)),
             "B&C&D" = 0, 
             "A&C&D"=length(subsetByOverlaps(multi_diffusion_rsfmri,multi_rsfmri_t1,invert=T)),
             "A&B&C" = length(subsetByOverlaps(multi_rsfmri_t1,multi_diffusion_rsfmri,invert=T)),
             "A&B&D" = length(subsetByOverlaps(multi_diffusion_t1,rsfmri.gr,invert=T)),
             "A&B&C&D" = sum(countOverlaps(multi_diffusion_rsfmri,multi_rsfmri_t1))))
plot(pd,quantities = T)
