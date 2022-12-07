#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)
library(ggplot2)
library(VennDiagram)
library(DescTools)
library(dplyr)
library(eulerr)

options(bitmapType="cairo")
Path <- "/path/"
setwd(Path)
signif_thr = 0.05/18877

########################################### 
# point towards loci txt and calc overlap #
###########################################

T1 <- fread(paste0(Path, "results/T1/magma.genes.out"))
rsfmri <- fread(paste0(Path, "results/rsfmri/magma.genes.out"))
diffusion <- fread(paste0(Path, "results/diffusion/magma.genes.out"))
multi <- fread(paste0(Path, "results/multimodal/magma.genes.out"))

T1 = T1[T1$P<signif_thr,1]
rsfmri = rsfmri[rsfmri$P<signif_thr,1]
diffusion = diffusion[diffusion$P<signif_thr,1]
multi = multi[multi$P<signif_thr,1]

############################
## Modality specific genes
############################

## T1
T1_diffusion <- inner_join(T1,diffusion)
T1_rsfmri <- inner_join(T1,rsfmri)
T1_multi = inner_join(T1,multi)

t1_overlap = unique(rbind(T1_diffusion,T1_rsfmri,T1_multi))
t1_unique =  T1 %>% anti_join(t1_overlap)

t1_overlap = unique(rbind(T1_diffusion,T1_rsfmri))
t1_specific = T1 %>% anti_join(t1_overlap)

## rsfmri
rsfmri_diffusion <- inner_join(rsfmri,diffusion)
rsfmri_t1 <- inner_join(T1,rsfmri)
rsfmri_multi = inner_join(rsfmri,multi)

rsfmri_overlap = unique(rbind(rsfmri_diffusion,rsfmri_t1,rsfmri_multi))
rsfmri_unique =  rsfmri %>% anti_join(rsfmri_overlap)

rsfmri_overlap = unique(rbind(rsfmri_diffusion,rsfmri_t1))
rsfmri_specific = rsfmri %>% anti_join(rsfmri_overlap)

## diffusion
diffusion_t1 <- inner_join(T1,diffusion)
diffusion_rsfmri <- inner_join(diffusion,rsfmri)
diffusion_multi = inner_join(diffusion,multi)

diffusion_overlap = unique(rbind(diffusion_t1,diffusion_rsfmri,diffusion_multi))
diffusion_unique =  diffusion %>% anti_join(diffusion_overlap)

diffusion_overlap = unique(rbind(diffusion_t1,diffusion_rsfmri))
diffusion_specific = diffusion %>% anti_join(diffusion_overlap)

############################
## Cross-modality loci
############################

## In >2 modalities
crossmodality = unique(rbind(T1_diffusion,T1_rsfmri,diffusion_rsfmri))

## Unique for multimodal
multi_overlap = unique(rbind(diffusion_multi,rsfmri_multi,T1_multi))
multi_unique = multi %>% anti_join(multi_overlap)

## Merge
multimodal = rbind(crossmodality,multi_unique)

############################
## Save loci
############################
write.table(t1_specific,file=paste0(Path,"results/T1/T1_specific_genes.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(rsfmri_specific,file=paste0(Path,"results/rsfmri/rsfmri_specific_genes.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(diffusion_specific,file=paste0(Path,"results/diffusion/diffusion_specific_genes.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
write.table(multimodal,file=paste0(Path,"results/multimodal/crossmodality_genes.txt"),quote = F,row.names = F,col.names = T,sep = "\t")


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
             "A&B" = nrow(anti_join(T1_multi,union(rsfmri,diffusion))),
             "A&C" = nrow(anti_join(rsfmri_multi,union(T1,diffusion))),
             "A&D" = nrow(anti_join(diffusion_multi,union(T1,rsfmri))),
             "B&C" = nrow(anti_join(T1_rsfmri,union(multi,diffusion))),
             "B&D" = nrow(anti_join(T1_diffusion,union(multi,rsfmri))),
             "C&D" = nrow(anti_join(diffusion_rsfmri,union(T1,diffusion))),
             "B&C&D" = nrow(anti_join(inner_join(inner_join(rsfmri,diffusion),T1),multi)), 
             "A&C&D"= nrow(anti_join(inner_join(inner_join(rsfmri,diffusion),multi),T1)),
             "A&B&C" = nrow(anti_join(inner_join(inner_join(rsfmri,multi),T1),diffusion)),
             "A&B&D" = nrow(anti_join(inner_join(inner_join(multi,diffusion),T1),rsfmri)),
             "A&B&C&D" = nrow(inner_join(inner_join(inner_join(rsfmri,multi),diffusion),T1))))
plot(pd,quantities = T)
