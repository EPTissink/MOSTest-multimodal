#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)
library(readr)

Path = "/cluster/projects/p33/users/elizabept/multimodal/discovery/phenos"
GenPath <- "/cluster/projects/p33/groups/imaging/ukbio/genetics"

Batch="42k"
Date="160320"

######################
# Merge all measures #
######################

daniel_edge = fread("/cluster/projects/p33/users/dtroelfs/NORMENT/network_genetics/roi_analysis/nets/schaefer_17networks_edge_features.txt",data.table = F)
daniel_node = fread("/cluster/projects/p33/users/dtroelfs/NORMENT/network_genetics/roi_analysis/nets/schaefer_17networks_node_features.txt",data.table = F)
DM = merge(daniel_edge,daniel_node,by="eid")
DM$FID <- DM$eid -> DM$IID

ICA <- fread(paste0(Path,"/results/QCed_norm/UKB42k_rsfMRI_qnorm_ecdf_070322.txt"),data.table=F)

DM = DM[DM$IID %in% ICA$IID,]
DM$eid = NULL

##########################
# Inverse-normal transformation #
##########################

##### Apply inverse-normal transformation
Features <- colnames(DM)[-c(154:155)]
for(f in Features) {
  DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
}

#### Save File
write.table(DM[,c("FID","IID",Features)],file=paste0(Path,"/results/QCed_norm/UKB",Batch,"_schaefer_qnorm_ecdf_230222.txt"),sep = "\t",quote=FALSE,row.names=FALSE)