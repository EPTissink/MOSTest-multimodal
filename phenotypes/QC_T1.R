#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)

Path = "/cluster/projects/p33/users/elizabept/multimodal/discovery/phenos/data"
GenPath <- "/cluster/projects/p33/groups/imaging/ukbio/genetics"

Batch="42k"
Date="160320"

######################
# Merge all measures #
######################

DM <- fread(paste0(Path,"/",Batch,"_Cov_phenos_T1_",Date,".txt"),data.table=F)
Dict <- fread(paste0(Path,"/FSdict_ver53_230519.txt"),data.table=F)
colnames(DM) <- gsub("-","_",colnames(DM))

##########################
# Apply inverse-normal transformation #
##########################

Subcort <- Dict$ROI_name[grep("SubcorticalVolume\\b",Dict$ROI_type)]
Area <- Dict$ROI_name[grep("CorticalArea\\b",Dict$ROI_type)]
Thick <- Dict$ROI_name[grep("CorticalThickness\\b",Dict$ROI_type)]
Features <- c(Subcort,Area,Thick)
for(f in Features) {
  DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
}

#### Save File
write.table(DM[,c("FID","IID",Features)],file=paste0(Path,"/UKB",Batch,"_GWAS_qnorm_ecdf_230222.txt"),sep = "\t",quote=FALSE,row.names=FALSE)