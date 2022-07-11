source ~/.bashrc
#! /usr/bin/env Rscript
#loadr4

rm(list=ls())
library(data.table)

Path <- "/cluster/projects/p33/groups/imaging/ukbio"
Batch="42k"
BatchDate="160320"
Date="051119"

#################
# Prep DM GWAS  #
#################

############### Add pop comps and genetic batch to demog file
Demog <- fread(paste0(Path,"/stats/batch_",Batch,"/UKB",Batch,"_demogShort_",BatchDate,".txt"),data.table=F)
PopComp <- fread(paste0(Path,"/genetics/popComps/UKB500k_popComps_",Date,".txt"),data.table=F)[,c("ID",paste0("C",1:20))]
DM <- merge(Demog,PopComp, by = "ID")

############### Clean up and remove missings, non-caucasians
DM$Sex <- as.data.frame(model.matrix(~ Sex - 1, data=DM))[,2]
DM <- DM[!duplicated(DM$ID),]
DM <- droplevels(DM)

Cauc <- fread(paste0(Path,"/demographics/UKB500k_caucasians_",Date,".txt"),data.table=F)
DM <- DM[DM$ID %in% Cauc$V2,]

AvGen <- fread(paste0(Path,"/genetics/UKB",Batch,"_",BatchDate,"/UKB",Batch,"_QCed_",BatchDate,".fam"),data.table=F)
DM <- DM[DM$ID %in% AvGen$V2,]

############## Remove bad scans based on Euler numbers
Euler <- fread(paste0(Path,"/stats/batch_",Batch,"/allEuler_UKB.csv"),data.table=F)
Euler$ID = substring(Euler$subject,4)
Euler$subject = NULL
DM <- merge(DM,Euler, by = "ID")
DM <- DM[complete.cases(DM),]

Outliers <- c()
for(i in 1:length(unique(DM$Scanner))){
  current <- DM[DM$Scanner==unique(DM$Scanner)[i],]
  current$euler_lh_resid <- scale(lm(euler_lh ~ Sex + Age, data = current)$residuals)
  current$euler_rh_resid <- scale(lm(euler_rh ~ Sex + Age, data = current)$residuals)
  Outliers <- c(Outliers,current$ID[which(current$euler_lh_resid < -3 | current$euler_rh_resid < -3)])
}

DM <- DM[!(DM$ID %in% Outliers),]
DM$Euler <- DM$euler_lh + DM$euler_rh

### Exclude retracted consent
Withdrawn = fread("/tsd/p33/data/durable/s3-api/ukblake/participant_withdrawal/w27412_20220222.csv",data.table=F)
colnames(Withdrawn) = "ID"
DM = DM[!(DM$ID %in% Withdrawn),]

#### Add global area and thick
measures <- c("area","thickness")
for(i in 1:2){
  templ <- fread(paste0(Path,"/stats/batch_",Batch,"/regularFSstats/lh.",measures[i],".UKBB.txt"),data.table=F)
  templ$ID = substring(templ[[paste0("lh.aparc.", measures[i])]],4)
  tempr <- fread(paste0(Path,"/stats/batch_",Batch,"/regularFSstats/rh.",measures[i],".UKBB.txt"),data.table=F)
  tempr$ID = substring(tempr[[paste0("rh.aparc.", measures[i])]],4)
  temp <- merge(templ,tempr,by="ID")
  temp <- temp[!duplicated(temp$ID),]
  DM <- merge(DM,temp,by="ID")
} 

DM$MeanThickness<-(DM$lh_MeanThickness_thickness + DM$rh_MeanThickness_thickness)/2
DM$TotalArea<-DM$lh_WhiteSurfArea_area + DM$rh_WhiteSurfArea_area

#### Add eICV for subcortical volumes
SubCort = fread(paste0(Path,"/stats/batch_",Batch,"/regularFSstats/subcorticalstats.UKBB.txt"),data.table=F) 
SubCort$ID = substring(SubCort$Measure,4)
DM = merge(DM,SubCort,by="ID")

#### Add functional covariates
Cov = fread(paste0(Path,"/demographics/UKB",Batch,"_demogAll_",BatchDate,".txt"),data.table=F)
Cov$Motion = Cov[,grepl("25741-2.0",names(Cov))]
Cov$SNR = Cov[,grepl("25744-2.0",names(Cov))]
DM = merge(DM,Cov[,c("ID","Motion","SNR")],by="ID")

### Add genotypes array
Cov$Array = as.numeric(Cov[, grepl( "22000-0.0" , names(Cov) ) ])>1
DM = merge(DM,Cov[,c("ID","Array")], by="ID")
DM$IID <- DM$ID -> DM$FID

##### Write everything to one file
write.table(DM,file=paste0("/cluster/projects/p33/users/elizabept/multimodal/discovery/phenos/data/",Batch,"_Cov_phenosT1_",BatchDate,".txt"),sep="\t",quote=F,row.names=F)

##### Write Covariates to file
write.table(DM[,c("FID","IID","Age","Sex","Scanner","Array","Euler","MeanThickness","TotalArea","EstimatedTotalIntraCranialVol","Motion","SNR",paste0("C",1:20))],
file=paste0("/cluster/projects/p33/users/elizabept/multimodal/discovery/phenos/data/",Batch,"_COV_",BatchDate,".txt"), sep = "\t", quote=FALSE, row.names=FALSE)