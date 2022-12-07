## Prepare temporal variance and network connectivity for ABCD
## Network connectivity between 12 networks:
# auditory (AUD)
# cingulo-opercular (CON)
# cingulo-parietal (CP) aka ParietoOccip?
# default mode (DMN)
# dorsal attention (DAN)
# frontoparietal (FPC)
# retrosplenial (RSP) aka MedialParietal
# sensorimotor-hand (SM-hand)
# sensorimotor-mouth (SM-mouth)
# salience (SAL)
# ventral attention (VAN)
# visual (VIS)
# 13th "none" network is excluded

## Temporal variance:
# 333 parcels devided over same 12 networks (again "none" network is excluded)

library(data.table)
library(dplyr)

############
## Load data
############

# Load Phenotypes
cortical = fread("/path/ABCDStudyNDA_4.0/abcd_betnet02.txt",data.table = F)
cortical = cortical[c(1,grep("baseline",cortical$rsfmri_c_ngd_visitid)),] # only first visits
subcortical = fread("/path/ABCDStudyNDA_4.0/mrirscor02.txt",data.table = F)
subcortical = subcortical[c(1,grep("baseline",subcortical$rsfmri_cor_ngd_scs_visitid)),] # only first visits
temporal = fread("/path/ABCDStudyNDA_4.0/abcd_mrirstv02.txt",data.table = F)
temporal = temporal[c(1,grep("baseline",temporal$rsfmri_var_visitid)),] # only first visits

# Recommended inclusion
inclusion = fread("/path/ABCDStudyNDA_4.0/abcd_imgincl01.txt",data.table = F)
inclusion = inclusion[startsWith(inclusion$eventname,"baseline"),] # only first visits
inclusion = inclusion[inclusion$imgincl_rsfmri_include == 1,] # quality control

# Covariate information
scanner_info = fread("/path/ABCDStudyNDA_4.0/abcd_mri01.txt",data.table = F)
scanner_info = scanner_info[endsWith(scanner_info$mri_info_visitid,"baseline"),] # only first visits
batch_info = fread("/path/ABCD_genotype/ABCD_release3.0_.batch_info.txt",data.table = F)
QC = fread("/path/ABCDStudyNDA_4.0/abcd_auto_postqc01.txt",data.table = F)
QC = QC[startsWith(QC$eventname,"baseline"),] # only first visits
PCs = fread("/path/ABCD_101219_QCed.pruned.pca.eigenvec",data.table = F)

data = left_join(inclusion,scanner_info,by="subjectkey")
data = left_join(data,batch_info,by=c("subjectkey"="abcd.id_redcap"))
data = left_join(data,QC,by="subjectkey")
data = left_join(data,temporal,by="subjectkey")
data = left_join(data,PCs,by=c("subjectkey"="IID"))

data = droplevels(data)
data$FID <- sub('NDAR_', '', data$subjectkey) -> data$IID

############                
# Select edges
############
cortical_connectivity = cortical[,grep("network",cortical[1,])]
cortical_connectivity[,c(1,15,29,43,57,71,85,99,113,127,141,155,169)] = NULL # exclude within network connectivity
cortical_connectivity <- cortical_connectivity[!duplicated(as.list(cortical_connectivity[-1,]))] # select lower triangle of connectivity matrix
cortical_connectivity = cortical_connectivity[,grep("none",cortical_connectivity[1,],invert = T)] # exclude none network
cortical_connectivity = as.data.frame(sapply(cortical_connectivity, as.numeric))
cortical_connectivity = cbind(cortical$subjectkey,cortical_connectivity)
names(cortical_connectivity)[names(cortical_connectivity) == "cortical$subjectkey"] <- "subjectkey"

edges = cortical_connectivity
colnames(edges) = c("subjectkey",gsub("rsfmri_cor_ngd", "edge", colnames(edges[,-1])))
colnames(edges) = c("subjectkey",gsub("rsfmri_c", "edge", colnames(edges[,-1])))
  
#############
# Select nodes
#############

# Load file from parcel to network assignment
assignment = fread("/path/Gordon_parcels.txt",data.table = F)
parcels = as.data.frame(sapply(temporal[-1,grep("gp",colnames(temporal))], as.numeric))
assignment$network = as.factor(unlist(strsplit(assignment$V2,"_"))[grep("rh.R|lh.L|.label",unlist(strsplit(assignment$V2,"_")),invert = T)])

# Average temporal variance per network
for (network in levels(assignment$network)) {
  
  parcels_network = as.integer(subset(rownames(assignment),assignment$network == network))
  parcels$average = rowMeans(parcels[,parcels_network])
  names(parcels)[names(parcels) == 'average'] <- paste0("node_",network)

}

# Exclude "none" network
nodes = cbind(temporal$subjectkey[-1],parcels[,grep("node_",colnames(parcels))])
names(nodes)[names(nodes) == "temporal$subjectkey[-1]"] <- "subjectkey"
nodes$node_None = NULL

########################
## Transform normal phenotypes
########################

rsfmri = merge(nodes,edges,by="subjectkey")
rsfmri = rsfmri[rsfmri$subjectkey %in% data$subjectkey,]

# Ranked-based inverse normal transformation
Features <- colnames(rsfmri[,-1])
for(f in Features) {
  rsfmri[,f] <- qnorm(ecdf(rsfmri[,f])(rsfmri[,f]) - 0.5/dim(rsfmri)[1])
}

# save
rsfmri$FID <- sub('NDAR_', '', rsfmri$subjectkey) -> rsfmri$IID
write.table(rsfmri[,c("FID","IID",Features)],file="/path/ABCD_9627_rsfMRI_QCed_qnorm_ecdf_040522.txt",
            sep = "\t",quote=FALSE,row.names=FALSE)

########################
## Covariate file
########################

# age, age2, sex, genotyping array, 20 genetic PCs, scanner site, software imaging, signal to noise, motion
data = data[!is.na(data$PC1),]
data$Sex = as.data.frame(model.matrix(~ sex.x - 1, data=data))[,2]
data$Age = poly(as.numeric(data$interview_age.x),2)
#data$Age = as.numeric(data$interview_age.x)
#data$Software = paste0("Software_",as.numeric(factor(data$mri_info_softwareversion)))
  
base_cov = c("Age","Sex","mri_info_deviceserialnumber","BATCH",paste0("PC",c(1:20)))
fmri_cov = c("rsfmri_var_meanmotion","apqc_fmri_intra_tsnr")

cov = data[,c("FID","IID",base_cov,fmri_cov)]
colnames(cov) = c("FID","IID","Age","Sex","Scanner","Batch",paste0("PC",c(1:20)),"Motion","SNR")

# Save
write.table(cov,file="/path/ABCD_cov_rsfMRI_130522.txt",
            sep = "\t",quote=FALSE,row.names=FALSE)
