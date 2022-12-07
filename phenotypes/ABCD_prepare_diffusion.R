library(data.table)
library(dplyr)

############
## Load data
############

# Load Phenotypes
PCs = fread("/path/ABCD11904_diffusion_QCed_230322.txt",data.table = F)
PCs$IID = PCs$FID = sapply(strsplit(PCs$Row,"_"), `[`, 1)
PCs$Row = paste0("NDAR_",PCs$IID)

# Recommended inclusion
inclusion = fread("/path/ABCDStudyNDA_4.0/abcd_imgincl01.txt",data.table = F)
inclusion = inclusion[startsWith(inclusion$eventname,"baseline"),] # only first visits
inclusion = inclusion[inclusion$imgincl_dmri_include == 1,] # quality control

PCs = PCs[PCs$Row %in% inclusion$subjectkey,]


#############
# Select PCs
#############

N0_PCs = PCs[,c("FID","IID",paste0("N0_PCs",1:65))] # N0 = 65
NF_PCs = PCs[,c("FID","IID",paste0("NF_PCs",1:87))] # NF = 87
ND_PCs = PCs[,c("FID","IID",paste0("ND_PCs",1:124))] # ND = 124

##########
## Save
##########

write.table(N0_PCs,file="/path/QCed/ABCD_10406_N0_QCed_040522.txt",
            sep = "\t",quote=FALSE,row.names=FALSE)

write.table(NF_PCs,file="/path/QCed/ABCD_10406_NF_QCed_040522.txt",
            sep = "\t",quote=FALSE,row.names=FALSE)

write.table(ND_PCs,file="/path/QCed/ABCD_10406_ND_QCed_040522.txt",
            sep = "\t",quote=FALSE,row.names=FALSE)

