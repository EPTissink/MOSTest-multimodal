## Setup
setwd("/cluster/p/p33/cluster/users/elizabept/multimodal/discovery/phenos/data/diffusion")
library(ggplot2)
library(data.table)
library(nFactors)

## Load data
PC = c(c(1:1000),c(1:1000),c(1:1000))
metric = c(rep("N0", 1000), rep("NF", 1000), rep("ND", 1000))
N0 = fread("N0_Eigenvalues.txt",data.table = F,header=T)
NF = fread("NF_Eigenvalues.txt")
ND = fread("ND_Eigenvalues.txt")

df = cbind(PC,metric,rbind(N0,NF,ND))

## Plot Eigenvalues / Scree plot
qplot(df$PC, df$Var1,colour=metric) + 
  xlab("Principal Component") + 
  ylab("Eigen value") +
  theme_classic()

## Determine elbow
nScree(x=as.vector(unlist(N0)),model = "components",aparallel = NULL) # optimal coordinates = 65
nScree(x=as.vector(unlist(NF)),model = "components",aparallel = NULL) # optimal coordinates = 87
nScree(x=as.vector(unlist(ND)),model = "components",aparallel = NULL) # optimal coordinates = 124


