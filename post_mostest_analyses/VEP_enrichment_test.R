library(data.table)
library(dplyr)
library(ggplot2)

## Load functional annotation output from VEP (website)
crossmodal_vep = fread("/path/crossmodality_lead.ids.VEP_rs_consequence.txt",data.table = F,header=F)
smri_vep = fread("/path/T1_specific_lead.ids.VEP_rs_consequence.txt",data.table = F,header=F)
dmri_vep = fread("/path/diffusion_specific_lead.ids.VEP_rs_consequence.txt",data.table = F,header=F)
fmri_vep = fread("/path/rsfmri_specific_lead.ids.VEP_rs_consequence.txt",data.table = F,header=F)

# Obtain frequencies/counts (and proportions)
data_cross = as.data.frame(table(crossmodal_vep$V2)/dim(crossmodal_vep)[1])
data_cross$count = table(crossmodal_vep$V2)
data_cross$modality = "Cross-modality lead SNPs"
data_smri = as.data.frame(table(smri_vep$V2)/dim(smri_vep)[1])
data_smri$count = table(smri_vep$V2)
data_smri$modality = "sMRI-specific lead SNPs"
data_dmri = as.data.frame(table(dmri_vep$V2)/dim(dmri_vep)[1])
data_dmri$count = table(dmri_vep$V2)
data_dmri$modality = "dMRI-specific lead SNPs"
data_fmri = as.data.frame(table(fmri_vep$V2)/dim(fmri_vep)[1])
data_fmri$count = table(fmri_vep$V2)
data_fmri$modality = "fMRI-specific lead SNPs"

# Create average reference out of 10 randomly pruned sets of 10K SNPs
for (n in 1:10) {
  ref_vep = fread(paste0("/path/random_lead.id_",n,".rand_10K.VEP_rs_consequence.txt"),data.table = F,header = F)
  data = as.data.frame(table(ref_vep$V2))
  if (n == 1) {data_ref = data}
  else data_ref = merge(data_ref,data,by="Var1")
}
data_ref$refcount = round(rowMeans(data_ref[,2:11]))

## Enrichment test function
# This function takes input "x", which is one row with two columns from the ref.count dataframe
# The values in x are (in order):
#    ref.count: How many times did we count this annotation for the reference panel?
#    count: How many times did we count this annotation for our selected variants?
# This function also uses two variables that are defined outside the function
#    N: How many times total have we counted any annotation for the reference data?
#    n: How many times total have we counted any annotation for our selected variants?
# It uses these values to set up the following table for a Fisher exact test:
#             | selected | nonselected | allvariants
#    curanno  |   x[2]   |             |   x[1]
#    othanno  |          |             |
#    allanno  |    n     |             |    N

calcFisherP <- function(x) {
  x <- as.numeric(x)
  curanno_selected <- x[2]
  othanno_selected <- n - x[2]
  curanno_nonselected <- x[1] - x[2]
  othanno_nonselected <- N - n - curanno_nonselected
  data <- 
    matrix(c(curanno_selected,othanno_selected,curanno_nonselected,othanno_nonselected), 
           ncol=2)
  return(fisher.test(data)$p.value)
}
calcFisherOR <- function(x) {
  x <- as.numeric(x)
  curanno_selected <- x[2]
  othanno_selected <- n - x[2]
  curanno_nonselected <- x[1] - x[2]
  othanno_nonselected <- N - n - curanno_nonselected
  data <- 
    matrix(c(curanno_selected,othanno_selected,curanno_nonselected,othanno_nonselected), 
           ncol=2)
  return(fisher.test(data)$estimate)
}

## Apply function
out = c()
dflist <- list(data_smri, data_dmri, data_fmri, data_cross)
for (i in dflist) {
  data = merge(i,data_ref[,c("Var1","refcount")],by="Var1")
  n = as.numeric(sum(i$count))
  N = as.numeric(sum(data_ref$refcount))
  data$fisher.P <- apply(data[,c("refcount","count")], 1, calcFisherP)
  data$OR <- apply(data[,c("refcount","count")], 1, calcFisherOR)
  out = rbind(out,data)
}

write.table(out, "/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/VEP.txt",
            quote = F,row.names = F,col.names = T,sep = "\t")

## Prepare for plotting
out$Var1 = gsub("_"," ",out$Var1)
out$Bonf = "P < 0.05"
out[out$fisher.P>0.05,"Bonf"] = "P > 0.05"
out[out$fisher.P<(0.05/nrow(out)),"Bonf"] = "P < 0.05/41"
out$Bonf = factor(out$Bonf,levels=c("P > 0.05","P < 0.05","P < 0.05/41"))
out$modality = factor(out$modality,levels=c("Cross-modality lead SNPs",
                                            "sMRI-specific lead SNPs",
                                            "dMRI-specific lead SNPs",
                                            "fMRI-specific lead SNPs"))
## Plot
library(wesanderson)
R2 <- wes_palettes$Royal1

ggplot(out, aes(log2(OR), reorder(Var1, OR), fill = factor(Bonf))) +
  facet_grid(cols=vars(modality)) +
  geom_bar(stat="identity",position="dodge",width = 0.7)  +
  scale_fill_manual(values=wes_palette(n=3, name="Royal1")[c(3,1,2)]) +
  ylab(label = element_blank()) +
  theme_light()+
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    legend.title = element_blank(),
    legend.position = "right") + 
  #scale_alpha_discrete(range = c(0.3, 1)) +
  labs(x ="Enrichment (log2(OR))")