# Set up environment
setwd("/path/h2")
library(data.table)
library(ggplot2)
library(ggridges)

options(bitmapType="cairo")

# Load data
data = fread("summary_schaefer.txt")

# Threshold diffusion elbow
NF_excl = paste0("NF_PCs",c(88:1000))
N0_excl = paste0("N0_PCs",c(66:1000))
ND_excl = paste0("ND_PCs",c(125:1000))
data = data[!data$pheno %in% NF_excl,]
data = data[!data$pheno %in% N0_excl,]
data = data[!data$pheno %in% ND_excl,]
data = data[!grepl('^edge|^node', data$pheno),]

# set Ratio < 0 (usually indicates GC correction) to NA
data[data$Ratio == "<",c("Ratio","Ratio_se")] = NA 

# Create modalities
Yeo = c("VisCent","VisPeri","SomMotA","SomMotB","DorsAttnA","DorsAttnB",
        "SalVentAttnA","SalVentAttnB","LimbicB","LimbicA","ContA",
        "ContB","ContC","DefaultA","DefaultB","DefaultC","TempPar")
data$modality = "sMRI"
data$modality[grepl('conn', data$pheno)] = "fMRI"
data$modality[grepl('PCs',data$pheno)] = "dMRI"  
data$modality[data$pheno %in% Yeo] = "fMRI"

# Create submodalities
data$submodality = "Subcortical volume"
data$submodality[grepl('conn', data$pheno)] = "Network functional connectivity"
data$submodality[data$pheno %in% Yeo] = "Network temporal variance"
data$submodality[grepl('NF',data$pheno)] = "Free water (NF)"
data$submodality[grepl('ND',data$pheno)] = "Axons (ND)"
data$submodality[grepl('N0',data$pheno)] = "Cell bodies (N0)"
data$submodality[grepl('thickness',data$pheno)] = "Cortical thickness"
data$submodality[grepl('area',data$pheno)] = "Cortical surface area"
data$submodality = factor(data$submodality,levels=c("Subcortical volume",
                                                    "Cortical thickness",
                                                    "Cortical surface area",
                                                    "Network temporal variance",
                                                    "Network functional connectivity",
                                                    "Free water (NF)",
                                                    "Axons (ND)",
                                                    "Cell bodies (N0)"))

# Remove insignificant phenotypes
data = data[data$h2_z>1.96,]

#
library(tidyverse)
dataMedian <- summarise(group_by(data, submodality), MD = median(h2))
dataIQR <- summarise(group_by(data, submodality), MD = IQR(h2))

# Density plot
ggplot(data, aes(x = h2,fill=modality)) +
  geom_density(position = 'identity', alpha = 0.5, adjust = 1) +
  labs(x = 'SNP-based heritability (h2)', y = 'Density')+
  labs(fill = 'Modality')+
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    text = element_text(family = 'Times'),
    legend.position = 'right'
  )



# Plotting individual phenotype h2 does not really work
ggplot(data, aes(fill = modality)) +
  geom_errorbar(data, mapping=aes(x=phenotype, ymin=(h2-h2_se), ymax=(h2+h2_se)), width=0.2, size=1) + 
  geom_point(data, mapping=aes(x=phenotype, y=h2), size=4, shape=21) +
  labs(x = '', y = 'SNP-based heritability')+
  scale_colour_brewer(palette = 'Dark2') +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    text = element_text(family = 'Times')
  )

# Ridge plot
ggplot(data, aes(x = h2, y = submodality , fill = modality)) +
  geom_density_ridges(position = 'identity', alpha = 0.8) +
  labs(x = 'SNP-based heritability (h2)')+
  labs(fill = 'Modality')+
  scale_fill_manual(values=c("#FDC0BE","#ADD8E6","#CDEECC"))  +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = 'right'
  )
