library(data.table)
require(plyr)
library(ggplot2)
library(RColorBrewer)
rm(list=ls())
setwd("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Genes/Celltypes/")


# Bhaduri et al
Bhaduri_Cross_modality = fread("Bhaduri_FetalBrain/MOSTest_Bhaduri_FisherExact.txt")
colnames(Bhaduri_Cross_modality)[1] = "CellType_2"
Bhaduri_sMRI_specific = fread("Bhaduri_FetalBrain/T1_Bhaduri_FisherExact.txt")
colnames(Bhaduri_sMRI_specific)[1] = "CellType_2"
Bhaduri_dMRI_specific = fread("Bhaduri_FetalBrain/diffusion_Bhaduri_FisherExact.txt")
colnames(Bhaduri_dMRI_specific)[1] = "CellType_2"
Bhaduri_fMRI_specific = fread("Bhaduri_FetalBrain/rsfmri_Bhaduri_FisherExact.txt")
colnames(Bhaduri_fMRI_specific)[1] = "CellType_2"

## rbind
dfnames <- c("Bhaduri_Cross_modality", "Bhaduri_sMRI_specific","Bhaduri_dMRI_specific","Bhaduri_fMRI_specific")
grouped = do.call(rbind, lapply(dfnames, function(x) cbind(get(x), Source=x)))

## Filter/QC
grouped_filtered = grouped[grouped$n_overlap>2,]

## Significance
grouped_threshold = 0.05/nrow(grouped_filtered)

##Plot
ggplot(data = grouped_filtered, mapping = aes(x = CellType_2,y = Source,fill = ifelse(n_overlap > 2, OR, NA))) +
  geom_tile() +
  geom_text(aes(label=ifelse(p.value < 0.05/nrow(grouped_filtered), "*",NA)), colour = "black", size=6, check_overlap = TRUE) +
  scale_fill_gradient(name = "OR",
                      low = "#edf8b1",
                      high = "#2c7fb8",
                      na.value="white") +
  theme_classic() +
  theme(strip.placement = "outside", # Move depth boxes to bottom of plot
        strip.text.x = element_text(size = 13),
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.y=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1,size=12),
        axis.title.x = element_blank(),
        #strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        legend.position = "top")
