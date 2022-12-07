## Manhattan plots for MOSTest project

#install.packages("qqman")
library(ggplot2)
library(data.table)
library(tidyverse)
path = "/path/"
  
# Load results
fMRI = fread(paste0(path,"rsfmri/rsfmri_eig0.most.orig.csv.gz"))
sMRI = fread(paste0(path,"T1/T1_eig16.most.orig.csv.gz"))
dMRI = fread(paste0(path,"diffusion/diffusion_eig10.most.orig.csv.gz"))
multi = fread(paste0(path,"multimodal/multimodal_mostest_eig26.most.orig.clump.indep.csv"))

gwas_data = fMRI[fMRI$PVAL<0.01,]
colnames(gwas_data) = c("chr","SNP","bp","A1","A2","N","p","Z")

## Prepare for plotting
data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)
gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)
axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))
ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)
sig <- 5e-8

## Plot
manhattan = ggplot(gwas_data, aes(x = bp_cum, y = -log10(p), 
                                  color = as_factor(chr), size = -log10(p))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave("/path/manhattan_fmri.png", manhattan)
