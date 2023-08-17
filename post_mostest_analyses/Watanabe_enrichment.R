library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

###### LOAD DATA #######
### 1. Use Watanabe et al 2019 as reference lead SNPs of 480 traits
ref_data <- fread("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/Watanabe_Table2.txt", header = T,dec = ",")[,1:3]
ref_data$LeadSNPs_count = round(ref_data$LeadSNPs_prop*43492) # Total amount of lead SNPs mentioned in the paper
ref_data$Genome_count = round(ref_data$Genome_prop*1740179) # Total amount of SNPs mentioned in the paper
colnames(ref_data)[1] = "Category"

#### 2. Cross-modality and single-modality lead SNPs with ANNOVAR from FUMA
crossmodal = fread("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/bonf3/crossmodality_loci_Bonf3.txt")
colnames(crossmodal)[3] = "LEADSNP"
smri = fread("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/bonf3/T1_specific_loci_Bonf3.txt")
dmri = fread("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/bonf3/diffusion_specific_loci_Bonf3.txt")
fmri = fread("/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/Pleiotropy/Loci/bonf3/rsfmri_specific_loci_Bonf3.txt")

# Add uniqueID
rs_to_uniq = fread("/Users/elleke/Downloads/multimodal_eig16.most.orig_uniqID.csv.gz")
colnames(rs_to_uniq)[9] = "uniqID"
crossmodal_merged = merge(crossmodal,rs_to_uniq[,c(1,2,3,9)],by.x="LEADSNP",by.y="SNP")
smri_merged = merge(smri,rs_to_uniq[,c(1,2,3,9)],by.x="LEADSNP",by.y="SNP")
dmri_merged = merge(dmri,rs_to_uniq[,c(1,2,3,9)],by.x="LEADSNP",by.y="SNP")
fmri_merged = merge(fmri,rs_to_uniq[,c(1,2,3,9)],by.x="LEADSNP",by.y="SNP")

# Retrieve ANNOVAR annotations from FUMA jobs
annovar <- fread("/Users/elleke/Downloads/FUMA_job273994/annov.txt")
crossmodal_merged = unique(merge(crossmodal_merged,annovar[,c(1,3)],by=c("uniqID")))
smri_merged = unique(merge(smri_merged,annovar[,c(1,3)],by=c("uniqID")))
dmri_merged = unique(merge(dmri_merged,annovar[,c(1,3)],by=c("uniqID")))
fmri_merged = unique(merge(fmri_merged,annovar[,c(1,3)],by=c("uniqID")))

# Merge data
crossmodal_merged$Dataset <- "Cross-modality lead SNPs"
smri_merged$Dataset <- "sMRI-modality lead SNPs"
dmri_merged$Dataset <- "dMRI-modality lead SNPs"
fmri_merged$Dataset <- "fMRI-modality lead SNPs"
main_df <- rbind(crossmodal_merged, smri_merged, dmri_merged, fmri_merged)[,c(2,9,10)]
colnames(main_df) = c("SNP","Category","Dataset")

# Get counts per modality per VEP category
category_counts_per_dataset <- main_df %>%
  dplyr::group_by(Category, Dataset) %>%
  dplyr::summarise(SNPs = dplyr::n()) %>%
  tidyr::pivot_wider(names_from = Dataset, values_from = SNPs, values_fill = 0)
category_counts_per_dataset <- rbind(category_counts_per_dataset, colSums(category_counts_per_dataset[, -1]))
category_counts_per_dataset[10,1] = "Total"

# Left join category_counts_per_dataset with average_ref_counts
combined_counts <- left_join(ref_data[,c("Category","LeadSNPs_count")], category_counts_per_dataset, by = "Category") 
combined_counts[is.na(combined_counts)] <- 0

#### Test for enrichment ###
# We use
#    AverageCount: How many times did we count this annotation for the reference panel?
#    count: How many times did we count this annotation for our selected variants?
#    N: How many times total have we counted any annotation for the reference data?
#    n: How many times total have we counted any annotation for our selected variants?
#             | leadsnps | refsnps
#    curanno  |   x[2]   |    x[1]     |   
#    othanno  |          |             |
#    allanno  |    n     |     N       |    

N <- combined_counts$LeadSNPs_count[12]

# Loop over the columns
result_list <- list()
for (col_name in colnames(combined_counts)[3:6]) { 
  n <- combined_counts[[col_name]][12]
  
  # Apply the calcFisherP function to each row
  result <- combined_counts %>%
    rowwise() %>%
    mutate(modality = col_name,
           obs_prop = .data[[col_name]]/n, 
           ref_prop = LeadSNPs_count/N,
           OR = fisher.test(matrix(c(.data[[col_name]],n-.data[[col_name]],LeadSNPs_count,N-LeadSNPs_count),ncol=2))$estimate,
           P = fisher.test(matrix(c(.data[[col_name]],n-.data[[col_name]],LeadSNPs_count,N-LeadSNPs_count),ncol=2))$p.value,
           CI_lower = fisher.test(matrix(c(.data[[col_name]],n-.data[[col_name]],LeadSNPs_count,N-LeadSNPs_count),ncol=2))$conf.int[1],
           CI_upper = fisher.test(matrix(c(.data[[col_name]],n-.data[[col_name]],LeadSNPs_count,N-LeadSNPs_count),ncol=2))$conf.int[2]
    )
  
  result_list[[col_name]] <- result
}

# Combine the results into a single data frame
enrichment_results <- bind_rows(result_list)
enrichment_results = enrichment_results[enrichment_results$Category!="Total",]



#### PLOTTING #####
## Prepare for plotting
enrichment_results$Category = gsub("_"," ",enrichment_results$Category)
enrichment_results$Bonf = "P < 0.05"
enrichment_results[enrichment_results$P>0.05,"Bonf"] = "P > 0.05"
enrichment_results[enrichment_results$P<(0.05/nrow(enrichment_results)),"Bonf"] = "P < 0.05/44"
enrichment_results$Bonf = factor(enrichment_results$Bonf,levels=c("P > 0.05","P < 0.05","P < 0.05/44"))

enrichment_results$enriched = "OR = 1"
enrichment_results[enrichment_results$Bonf == "P < 0.05/44" & enrichment_results$OR > 1,"enriched"] = "OR > 1"
enrichment_results[enrichment_results$Bonf == "P < 0.05/44" & enrichment_results$OR < 1,"enriched"] = "OR < 1"

enrichment_results$modality = factor(enrichment_results$modality,levels=c("Cross-modality lead SNPs",
                                                                          "sMRI-modality lead SNPs",
                                                                          "dMRI-modality lead SNPs",
                                                                          "fMRI-modality lead SNPs"))
# Plot
ggplot(enrichment_results, aes(x=obs_prop,y=factor(modality),fill=factor(Category),color=factor(enriched),linetype=factor(enriched))) +
  geom_bar(position="fill",stat="identity") +
  scale_fill_manual(values=c("#8b4513","#ffd700","#4682b4","#eee8aa","#2f4f4f","#228b22","#cd853f","#00008b","#1e90ff",
                             "#ffff54", "#ff69b4")) +
  scale_color_manual(values=c("black","#FFFFFF00","black")) +
  scale_linetype_manual(values = c("dotted", "solid","solid"),) +
  theme(
    aspect.ratio = 4/3,
    axis.text.y=element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(size = 0),
  ) +  
  theme_bw() +
  labs(y="",x="Proportion of lead SNPs")+
  coord_flip()

