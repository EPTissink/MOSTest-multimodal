# Load required library
library(data.table)
library(ggplot2)
library(gridExtra)

# Create an empty list to store the data.frames
data_list <- list()

# Get a list of all directories in the specified path
path_to_directories <- "/Users/elleke/SURFdrive/Projects/MOSTest/Results/Discovery/cFDR/pleioprs"
dirs <- list.dirs(path_to_directories, recursive = FALSE)

# Loop through each directory
for (dir in dirs) {
  # Check if the directory name contains "FUMA_gene2func"
  if (grepl("FUMA_gene2func", dir)) {
    # Extract disorder and condition from the directory name
    dir_parts <- strsplit(basename(dir), "_")[[1]]
    disorder <- dir_parts[3]
    condition <- dir_parts[4]
    
    # Create the file path for gtex_v8_ts_DEG.txt
    file_path <- file.path(dir, "gtex_v8_ts_DEG.txt")
    
    # Read the data from the file using fread (assuming the file is tab-separated)
    if (file.exists(file_path)) {
      data <- fread(file_path, sep = "\t", header = TRUE, data.table = FALSE)
      
      # Add columns for disorder and condition
      data$disorder <- disorder
      data$condition <- condition
      
      # Append the data.frame to the list
      data_list[[length(data_list) + 1]] <- data
    }
  }
}

# Combine all data.frames into one data.table
combined_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

# Preprocess data to create "Brain tissue" and "Body tissue" factor column
combined_data$Tissue <- ifelse(grepl("^Brain", combined_data$GeneSet), "Brain tissue", "Body tissue")
combined_data$GeneSet <- gsub("_", " ", combined_data$GeneSet)
combined_data = combined_data[combined_data$Category == "DEG.up",]
df_filtered <- combined_data[order(-log10(combined_data$p), combined_data$Tissue, combined_data$GeneSet,decreasing = F), ]
df_filtered$condition = factor(df_filtered$condition,levels=c("original","condfdr"))

ggplot(combined_data, aes(x = reorder(GeneSet, log10(p)), y = -log10(p), fill = disorder)) +
  geom_bar(stat = "identity", position="dodge", width = 1) +
  facet_grid(factor(condition, levels=c("originalGWAS", "condFDR"))~factor(Tissue, levels=c("Brain tissue", "Body tissue")), scales = "free_x", space = "free") +
  labs( y = "-log10(p)", fill = "Disorder") +
  scale_fill_manual(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position = "top") +
  geom_hline(yintercept = -log10(0.05/540), linetype = "dashed", color = "red")
