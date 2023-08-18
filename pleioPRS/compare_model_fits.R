#
library(data.table)
data = fread("/cluster/projects/p33/users/alexeas/elleke/pleioprs/src/multimodal_PRS_R2.SCZ_BIP_DEP_ASD_ADHD.200_two_thirds_samples.txt")


# Filter the data for pleioPGS and originalPGS separately
pleioPGS_data <- subset(data, PRS == "Disorder | Multimodal")
originalPGS_data <- subset(data, PRS == "Disorder")

# Create an empty vector to store the p-values
p_values <- as.data.frame(cbind(levels(factor(data$TRAIT)),c(0,0,0,0,0)))
t_values <- as.data.frame(cbind(levels(factor(data$TRAIT)),c(0,0,0,0,0)))

# Perform the Wilcoxon signed-rank test for each disorder type
for (disorder_type in levels(factor(data$TRAIT))) {
  # Subset the data for the current disorder type
  pleioPGS_r2 <- subset(pleioPGS_data, TRAIT == disorder_type)$R2
  originalPGS_r2 <- subset(originalPGS_data, TRAIT == disorder_type)$R2
  
  # Perform paired t-test
  test_result <- t.test(pleioPGS_r2,originalPGS_r2, paired = TRUE)
  
  # Store the p-value in the vector
  p_values[p_values$V1==disorder_type,2] <- test_result$p.value
  t_values[t_values$V1==disorder_type,2] <- test_result$statistic
}


