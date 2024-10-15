# This script performs a cumulative overlap analysis between wild-type (WT) and Kabuki Syndrome (KS) 
# genomic datasets across multiple simulations, focusing on three different methylation states (H3K4me1, H3K4me2, and H3K4me3).
# 
# Key functionalities include:
# - Accumulating overlap counts across multiple simulations
# - Calculating differences in methylation patterns between WT and KS datasets
# - Generating null distributions to assess the significance of observed differences
# - Identifying outliers and regions of interest (protected or with the highest differences)
# - Plotting cumulative counts and comparing observed vs. null distributions using ggplot2
# 
# The script also includes functionality to save the results as CSV files and visualize the results 
# through histograms and combined plots for easier interpretation.

# Author: Mary Yuhanna


# Load required libraries
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(cowplot)


## Accumulating overlap counts across multiple simulations
# Function to calculate overlaps and accumulate counts in a matrix across simulations

perform_cumulative_overlap_analysis <- function(simulation_list, granges_true_data, n_simulations = 100) {
  n_jung_regions <- length(granges_true_data)
  cumulative_matrix <- matrix(0, nrow = n_jung_regions, ncol = 1) 
  
  for (i in seq_len(n_simulations)) {
    simulation_gr <- simulation_list[[i]]
    overlap_vector <- overlapsAny(granges_true_data, simulation_gr)
    cumulative_matrix[, 1] <- cumulative_matrix[, 1] + as.integer(overlap_vector)
  }
  
  return(cumulative_matrix)
}

# Generate cumulative matrices for wild-type (WT) and Kabuki syndrome (KS) datasets
H3K4me1_wild_cumulative <- perform_cumulative_overlap_analysis(H3K4me1_layer_new, granges_true_data)
H3K4me2_wild_cumulative <- perform_cumulative_overlap_analysis(H3K4me2_layer_new, granges_true_data)
H3K4me3_wild_cumulative <- perform_cumulative_overlap_analysis(H3K4me3_layer_new, granges_true_data)

H3K4me1_ks_cumulative <- perform_cumulative_overlap_analysis(H3K4me1_layer_ks, granges_true_data)
H3K4me2_ks_cumulative <- perform_cumulative_overlap_analysis(H3K4me2_layer_ks, granges_true_data)
H3K4me3_ks_cumulative <- perform_cumulative_overlap_analysis(H3K4me3_layer_ks, granges_true_data)

# Calculate differences between WT and KS for H3K4me1, H3K4me2, and H3K4me3
diffs_H3K4me1 <- H3K4me1_wild_cumulative - H3K4me1_ks_cumulative
diffs_H3K4me2 <- H3K4me2_wild_cumulative - H3K4me2_ks_cumulative
diffs_H3K4me3 <- H3K4me3_wild_cumulative - H3K4me3_ks_cumulative

# Calculate mean fractions for WT and KS in each methylation state
mean_fraction_wt_H3K4me1 <- sum(H3K4me1_wild_cumulative) / (length(H3K4me1_wild_cumulative) * 100)
mean_fraction_ks_H3K4me1 <- sum(H3K4me1_ks_cumulative) / (length(H3K4me1_ks_cumulative) * 100)

mean_fraction_wt_H3K4me2 <- sum(H3K4me2_wild_cumulative) / (length(H3K4me2_wild_cumulative) * 100)
mean_fraction_ks_H3K4me2 <- sum(H3K4me2_ks_cumulative) / (length(H3K4me2_ks_cumulative) * 100)

mean_fraction_wt_H3K4me3 <- sum(H3K4me3_wild_cumulative) / (length(H3K4_wild_cumulative) * 100)
mean_fraction_ks_H3K4me3 <- sum(H3K4me3_ks_cumulative) / (length(H3K4me3_ks_cumulative) * 100)

# Generate counts tables for WT and KS based on mean fractions
numberOfRegions <- length(H3K4me1_wild_cumulative)

countsTable_H3K4me1 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_H3K4me1),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_H3K4me1)
)

countsTable_H3K4me2 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_H3K4me2),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_H3K4me2)
)

countsTable_H3K4me3 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_H3K4me3),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_H3K4me3)
)




## Generating null distributions to assess the significance of observed differences
# Generate null distributions for WT and KS

nullSize <- 40000

null_diffs_H3K4me1 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_H3K4me1) -
  rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_H3K4me1)

null_diffs_H3K4me2 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_H3K4me2) -
  rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_H3K4me2)

null_diffs_H3K4me3 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_H3K4me3) -
  rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_H3K4me3)

# Plot histograms of observed differences and null distributions
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

# H3K4me1
hist(null_diffs_H3K4me1, breaks = 50, col = "lightblue", main = "Null Distribution (H3K4me1)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("a)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

hist(diffs_H3K4me1, breaks = 50, col = "lightcoral", main = "Observed Differences (H3K4me1)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("b)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# H3K4me2
hist(null_diffs_H3K4me2, breaks = 50, col = "lightblue", main = "Null Distribution (H3K4me2)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("c)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

hist(diffs_H3K4me2, breaks = 50, col = "lightcoral", main = "Observed Differences (H3K4me2)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("d)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# H3K4me3
hist(null_diffs_H3K4me3, breaks = 50, col = "lightblue", main = "Null Distribution (H3K4me3)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("e)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

hist(diffs_H3K4me3, breaks = 50, col = "lightcoral", main = "Observed Differences (H3K4me3)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("f)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# Identify outliers based on the null distribution quantiles
cutoffs_H3K4me1 <- quantile(null_diffs_H3K4me1, probs = c(0.025, 0.975))
cutoffs_H3K4me2 <- quantile(null_diffs_H3K4me2, probs = c(0.025, 0.975))
cutoffs_H3K4me3 <- quantile(null_diffs_H3K4me3, probs = c(0.025, 0.975))

outliers_H3K4me1 <- which(diffs_H3K4me1 <= cutoffs_H3K4me1[1] | diffs_H3K4me1 >= cutoffs_H3K4me1[2])
outliers_H3K4me2 <- which(diffs_H3K4me2 <= cutoffs_H3K4me2[1] | diffs_H3K4me2 >= cutoffs_H3K4me2[2])
outliers_H3K4me3 <- which(diffs_H3K4me3 <= cutoffs_H3K4me3[1] | diffs_H3K4me3 >= cutoffs_H3K4me3[2])

cat("Number of outliers for H3K4me1:", length(outliers_H3K4me1), "\n")
cat("Number of outliers for H3K4me2:", length(outliers_H3K4me2), "\n")
cat("Number of outliers for H3K4me3:", length(outliers_H3K4me3), "\n")

# Save protected methylation regions to CSV
protected_methylation_H3K4me2 <- which(diffs_H3K4me2 < 0)
protected_granges_H3K4me2 <- granges_true_data[protected_methylation_H3K4me2]
write.csv(as.data.frame(protected_granges_H3K4me2), "protected_methylation_regions_H3K4me2.csv", row.names = FALSE)

# Save regions with max/min differences in H3K4me2
most_protected_index_H3K4me2 <- which(diffs_H3K4me2 == min(diffs_H3K4me2))
highest_diff_index_H3K4me2 <- which(diffs_H3K4me2 == max(diffs_H3K4me2))

most_protected_region_H3K4me2 <- granges_true_data[most_protected_index_H3K4me2]
highest_diff_region_H3K4me2 <- granges_true_data[highest_diff_index_H3K4me2]

write.csv(as.data.frame(most_protected_region_H3K4me2), "most_protected_region_H3K4me2.csv", row.names = FALSE)
write.csv(as.data.frame(highest_diff_region_H3K4me2), "highest_diff_region_H3K4me2.csv", row.names = FALSE)

# Plot cumulative counts for H3K4me1, H3K4me2, and H3K4me3 using ggplot2
plot_H3K4me1 <- ggplot(plot_data_H3K4me1, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Methylation Counts (H3K4me1)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -0.1, vjust = 1.5, size = 6)

# Plot for H3K4me2
plot_H3K4me2 <- ggplot(plot_data_H3K4me2, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Methylation Counts (H3K4me2)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "b)", hjust = -0.1, vjust = 1.5, size = 6)

# Plot for H3K4me3
plot_H3K4me3 <- ggplot(plot_data_H3K4me3, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Methylation Counts (H3K4me3)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "c)", hjust = -0.1, vjust = 1.5, size = 6)

# Combine the plots for H3K4me1, H3K4me2, and H3K4me3 using cowplot
combined_plot <- plot_grid(plot_H3K4me1, plot_H3K4me2, plot_H3K4me3, ncol = 1)

# Display the combined plot
combined_plot




##Identifying outliers and regions of interest (protected or with the highest differences)
# Identify regions where the difference between WT and KS in H3K4me2 is below 0

protected_methylation_H3K4me2 <- which(diffs_H3K4me2 < 0)

# Get the corresponding genomic regions (GRanges) for protected methylation areas in H3K4me2
protected_granges_H3K4me2 <- granges_true_data[protected_methylation_H3K4me2]

# Print the protected regions
print(protected_granges_H3K4me2)


# Alternatively, save as a CSV file for simpler viewing
protected_methylation_df <- as.data.frame(protected_granges_H3K4me2)
write.csv(protected_methylation_df, "/Users/username/desktop/protected_methylation_regions_H3K4me2.csv", row.names = FALSE)

# Find the minimum difference in H3K4me2 (most negative value)
min_diff_H3K4me2 <- min(diffs_H3K4me2)

# Find the index of the region with the most negative difference
most_protected_index_H3K4me2 <- which(diffs_H3K4me2 == min_diff_H3K4me2)

# Get the corresponding GRanges for the most protected region in H3K4me2
most_protected_region_H3K4me2 <- granges_true_data[most_protected_index_H3K4me2]

# Print the most protected region
print(most_protected_region_H3K4me2)

# Convert to a data frame and save as CSV
most_protected_region_df <- as.data.frame(most_protected_region_H3K4me2)
write.csv(most_protected_region_df, "/Users/username/desktop/most_protected_region_H3K4me2.csv", row.names = FALSE)

# Find the maximum difference in H3K4me2 (largest positive value)
max_diff_H3K4me2 <- max(diffs_H3K4me2)

# Find the index/indices of the region(s) with the highest difference in H3K4me2
highest_diff_index_H3K4me2 <- which(diffs_H3K4me2 == max_diff_H3K4me2)

# Get the corresponding GRanges for the region(s) with the highest difference
highest_diff_region_H3K4me2 <- granges_true_data[highest_diff_index_H3K4me2]

# Print the regions with the highest difference
print(highest_diff_region_H3K4me2)

# Optionally convert to a data frame and save as CSV
highest_diff_region_H3K4me2_df <- as.data.frame(highest_diff_region_H3K4me2)
write.csv(highest_diff_region_H3K4me2_df, "/Users/username/desktop/highest_diff_region_H3K4me2.csv", row.names = FALSE)




##Plotting cumulative counts and comparing observed vs. null distributions
# Overlay of observed and null distributions for H3K4me1, H3K4me2, and H3K4me3

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

# Custom function for plotting observed vs. null distributions
plot_distribution_comparison <- function(observed, null_dist, title, cutoffs) {
  hist(null_dist, breaks = 50, col = "lightblue", main = title, xlab = "Difference (WT - KS)", ylab = "Frequency")
  hist(observed, breaks = 50, col = rgb(1, 0, 0, alpha = 0.5), add = TRUE)
  abline(v = cutoffs, col = "red", lty = 2, lwd = 2)
}

# Plot overlay of observed and null distributions for H3K4me1
plot_distribution_comparison(diffs_H3K4me1, null_diffs_H3K4me1, "Overlay of Observed and Null Distributions (H3K4me1)", cutoffs_H3K4me1)
mtext("a)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# Plot overlay of observed and null distributions for H3K4me2
plot_distribution_comparison(diffs_H3K4me2, null_diffs_H3K4me2, "Overlay of Observed and Null Distributions (H3K4me2)", cutoffs_H3K4me2)
mtext("b)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# Plot overlay of observed and null distributions for H3K4me3
plot_distribution_comparison(diffs_H3K4me3, null_diffs_H3K4me3, "Overlay of Observed and Null Distributions (H3K4me3)", cutoffs_H3K4me3)
mtext("c)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)









