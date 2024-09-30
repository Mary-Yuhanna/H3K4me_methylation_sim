
# attempt to use permutation testing on genome-wide scores
#  set of pre-defined regions
# each region has a score (here a count) derived from two groups.
#  are any specific regions outliers based on global averages? 
#  Expect most regions to be similar or zero in both sets. 
#  ? test on numbers or test on difference in numbers?




# Strategy, calculate the difference in counts for a set of regions between two grouped samples
# can use the binomial function rbinom to generate count for each region for each treatment.
# include expected difference in count totals between treatments.
mean.WT <- 0.3    # based on observation that the regions are marked on average in 30% of wild-type samples/simulations
mean.KO <- 0.15   # based on observation that the regions are marked on average in 15% of wild-type samples/simulations


numberOfRegions <- 2000
# your data might be something like this:-
countsTable <- data.frame(count.WT = rbinom(n=numberOfRegions, size=100, prob = mean.WT),
                         count.KO =  rbinom(n=numberOfRegions, size=100, prob = mean.KO)   ) 


#  if the region scores are Random, with respect to the position, then differences will 
# follow a null distribution.
# let's generate the null distribution based on the difference of pairs of counts.
nullSize <- 20000    # make this a number ~10 larger than your number of regions
randScores.WT <-  rbinom(n=nullSize, size=100, prob = mean.WT)
randScores.KO <-  rbinom(n=nullSize, size=100, prob = mean.KO)

diffsNull <- randScores.WT- countsTable$count.KO
hist(diffsNull, breaks=50)   # in this case fits a normal distribution.  your data may not..

# can use the quantile function find what "diff" value cut's off a given % of values. 
# N.B. here for a one-tailed test where we only care about very high diffs
quantile(diffsNull, probs=c( 0.95, 0.99, 0.999))

# here a two-tailed set of values. where the "significant" region falls at both ends of the distribution
quantile(diffsNull, probs=c( 0.0005, 0.005, 0.025, 0.975, 0.995, 0.9995))

# output should look something like
# 0.05%   0.5%   2.5%  97.5%  99.5% 99.95% 
#   -5     -1      3     27     30     35

# a 5% alpha value, you would look for any regions showing a difference less than 3 or greater than 27
# at 1%, you would use -1 and 30, respectively


print(diffsData)

diffsData <- countsTable$count.WT - countsTable$count.KO
hist(diffsData, breaks=50)  # check this looks something like the null distribution. It might not.
# approx 5% of values in the real data should fall outside of these cut-offs.

length(diffsData) * 0.025
sum(diffsData <= 3)
sum(diffsData >= 27)

# if you have very many regions, you may want to go down to very low p-values
quantile(diffsNull, probs=c( 0.000005, 0.000005, 0.00005, 0.0005))  
# if any of the numbers are the same, you don't have any power from your null to go lower.









# Defining the number of regions for WT and KO - these numbers I have found with a different line of code. 
numberOfRegions_WT <- 2083168 + 450274 + 97494
numberOfRegions_KO <- 1053759 + 107112 + 11068

# I have defined a function to perform the analysis for a single state with modified plot titles and labels
analyze_state_with_labels <- function(mean_WT, mean_KO, me_layer_new, me_layer_ks, state_name, trend_desc, label) {
  # Calculating the widths (lengths) of methylated regions for WT and KO
  countsTable <- data.frame(
    count.WT = sapply(me_layer_new, function(x) sum(width(x))),
    count.KO = sapply(me_layer_ks, function(x) sum(width(x)))
  )
  
  # Defining the size of the null distribution
  nullSize <- 10000  # Increased null size for better accuracy - 10,000 sims
  
  # Generate random scores for the null distribution based on binomial distributions
  randScores.WT <- rbinom(n = nullSize, size = numberOfRegions_WT, prob = mean_WT)
  randScores.KO <- rbinom(n = nullSize, size = numberOfRegions_KO, prob = mean_KO)
  
  # Calculate the difference in the null distribution
  diffsNull <- randScores.WT - randScores.KO
  
  # Standardise the null distribution
  diffsNull <- scale(diffsNull)
  
  # Calculate the actual differences in the data
  diffsData <- countsTable$count.WT - countsTable$count.KO
  
  # Standardise the actual differences using the same scale as the null distribution
  diffsData <- (diffsData - mean(diffsData)) / sd(diffsData)
  
  # Calculate p-values for each observed difference - 
  # For each observed difference, the p-value is calculated as the proportion of null differences that are greater than or equal to the observed difference.
  p_values <- sapply(diffsData, function(obs_diff) {
    p_value <- mean(diffsNull >= obs_diff)
    return(p_value)
  })
  
  # Plot the histogram of the null distribution differences
  hist(diffsNull, breaks = 50, main = paste("Null Distribution of", state_name), 
       xlab = "Null Data", ylab = "Frequency", col = "lightgray")
  mtext(label[1], side = 3, line = 1, adj = -0.1, cex = 1.2)  # Adjusted cex to 1.2 for smaller labels
  
  # Plot the histogram of the actual data differences
  hist(diffsData, breaks = 50, main = paste("Positively Skewed Distribution of", state_name), 
       xlab = "Observed Data (WT - KO)", ylab = "Frequency", col = "lightgray")
  mtext(label[2], side = 3, line = 1, adj = -0.1, cex = 1.2)  # Adjusted cex to 1.2 for smaller labels
  
  # Print the results
  cat("Null distribution mean and sd for", state_name, ":", mean(diffsNull), sd(diffsNull), "\n")
  cat("P-values for observed differences in", state_name, ":", p_values, "\n")
  
  # Return both counts table and p-values for further inspection
  return(list(countsTable = countsTable, p_values = p_values))
}

# Set up a 3x2 plotting layout
par(mfrow = c(3, 2))

# Run the analysis for each state and generate the plots with labels using the function I created earlier. 
results_ME1 <- analyze_state_with_labels(mean_WT_ME1, mean_KO_ME1, me1_layer_new, me1_layer_ks, "ME1", "Trend Desc ME1", c("a)", "b)"))
results_ME2 <- analyze_state_with_labels(mean_WT_ME2, mean_KO_ME2, me2_layer_new, me2_layer_ks, "ME2", "Trend Desc ME2", c("c)", "d)"))
results_ME3 <- analyze_state_with_labels(mean_WT_ME3, mean_KO_ME3, me3_layer_new, me3_layer_ks, "ME3", "Trend Desc ME3", c("e)", "f)"))



print(me1_layer_ks)


#TAKE JUNG REIGONS - 58000, OVERLAPSANY FUNCTION ON ME1_LAYER_KO - WILL GIVE A VECTOR OF TRUE OR FALSE AND DO IT FOR ALL 100 REPLICATIONS, WHICH WILL GIVE 100 VECTORS OF TRUE/ FALSE, PUT THIS IN A MATRIX 
# OR CONVERT TRUES TO 1S AND FALSE FOR 0S, 
# TAKE EACH SIMULATION AND COMPARE TO THE JUNG DATA, SO FOR EXAMPLE, IN ME1_LAYER_NEW THERE ARE 100 SIMULATIONS, WE TAKE ONE SIMULATION AND THEIR GRANGES LIST AND COMPARE AGAINST JUNG. WE MAKE THE 58,000 MATRIX, THEN DO THE SAME FOR THE SECOND AND THIRD SIMULATION ETC.. EACH TIME THE NEW MATRIX IS MADE, WE ADD TO THE PREVIOUS MATRIX. IN THE END WE SHOULD    ARE SPECIFIC REIGONS DIFFERENT BETWEEN KNOCKOUT AND CONTROL 
# BOTTOM 1% 







library(GenomicRanges)
library(ggplot2)

# Function to calculate overlaps and accumulate counts in a matrix
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

# Generate cumulative matrices for wild-type and Kabuki syndrome datasets

# Wild type 
me1_wild_cumulative <- perform_cumulative_overlap_analysis(me1_layer_new, granges_true_data)
me2_wild_cumulative <- perform_cumulative_overlap_analysis(me2_layer_new, granges_true_data)
me3_wild_cumulative <- perform_cumulative_overlap_analysis(me3_layer_new, granges_true_data)

# Kabuki syndrome 
me1_ks_cumulative <- perform_cumulative_overlap_analysis(me1_layer_ks, granges_true_data)
me2_ks_cumulative <- perform_cumulative_overlap_analysis(me2_layer_ks, granges_true_data)
me3_ks_cumulative <- perform_cumulative_overlap_analysis(me3_layer_ks, granges_true_data)

# Calculate the difference in cumulative overlaps between WT and KS for each methylation state
diffs_me1 <- me1_wild_cumulative - me1_ks_cumulative
diffs_me2 <- me2_wild_cumulative - me2_ks_cumulative
diffs_me3 <- me3_wild_cumulative - me3_ks_cumulative





# Step 1: Calculate mean fractions for WT and KS for ME1, ME2, and ME3
mean_fraction_wt_me1 <- sum(me1_wild_cumulative) / (length(me1_wild_cumulative) * 100)
mean_fraction_ks_me1 <- sum(me1_ks_cumulative) / (length(me1_ks_cumulative) * 100)

mean_fraction_wt_me2 <- sum(me2_wild_cumulative) / (length(me2_wild_cumulative) * 100)
mean_fraction_ks_me2 <- sum(me2_ks_cumulative) / (length(me2_ks_cumulative) * 100)

mean_fraction_wt_me3 <- sum(me3_wild_cumulative) / (length(me3_wild_cumulative) * 100)
mean_fraction_ks_me3 <- sum(me3_ks_cumulative) / (length(me3_ks_cumulative) * 100)

print(numberOfRegions)

# Step 2: Generate counts and differences for ME1, ME2, and ME3
numberOfRegions <- length(me1_wild_cumulative) 


print(countsTable_me1)


# Generate the counts table based on the calculated mean fractions
countsTable_me1 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_me1),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_me1)
)

countsTable_me2 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_me2),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_me2)
)

countsTable_me3 <- data.frame(
  count.WT = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_wt_me3),
  count.KS = rbinom(n = numberOfRegions, size = 100, prob = mean_fraction_ks_me3)
)

# Step 3: Generate the null distributions
nullSize <- 40000  


# For ME1
randScores_WT_me1 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_me1)
randScores_KS_me1 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_me1)
null_diffs_me1 <- randScores_WT_me1 - randScores_KS_me1

# For ME2
randScores_WT_me2 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_me2)
randScores_KS_me2 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_me2)
null_diffs_me2 <- randScores_WT_me2 - randScores_KS_me2

# For ME3
randScores_WT_me3 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_wt_me3)
randScores_KS_me3 <- rbinom(n = nullSize, size = 100, prob = mean_fraction_ks_me3)
null_diffs_me3 <- randScores_WT_me3 - randScores_KS_me3

# Step 4: Calculate observed differences in your actual data
diffs_me1 <- me1_wild_cumulative - me1_ks_cumulative
diffs_me2 <- me2_wild_cumulative - me2_ks_cumulative
diffs_me3 <- me3_wild_cumulative - me3_ks_cumulative




# Step 5: Plot histograms of the observed differences and null distributions for comparison
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))  # Set up a 3x2 grid for plots and adjust margins

# ME1 Null Distribution
hist(null_diffs_me1, breaks = 50, col = "lightblue", main = "Null Distribution (ME1)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("a)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# ME1 Observed Differences
hist(diffs_me1, breaks = 50, col = "lightcoral", main = "Observed Differences (ME1)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("b)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# ME2 Null Distribution
hist(null_diffs_me2, breaks = 50, col = "lightblue", main = "Null Distribution (ME2)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("c)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# ME2 Observed Differences
hist(diffs_me2, breaks = 50, col = "lightcoral", main = "Observed Differences (ME2)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("d)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# ME3 Null Distribution
hist(null_diffs_me3, breaks = 50, col = "lightblue", main = "Null Distribution (ME3)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("e)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# ME3 Observed Differences
hist(diffs_me3, breaks = 50, col = "lightcoral", main = "Observed Differences (ME3)", 
     xlab = "Difference (WT - KS)", ylab = "Frequency")
mtext("f)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)

# Step 6: Identify outliers based on quantiles of the null distribution
cutoffs_me1 <- quantile(null_diffs_me1, probs = c(0.025, 0.975))
cutoffs_me2 <- quantile(null_diffs_me2, probs = c(0.025, 0.975))
cutoffs_me3 <- quantile(null_diffs_me3, probs = c(0.025, 0.975))

outliers_me1 <- which(diffs_me1 <= cutoffs_me1[1] | diffs_me1 >= cutoffs_me1[2])
outliers_me2 <- which(diffs_me2 <= cutoffs_me2[1] | diffs_me2 >= cutoffs_me2[2])
outliers_me3 <- which(diffs_me3 <= cutoffs_me3[1] | diffs_me3 >= cutoffs_me3[2])

cat("Number of outliers for ME1 (revised):", length(outliers_me1), "\n")
cat("Number of outliers for ME2 (revised):", length(outliers_me2), "\n")
cat("Number of outliers for ME3 (revised):", length(outliers_me3), "\n")

# Overlay of Observed and Null Distributions
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
plot_distribution_comparison(diffs_me1, null_diffs_me1, "Overlay of Observed and Null Distributions (ME1)", cutoffs_me1)
mtext("a)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)
plot_distribution_comparison(diffs_me2, null_diffs_me2, "Overlay of Observed and Null Distributions (ME2)", cutoffs_me2)
mtext("b)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)
plot_distribution_comparison(diffs_me3, null_diffs_me3, "Overlay of Observed and Null Distributions (ME3)", cutoffs_me3)
mtext("c)", side = 3, line = 0.5, adj = 0, cex = 1.2, font = 2)









# Identify regions where the difference between WT and KS in ME2 is below 0
protected_methylation_me2 <- which(diffs_me2 < 0)
# Get the corresponding genomic regions (GRanges) for protected methylation areas in ME2
protected_granges_me2 <- granges_true_data[protected_methylation_me2]

# Print or view the protected regions
print(protected_granges_me2)

# Save the protected methylation regions to a file (e.g., BED format or CSV)
library(rtracklayer)


# Alternatively, save as a CSV file for simpler viewing
protected_methylation_df <- as.data.frame(protected_granges_me2)
write.csv(protected_methylation_df, "/Users/maryyuhanna/desktop/protected_methylation_regions_me2.csv", row.names = FALSE)




# Find the minimum difference in ME2 (most negative value)
min_diff_me2 <- min(diffs_me2)

# Find the index of the region with the most negative difference
most_protected_index_me2 <- which(diffs_me2 == min_diff_me2)
# Get the corresponding GRanges for the most protected region in ME2
most_protected_region_me2 <- granges_true_data[most_protected_index_me2]

# Print the most protected region
print(most_protected_region_me2)

# Convert to a data frame and save as CSV
most_protected_region_df <- as.data.frame(most_protected_region_me2)
write.csv(most_protected_region_df, "/Users/maryyuhanna/desktop/most_protected_region_me2.csv", row.names = FALSE)






# Find the maximum difference in ME2 (largest positive value)
max_diff_me2 <- max(diffs_me2)

# Find the index/indices of the region(s) with the highest difference in ME2
highest_diff_index_me2 <- which(diffs_me2 == max_diff_me2)

# Get the corresponding GRanges for the region(s) with the highest difference
highest_diff_region_me2 <- granges_true_data[highest_diff_index_me2]

# Print the region(s) with the highest difference
print(highest_diff_region_me2)



# Optionally convert to data frame and save as CSV
highest_diff_region_me2_df <- as.data.frame(highest_diff_region_me2)
write.csv(highest_diff_region_me2_df, "/Users/maryyuhanna/desktop/highest_diff_region_me2.csv", row.names = FALSE)







































install.packages("cowplot")




# Ensure the necessary libraries are loaded
library(ggplot2)
library(cowplot)

# Plotting the cumulative counts for ME1, ME2, and ME3

# ME1
plot_me1 <- ggplot(plot_data_me1, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Mehtylation Counts (ME1)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -0.1, vjust = 1.5, size = 6)

# ME2
plot_me2 <- ggplot(plot_data_me2, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Mehtylation Counts (ME2)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "b)", hjust = -0.1, vjust = 1.5, size = 6)

# ME3
plot_me3 <- ggplot(plot_data_me3, aes(x = Counts, fill = Condition)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Raw Mehtylation Counts (ME3)", x = "Counts", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -Inf, y = Inf, label = "c)", hjust = -0.1, vjust = 1.5, size = 6)

# Combine the plots using cowplot
combined_plot <- plot_grid(plot_me1, plot_me2, plot_me3, ncol = 1)

# Display the combined plot
combined_plot
