# Genomic Simulation Results Analysis
#
# This script reads and combines two CSV files containing simulation results from wild-type and KMT2D knockout experiments.
# It performs two-way ANOVA on several metrics (Specificity, Sensitivity, Accuracy, and Percentage) to analyze the 
# differences between genomic layers and simulation types (wild-type vs knockout).
# 
# Key functionalities:
# - Combine simulation results from wild-type and KMT2D knockout into a single dataframe.
# - Perform two-way ANOVA to evaluate the interaction between genomic layers and simulation types.
# 
# Author: Mary Yuhanna

# Reading the two CSV files
wild_type <- read.csv("/Users/username/desktop/simulation_results_ks_100.csv")
ks_sim <- read.csv("/Users/username/desktop/simulation_results_ksKMT2D_100.csv")

# Combine the two dataframes if they have the same structure
combined_results_df <- rbind(wild_type, ks_sim)

# Check the structure of the combined dataframe
str(combined_results_df)

# View the first few rows
head(combined_results_df)

# Add a 'Simulation' column to identify the simulation type
wild_type$Simulation <- "WildType"
ks_sim$Simulation <- "KMT2D_Knockout"

# Combine the dataframes again with the new 'Simulation' column
combined_results_df <- rbind(wild_type, ks_sim)

# Check the structure and first few rows of the combined dataframe
str(combined_results_df)
head(combined_results_df)

# Perform two-way ANOVA for various metrics
# ANOVA for Specificity
specificity_anova <- lm(Specificity ~ Layer * Simulation, data = combined_results_df)
summary(specificity_anova)

# ANOVA for Sensitivity
sensitivity_anova <- lm(Sensitivity ~ Layer * Simulation, data = combined_results_df)
summary(sensitivity_anova)

# ANOVA for Accuracy
accuracy_anova <- lm(Accuracy ~ Layer * Simulation, data = combined_results_df)
summary(accuracy_anova)

# ANOVA for Percentage
percentage_anova <- lm(Percentage ~ Layer * Simulation, data = combined_results_df)
summary(percentage_anova)

