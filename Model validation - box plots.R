# Visualisation of Genomic Simulation Results
#
# This script generates a series of boxplots to compare the results of wild-type (WT) and KMT2D knockout (KS) genomic simulations.
# The plots include comparisons of Percentage, Accuracy, Sensitivity, and Specificity across various H3K4me methylation states.
#
# Key functionalities:
# - Combine simulation results from WT and KS into a single dataframe.
# - Create and customise boxplots to visualise differences between WT and KS for each metric.
# - Arrange the plots into a grid for visual comparison.
# - Save the combined plots as a PNG file.
#
# Author: Mary Yuhanna
# ----------------------------------------------------

# Load necessary libraries
install.packages("ggplot2")
library(ggplot2)
library(gridExtra)
library(grid)

# ----------------------------------------------------
# Data Preparation: Combine WT and KS simulation results
# ----------------------------------------------------

# Add a column to distinguish between the simulations
results_df$Simulation <- "Wild-Type"
results_df_ks$Simulation <- "KMT2D Knockout"

# Combine the two dataframes using rbind
combined_results_df <- rbind(results_df, results_df_ks)

# ----------------------------------------------------
# Visualisation: Create boxplots comparing WT and KS for multiple metrics
# ----------------------------------------------------

# Define custom colors for the simulations
custom_colors <- c("Wild-Type" = "#1f77b4", "KMT2D Knockout" = "#d62728")

# Define the a), b), c), d) annotations for the plots
label_a <- "a)"
label_b <- "b)"
label_c <- "c)"
label_d <- "d)"

# Plot a: Percentage comparison
plot_a <- ggplot(combined_results_df, aes(x = Layer, y = Percentage, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Reduction in Global H3K4me Marks", 
       x = "H3K4me Methylation State", y = "Percentage") +
  theme_minimal() +
  ylim(0, 0.15) +
  annotate("text", x = 0.5, y = 0.14, label = label_a, size = 5, hjust = 0, vjust = 0)

# Plot b: Accuracy comparison
plot_b <- ggplot(combined_results_df, aes(x = Layer, y = Accuracy, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Decrease in Accuracy for ME1 State", 
       x = "H3K4me Methylation State", y = "Accuracy") +
  theme_minimal() +
  annotate("text", x = 0.5, y = 0.99, label = label_b, size = 5, hjust = 0, vjust = 0)

# Plot c: Sensitivity comparison
plot_c <- ggplot(combined_results_df, aes(x = Layer, y = Sensitivity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Sensitivity in KMT2D Knockout", 
       x = "H3K4me Methylation State", y = "Sensitivity") +
  theme_minimal()

# Plot d: Specificity comparison
plot_d <- ggplot(combined_results_df, aes(x = Layer, y = Specificity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Specificity in Wild-Type", 
       x = "H3K4me Methylation State", y = "Specificity") +
  theme_minimal()

# ----------------------------------------------------
# Plot Arrangement: Arrange the boxplots in a 2x2 grid
# ----------------------------------------------------

# Arrange the plots in a 2x2 grid
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)

# Save the combined plot as a PNG file
ggsave(filename = "/Users/username/desktop/combined_plots.png", 
       plot = grid.arrange(plot_a, plot_b, plot_c, plot_d, nrow = 2), 
       width = 12, height = 10, units = "in", dpi = 300)

