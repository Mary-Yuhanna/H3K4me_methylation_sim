install.packages("ggplot2")
library(ggplot2)


# Combined box plot for Percentage of ME1, ME2, and ME3
ggplot(results_df, aes(x = Layer, y = Percentage, fill = Layer)) +
  geom_boxplot() +
  labs(title = "Proportion of H3K4me marks in the hg38 genome for the Wild-Type simulation", x = "Layer", y = "Percentage") +
  theme_minimal() +
  ylim(0, 0.15)



# Box plot for Accuracy of ME1, ME2, and ME3
ggplot(results_df, aes(x = Layer, y = Accuracy, fill = Layer)) +
  geom_boxplot() +
  labs(title = "Model Accuracy for the prediction of of H3K4me marks in the hg38 genome for the Wild-Type simulation", x = "Layer", y = "Accuracy") +
  theme_minimal() +
  ylim(0, 1)

# Box plot for Sensitivity of ME1, ME2, and ME3
ggplot(results_df, aes(x = Layer, y = Sensitivity, fill = Layer)) +
  geom_boxplot() +
  labs(title = "Model Sensivity for the prediction of of H3K4me marks in the hg38 genome for the Wild-Type simulation", x = "Layer", y = "Sensitivity") +
  theme_minimal()

# Box plot for Specificity of ME1, ME2, and ME3
ggplot(results_df, aes(x = Layer, y = Specificity, fill = Layer)) +
  geom_boxplot() +
  labs(title = "Model Specificity for the prediction of of H3K4me marks in the hg38 genome for the Wild-Type simulation", x = "Layer", y = "Specificity") +
  theme_minimal() +
  ylim(0.9, 1)




























# Load necessary libraries
library(ggplot2)

# Add a column to distinguish between the simulations
results_df$Simulation <- "Wild-Type"
results_df_ks$Simulation <- "KMT2D Knockout"

# Combine the two dataframes
combined_results_df <- rbind(results_df, results_df_ks)

# Define custom colors
custom_colors <- c("Wild-Type" = "#1f77b4", "KMT2D Knockout" = "#d62728")

# Combined box plot for Percentage of ME1, ME2, and ME3
ggplot(combined_results_df, aes(x = Layer, y = Percentage, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Reduction in Global H3K4me marks in KMT2D Knockout Simulation", 
       x = "H3K4me Methylation State", y = "Percentage") +
  theme_minimal() +
  ylim(0, 0.15)

# Box plot for Accuracy of ME1, ME2, and ME3
ggplot(combined_results_df, aes(x = Layer, y = Accuracy, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "High Accuracy in Predicting Global H3K4me Marks with a Disproportionate Decrease in Accuracy for ME1 State, Most Pronounced in Wild-Type", 
       x = "H3K4me Methylation State", y = "Accuracy") +
  theme_minimal() 


# Box plot for Sensitivity of ME1, ME2, and ME3
ggplot(combined_results_df, aes(x = Layer, y = Sensitivity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Sensitivity in KMT2D Knockout for the Prediction of Global H3K4me marks", 
       x = "H3K4me Methylation State", y = "Sensitivity") +
  theme_minimal()

# Box plot for Specificity of ME1, ME2, and ME3
ggplot(combined_results_df, aes(x = Layer, y = Specificity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Specificity in Wild-Type for the Prediction of Global H3K4me marks", 
       x = "H3K4me Methylation State", y = "Specificity") +
  theme_minimal() 


















# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(grid) 

# Add a column to distinguish between the simulations
results_df$Simulation <- "Wild-Type"
results_df_ks$Simulation <- "KMT2D Knockout"

# Combine the two dataframes
combined_results_df <- rbind(results_df, results_df_ks)

# Define custom colors
custom_colors <- c("Wild-Type" = "#1f77b4", "KMT2D Knockout" = "#d62728")

# Define the a), b), c), d) annotations
label_a <- "a)"
label_b <- "b)"
label_c <- "c)"
label_d <- "d)"

# Plot a)
plot_a <- ggplot(combined_results_df, aes(x = Layer, y = Percentage, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Reduction in Global H3K4me Marks", 
       x = "H3K4me Methylation State", y = "Percentage") +
  theme_minimal() +
  ylim(0, 0.15) +
  annotate("text", x = 0.5, y = 0.14, label = label_a, size = 5, hjust = 0, vjust = 0)

# Plot b)
plot_b <- ggplot(combined_results_df, aes(x = Layer, y = Accuracy, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Decrease in Accuracy for ME1 State", 
       x = "H3K4me Methylation State", y = "Accuracy") +
  theme_minimal() +
  annotate("text", x = 0.5, y = 0.99, label = label_b, size = 5, hjust = 0, vjust = 0)

print(plot_b)


# Plot c)
plot_c <- ggplot(combined_results_df, aes(x = Layer, y = Sensitivity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Sensitivity in KMT2D Knockout", 
       x = "H3K4me Methylation State", y = "Sensitivity") +
  theme_minimal()

print(plot_c)


# Plot d)
plot_d <- ggplot(combined_results_df, aes(x = Layer, y = Specificity, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Higher Specificity in Wild-Type", 
       x = "H3K4me Methylation State", y = "Specificity") +
  theme_minimal()



# Arrange the plots in a 2x2 grid
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)


combined_plot <- grid.arrange(
  plot_a, plot_b, plot_c, plot_d, 
  nrow = 2
)



ggsave(filename = "/Users/maryyuhanna/desktop/combined_plots.png", plot = combined_plot, width = 12, height = 10, units = "in", dpi = 300)






# Plot b) with adjusted title size and spacing
plot_b <- ggplot(combined_results_df, aes(x = Layer, y = Accuracy, fill = Simulation, color = Simulation)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Disproportionate Decrease in Accuracy for ME1 State", 
       x = "H3K4me Methylation State", y = "Accuracy") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, margin = margin(b = 20))) +  # Adjust title size and margin
  annotate("text", x = 0.5, y = 0.99, label = label_b, size = 5, hjust = 0, vjust = 0)

# Recreate the combined plot layout
combined_plot <- grid.arrange(
  plot_a, plot_b, plot_c, plot_d, 
  nrow = 2
)

# Save the combined plot without the title cutting off
ggsave(filename = "/Users/maryyuhanna/desktop/combined_plots_adjusted_10.png", plot = combined_plot, width = 12, height = 10, units = "in", dpi = 300)


