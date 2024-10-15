# Genomic Layer Results Analysis
#
# This script loads and processes genomic layer results (H3K4me1, H3K4me2, H3K4me3) for both wild-type and KMT2D knockout simulations. 
# It calculates key metrics (Accuracy, Sensitivity, Specificity, Percentage) based on comparisons between simulation data 
# and Jung et al. data. The results are stored in data frames and prepared for further analysis.
#
# Key functionalities:
# - Load GRanges objects for different genomic layers (ME1, ME2, ME3) from RData files.
# - Convert Jung et al.'s data into a GRanges object for comparison.
# - Calculate confusion matrices for genomic ranges to assess overlaps between simulation data and true data.
# - Compute key performance metrics (Accuracy, Sensitivity, Specificity, Percentage).
# - Combine results from multiple layers into data frames for further analysis.
#
# Author: Mary Yuhanna
# ----------------------------------------------------

# Load necessary libraries
library(GenomicRanges)
library(dplyr)

# Set working directory (Update this path for your environment)
setwd("/Users/username/desktop/genomewide-prelim")

# ----------------------------------------------------
# Function to load GRanges objects from RData files matching a pattern
# ----------------------------------------------------
load_layer_results <- function(i) {
  thisLayerFiles <- dir(pattern = paste0(i, "_layer.RData"))
  layer_results <- GRangesList()
  
  for (thisFile in thisLayerFiles) {
    load(thisFile)
    if (exists("me1_layer")) {
      thisList <- GRangesList(me1_layer)
      rm(me1_layer)
    } else if (exists("me2_layer")) {
      thisList <- GRangesList(me2_layer)
      rm(me2_layer)
    } else if (exists("me3_layer")) {
      thisList <- GRangesList(me3_layer)
      rm(me3_layer)
    } else {
      stop("No expected layer object found in the loaded file.")
    }
    names(thisList) <- thisFile
    layer_results <- c(layer_results, thisList)
  }
  
  return(layer_results)
}

# ----------------------------------------------------
# Load ME1, ME2, and ME3 layer results
# ----------------------------------------------------
me1_layer_new <- load_layer_results(1)
me2_layer_new <- load_layer_results(2)
me3_layer_new <- load_layer_results(3)

# Load Jung Data and convert it to a GRanges object
true_data <- read.csv("/Users/maryyuhanna/desktop/Files for Kabuki/H3K4me2_results_KS.share copy.csv")
granges_true_data <- GRanges(seqnames = Rle(true_data$Chr),
                             ranges = IRanges(start = true_data$Start, end = true_data$End))

# ----------------------------------------------------
# Function to calculate the confusion matrix for genomic ranges
# ----------------------------------------------------
calcConfMat.GR <- function(query, subject, genomeSize) {
  try({
    strand(query) <- '*'
    strand(subject) <- '*'
    query <- reduce(query)
    subject <- reduce(subject)
    
    TP <- sum(width(GenomicRanges::intersect(query, subject)))
    FP <- (sum(width(query))) - TP
    FN <- (sum(width(subject))) - TP
    TN <- genomeSize - (TP + FP + FN)
    
    if (any(c(TP, FN, FP, TN) < 0)) {
      stop("Negative values in confusion matrix")
    }
    
    outMat <- matrix(c(TP, FN, FP, TN), ncol = 2)
    dimnames(outMat) <- list(subject = c("TRUE", "FALSE"), query = c("TRUE", "FALSE"))
    
    return(outMat)
  }, silent = TRUE)
}

# ----------------------------------------------------
# Function to calculate and store performance metrics
# ----------------------------------------------------
calculate_and_store_metrics <- function(layer_results, label) {
  local_results_df <- data.frame(Layer = character(), Percentage = numeric(), Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric(), stringsAsFactors = FALSE)
  
  for (layer_name in names(layer_results)) {
    print(paste("Processing layer:", layer_name))
    layer <- layer_results[[layer_name]]
    
    # Debugging: Check the structure of the layer
    print("Layer structure:")
    print(layer)
    
    try({
      conf_matrix <- calcConfMat.GR(query = layer, subject = granges_true_data, genomeSize = genome_size)
      print("Confusion matrix:")
      print(conf_matrix)
      
      TP <- conf_matrix[1, 1]
      FN <- conf_matrix[1, 2]
      FP <- conf_matrix[2, 1]
      TN <- conf_matrix[2, 2]
      
      # Accuracy
      acc <- function(TP, TN, FP, FN) {
        return((TP + TN) / (TP + FP + FN + TN))
      }
      
      # Sensitivity (TPR)
      tpr <- function(TP, TN, FP, FN) {
        return(TP / (TP + FN))
      }
      
      # Specificity (TNR)
      tnr <- function(TP, TN, FP, FN) {
        return(TN / (TN + FP))
      }
      
      accuracy <- acc(TP, TN, FP, FN)
      sensitivity <- tpr(TP, TN, FP, FN)
      specificity <- tnr(TP, TN, FP, FN)
      
      percentage <- sum(width(layer)) / genome_size
      
      print(paste("Accuracy:", accuracy))
      print(paste("Sensitivity:", sensitivity))
      print(paste("Specificity:", specificity))
      print(paste("Percentage:", percentage))
      
      local_results_df <- rbind(local_results_df, data.frame(Layer = label, Percentage = percentage, Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity))
    }, silent = FALSE)
  }
  
  return(local_results_df)
}

# ----------------------------------------------------
# Initialize data frame and calculate metrics for ME1, ME2, and ME3 layers
# ----------------------------------------------------
results_df <- data.frame(Layer = character(), Percentage = numeric(), Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric(), stringsAsFactors = FALSE)
genome_size <- 3.2e9  # Define the total genome size

# Calculate and store metrics for ME1, ME2, and ME3
results_df <- rbind(results_df, calculate_and_store_metrics(me1_layer_new, "ME1"))
results_df <- rbind(results_df, calculate_and_store_metrics(me2_layer_new, "ME2"))
results_df <- rbind(results_df, calculate_and_store_metrics(me3_layer_new, "ME3"))

# End of script for wild-type data
# ----------------------------------------------------

# ----------------------------------------------------
# Repeat the same for KMT2D knockout data
# ----------------------------------------------------

# Set working directory for KS data (Update path for your environment)
setwd("/Users/username/desktop/genomewide-prelim-ks")

# Load ME1, ME2, and ME3 layer results for KMT2D knockout
me1_layer_ks <- load_layer_results_ks(1)
me2_layer_ks <- load_layer_results_ks(2)
me3_layer_ks <- load_layer_results_ks(3)

# Initialize data frame for KS data
results_df_ks <- data.frame(Layer = character(), Percentage = numeric(), Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric(), stringsAsFactors = FALSE)

# Calculate and store metrics for ME1, ME2, and ME3 (KS)
results_df_ks <- rbind(results_df_ks, calculate_and_store_metrics(me1_layer_ks, "ME1"))
results_df_ks <- rbind(results_df_ks, calculate_and_store_metrics(me2_layer_ks, "ME2"))
results_df_ks <- rbind(results_df_ks, calculate_and_store_metrics(me3_layer_ks, "ME3"))
