# Performance Metrics for Model Evaluation
# 
# This script provides various performance metrics commonly used in classification tasks 
# such as accuracy, sensitivity, specificity, precision, and more. These functions can be 
# used to evaluate classification models, particularly for genomic data or other types of 
# bioinformatics-related analyses. The script includes functions for calculating confusion 
# matrix-based measures like Matthews correlation coefficient (MCC), false positive rate 
# (FPR), and false discovery rate (FDR), among others.
#
# Key functionalities:
# - Calculates a variety of evaluation metrics (accuracy, sensitivity, specificity, etc.)
# - Implements Matthews correlation coefficient, suitable for balanced datasets
# - Defines a custom scoring function 'score.hits' for measuring performance in overlap 
#   analysis between two IRanges objects (query and target)
# - Supports further extension for genomic range comparisons
#
# Author: Dr Dave Gerrard

## Function to calculate accuracy
acc <- function(TP, TN, FP, FN) {
  return((TP + TN) / (TP + FP + FN + TN))
}

# Example: acc(50, 50, 3, 20)

## Function to calculate Matthews correlation coefficient
matthews.cor <- function(TP, TN, FP, FN) {
  return((TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
}

# Example: matthews.cor(50, 50, 3, 20)

## Function to calculate True Positive Rate (Sensitivity)
tpr <- function(TP, TN, FP, FN) {
  return(TP / (TP + FN))
}

# Example: tpr(50, 50, 3, 20)

## Function to calculate True Negative Rate (Specificity)
tnr <- function(TP, TN, FP, FN) {
  return(TN / (TN + FP))
}

# Example: tnr(50, 50, 3, 20)

## Function to calculate Positive Predictive Value (Precision)
ppv <- function(TP, TN, FP, FN) {
  return(TP / (TP + FP))
}

# Example: ppv(50, 50, 3, 20)

## Function to calculate Negative Predictive Value
npv <- function(TP, TN, FP, FN) {
  return(TN / (TN + FN))
}

# Example: npv(50, 50, 3, 20)

## Function to calculate False Positive Rate
fpr <- function(TP, TN, FP, FN) {
  return(FP / (FP + TN))
}

# Example: fpr(50, 50, 3, 20)

## Function to calculate False Negative Rate
fnr <- function(TP, TN, FP, FN) {
  return(FN / (TP + FN))
}

# Example: fnr(50, 50, 3, 20)

## Function to calculate False Discovery Rate
fdr <- function(TP, TN, FP, FN) {
  return(FP / (TP + FP))
}

# Example: fdr(50, 50, 3, 20)

## Custom scoring function for genomic range comparisons
# This function calculates performance metrics (e.g., accuracy, FDR) for overlap analysis
# between query and target genomic ranges using IRanges objects.
score.hits <- function(query, target, method = acc) {
  FP <- sum(width(setdiff(query, target)))
  TP <- sum(width(intersect(query, target)))
  FN <- sum(width(setdiff(target, query)))
  TN <- sum(width(intersect(gaps(query), gaps(target))))
  
  score <- do.call(method, args = list(TP = TP, TN = TN, FP = FP, FN = FN))
  return(score)
}

# Example usage:
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR)  # Default method is accuracy
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method = fdr)
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method = matthews.cor)

## Calculate specificity from a confusion matrix
specFromConfMat <- function(confMatrix) {
  return(confMatrix[2, 2] / (confMatrix[2, 2] + confMatrix[2, 1]))
}

## Calculate sensitivity from a confusion matrix
sensFromConfMat <- function(confMatrix) {
  return(confMatrix[1, 1] / (confMatrix[1, 1] + confMatrix[1, 2]))
}
