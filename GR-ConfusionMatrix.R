# calcConfMat.GR
# This function calculates a confusion matrix for two sets of genomic ranges (GRanges objects),
# comparing the overlap between the "query" and "subject" sets. The original code was written 
# by Dr. Dave Gerrard at the University of Manchester in collaboration with Jeremy George, an 
# undergraduate project student (2022-23).
# 
# The function takes two sets of GRanges objects and computes the number of true positives (TP),
# false positives (FP), false negatives (FN), and true negatives (TN) based on the overlap
# between the two sets and the provided genome size.
#
# Parameters:
#   - query: GRanges object representing the predicted data (e.g., predictions)
#   - subject: GRanges object representing the true data (e.g., actual data)
#   - maxgap: Optional parameter (currently unused) for maximum gap allowed between ranges.
#   - minoverlap: Optional parameter (currently unused) for minimum overlap between ranges.
#   - genomeSize: The total size of the genome, which should be the sum of the lengths of all
#                 sequences in the genome, used to calculate the number of true negatives.
#
# Returns:
#   - A 2x2 confusion matrix with TP, FP, FN, and TN counts for the query and subject.
#   - Useful for statistical tests such as chi-square or for calculating performance metrics.
#
# Example usage:
#   calcConfMat.GR(query=gr1, subject=gr2, genomeSize=30)

calcConfMat.GR <- function(query, subject, maxgap = -1L, minoverlap = 0L, genomeSize) {
  
  # Ensure that strand information is ignored for comparison purposes
  strand(query) <- '*'
  strand(subject) <- '*'
  
  # Reduce each GRanges object to a non-overlapping set
  query <- reduce(query)
  subject <- reduce(subject)
  
  # Calculate true positives (TP) - overlapping regions between query and subject
  TP <- sum(width(intersect(query, subject)))
  
  # Calculate false positives (FP) - regions in query but not in subject
  FP <- sum(width(query)) - TP
  
  # Calculate false negatives (FN) - regions in subject but not in query
  FN <- sum(width(subject)) - TP
  
  # Calculate true negatives (TN) - remaining non-overlapping genome
  TN <- genomeSize - (TP + FP + FN)
  
  # Create the confusion matrix
  outMat <- matrix(c(TP, FN, FP, TN), ncol = 2)
  dimnames(outMat) <- list(subject = c("TRUE", "FALSE"), query = c("TRUE", "FALSE"))
  
  # Ensure there are no negative values, which indicates an incorrect genome size
  stopifnot(!any(outMat < 0))
  
  return(outMat)
}

# Example datasets for testing the function (from the GenomicRanges package)
# Usage of the function:
# calcConfMat.GR(query=gr, subject=gr1, genomeSize=40)

# Calculate the true positive rate (TPR) from the confusion matrix
tprFromConfMat <- function(confMatrix) {
  # TPR = TP / (TP + FN)
  return(confMatrix[1,1] / (confMatrix[1,1] + confMatrix[1,2]))
}

# Example usage:
# confMat <- calcConfMat.GR(gr1, gr3, genomeSize = 20)
# tpr <- tprFromConfMat(confMat)
# print(tpr)
