# Genomic Ranges Filtering Function for GenomicLayers
#
# This script defines a function to filter genomic ranges (GRanges) based on their size, 
# allowing the user to remove ranges that fall outside specified size thresholds. 
# This is particularly useful in genomic simulations and analyses where certain ranges 
# (e.g., nucleosome positions) must adhere to specific size constraints.
#
# Key functionalities:
# - Filters out genomic ranges (GRanges) that are smaller than a minimum size or larger than a maximum size.
# - Provides flexibility to set both minimum and maximum size constraints.
# - Includes an optional verbose mode to display information about the number of removed ranges.
# 
# Author: Dr Dave Gerrard

## Function to remove GRanges by size
# Args:
#   x: A GRanges object to filter.
#   maxSize: Maximum allowable size for the ranges. Default is -1 (no maximum limit).
#   minSize: Minimum allowable size for the ranges. Default is 1.
#   verbose: Boolean flag to print details about removed ranges. Default is FALSE.
# Returns:
#   Filtered GRanges object with ranges within the specified size constraints.
removeGRangesBySize <- function(x, maxSize = -1, minSize = 1, verbose = FALSE) {
  require(GenomicRanges)
  
  # Ensure the input is a GRanges object
  stopifnot(class(x) == "GRanges")
  
  # Apply size filters
  aboveMin <- width(x) >= minSize
  if(maxSize == -1) {
    withinMax <- rep(TRUE, length(x))  # No upper limit
  } else {
    withinMax <- width(x) <= maxSize
  }
  
  # Keep ranges that are within both min and max size limits
  keepIndex <- withinMax & aboveMin
  if(verbose) {
    print(paste("Removing features:", sum(!keepIndex), "/", length(x)))
  }
  
  # Return the filtered GRanges object
  return(x[keepIndex])
}

# Example usage:
# hist(width(removeGRangesBySize(x = nucLayerSet_2$layerSet[["H3K4me1"]], minSize = 147)), breaks = 50)
# range(width(removeGRangesBySize(x = nucLayerSet_2$layerSet[["H3K4me1"]], minSize = 147)))
