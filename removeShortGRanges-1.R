# Genomic Ranges Filtering by Minimum Size
#
# This script defines a function to remove genomic ranges (GRanges) that are smaller than a specified minimum size.
# It is commonly used in genomic analyses where short ranges (e.g., small segments or noise) need to be filtered out.
# 
# Key functionalities:
# - Filters out ranges from a GRanges object that do not meet the specified minimum size requirement.
# - Provides a simple mechanism for cleaning up data in genomic simulations and analyses.
# 
# Author: Dr Dave Gerrard

## Function to remove GRanges shorter than a given minimum size
# Args:
#   x: A GRanges object to filter.
#   minSize: Minimum allowable size for the ranges. Default is 1.
# Returns:
#   Filtered GRanges object with ranges that are at least minSize in width.
removeShortGRanges <- function(x, minSize = 1) {
  require(GenomicRanges)
  
  # Ensure the input is a GRanges object
  stopifnot(class(x) == "GRanges")
  
  # Apply the minimum size filter
  keepIndex <- width(x) >= minSize
  
  # Return the filtered GRanges object
  return(x[keepIndex])
}

# Example usage:
# hist(width(removeShortGRanges(x = nucLayerSet_2$layerSet[["H3K4me1"]], minSize = 147)), breaks = 50)
# range(width(removeShortGRanges(x = nucLayerSet_2$layerSet[["H3K4me1"]], minSize = 147)))
