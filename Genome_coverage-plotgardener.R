# Genomic Coverage and Signal Plotting for H3K4me2 Wild-Type and Kabuki Syndrome
#
# This script loads genomic data for both wild-type and Kabuki syndrome simulations (ME2 layer).
# It calculates coverage for specific chromosomes (chr20 and chr3) and compares the results with Jung et al.'s data.
# The data is then visualised using plotgardener to generate chromosome-level signal plots.
#
# Key functionalities:
# - Load genomic data for wild-type and Kabuki syndrome for specific chromosomes.
# - Calculate coverage for each dataset using GenomicRanges.
# - Plot signals and genomic ranges for comparison on chromosome 20 and chromosome 3.
# - Annotate the plots for visual clarity.
#
# Author: Mary Yuhanna
# ----------------------------------------------------

# Load necessary libraries
library(GenomicRanges)
library(plotgardener)
library(readr)  # For reading CSV files

# ----------------------------------------------------
# Load and process Jung Data
# ----------------------------------------------------
file_path <- "/Users/maryyuhanna/desktop/Files for Kabuki/H3K4me2_results_KS.share copy.csv"
jung_data <- read_csv(file_path)

# Filter for chromosome 20 and convert to GRanges object
jung_data_chr20 <- jung_data[jung_data$Chr == "chr20", ]  
granges_jung_20 <- GRanges(
  seqnames = Rle(jung_data_chr20$Chr),        
  ranges = IRanges(start = jung_data_chr20$Start, end = jung_data_chr20$End)
)

# Filter for chromosome 3 and convert to GRanges object
jung_data_chr3 <- jung_data[jung_data$Chr == "chr3", ]  
granges_jung_3 <- GRanges(
  seqnames = Rle(jung_data_chr3$Chr),        
  ranges = IRanges(start = jung_data_chr3$Start, end = jung_data_chr3$End)
)

# Calculate coverage for Jung Data for chromosome 20 and 3
scoresGR_jung_20 <- as(coverage(granges_jung_20), "GRanges")
scoresGR_jung_3 <- as(coverage(granges_jung_3), "GRanges")

# ----------------------------------------------------
# Calculate coverage for Wild-Type and Kabuki Syndrome data
# ----------------------------------------------------

# Ensure me2_layer_new (wild-type) is defined
if (!exists("me2_layer_new")) {
  stop("me2_layer_new is not defined. Please ensure you load or calculate this object.")
}

# Calculate coverage for wild-type
scoresGR_wt <- as(coverage(me2_layer_new), "GRanges")

# Subset for chromosome 20 and chromosome 3 for wild-type
scoresGR_chr20_wt <- scoresGR_wt[seqnames(scoresGR_wt) == "chr20", ]
scoresGR_chr3_wt <- scoresGR_wt[seqnames(scoresGR_wt) == "chr3", ]

# Ensure me2_layer_ks (Kabuki syndrome) is defined
if (!exists("me2_layer_ks")) {
  stop("me2_layer_ks is not defined. Please ensure you load or calculate this object.")
}

# Calculate coverage for Kabuki syndrome
scoresGR_ks <- as(coverage(me2_layer_ks), "GRanges")

# Subset for chromosome 20 and chromosome 3 for Kabuki syndrome
scoresGR_chr20_ks <- scoresGR_ks[seqnames(scoresGR_ks) == "chr20", ]
scoresGR_chr3_ks <- scoresGR_ks[seqnames(scoresGR_ks) == "chr3", ]

# ----------------------------------------------------
# Set up plotgardener page and plot signals
# ----------------------------------------------------

# Create a plot page
pageCreate(width = 10, height = 10)  # Set page dimensions

# Coordinates for chromosome 20 and chromosome 3
chrom_20 <- "chr20"
chromstart_20 <- 57495133
chromend_20 <- 57495600

chrom_3 <- "chr3"
chromstart_3 <- 71488262
chromend_3 <- 71494719

# ---- PANEL A: Plot for chromosome 20 (Wild-type and Kabuki Syndrome) ----

# Wild-type coverage for chromosome 20
plotSignal(
  data = scoresGR_chr20_wt,
  chrom = chrom_20,
  chromstart = chromstart_20,
  chromend = chromend_20, 
  assembly = "hg38",
  x = 1, y = 2, width = 4, height = 1.5,  
  linecolor = "blue"
)
plotText(
  label = "Panel A: ME2 Wild-type (Chromosome 20)",
  x = 1, y = 1.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Kabuki syndrome coverage for chromosome 20
plotSignal(
  data = scoresGR_chr20_ks,
  chrom = chrom_20,
  chromstart = chromstart_20,
  chromend = chromend_20, 
  assembly = "hg38",
  x = 1, y = 4, width = 4, height = 1.5,  
  linecolor = "green"
)
plotText(
  label = "ME2 Kabuki syndrome (Chromosome 20)",
  x = 1, y = 3.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Plot Jung Data ranges for chromosome 20
plotRanges(
  data = granges_jung_20,
  chrom = chrom_20,
  chromstart = chromstart_20,
  chromend = chromend_20, 
  assembly = "hg38",
  x = 1, y = 5, width = 4, height = 1.5,  
  fill = "black",
  linecolor = "black"
)
plotText(
  label = "Jung Data Ranges (Chromosome 20)",
  x = 1, y = 5.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Add genome label for chromosome 20
plotGenomeLabel(
  chrom = chrom_20,
  chromstart = chromstart_20,
  chromend = chromend_20, 
  assembly = "hg38",
  x = 1, y = 7.5, length = 4
)

# ---- PANEL B: Plot for chromosome 3 (Wild-type and Kabuki Syndrome) ----

# Wild-type coverage for chromosome 3
plotSignal(
  data = scoresGR_chr3_wt,
  chrom = chrom_3,
  chromstart = chromstart_3,
  chromend = chromend_3, 
  assembly = "hg38",
  x = 6, y = 2, width = 4, height = 1.5,  
  linecolor = "blue"
)
plotText(
  label = "Panel B: ME2 Wild-type (Chromosome 3)",
  x = 6, y = 1.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Kabuki syndrome coverage for chromosome 3
plotSignal(
  data = scoresGR_chr3_ks,
  chrom = chrom_3,
  chromstart = chromstart_3,
  chromend = chromend_3, 
  assembly = "hg38",
  x = 6, y = 4, width = 4, height = 1.5,  
  linecolor = "green"
)
plotText(
  label = "ME2 Kabuki syndrome (Chromosome 3)",
  x = 6, y = 3.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Plot Jung Data ranges for chromosome 3
plotRanges(
  data = granges_jung_3,
  chrom = chrom_3,
  chromstart = chromstart_3,
  chromend = chromend_3, 
  assembly = "hg38",
  x = 6, y = 5, width = 4, height = 1.5,  
  fill = "black",
  linecolor = "black"
)
plotText(
  label = "Jung Data Ranges (Chromosome 3)",
  x = 6, y = 5.8, just = "left", fontsize = 10,
  fontface = "bold"
)

# Add genome label for chromosome 3
plotGenomeLabel(
  chrom = chrom_3,
  chromstart = chromstart_3,
  chromend = chromend_3, 
  assembly = "hg38",
  x = 6, y = 7.5, length = 4
)

# Close the plot
pageFinish()
