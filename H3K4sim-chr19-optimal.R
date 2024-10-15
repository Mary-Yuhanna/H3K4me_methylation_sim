# This script utilises the GenomicLayers package, originally developed by Dr. Dave Gerrard at the University of Manchester.
# The original skeleton code was focused on general genomic data analysis, including nucleosome positioning 
# and the methylation/demethylation processes of H3K4 in nucleosomes.
# 
# In collaboration with Dr. Gerrard, I, Mary Yuhanna, have optimised and extensively modified this script 
# for specific application to the human genome (hg38), particularly focusing on chromosome 19. 
# These modifications include adjusting binding factors, updating simulation parameters, and refining the code 
# for optimal performance with human genomic data. I have also removed original components that were not applicable 
# and replaced them with new, functional parts to enhance the overall efficiency of the simulation.


# Load necessary libraries
library(getopt)
library(GenomicLayers)
library(BiocGenerics)
library(Biostrings)
library(devtools)

# Define command-line options
spec <- matrix(c(
  'verbose',   'v', 2, "integer",
  'help',      'h', 0, "logical",
  'eraser',    'e', 1, "integer",
  'out.file',  'o', 2, "character",
  'writer',    'w', 2, "integer"
), byrow = TRUE, ncol = 4)

# Parse command-line options
opt <- getopt(spec)

# Set default values for options if they are not provided
if (is.null(opt$writer))  { opt$writer  <- 24000 }
if (is.null(opt$eraser))  { opt$eraser  <- 48000 }
if (is.null(opt$verbose)) { opt$verbose <- FALSE }

# Output file
outFile <- opt$out.file
if (is.null(opt$out.file)) { outFile <- "outputfile" }

# Print the values for verification
print(paste("Writer:", opt$writer))
print(paste("Eraser:", opt$eraser))

# Assign the values to variables
thisWriterValue <- opt$writer
thisEraserValue <- opt$eraser

# LOAD LIBRARIES ----------------------------------------------------------

# Define the path to the GenomicLayers package
genomic_layers_path <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/GenomicLayers"

# Source additional scripts (update paths accordingly)
# source("/path/to/removeShortGRanges.R")
# source("/path/to/removeGRangesBySize.R")
# source("/path/to/plotLayers.R")

# Load genome data
library(BSgenome.Hsapiens.UCSC.hg38)
genome1 <- BSgenome.Hsapiens.UCSC.hg38
sequences_to_keep <- "chr19"  # Focus on chromosome 19
genome <- keepBSgenomeSequences(genome1, sequences_to_keep)

# CREATE METHYLATING BINDING FACTORS --------------------------------------

# This section simulates the process of DNA methylation and demethylation at nucleosome positions.
# Binding factors are created to recognise specific patterns in the genome and apply or remove methylation marks.
# The simulation models methylation states from me0 to me3 and transitions between these states.

nucleosomeWidth <- 147  

allNucleosomeLayers <- c(
  "sampled_meUp", "sampled_meDown", "nucleosome", "H3K4me_promotion",
  "H3K4me_any", "H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"
)

# Create a DNA-binding factor used to initiate methylation at CG sites 
bf_TF <- createBindingFactor.DNA_consensus(
  name = "bf_TF",
  patternString = "CG",
  mod.layers = "tf_active",
  mod.marks = 1
)
bf_TF$profile$LAYER.1 <- NULL  # Remove unnecessary layer

# Remove bf_TF to test if the simulation works without it
# This allows checking if the pattern is independent of DNA sequence

# Select nucleosomes at random to be methylated
bf_meUp_sampler <- createBindingFactor.layer_region(
  name = "bf_meUp_sampler", 
  profile.layers = "nucleosome", 
  patternLength = nucleosomeWidth,
  profile.marks = 1,
  mod.layers = "sampled_meUp",
  mod.marks = 1,
  stateWidth = nucleosomeWidth
)

# Create bf_meUp binding factor to span 2 nucleosomes plus pattern length (300bp)
# This BF marks the areas as ready for H3K4 methylation
bf_meUp <- createBindingFactor.layer_region(
  name = "bf_meUp",
  profile.layers = c("tf_active", "sampled_meUp"), 
  patternLength = 2,
  profile.marks = c(1, 1),
  mod.layers = "H3K4me_promotion",
  mod.marks = 1,
  stateWidth = 2 + (2 * nucleosomeWidth)
)

# Define binding factors that change the methylation states
# Profile marks: 1 means marked, 0 means unmarked
# For example, bf_promotion_0_1 looks for nucleosomes marked for H3K4me_promotion but unmarked for me1, me2, me3

# Transition from me0 to me1
bf_promotion_0_1 <- createBindingFactor.layer_region(
  name = "bf_promotion_0_1", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 1, 0, 0, 0), 
  mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 1, 0, 0, 1)
)

# Transition from me1 to me2
bf_promotion_1_2 <- createBindingFactor.layer_region(
  name = "bf_promotion_1_2", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 1, 1, 0, 0), 
  mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 1, 0, 1)
)

# Transition from me2 to me3
bf_promotion_2_3 <- createBindingFactor.layer_region(
  name = "bf_promotion_2_3", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 1, 0, 1, 0), 
  mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 0, 1, 1)
)

# Handle me3 state without promotion
bf_promotion_3_3 <- createBindingFactor.layer_region(
  name = "bf_promotion_3_3", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 1, 0, 0, 1), 
  mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 0, 1, 1)
)

# CREATE DEMETHYLATING BINDING FACTORS ------------------------------------

# Create binding factors for demethylation processes
# A core BF marks regions for demethylation, and accessory BFs perform state transitions

# Mark regions for demethylation
bf_meDown <- createBindingFactor.layer_region(
  name = "bf_meDown",
  profile.layers = "nucleosome", 
  patternLength = nucleosomeWidth,
  profile.marks = 1,
  mod.layers = "H3K4me_demethylate",
  mod.marks = 1,
  stateWidth = nucleosomeWidth
)

# Transition from me3 to me2
bf_demethylate_3_2 <- createBindingFactor.layer_region(
  name = "bf_demethylate_3_2", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 0, 0, 1), 
  mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 1, 0, 1)
)

# Transition from me2 to me1
bf_demethylate_2_1 <- createBindingFactor.layer_region(
  name = "bf_demethylate_2_1", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 0, 1, 0), 
  mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 1, 0, 0, 1)
)

# Transition from me1 to me0
bf_demethylate_1_0 <- createBindingFactor.layer_region(
  name = "bf_demethylate_1_0", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 1, 0, 0), 
  mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 0, 0, 0)
)

# Handle me0 state without demethylation
bf_demethylate_0_0 <- createBindingFactor.layer_region(
  name = "bf_demethylate_0_0", 
  patternLength = nucleosomeWidth, 
  profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
  profile.marks = c(1, 0, 0, 0), 
  mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
  mod.marks = c(0, 0, 0, 0, 0)
)

# Combine methylation binding factors into a list
bfList_H3K4_methylate <- list(
  bf_TF = bf_TF,
  bf_meUp_sampler = bf_meUp_sampler,
  bf_meUp = bf_meUp,
  bf_promotion_0_1 = bf_promotion_0_1, 
  bf_promotion_1_2 = bf_promotion_1_2,
  bf_promotion_2_3 = bf_promotion_2_3,
  bf_promotion_3_3 = bf_promotion_3_3
)

# Assign abundances to the binding factors
saturationAbundance <- 1e10
bf_methylate_abu <- c(thisWriterValue, rep(saturationAbundance, 6))  
names(bf_methylate_abu) <- names(bfList_H3K4_methylate)

# CREATE LAYERSET USING THE GENOME ----------------------------------------

# Set up the layerset with the genome
scLayerSetNuc <- createLayerSet.BSgenome(
  genome = genome, 
  layer.names = c(
    "sampled_meUp", "sampled_meDown", "nucleosome", "tf_active",
    "H3K4me_promotion", "H3K4me_any", "H3K4me_demethylate",
    "H3K4me1", "H3K4me2", "H3K4me3"
  ), 
  n.layers = 11,
  verbose = TRUE
)

# Test some of the binding factors
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_TF)
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_meUp)

# Deploy demethylation agents
bfList_H3K4_demethyl <- list(
  bf_meDown = bf_meDown,
  bf_demethylate_3_2 = bf_demethylate_3_2,
  bf_demethylate_2_1 = bf_demethylate_2_1,
  bf_demethylate_1_0 = bf_demethylate_1_0,
  bf_demethylate_0 = bf_demethylate_0_0
)

# Assign abundances to demethylation binding factors
bf_demethylate_abu <- c(thisEraserValue, rep(saturationAbundance, 4))
names(bf_demethylate_abu) <- names(bfList_H3K4_demethyl)

# Test demethylation binding factors
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_meDown)
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_demethylate_3_2)

# Combine all binding factors and abundances
bfLIst_both <- c(bfList_H3K4_methylate, bfList_H3K4_demethyl)
bfAbund_both <- c(bf_methylate_abu, bf_demethylate_abu)
layersToClean <- c(
  "sampled_meUp", "sampled_meDown", "H3K4me_promotion", "H3K4me1",
  "H3K4me2", "H3K4me3", "H3K4me_any", "H3K4me_demethylate"
)

# SECTION LOOP 1 - POSITION THE NUCLEOSOMES --------------------------------

# Position nucleosomes randomly, ensuring they do not overlap and are not closer than 10bp
scLayerSetNuc$layerSet[["nucleosome"]] <- randGrangesBigGenome(
  genome = genome, 
  sizeFunc = function(value, n) rep(x = value, times = n), 
  gapFunc = function(n, value) rpois(n = n, lambda = value), 
  argsSizeFunc = list(value = 147, n = 5),
  argsGapFunc = list(value = 40, n = 7)
)

# Save a copy of the layerset for future use
fileOut <- paste0(outFile, "scLayerSetNuc.object.Rdata")
save(scLayerSetNuc, file = fileOut) 

# SECTION LOOP 2 - APPLY METHYLATION AND DEMETHYLATION ---------------------

# Initialize the layerset for the simulation
scLayerSetBoth <- scLayerSetNuc

# Prepare binding factors and abundances for simulation
clean <- TRUE
bfAbund_both_eFirst <- c(bf_demethylate_abu, bf_methylate_abu) 
bfLIst_both_eFirst <- bfLIst_both[names(bfAbund_both_eFirst)] 

# Adjust abundances based on command-line inputs
bfAbund_both_eFirst["bf_meUp_sampler"] <- thisWriterValue
bfAbund_both_eFirst["bf_meDown"] <- thisEraserValue

# Run the simulation for a specified number of iterations
execution_time <- system.time({
  for (i in 1:100) {  # Adjust the number of iterations as needed
    print(i)
    scLayerSetBoth <- runLayerBinding.BSgenome(
      layerList = scLayerSetBoth,
      factorSet = bfLIst_both_eFirst,
      bf.abundances = bfAbund_both_eFirst,
      verbose = TRUE,
      collect.stats = TRUE,
      keep.stats = TRUE
    )
    if (clean) {
      # Clean up layers by removing short ranges and ensuring proper sizes
      for (thisLayer in layersToClean) {
        scLayerSetBoth$layerSet[[thisLayer]] <- removeShortGRanges(
          x = scLayerSetBoth$layerSet[[thisLayer]],
          minSize = nucleosomeWidth
        )
      }
      scLayerSetBoth$layerSet[["nucleosome"]] <- removeGRangesBySize(
        x = scLayerSetBoth$layerSet[["nucleosome"]],
        verbose = TRUE,
        minSize = nucleosomeWidth,
        maxSize = nucleosomeWidth
      )
      scLayerSetBoth$layerSet[["sampled_meUp"]] <- removeGRangesBySize(
        x = scLayerSetBoth$layerSet[["sampled_meUp"]],
        verbose = TRUE,
        maxSize = 0
      )
    }
  }
})

# Export the output from the simulation
fileOut <- paste0(outFile, ".history.tsv")
write.table(scLayerSetBoth$history, file = fileOut, quote = FALSE, sep = "\t")

# Generate a plot of the simulation results
fileOut <- paste0(outFile, ".nBlocks.png")
png(filename = fileOut, width = 2400, height = 1600, res = 300, type = "cairo")
plot(
  scLayerSetBoth$history$nBlocks.H3K4me3, pch = 6, ylab = "Nucleosomes", 
  xlab = "Iteration", main = "Simulation Results", ylim = c(0, 350000)
)
points(scLayerSetBoth$history$nBlocks.H3K4me2, pch = 4)
points(scLayerSetBoth$history$nBlocks.H3K4me1, pch = 1)
legend("topleft", legend = c("H3K4me3", "H3K4me2", "H3K4me1"), pch = c(6, 4, 1))
dev.off()

print(execution_time)

# Data extraction from simulation for later analysis
me1_layer <- scLayerSetBoth$layerSet[["H3K4me1"]]
me2_layer <- scLayerSetBoth$layerSet[["H3K4me2"]]
me3_layer <- scLayerSetBoth$layerSet[["H3K4me3"]]
history_layer <- scLayerSetBoth$history

# Save the methylation layers for future use
save(me1_layer, file = paste0(outFile, "_me1_layer.RData"))
save(me2_layer, file = paste0(outFile, "_me2_layer.RData"))
save(me3_layer, file = paste0(outFile, "_me3_layer.RData"))

# Optional: Calculate and print the proportion of the genome covered by each methylation state
print(sum(width(me1_layer)))
print(sum(width(me1_layer)) / 3.2e9)  # Assuming genome size of 3.2 Gb

print(sum(width(me2_layer)))
print(sum(width(me2_layer)) / 3.2e9)

print(sum(width(me3_layer)))
print(sum(width(me3_layer)) / 3.2e9)

# Clean up large objects if no longer needed
# unlink(paste0(outputDir, "scLayerSetBoth.object.Rdata"))
