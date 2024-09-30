# Load the getopt library
library(getopt)

# Define command-line options
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'eraser'  , 'e', 1, "integer",
  'out.file' , 'o', 2, "character",
  'writer' , 'w', 2, "integer"
), byrow=TRUE, ncol=4)

# Parse command-line options
opt <- getopt(spec)

# Set default values for options if they are not provided
if ( is.null(opt$writer) ) { opt$writer <- 24000 }
if ( is.null(opt$eraser) ) { opt$eraser <- 48000 }
if ( is.null(opt$verbose)) { opt$verbose <- FALSE }

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
#MARY'S NOTES - as described, this loads the GenomicLayers library and sources three R scripts that are essential for the simulation: removeShortGRanges, removeGRangesBySize and plotlayers


library(GenomicLayers)  # to use installed version
library(BiocGenerics)
library(Biostrings)
library(GenomicLayers)
library(devtools)

#defining the path
genomic_layers_path <- "/Library/Frameworks/R.framework/Versions/4.2/Resources/library/GenomicLayers"


source("/scratch/j52505my/startingResources/removeShortGRanges.R")
source("/scratch/j52505my/startingResources/removeGRangesBySize.R")
source("/scratch/j52505my/startingResources/plotLayers.R")

require(Biostrings)
require(BSgenome.Scerevisiae.UCSC.sacCer3)

#MARY'S NOTES - I believe this sets the genome as the Saccrohmyces genome...
library(BSgenome.Hsapiens.UCSC.hg38)
genome1 <- BSgenome.Hsapiens.UCSC.hg38
seqinfo(genome1)
seqnames(genome1)
sequences_to_keep <- "chr19"   # nuclear chroms only
genome <- keepBSgenomeSequences(genome1, sequences_to_keep)
seqinfo(genome)


# CREATE METHYLATING BINDING FACTORS --------------------------------------

#MARY'S NOTES - I believe this code is to simulate the process of DNA methylation and demethylation at nucleosome positions in the genome. We create various binding factors that recognise specific patterns in the genome and apply or remove methylation marks. The simulation aims to model the methylation states from me0 to me3 and the transitions between these states.

# Modelling allowing one factor to remove sequential methylation 

# Requires accessory BFs that variously match me0, me1, me2 or me3 and then set the next lower layer.
#  see sequential_H3K4me.R for the opposing writer
# this script goes on to combine the two

# One core bf recognises the main pattern but instead of making the change, it marks that part of the genome
# as "me_promotion".  Then four separate bfs operate to deliver the four possible state changes:-
#  me0 -> me1  
#  me1 -> me2
#  me2 -> me3
#  me3 -> me4
# as they operate, they set all other layers to zero (unmarked). 

nucleosomeWidth <- 147  

allNucleosomeLayers <- c("sampled_meUp", "sampled_meDown", "nucleosome","H3K4me_promotion", "H3K4me_any", "H3K4me_demethylate","H3K4me1", "H3K4me2", "H3K4me3")

#MARY'S NOTES - this part of the code is commented out - possibly to see if the simulation focuses on nucleosome position rather than DNA sequence? Not too sure may need to confirm with Dave. 
##  a DNA-binding factor used to initiate methylation at certain sites. 
bf_TF <-  createBindingFactor.DNA_consensus(name="bf_TF", patternString = "CG",
                                            mod.layers="tf_active", mod.marks = 1)
bf_TF$profile$LAYER.1 <- NULL

# in this script I am removing the bf_TF to check that this will work without it. 
# IF i pre-specify a number of nucleosomes as TF_ACTIVE
# I want to check the pattern is independent of DNA sequence before
# I then measure the effect of DNA sequence by re-introducing TF_BF or a CpG element. 
#MARY'S NOTES - I guess this means we can compare the effect of specific DNA sequences on the methylation process.


#MARY'S NOTES - creating a binding factor that selects nucleosomes to be methylated at random. I think this mimics the natural biological process of priming nucleosomes? 
# selects nucleosomes at random to be methylated.  
bf_meUp_sampler <- createBindingFactor.layer_region(name="bf_meUp_sampler", 
                                                    profile.layers="nucleosome", 
                                                    patternLength = nucleosomeWidth,
                                                    profile.marks = 1,
                                                    mod.layers="sampled_meUp",
                                                    mod.marks=1,
                                                    stateWidth= nucleosomeWidth)



#MARY'S NOTES - creating bf_meUp binding factor to incorporate the BF spanning 2 nucleosomes plus pattern length, which affects a larger genomic context. In this case 300bp. This BF marks the areas as H3K4me ready.  
# in this script, where I am trying to tie upgrades to nucleosome positions, 
#  this meUp bf needs to widen the effect of tfs from 147bp (previous sims) to (147bp * 2) + patternlength
# Any nucleosome fully within this region could then be a candidate. 
# the other bfs should continue to use nucleosomeWidth but now require "nucleosome" to be marked
bf_meUp <- createBindingFactor.layer_region(name="bf_meUp", profile.layers=c("tf_active", "sampled_meUp"), 
                                            patternLength = 2,
                                            profile.marks = c(1,1),
                                            mod.layers="H3K4me_promotion",
                                            mod.marks=1,
                                            stateWidth= 2 + (2*nucleosomeWidth))



#MARY'S NOTES - these sections define more BFs that change the methylation states. The profile marks 1 means marked and 0 means unmarked. For example, in bf_0_1 the profile mark it is looking for is (1,1,0,0,0), marked for nucleosome, marked for H3K4me promotation, but unmarked for any methylation states (me1,2,3). Once the criterion is satisfied, this modifies the marks to (0,1,0,0,1) which is marked for me1, unmarked for me2 and me3 and sets H3K4any to marked.   
# this accessory BF acts on the hits identified by bf_meUp and looks for 
# those which have zero methylation. These become marked at H3K4me1
bf_promotion_0_1 <- createBindingFactor.layer_region(name="bf_promotion_0_1", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1, 0,0,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 1,0,0,1))
# this accessory bf promotes me1 to me2 where "me_promotion" is marked
bf_promotion_1_2 <- createBindingFactor.layer_region(name="bf_promotion_1_2", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c( 1,1, 1,0,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,1,0,1))
# this accessory bf promotes me1 to me2 where "me_promotion" is marked
bf_promotion_2_3 <- createBindingFactor.layer_region(name="bf_promotion_2_3", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1, 0,1,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,0,1,1))
# need a NULL option if me3 is an allowable pattern - this one does nothing but remove the promotion.
bf_promotion_3_3 <- createBindingFactor.layer_region(name="bf_promotion_3_3", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1, 1, 0,0,1), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,0,1,1))


# CREATE DE-METHYLATING BINDING FACTORS -----------------------------------


# One core bf recognises the main pattern but instead of making the change, it marks that part of the genome
# as "H3K4me_demethylate".  Then four separate bfs operate to deliver the four possible state changes:-
#  me0 <- me1  
#  me1 <- me2
#  me2 <- me3
#  me3 <- me3
# as they operate, they set all other layers to zero (unmarked). 

# does not automatically "home in" on any methylation.  
#   Also samples unmethylated nucleosomes - MARYS NOTES - why does it also sample unmethylated nucleosomes, if the purpose is to demethylate? Surely we don't need should start with methylated nucleosomes? 
bf_meDown <- createBindingFactor.layer_region(name="bf_meDown", profile.layers="nucleosome", 
                                              patternLength = nucleosomeWidth,
                                              profile.marks = 1,
                                              mod.layers="H3K4me_demethylate",
                                              mod.marks=1,
                                              stateWidth=nucleosomeWidth)

# this accessory BF acts on the hits identified by bf_meDown and looks for 
# this accessory bf demethylates me3 to me2 where "H3K4me_demethylate" is marked
bf_demethylate_3_2 <- createBindingFactor.layer_region(name="bf_demethylate_3_2", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,0,1), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 0,1,0,1))
# this accessory bf demethylates me2 to me1 where "H3K4me_demethylate" is marked
bf_demethylate_2_1 <- createBindingFactor.layer_region(name="bf_demethylate_2_1", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,1,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 1,0,0,1))

# this accessory bf demethylates me1 to me0 where "H3K4me_demethylate" is marked
bf_demethylate_1_0 <- createBindingFactor.layer_region(name="bf_demethylate_1_0", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 1,0,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c( 0,0,0,0,0))

#  this one should do nothing but set H3K4me_any to zero/empty
bf_demethylate_0_0 <- createBindingFactor.layer_region(name="bf_demethylate_0_0", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,0,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 0,0,0,0))

# all of the above need to be added along with bf_meUp in quantities to saturate the genome (every promotion site).

#MARY'S NOTES - this creates a list of BFs associated with the methylation process. We then assign a numerical value for the BFs saturation to ensure full genome saturation. 
bfList_H3K4_methylate <- list(bf_TF=bf_TF, bf_meUp_sampler=bf_meUp_sampler, bf_meUp=bf_meUp, bf_promotion_0_1=bf_promotion_0_1, 
                              bf_promotion_1_2=bf_promotion_1_2, bf_promotion_2_3=bf_promotion_2_3,
                              bf_promotion_3_3=bf_promotion_3_3)


# Assign abundances to the binding factors. MARY'S NOTES - this sets meUp sampler initially to 2400 and the remaining five BFs to 1,000,000
saturationAbundance <- 10000000000
#changed abundance from 24000 to 48000 to increase sensitvity
bf_methylate_abu <- c( thisWriterValue, rep(saturationAbundance, 6) )  
names(bf_methylate_abu) <- names(bfList_H3K4_methylate)



# CREATE LAYERSET USING THE GENOME ----------------------------------------

#MARY'S NOTES - This creates the layers with the previously specified genome and ensures that we have some output for debugging purposes. Some factors are commented out which ensures proper function of the BFs.
# Now set up the layerset, with the genome we are using.  
scLayerSetNuc <- createLayerSet.BSgenome(genome=genome, 
                                         layer.names=c("sampled_meUp", "sampled_meDown",
                                                       "nucleosome", "tf_active","H3K4me_promotion", 
                                                       "H3K4me_any", "H3K4me_demethylate",
                                                       "H3K4me1", "H3K4me2", "H3K4me3"), 
                                         n.layers=11,verbose=TRUE)

# test some of the factors.   
#matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf_nucleosome)
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_TF)  # ~ 5000 hits
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_meUp)
#matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf_promotion_0_1)

# now deploy the demethylation agents. 
bfList_H3K4_demethyl <- list(bf_meDown=bf_meDown, bf_demethylate_3_2=bf_demethylate_3_2, bf_demethylate_2_1=bf_demethylate_2_1, bf_demethylate_1_0=bf_demethylate_1_0, bf_demethylate_0= bf_demethylate_0_0)
# I suspect bf_demethylate_0_0 is not needed
# however, it might be worth having something like this to sample randomly (see also recent "sampler" type bf)

#MARY'S NOTES - Once the list of demethylation BFs has been created, we set their abundances accordingly. The first value of 0 acts as an eraser and the subsequent BFs are set to higher values to ensure saturation of the genome. The commented out section is for validation and debugging purposes. 
bf_demethylate_abu <- c( thisEraserValue, rep(saturationAbundance, 4))  # setting the eraser to Zero for the burn-in.changed to 24000, which makes the ratio of methylators to demethylators 2:1
names(bf_demethylate_abu) <- names(bfList_H3K4_demethyl)

matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_meDown)    # should be many
matchBindingFactor.BSgenome(layerSet = scLayerSetNuc, bindingFactor = bf_demethylate_3_2)  # should be none yet

##simulation with both methyl writers and erasers. 

bfLIst_both <- c(bfList_H3K4_methylate, bfList_H3K4_demethyl)
bfAbund_both <- c(bf_methylate_abu, bf_demethylate_abu)  # there are no shared names so this is OK. Otherwise write out in full
layersToClean <- c("sampled_meUp", "sampled_meDown", 'H3K4me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3',  'H3K4me_any', 'H3K4me_demethylate') 


# SECTION LOOP 1 - POSITION THE NUCLEOSOMES  --------------------------------------------

#MARY'S NOTES - We are using randomness to position the nucelosomes to ensure they do not overlap. This initial layerset is then saved to scLayerSetNuc.object.Rdata for future reruns. 
# Nucleosome should be non-overlapping and no closer that 10bp. 

# new (much faster) version
scLayerSetNuc$layerSet[["nucleosome"]] <- randGrangesBigGenome(genome=genome, 
                                                               sizeFunc = function(value, n) rep(x=value, times=n), 
                                                               gapFunc = function(n, value) rpois(n = n , lambda = value), 
                                                               argsSizeFunc = list(value=147, n=5),
                                                               argsGapFunc = list(value=40, n=7))

#save a copy  layerSet used as the start for all these in case we need to re-run
fileOut <- paste0(outFile,"scLayerSetNuc.object.Rdata")
save(scLayerSetNuc, file = fileOut) 

# SECTION LOOP 2 - Apply the methylation and de-methylation with  --------------------------------------------

#This loop chooses which nucleosomes are going to be modified 


scLayerSetBoth <- scLayerSetNuc   # in case we want to keep a clean copy

# for this script, dicate a subset of nucleosomes to be modified and are marked as tf_active
# this is to check the model is NOT dependent on sequence. 
# going to leave this as quite rare for this script 
#  5% of all nucleosome = 3240
# take this section out as we need to see if the simulation is dependent on sequence. 
#randNucs <- sample(1:length(scLayerSetNuc$layerSet$nucleosome), 3240, replace = FALSE)
#scLayerSetBoth$layerSet$tf_active <- scLayerSetNuc$layerSet$nucleosome[randNucs]
#scLayerSetBoth$layerSet$tf_active <- resize(scLayerSetBoth$layerSet$tf_active, width=6, fix = "center" )  # set to minimum pattern size for meUP to avoid clashing.
#bfAbund_both["bf_meDown"] <- 0   # use this in a dry run to determine unopposed kinetics. MARY'S NOTES - (To understand how the simulation) behaves without the demethylation process. 

# Now run through a variety of abundance levels and runLayerBindng from scLayerSetBoth

# here a for loop to test if running the demethylase first, appreciably increases
# the final levels of me1,2,3
# using a differently ordered set of bf abundances. 
clean <- TRUE
abundSpread <- seq(24000, 2000, by=-2000)   
bfAbund_both_eFirst <- c(bf_demethylate_abu, bf_methylate_abu) 
bfLIst_both_eFirst <- bfLIst_both[names(bfAbund_both_eFirst)]   # this determines the order the binding factors are used in runLayerBinding

# a for loop where the value of the writer "bf_meUp" is varied but it's abundance is controlled by "bf_meUp_sampler"

#for(thisAbund in abundSpread)  {
#thisAbund <- 2400   # use this to try out just one loop
scLayerSetBothAbund <- scLayerSetBoth

##changed abundance from 24000 to 48000 to increase sensitvity
bfAbund_both_eFirst["bf_meUp_sampler"] <- thisWriterValue
bfAbund_both_eFirst["bf_meDown"] <- thisEraserValue
#subSimName <- paste0("abundSpread.w",bfAbund_both_eFirst["bf_meUp_sampler"] ,".e",bfAbund_both_eFirst["bf_meDown"])
#print(subSimName)

# Measure the time taken for the execution of the inner loop
execution_time <- system.time({
  for(i in 1:100) { # MARY'S NOTES - running the simulation 100 times 
    print(i)

#bfAbund_both_eFirst["bf_meUp_sampler"] <- thisWriterValue
    #bfAbund_both_eFirst["bf_meDown"] <- thisEraserValue

    scLayerSetBothAbund <- runLayerBinding.BSgenome(layerList = scLayerSetBothAbund, factorSet = bfLIst_both_eFirst, bf.abundances = bfAbund_both_eFirst, verbose = T, collect.stats = T, keep.stats = T) # MARY'S NOTES - this applied the BFs to layerset. 
    if(clean) {
      # TIDY UP BY REMOVING HISTONE MODS < 147BP..  # I'M NOT SURE THIS IS A GOOD IDEA.
      for(thisLayer in layersToClean) {
        scLayerSetBothAbund$layerSet[[thisLayer]] <- removeShortGRanges(x=scLayerSetBothAbund$layerSet[[thisLayer]], minSize = nucleosomeWidth) # MARY'S NOTES - This removes histone modifications shorter than the nucleosome width?
      }    
      scLayerSetBothAbund$layerSet[["nucleosome"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["nucleosome"]], verbose=T, minSize = nucleosomeWidth, maxSize = nucleosomeWidth)  # remove anything not exactly nucleosomeWidth
      scLayerSetBothAbund$layerSet[["sampled_meUp"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["sampled_meUp"]], verbose=T, maxSize = 0)  # clean out the sampler each iteration
    }
  }
})


# export the output from this run.
fileOut <- paste0(outFile, ".history.tsv")
write.table(scLayerSetBothAbund$history, file= fileOut, quote=F, sep="\t")

fileOut <- paste0(outFile, ".nBlocks.png")
png(filename = fileOut, width=2400, height=1600, res=300, type= "cairo")
plot(scLayerSetBothAbund$history$nBlocks.H3K4me3, pch=6, ylab="nucleosomes", 
     xlab="iteration", main=subSimName, ylim=c(0,350000))
points(scLayerSetBothAbund$history$nBlocks.H3K4me2, pch=4)
points(scLayerSetBothAbund$history$nBlocks.H3K4me1, pch=1)
legend("topleft", legend=c("H3K4me3", "H3K4me2", "H3K4me1"), pch=c(6, 4, 1))
dev.off()

print(execution_time)


#Data extraction from human Chr19 simulation for later overlap analysis

me1_layer <- scLayerSetBothAbund$layerSet[["H3K4me1"]]
me2_layer <- scLayerSetBothAbund$layerSet[["H3K4me2"]]
me3_layer <- scLayerSetBothAbund$layerSet[["H3K4me3"]]
history_layer <- scLayerSetBothAbund$history

print(me1_layer)
print(me2_layer)
print(me3_layer)
print(head(history_layer))

save(me1_layer, file = paste0(outFile, "me1_layer.RData"))
save(me2_layer, file = paste0(outFile, "me2_layer.RData"))
save(me3_layer, file = paste0(outFile, "me3_layer.RData"))
#write.table(history_layer, file = paste0(outputDir, "history_layer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


sum(width(me1_layer))
sum(width(me1_layer)) / 3.2e9

sum(width(me2_layer))
sum(width(me2_layer)) / 3.2e9

sum(width(me3_layer))
sum(width(me3_layer)) / 3.2e9




#unlink(paste0(outputDir,"scLayerSetBoth.object.Rdata"))  # large object that is probably no longer needed.
