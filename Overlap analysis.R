#Data extraction from human Chr19 simulation for later overlap analysis

me1_layer <- scLayerSetBothAbund$layerSet[["H3K4me1"]]
me2_layer <- scLayerSetBothAbund$layerSet[["H3K4me2"]]
me3_layer <- scLayerSetBothAbund$layerSet[["H3K4me3"]]
history_layer <- scLayerSetBothAbund$history

print(me1_layer)
print(me2_layer)
print(me3_layer)
print(head(history_layer))

save(me1_layer, file = paste0(outputDir, "me1_layer.RData"))
save(me2_layer, file = paste0(outputDir, "me2_layer.RData"))
save(me3_layer, file = paste0(outputDir, "me3_layer.RData"))
write.table(history_layer, file = paste0(outputDir, "history_layer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


#loading in Jung et al's me2 data and converting them into genomic ranges. 

jung_data <- read.csv("/Users/maryyuhanna/desktop/Files for Kabuki/H3K4me2_results_KS.share.csv")

jung_chr19 <- subset(jung_data, Chr == "chr19")

# Convert the filtered data frame to a GRanges object
granges_chr19_jung <- GRanges(seqnames = Rle(jung_chr19$Chr),
                         ranges = IRanges(start = jung_chr19$Start, end = jung_chr19$End))
                         
sum(width(granges_chr19_jung))
sum(width(reduce(granges_chr19_jung)))
print(granges_chr19_jung)


#OVERLAP ANALYSIS

#providing genome size
library(BSgenome.Hsapiens.UCSC.hg38)
genome_10 <- BSgenome.Hsapiens.UCSC.hg38
genome_size_chr19 <- sum(seqlengths(genome_10)["chr19"])

#caluclating the confusion matrix

# Function to calculate the confusion matrix for genomic ranges
calcConfMat.GR <- function(query, subject, maxgap = -1L, minoverlap = 0L, genomeSize)  {
  strand(query) <- '*'
  strand(subject) <- '*'
  query <- reduce(query)
  subject <- reduce(subject)
  
  TP <- sum(width(intersect(query, subject)))
  FP <- (sum(width(query))) - TP
  FN <- (sum(width(subject))) - TP
  TN <- genomeSize - (TP + FP + FN)
  
  outMat <- matrix(c(TP, FN, FP, TN), ncol = 2)
  dimnames(outMat) <- list(subject = c("TRUE", "FALSE"), query = c("TRUE", "FALSE"))
  
  stopifnot(!any(outMat < 0))
  
  return(outMat)
}

# Calculate the confusion matrix
conf_matrix <- calcConfMat.GR(query = me2_layer, subject = granges_chr19_jung, genomeSize = genome_size_chr19)
print(conf_matrix)


# Extracting values from the confusion matrix to use for further analysis
TP <- conf_matrix[1, 1]
FN <- conf_matrix[1, 2]
FP <- conf_matrix[2, 1]
TN <- conf_matrix[2, 2]

# Accuracy
acc <- function(TP, TN, FP, FN) { 
  
  return( (TP + TN) / (TP +FP +FN +TN))
  
}

accuracy <- acc(TP, TN, FP, FN)
print(paste("Accuracy:", accuracy))


# True positive rate (sensitivity)
tpr <- function(TP, TN, FP, FN) {   #  true positive rate, TPR a.k.a.sensitivity
  return(TP/(TP+FN))
}

sensitivity <- tpr(TP, TN, FP, FN)
print(paste("Sensitivity (TPR):", sensitivity))

# True negative rate (specificity)
tnr <- function(TP, TN, FP, FN) { #  true negative rate, TNR a.k.a. specificity
  return(TN/(TN+FP))
}

specificity <- tnr(TP, TN, FP, FN)
print(paste("Specificity (TNR):", specificity))

# Positive predictive value (precision)
ppv  <- function(TP, TN, FP, FN) { # positive predictive value, PPV a.k.a. precision
  return(TP/(TP+FP))
}

precision <- ppv(TP, TN, FP, FN)
print(paste("Precision (PPV):", precision))

# Negative predictive value
npv <- function(TP, TN, FP, FN) { # negative predictive value, NPV
  return(TN/(TN+FN))
}

negative_predictive_value <- npv(TP, TN, FP, FN)
print(paste("Negative Predictive Value (NPV):", negative_predictive_value))

# False positive rate
fpr <- function(TP, TN, FP, FN) { #false positive rate
  return(FP/(FP+TN))
}

false_positive_rate <- fpr(TP, TN, FP, FN)
print(paste("False Positive Rate (FPR):", false_positive_rate))

# False negative rate
fnr <- function(TP, TN, FP, FN) { #false negative rate
  return(FN/(TP+FN))
}

false_negative_rate <- fnr(TP, TN, FP, FN)
print(paste("False Negative Rate (FNR):", false_negative_rate))

# False discovery rate
fdr <- function(TP, TN, FP, FN) { #false discovery rate
  return(FP/(TP+FP))
}

false_discovery_rate <- fdr(TP, TN, FP, FN)
print(paste("False Discovery Rate (FDR):", false_discovery_rate))









