library(GenomicLayers)
library(GenomicRanges)
library(MotifDb)
library(Biostrings)
library(seqLogo)

## DNMT3 adds DNA methylation (adds to two cpgs) - cant bind when H3K4me3 isnt there?
dnmt3.bf <- createBindingFactor.DNA_consensus(name='dnmt3.bf', 
                                              patternString = 'CG',
                                              profile.layers = c('H3K4me3', 'H3K4me2'), 
                                              profile.marks = c(0,0), 
                                              mod.layers = 'DNA_methyl', 
                                              mod.marks = 1)

## MLL
HNF4.motif <- query(MotifDb, andStrings = c('hnf4', 'musculus'),
                    orStrings = c('hocomoco', 'jaspar', 'jolma', 'uniprobe'),
                    notStrings = c('hnf4g'))

seqLogo(HNF4.motif[[1]]) # 
seqLogo(HNF4.motif[[2]])
seqLogo(HNF4.motif[[3]])
seqLogo(HNF4.motif[[4]])
seqLogo(HNF4.motif[[5]])

pwm.hnf4 <- HNF4.motif[[1]]
hnf4a.bf <- createBindingFactor.DNA_motif(name = 'hnf4a.bf', 
                                          pwm = pwm.hnf4, 
                                          min.score = '80%', 
                                          with.score = FALSE, 
                                          mod.layers =  'mll_tfs', 
                                          mod.marks = 1 ) 

mll.bf <- createBindingFactor.layer_region(name = 'mll.bf',
                                           patternLength = 2,
                                           patternString = 'N',
                                           stateWidth = 147, 
                                           profile.layers = "mll_tfs",
                                           profile.marks = 1,       
                                           mod.layers = "H3K4me1", 
                                           mod.marks = 1)

## COMPASS which binds CpGs and adds H3K4me3 - binds when no DNA methylation
cfp1.motif2 <- query(MotifDb, andStrings = c('cxxc', 'musculus'))
length(cfp1.motif2)
pwm.cfp12 <- cfp1.motif2[[1]]
seqLogo(pwm.cfp12)
Cfp1.bf <- createBindingFactor.DNA_motif(name = 'Cfp1.bf', 
                                      pwm = pwm.cfp12, 
                                      min.score = '80%', 
                                      with.score = FALSE, 
                                      mod.layers =  'active_tf', 
                                      mod.marks = 1 ) 

## add sequential H3K4me by set1
set1_meUP <- createBindingFactor.layer_region(name = 'set1_meUP', patternLength = 6, profile.layers = c('active_tf', 'DNA_methyl', 'H3K4me3'),
                                              profile.marks = c(1,0,0), mod.layers = 'me_promotion', 
                                              mod.marks = 1, stateWidth = 147)

set1.promotion.0_1 <- createBindingFactor.layer_region(name = 'set1.promotion.0_1', patternLength = 147, 
                                                      profile.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                      profile.marks = c(1,0,0,0), 
                                                      mod.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                      mod.marks = c(0,1,0,0,1))

set1.promotion.1_2 <- createBindingFactor.layer_region(name = 'set1.promotion.1_2', patternLength = 147, 
                                                       profile.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                       profile.marks = c(1,1,0,0), 
                                                       mod.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                       mod.marks = c(0,0,1,0,1))

set1.promotion.2_3 <- createBindingFactor.layer_region(name = 'set1.promotion.2_3', patternLength = 147, 
                                                       profile.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                       profile.marks = c(1,0,1,0), 
                                                       mod.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'),
                                                       mod.marks = c(0,0,0,1,1)) 

set1.promotion.3_3 <- createBindingFactor.layer_region(name = 'set1.promotion.3_3', patternLength = 147,  
                                                       profile.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                       profile.marks = c(1,0,0,1), 
                                                       mod.layers = c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'),
                                                       mod.marks = c(0,0,0,1,1))

## sequential H3K4 demethylation
kdm5.meDOWN <- createBindingFactor.layer_region(name = 'kdm5.meDOWN', profile.layers = 'H3K4me_any', 
                                                patternLength = 147, profile.marks = 1, 
                                                mod.layers = 'H3K4me_demethylation', mod.marks = 1, 
                                                stateWidth = 147)

kdm5.demethylation.3_2 <- createBindingFactor.layer_region(name = 'kdm5.demethylation.3_2', patternLength = 147, 
                                                           profile.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                           profile.marks = c(1,0,0,1), 
                                                           mod.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                           mod.marks = c(0,0,1,0,1))

kdm5.demethylation.2_1 <- createBindingFactor.layer_region(name = 'kdm5.demethylation.2_1', patternLength = 147, 
                                                           profile.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                           profile.marks = c(1,0,1,0), 
                                                           mod.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                           mod.marks = c(0,1,0,0,1))

kdm5.demethylation.1_0 <- createBindingFactor.layer_region(name = 'kdm5.demethylation.1_0', patternLength = 147, 
                                                           profile.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                           profile.marks = c(1,1,0,0), 
                                                           mod.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                           mod.marks = c(0,0,0,0,0))

kdm5.demethylation.0_0 <- createBindingFactor.layer_region(name = 'kdm5.demethylation.0_0', patternLength = 147, 
                                                           profile.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3'), 
                                                           profile.marks = c(1,0,0,0), 
                                                           mod.layers = c('H3K4me_demethylation', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any'), 
                                                           mod.marks = c(0,0,0,0,0))

## genome
library(BSgenome.Mmusculus.UCSC.mm10)
genome1 <- BSgenome.Mmusculus.UCSC.mm10
seqinfo(genome1)
seqnames(genome1)
sequences_to_keep <- "chr19"   # nuclear chroms only
genomeNuc <- keepBSgenomeSequences(genome1, sequences_to_keep)
seqinfo(genomeNuc)
nucLayerSet_2 <- createLayerSet.BSgenome(genome=genomeNuc, 
                                         layer.names=c('active_tf','me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'DNA_methyl', 'H3K4me_any', 'H3K4me_demethylation'), 
                                         n.layers=8, verbose = TRUE)


## adding methylation
bf.list <- list(Cfp1.bf=Cfp1.bf,
                dnmt3.bf=dnmt3.bf,
                set1_meUP=set1_meUP,
                set1.promotion.0_1=set1.promotion.0_1,
                set1.promotion.1_2=set1.promotion.1_2,
                set1.promotion.2_3=set1.promotion.2_3,
                set1.promotion.3_3=set1.promotion.3_3)

saturationAbundance <- 1000000
bf_abundances_2 <- c(50000,25000,100000, rep(saturationAbundance, 4))
names(bf_abundances_2) <- names(bf.list)
bf_abundances_2

matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = Cfp1.bf)
matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = dnmt3.bf)
matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = set1_meUP)
matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = set1.promotion.0_1)

names(nucLayerSet_2$layerSet)
layersToClean <- c('me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3',  'H3K4me_any', 'H3K4me_demethylation') 

for (i in 1:200) {
  print(i)
  nucLayerSet_2 <- runLayerBinding.BSgenome(layerList = nucLayerSet_2, 
                                            factorSet = bf.list, 
                                            bf.abundances = bf_abundances_2, 
                                            verbose = T, collect.stats = T, keep.stats = T)
  plotLayers(nucLayerSet_2, layerNames=c("H3K4me1", "H3K4me2","H3K4me3"), chrom="chr19", xlim=c(3762000,3772000))
    
  # TIDY UP BY REMOVING HISTONE MODS < 147BP..  # I'M NOT SURE THIS IS A GOOD IDEA.
  for(thisLayer in layersToClean) {
    nucLayerSet_2$layerSet[[thisLayer]]  <- removeShortGRanges(x=nucLayerSet_2$layerSet[[thisLayer]], minSize = 147)
  }
  
}

nucLayerSet_2$history

plot(nucLayerSet_2$history$Coverage.DNA_methyl,  pch=24,ylab="Coverage", xlab="iteration", ylim=c(0, 3e+07))
points(nucLayerSet_2$history$Coverage.H3K4me2, pch=22)
points(nucLayerSet_2$history$Coverage.H3K4me1, pch=23)
points(nucLayerSet_2$history$Coverage.H3K4me3, pch=21)
legend("topleft", legend=c("H3K4me3", "H3K4me2", "H3K4me1", "DNA_methyl"), pch=c(21, 22, 23, 24))

## adding demethylation
bf.list_demethyl <- list(kdm5.meDOWN=kdm5.meDOWN,
                         kdm5.demethylation.3_2=kdm5.demethylation.3_2,
                         kdm5.demethylation.2_1=kdm5.demethylation.2_1,
                         kdm5.demethylation.1_0=kdm5.demethylation.1_0,
                         kdm5.demethylation.0_0=kdm5.demethylation.0_0)

bf_abundances_3 <- c(25000, rep(saturationAbundance, 4))
names(bf_abundances_3) <- names(bf.list_demethyl)

matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = kdm5.meDOWN)
matchBindingFactor.BSgenome(layerSet = nucLayerSet_2, bindingFactor = kdm5.demethylation.3_2)

nucLayerSet_down <- nucLayerSet_2


for (i in 1:50) {
  print(i)
  nucLayerSet_down <- runLayerBinding.BSgenome(layerList = nucLayerSet_down, 
                                               factorSet = bf.list_demethyl, 
                                               bf.abundances = bf_abundances_3, 
                                               verbose = T, collect.stats = T, keep.stats = T)

  plotLayers(nucLayerSet_down, layerNames=c("H3K4me1", "H3K4me2","H3K4me3"), chrom="chr19", xlim=c(3762000,3772000))
  
  # TIDY UP BY REMOVING HISTONE MODS < 147BP..  # I'M NOT SURE THIS IS A GOOD IDEA.
  for(thisLayer in layersToClean) {
    nucLayerSet_down$layerSet[[thisLayer]]  <- removeShortGRanges(x=nucLayerSet_down$layerSet[[thisLayer]], minSize = 147)
  }
}

nucLayerSet_2$history
nucLayerSet_down$history

plot(nucLayerSet_down$history$Coverage.H3K4me1, ylab="Coverage", xlab="iteration")
points(nucLayerSet_down$history$Coverage.H3K4me2, pch=22)
points(nucLayerSet_down$history$Coverage.H3K4me3, pch=23)
legend("topleft", legend=c("H3K4me3", "H3K4me2", "H3K4me1"), pch=c(21, 22, 23))

## both H3K4 methylation and demethylation
bf.list_both <- c(bf.list, bf.list_demethyl)
bf_abundances_both <- c(bf_abundances_2, bf_abundances_3)

nucLayerSet_both <- nucLayerSet_2

for(i in 101:250) {
  print(i)
  nucLayerSet_both <- runLayerBinding.BSgenome(layerList = nucLayerSet_both, factorSet = bf.list_both, bf.abundances = bf_abundances_both, verbose = T, collect.stats = T, keep.stats = T)

  plotLayers(nucLayerSet_both, layerNames=c("H3K4me1", "H3K4me2","H3K4me3"), chrom="chr19", xlim=c(3762000,3772000))
  
  # TIDY UP BY REMOVING HISTONE MODS < 147BP..  # I'M NOT SURE THIS IS A GOOD IDEA.
  for(thisLayer in layersToClean) {
    nucLayerSet_both$layerSet[[thisLayer]]  <- removeShortGRanges(x=nucLayerSet_both$layerSet[[thisLayer]], minSize = 147)
  }
}

par(mar=c(5,7,3,2))
plot(nucLayerSet_both$history$Coverage.DNA_methyl,  pch=24,ylab="", xlab="iteration", ylim=c(0, 3e+07), las=2)
title(ylab="Coverage", line=5)
points(nucLayerSet_both$history$Coverage.H3K4me2, pch=22)
points(nucLayerSet_both$history$Coverage.H3K4me1, pch=23)
points(nucLayerSet_both$history$Coverage.H3K4me3, pch=21)
legend("topleft", legend=c("H3K4me3", "H3K4me2", "H3K4me1", "DNA_methyl"), pch=c(21, 22, 23, 24))
abundanceLegend <- c("Cfp1.bf", "dnmt3.bf",  "set1_meUP", "kdm5.meDOWN")
legend("top", legend=paste(abundanceLegend, bf_abundances_both[abundanceLegend]))


for(i in 151:200) {
  print(i)
  nucLayerSet_both <- runLayerBinding.BSgenome(layerList = nucLayerSet_both, factorSet = bf.list_both, bf.abundances = bf_abundances_both, verbose = T, collect.stats = T, keep.stats = T)
  plot(nucLayerSet_both$history$Coverage.H3K4me3, ylab="Coverage", xlab="iteration", main=paste("iteration", i))
  points(nucLayerSet_both$history$Coverage.H3K4me2, pch=22)
  points(nucLayerSet_both$history$Coverage.H3K4me1, pch=23)
  legend("topleft", legend=c("H3K4me3", "H3K4me2", "H3K4me1"), pch=c(21, 22, 23))
}

## MLL3/4 - binds with some level of DNA methyaltion and CpGs

## MLL3/4 adds H3K4me1 - doesnt care about DNA methylation (just assume for now)
HNF4.motif <- query(MotifDb, andStrings = c('hnf4', 'musculus'),
                    orStrings = c('hocomoco', 'jaspar', 'jolma', 'uniprobe'),
                    notStrings = c('hnf4g'))

seqLogo(HNF4.motif[[1]]) # 
seqLogo(HNF4.motif[[2]])
seqLogo(HNF4.motif[[3]])
seqLogo(HNF4.motif[[4]])
seqLogo(HNF4.motif[[5]])

pwm.hnf4 <- HNF4.motif[[1]]
hnf4a.bf <- createBindingFactor.DNA_motif(name = 'hnf4a.bf', 
                                          pwm = pwm.hnf4, 
                                          min.score = '80%', 
                                          with.score = FALSE, 
                                          mod.layers =  'active.tfs', 
                                          mod.marks = 1 ) 

mll.bf <- createBindingFactor.layer_region(name = 'mll.bf',
                                           patternLength = 2,
                                           patternString = 'N',
                                           stateWidth = 147, 
                                           profile.layers = "active.tfs",
                                           profile.marks = 1,       
                                           mod.layers = "H3K4me1", 
                                           mod.marks = 1)

## add planned changes 














## test: lots of dna methylation will damage promoters because set1 cant bind, adding histone marks to high density cg sites protects them from methylation
## use specificity or senisitivity, or overlap of tss and enhancers 
## list of layers used 
H3K4me3
DNA_methyl

H3K4me1
active.tfs
