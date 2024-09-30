#TERMINAL CODE - CSF3 START TO FINISH:
  
  #LOGIN
ssh -X j52505my@csf3.itservices.manchester.ac.uk
qrsh -l short #this requests an interactive session and you will see this [j52505my@node602 [csf3] ~]$
exit

#DOWNLOADING HUMAN GENOME:
  wget https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz


#LOADING R:
  module load apps/gcc/R/4.3.1
R
install.packages("devtools")
install.packages("getopt")
library(getopt)

qsub -cwd -V -b y -N chr19-optimal Rscript /scratch/j52505my/chr19-optimal.R -w 20000 -e 40000 -o /scratch/j52505my/output_folder/outfile.20000.40000.chr19.1.Rdata


#CREATION OF R SCRIPT FILE:
  vi chr19.R
  

#LOGIN NODE SUBMISSION OF JOBS - from login2 node:
  
  cd /scratch/j52505my
  cd /mnt/iusers01/fatpou01/bmh01/msc-bioinf-2023-2024/j52505my/scripts
  cd /scratch/$maryyuhanna/
  qsub -cwd -V -b y -N chr19-optimal Rscript $HOME/scripts/chr19-4.R    #[j52505my@login1 [csf3] ~]$
qstat #checks status of the job submission
ls #BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz  Hello-World.o5018047  R
#Hello-World.e5018047                      Hello-World.o5018085  scratch
#Hello-World.e5018085                      localscratch          scripts

head Hello-World.o5018085


qsub -cwd -V -b y -N chr19-optimal Rscript /scratch/j52505my/chr19-optimal.R    #[j52505my@login1 [csf3] ~]$


/scratch/j52505my/startingResources/removeShortGRanges.R
/scratch/j52505my/startingResources/removeGRangesBySize.R
/scratch/j52505my/startingResources/plotLayers.R
  




human19_SamplerNucNoSeqAbundSpread_ZeroStart_T

scp j52505my@csf3-server:/scratch/j52505my/output_folder/human19_SamplerNucNoSeqAbundSpread_ZeroStart_T/outfile.20000.40000.chr19.1.nBlocks.png /Users/maryyuhanna/desktop/
  

#this works:
  
qsub -cwd -V -b y -N chr19.finaltest -pe smp.pe 2 Rscript /scratch/j52505my/chr19-optimal.R -w 150000 -e 250000 -o outfile.150000.250000.chr19.finaltest

qsub -cwd -V -b y -N genome-wide2 -pe smp.pe 2 Rscript /scratch/j52505my/genomewide-optimal.R -w 40000 -e 50000 -o outfile.40000.50000.genomewide2

qsub -b y -cwd -m e -M $USER -l short zip test_data.zip Rplots.pdf

  
#This is to save files to my desktop - you need to logout of CSF3 and do this on a normal terminal
scp j52505my@csf3.itservices.manchester.ac.uk:/scratch/j52505my/output_folder/human19_SamplerNucNoSeqAbundSpread_ZeroStart_T/outfile.1000000.2200000.genome.2.100.Rdata.nBlocks.png
 ~/Desktop/
  
  
  
