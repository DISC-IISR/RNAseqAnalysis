## ---- eval = FALSE---------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("Blessy/RNAseqAnalysis")

## ---- eval = FALSE---------------------------------------------------------
#  library("RNAseqAnalysis") # Loads the package
#  library(help="RNAseqAnalysis") # Lists package info
#  vignette("RNAseqAnalysis") # Opens vignette

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
#library("RNAseqAnalysis")
#setwd("~/RNAseqAnalysis")
#import_dataset()

## --------------------------------------------------------------------------
library(RNAseqAnalysis)
targetspath <- "targets.txt"
read.delim(targetspath, comment.char = "#")

## --------------------------------------------------------------------------
library(RNAseqAnalysis)
targetspath <- "targetsPE.txt"
read.delim(targetspath, comment.char = "#")

## --------------------------------------------------------------------------
library(RNAseqAnalysis)
targetspath <- "targetsPE.txt"
readLines(targetspath)

## --------------------------------------------------------------------------
library(RNAseqAnalysis)
library(systemPipeR)
targetspath <- "targetsPE.txt"
readComp(file=targetspath, format="vector", delim="-")

## --------------------------------------------------------------------------
library(RNAseqAnalysis)
parampath <- "tophat.param"
read.delim(parampath, comment.char = "#")

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
files <- list.files(path = "./dataset/", pattern = ".fastq.gz", all.files = FALSE, full.names = TRUE)
fqcR(groups=rep("None", 28), pairs=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14))

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
frqc()

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
filteR()

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library("biomartr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
library(BSgenome.Celegans.UCSC.ce11)
export(Celegans, "dataset/ce11.fasta", compress=FALSE)
getGFF(db = "ensemblgenomes", organism = "Caenorhabditis elegans", reference = TRUE, path = file.path("dataset"))
system("gunzip dataset/Caenorhabditis_elegans.WBcel235.39_ensemblgenomes.gff3.gz")
library("rtracklayer", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
gffRangedData <- import.gff("dataset/Caenorhabditis_elegans.WBcel235.39_ensemblgenomes.gff3")
myGranges<-as(gffRangedData, "GRanges")
seqlevelsStyle(myGranges)
seqlevelsStyle(myGranges) <- "UCSC"
myGranges
export(object = myGranges, con = "dataset/Caenorhabditis_elegans.gff3", format = "gff3", index = FALSE)
system("gunzip dataset/Caenorhabditis_elegans.gff3.bgz")
condition <- factor(c("N1","N1","N2","N2","M1","M1","M2","M2","M3","M3","M4","M4","U1","U1"))
countT(sysma = "./param/tophat.param", mytargets = "targetsPE.txt", file="dataset/Caenorhabditis_elegans.gff3", filter = c(Fold=2, FDR=100), con = condition)

## ---- results = "hide", warning = FALSE, message = FALSE, fig.show = 'hide'----
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
pairs <- c( "paired-end", "paired-end", "paired-end","paired-end","paired-end","paired-end","paired-end", "paired-end", "paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
count <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
condition <- factor(c("N1","N1","N2","N2","M1","M1","M2","M2","M3","M3","M4","M4","U1","U1"))
gene_exp(sysma="./param/tophat.param", mytargets="targetsPE.txt", x = pairs, y = count, con = condition)

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
annotate(filter=c(Fold=2, FDR=100), id_type="gene", CLSZ=2, cutoff=0.01)

## ---- results = "hide", warning = FALSE, message = FALSE-------------------
library(RNAseqAnalysis)
setwd("~/RNAseqAnalysis")
bamfiles <- c("extdata/SRR5253683.bam", "extdata/SRR5253684.bam", "extdata/SRR5253685.bam", "extdata/SRR5253686.bam", "extdata/SRR5253687.bam", "extdata/SRR5253688.bam", "extdata/SRR5253689.bam", "extdata/SRR5253690.bam", "extdata/SRR5253691.bam", "extdata/SRR5253692.bam", "extdata/SRR5253694.bam", "extdata/SRR5253695.bam","extdata/SRR5253696.bam", "extdata/SRR5253697.bam")
col <- c("red","yellow","green","blue", "black", "pink", "violet")
grp <- c(rep("N1", 2),rep("N2", 2),rep("M1", 2),rep("M2", 2),rep("M3", 2),rep("M4", 2),rep("U1", 2))
#grp <- c(rep("U1", 2),rep("M3", 2))
splice(pair = c("U1", "M3"), group = grp, colour = col, bamFiles=bamfiles)

## --------------------------------------------------------------------------
sessionInfo()

