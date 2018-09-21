#'@title alternative splicing analysis
#'
#'@description Its an Integrative function for the analysis of alternative splicing. it provides with Alternative splicing discovery using junctions and DE and DU estimation.
#'@import ASpli
#'@import GenomicFeatures
#'@import BiocParallel
#'@import AnnotationDbi
#'@param pairs	Vector of length two, either numeric or character, providing the pair of groups to be compared
#'@param group	Factorial vector with tags for each sample
#'@param colour Vector containing color for each condition
#'@param bamFiles The path to bam files
#'@export
#'@seealso \code{\link{ASpli}}
#'@examples
#'
#'\dontrun{
#' bamfiles <- c("extdata/A83.bam", "extdata/A84.bam", "extdata/A85.bam", "extdata/A86.bam", "extdata/A94.bam", "extdata/A95.bam")
#' bamfiles <- c("extdata/SRR5253683.bam", "extdata/SRR5253684.bam", "extdata/SRR5253685.bam",
#'  "extdata/SRR5253686.bam", "extdata/SRR5253687.bam", "extdata/SRR5253688.bam","extdata/SRR5253689.bam",
#'  "extdata/SRR5253690.bam","extdata/SRR5253691.bam","extdata/SRR5253692.bam","extdata/SRR5253694.bam",
#' "extdata/SRR5253695.bam","extdata/SRR5253696.bam","extdata/SRR5253697.bam")
#' col <- c("red","yellow","green")
#' col <- c("red","yellow","green","blue", "black", "pink", "violet")
#' grp <- c(rep("U1", 2),rep("M3", 2),rep("N1", 2))
#' grp <- c(rep("N1", 2),rep("N2", 2),rep("M1", 2),rep("M2", 2),rep("M3", 2),rep("M4", 2),rep("U1", 2))
#' splice(pairs = c("U1", "M3"), group = grp, colour = col, bamFiles=bamfiles)
#'}
#'
#'
#'\dontrun{
#' bamfiles <- c("extdata/A83.bam", "extdata/A84.bam", "extdata/A85.bam", "extdata/A86.bam", "extdata/A94.bam", "extdata/A95.bam")
#' bamfiles <- c("extdata/SRR5253683.bam", "extdata/SRR5253684.bam", "extdata/SRR5253685.bam",
#'  "extdata/SRR5253686.bam", "extdata/SRR5253687.bam", "extdata/SRR5253688.bam","extdata/SRR5253689.bam",
#'  "extdata/SRR5253690.bam","extdata/SRR5253691.bam","extdata/SRR5253692.bam","extdata/SRR5253694.bam",
#' "extdata/SRR5253695.bam","extdata/SRR5253696.bam","extdata/SRR5253697.bam")
#' col <- c("red","yellow","green")
#' col <- c("red","yellow","green","blue", "black", "pink", "violet")
#' grp <- c(rep("U1", 2),rep("M3", 2),rep("N1", 2))
#' grp <- c(rep("N1", 2),rep("N2", 2),rep("M1", 2),rep("M2", 2),rep("M3", 2),rep("M4", 2),rep("U1", 2))
#' splice(pairs = c("U1", "M3"), group = grp, colour = col, bamFiles=bamfiles)
#'}

splice <- function (pairs, group, colour, bamFiles)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  
  dir.create("results/splice_result")
  TxDb <- loadDb("./dataset/database.sqlite")
  
  features <- binGenome(TxDb)
  GenesCoord <- featuresg(features)
  BinsCoord <- featuresb(features)
  JunctionsCoord <- featuresj(features)
  
  if (interactive()) {
    n <- readline(prompt="Enter an integer: ")
    n <- as.integer(n)
    i <- 1
    bkg <- list()
    while (i<=n) {
      variable <- readline(prompt="Enter variable: ")
      bkg <- c(as.list(bkg), variable)
      i <- i + 1
    }
    bkg1 <- rep(bkg, each=2)
  } else {
    bkg1 <- rep(c("N1", "N2", "M1", "M2", "M3", "M4", "U1"), each=2)
  }
  
  targets <- data.frame(bam=bamFiles, condition=paste(bkg1, sep="."), row.names=sub("\\.bam$", "", bamFiles))
  bam <- loadBAM(targets)
  counts <- readCounts(features, bam, l = 100L, targets = targets, maxISize = 5000)
  GeneCounts <- countsg(counts)
  GeneRd <- rdsg(counts)
  BinCounts <- countsb(counts)
  BinRd <- rdsb(counts)
  JunctionCounts <- countsj(counts)
  e1iCounts <- countse1i(counts)
  ie2Counts <- countsie2(counts)
  
  if (interactive()) {
    pair <- bkg
  } else {
    pair <- c("N1", "N2", "M1", "M2", "M3", "M4", "U1")
  }
  as <- AsDiscover(counts, targets, features, bam, l=100L, pair=pair)
  irPIR <- irPIR(as)
  altPSI <- altPSI(as)
  esPSI <- esPSI(as)
  junctionsPIR <- junctionsPIR(as)
  junctionsPSI <- junctionsPSI(as)
  
  du <- DUreport(counts, targets, pairs, group)
  writeDU(du, output.dir="results/splice_result/example")
  genesde <- genesDE(du)
  binsdu <- binsDU(du)
  junctionsdu <- junctionsDU(du)
  
  writeCounts(counts, "results/splice_result/counts")
  writeRds(counts, "results/splice_result/rds")
  writeDU(du, output.dir="results/splice_result/du");
  writeAS(as=as, output.dir="results/splice_result/as");
  writeAll(counts=counts, du=du, as=as, output.dir="results/splice_result/all")
  
  targetsPlot <- data.frame(bam=bamFiles, sample=sub("\\.bam$", "", bamFiles), color=colour, stringsAsFactors=FALSE)
  bins <- binsDU(du)
  topTagsBins <- which(bins$bin.fdr <= 0.1 & abs(bins$logFC) >=0.58)
  auxdf<-bins[topTagsBins,]
  plotTopTags(auxdf, TxDb, targetsPlot, output.dir="results/splice_result/Plots")
  
}