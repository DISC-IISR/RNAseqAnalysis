#'@title quality check report
#'
#'@description takes fastq.gz file present within the data folder as input and provides the user with various graph related to the quality of the reads
#'@import Rqc
#'@import knitr
#'@importFrom utils browseURL
#'@param samples	It reads a random sample from files if this parameter is TRUE.
#'@param no Number of sequences to read from each input file. This represents sample size if 'sample' parameter is TRUE, if not represents the chunk size to read on each iteration. Default is read a sample of one million sequences from each input file.
#'@param groups group name for each input file.
#'@param tops number of top over-represented reads. Default is 10 reads.
#'@param pairs	combination of files for paired-end reads. By default, all input files are treated as single-end. For paired-end, please define a vector of numbers where two index with the same value represent a pair. Examples, single-end c(1,2,3,4) and paired-end c(1,1,2,2).
#'@param worker number of parallel workers.Default is 1.
#'@details Multiple input files present in the folder data of the workspace are processed in parallel using the function of Rqc. The result generated is saved in folder results/QC-report. All the images generated are stored in jpeg format and also a html report is also generated that can be opened in browser for analysis.
#'@export
#'@seealso \code{\link{rqc}}
#'@examples
#'
#'\dontrun{
#'fqcR()
#'fqcR(pair = c(1,1,2,2,3,3,4,4,5,5,6,6))
#'fqcR(pair = c(1,1,2,2,3,3,4,4),  worker = multicoreWorkers())
#'fqcR(pair = c(1,2,3,4))
#'}
#'

fqcR <- function (samples = TRUE, no = 1e+06, groups = NULL, tops = 10, pairs = NULL, worker=1)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  dir.create("results/QC-report")
  x <- list.files("dataset/", pattern =".fastq.gz", full.names = TRUE)
  p <- pairs
  y <- samples
  n1 <- no
  g <- groups
  t <- tops
  w <- worker

  qa <- rqcQA(x, sample = y, n = n1, group = g, top = t, pair = p, workers = w)

  print(perFileInformation(qa))

  rqcReadQualityBoxPlot(qa)
  ggsave("results/QC-report/read-mean-dist.jpeg")

  rqcReadQualityPlot(qa)
  ggsave("results/QC-report/average-quality-plo.jpeg")

  rqcCycleAverageQualityPlot(qa)
  ggsave("results/QC-report/cycle-average-quality-plot.jpeg")

  rqcReadFrequencyPlot(qa)
  ggsave("results/QC-report/readfrequency.jpeg")

  rqcFileHeatmap(qa[[1]])
  ggsave("results/QC-report/heatmap-reads.jpeg")

  rqcReadWidthPlot(qa)
  ggsave("results/QC-report/read-width-plot.jpeg")

  rqcCycleGCPlot(qa)
  ggsave("results/QC-report/cycle-gc-plo.jpeg")

  rqcCycleAverageQualityPcaPlot(qa)
  ggsave("results/QC-report/biplot.jpeg")


  files <- list.files("data/", pattern =".fastq.gz", full.names = TRUE)

  ## ----rqcQA---------------------------------------------------------------


  qa <- rqcQA(x, pair = y, workers=1)
  reportFile <- rqcReport(qa, outdir = "./results/QC-report")
  browseURL(reportFile)
}
