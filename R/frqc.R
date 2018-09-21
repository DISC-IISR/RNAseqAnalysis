#'@title qc
#'
#'@description reports basic quality and summary statistics on FASTQ files, including base and quality distribution by position, sequence length distribution, k-mers by position
#'@import qrqc
#'@import ShortRead
#'@export
#'@seealso \code{\link{filterFastq}}
#'@examples
#'
#'\dontrun{
#'frqc()
#'}
#'

frqc <- function () {
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  dir.create("results/QC-report")
  x <- 1:20

  for (val in x) {
    if (interactive())
      {
      fl <- file.choose()
      setwd("results/QC-report")
       dir_name <- readline(prompt="enter a output_dir name for the selected sample - ")
      dir.create(dir_name)
      setwd(dir_name)
    }
    else {

      fl <- "./dataset/SRR5253683_1.fastq.gz"

      }
    s.fastq <- readSeqFile(fl)
    print(s.fastq)
    qualPlot(s.fastq)
    ggsave("A_plot_of_quality_by_base_position.jpeg")
    basePlot(s.fastq)
    ggsave("Base_frequencies_by_position_in_sequence.jpeg")
    basePlot(s.fastq, bases=c("G", "C", "A", "T"), geom="bar", type="proportion")
    ggsave("Base_proportions_by_position_in_sequence.jpeg")
    gcPlot(s.fastq) + geom_hline(yintercept=0.5, color="purple")
    ggsave("GC_content_by_position.jpeg")
    kmerKLPlot(s.fastq)
    ggsave("K-L_terms_for_a_subset_of_top_k-mers.jpeg")

   if (interactive()) {
     setwd("..")
     setwd("..")
     setwd("..")
    n <- readline(prompt="Enter 0 TO STOP EXECUTION/ 1 TO CONTINUE: ")
  }
  else {
    n = 0
    setwd("..")
  }

   if (n == 0)
    {
      break
    }
  }
  return()
}
