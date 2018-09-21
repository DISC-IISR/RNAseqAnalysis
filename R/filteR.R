#'@title Trimming and filtering
#'
#'@description trim and filter raw reads containing N, trim as soon as 2 of 5 nucleotides has quality encoding less than "4" (phred score 20) and drop reads that are less than 36nt
#'@import ShortRead
#'@param a a character(1) with nchar(a) == 1L giving the letter at or below which a nucleotide is marked as failing. The default is "4" (phred score 20).
#'@param b drop reads that are less than the number of nucleotides given in b. Default is 36nt
#'@details Filters and trims each user selected raw data file (in .fastq.gz format), based on three parameters, raw reads containing N are removed, reagins having a phred score less than the value defined by user or less tha a phred score of 20 (by default) is trimmed out. As well as sequence reads that have no of bases less than the value given by user or less than 36nt (default value) is removed.
#'@export
#'@seealso \code{\link{filterFastq}}
#'@examples
#'
#'\dontrun{
#'filteR()
#'filteR(a = "4", b = 38)
#'}
#'

filteR <- function (a = "4", b = 36) {

  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  setwd("results")
  dir.create("filtered_data")
  setwd("filtered_data")
  x <- 1:20
  for (val in x) {
    if (interactive()) {
      fl <- file.choose()
    } else {
      fl <- "/home/disc/new_test/inst/dataset/SRR5253683_1.fastq.gz"
    }

    if (interactive()) {
      destination <- readline(prompt="Enter name for output file ")
    } else {
      destination <- "SRR5253683_1.fastq.gz"
    }


      stream <- open(FastqStreamer(fl))
      on.exit(close(stream))
      ft <- readFastq(fl)
      print(ft)

      repeat {
        ## input chunk
        fq <- yield(stream)
        if (length(fq) == 0)
          break
        ## An Introduction to ShortRead 5
        ## trim and filter, e.g., reads cannot contain 'N'...
        fq <- fq[nFilter()(fq)] # see ?srFilter for pre-defined filters
        ## trim as soon as 2 of 5 nucleotides has quality encoding less
        ## than "4" (phred score 20)
        fq <- trimTailw(fq, 2, a = a, 2)
        ## drop reads that are less than 36nt
        fq <- fq[width(fq) >= b]

        ## append to destination
        writeFastq(fq, destination, "a")
      }

      #fs <- readFastq(fq)
     # print(fs)

      if (interactive()) {
        n <- readline(prompt="Enter 0 TO STOP EXECUTION/ 1 TO CONTINUE: ")
      } else {
        n <- 0
      }
    if (n == 0)
    {
      break
    }
      
  }
  setwd("..")
}

