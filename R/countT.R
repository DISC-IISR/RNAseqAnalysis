#'@title generate Count Table
#'
#'@description the function generates count table using the raw reads and genome sequence. it also provides with a overall differential expression analysis plot of the reads. The count table generated can be used for the downstream analysis using tools like \code{edgeR} and \code{DESeq}
#'@import GenomicFeatures
#'@import GenomicAlignments
#'@import ggplot2
#'@import BiocParallel
#'@import grDevices
#'@import AnnotationDbi
#'@importFrom stats cor
#'@importFrom stats hclust
#'@importFrom stats dist
#'@importFrom utils write.table
#'@importFrom systemPipeR systemArgs
#'@importFrom systemPipeR runCommandline
#'@importFrom systemPipeR outpaths
#'@importFrom systemPipeR filterDEGs
#'@importFrom systemPipeR run_edgeR
#'@importFrom systemPipeR readComp
#'@importFrom systemPipeR returnRPKM
#'@importFrom systemPipeR alignStats
#'@importFrom systemPipeR overLapper
#'@importFrom systemPipeR vennPlot
#'@importFrom easyRNASeq BamFileList
#'@importFrom easyRNASeq assays
#'@importFrom ape plot.phylo
#'@importFrom ape as.phylo
#'@param sysma path to 'param' file; file structure follows a simple name/value syntax that converted into JSON format; for details about the file structure see sample files provided by package. Assign NULL to run the pipeline without 'param' file. This can be useful for running partial workflows, e.g. with pregenerated BAM files.
#'@param mytargets path to targets file
#'@param type type="SYSargs" returns SYSargs, type="json" returns param file content in JSON format (requires rjson library)
#'@param file Input GFF3 or GTF file. Can be a path to a file, or an URL, or a connection object, or a GFF3File or GTFFile object.
#'@param format Format of the input file. Accepted values are: "auto" (the default) for auto-detection of the format, "gff3", or "gtf". Use "gff3" or "gtf" only if auto-detection failed.
#'@param dataSource A single string describing the origin of the data file. Please be as specific as possible.
#'@param organism What is the Genus and species of this organism. Please use proper scientific nomenclature for example: "Homo sapiens" or "Canis familiaris" and not "human" or "my fuzzy buddy". If properly written, this information may be used by the software to help you out later.
#'@param taxonomyId By default this value is NA and the organism provided will be used to look up the correct value for this. But you can use this argument to override that and supply your own taxonomy id here (which will be separately validated). Since providing a valid taxonomy id will not require us to look up one based on your organism: this is one way that you can loosen the restrictions about what is and isn't a valid value for the organism.
#'@param circ_seqs A character vector to list out which chromosomes should be marked as circular.
#'@param chrominfo Data frame containing information about the chromosomes. Will be passed to the internal call to makeTxDb. See ?makeTxDb for more information. Alternatively, can be a Seqinfo object.
#'@param miRBaseBuild Specify the string for the appropriate build Information from mirbase.db to use for microRNAs. This can be learned by calling supportedMiRBaseBuildValues. By default, this value will be set to NA, which will inactivate the microRNAs accessor.
#'@param metadata A 2-column data frame containing meta information to be included in the TxDb object. See ?makeTxDb for more information about the format of metadata.
#'@param dbxrefTag If not missing, the values in the Dbxref attribute with the specified tag (like “GeneID”) are used for the feature names.
#'@param filter	Named vector with filter cutoffs of format c(Fold=2, FDR=1) where Fold refers to the fold change cutoff (unlogged) and FDR to the p-value cutoff.
#'@param con a formula which expresses how the counts for each gene depend on the variables in colData
#'@export
#'@seealso \code{\link{systemArgs}}
#'@examples
#'
#'\dontrun{
#'library(BSgenome.Celegans.UCSC.ce11)
#'export(Celegans, "ce11.fasta", compress=FALSE)
#'library("biomartr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
#'getGFF(db = "ensembl", organism = "Caenorhabditis elegans", reference = TRUE, path = file.path("inst", "dataset"))
#'system("gunzip inst/dataset/Caenorhabditis_elegans.WBcel235.91_ensembl.gff3.gz")
#'condition <- factor(c("N1","N1","N2","N2","M1","M1","M2","M2","M3","M3","M4","M4","U1","U1"))
#'countT(sysma = "./param/tophat.param", mytargets = "targetsPE.txt", file="dataset/WG.gff3", filter=c(Fold=2, FDR=100),  con = condition)
#'}
#'
countT <- function (sysma, mytargets, type = "SYSargs", file, format=c("auto", "gff3", "gtf"), dataSource=NA, organism=NA, taxonomyId=NA, circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL, miRBaseBuild=NA, metadata=NULL, dbxrefTag, filter, con)
{
 
  args <- systemArgs(sysma, mytargets, type)
  if (interactive()) {
    
    prompt<- "enter 'bowtie2-build' path/to/fasta_file path/to/fasta_file  "
    system(readline(prompt))
  } else {
    system("bowtie2-build ./dataset/ce11.fasta ./dataset/ce11.fasta")
  }
  bampaths <- runCommandline(args=args)

  read_statsDF <- alignStats(args)
  write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")


  txdb <- makeTxDbFromGFF(file, format, dataSource, organism, taxonomyId, circ_seqs, chrominfo, miRBaseBuild, metadata, dbxrefTag)
  saveDb(txdb, file="dataset/database.sqlite") # needs to be done only once
  txdb <- loadDb("./dataset/database.sqlite")
  if (interactive()) {
    eByg <- exonsBy(txdb, by=c(readline(prompt="provide one of 'gene', 'exon', 'cds' or 'tx' to determines the grouping.")))
  } else {
    eByg <- exonsBy(txdb, by=c("tx"))
  }

  bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
  if (interactive()) {
    counteByg <- summarizeOverlaps(eByg, bfl, mode="Union", if(readline(prompt="enter FALSE for strand specific RNA /TRUE for all other")=="TRUE"){ignore.strand=TRUE} else {ignore.strand=FALSE}, inter.feature=FALSE, if(readline(prompt="enter TRUE for SINGLE END RNA /FALSE for PAIREND ")=="TRUE"){singleEnd=TRUE} else {singleEnd=FALSE})
  } else {
    counteByg <- summarizeOverlaps(eByg, bfl, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=FALSE)
  }
  # Note: for strand-specific RNA-Seq set 'ignore.strand=FALSE' and for PE data set 'singleEnd=FALSE'

  countDFeByg <- assays(counteByg)$counts
  countDFeByg <- as.data.frame(countDFeByg)
  countDFeByg[is.na(countDFeByg)] <- 0
  write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
  
  rpkmDFeByg <- apply(countDFeByg, 2, function(bfl) returnRPKM(counts=bfl, ranges=eByg))
  rpkmDFeByg <- as.data.frame(rpkmDFeByg)
  rpkmDFeByg[is.na(rpkmDFeByg)] <- 0
  write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
  

  countDF <- read.table("results/countDFeByg.xls")

  cmp <- readComp(args, format="matrix", delim="-")

  #degseqDF <- run_DESeq2(countDF=countDF, targets=targetsin(args), cmp=cmp[[1]], independent=FALSE)
  #write.table(degseqDF, "./results/DESeq2comp.xls", quote=FALSE, sep="\t", col.names = NA)
  #DEG_list2 <- filterDEGs(degDF=degseqDF, filter)
  edgeDF <- run_edgeR(countDF=countDF, targets=targetsin(args), cmp=cmp[[1]], independent=FALSE, mdsplot="")

  write.table(edgeDF, "results/edgeRcomp.xls", quote=FALSE, sep="\t", col.names = NA)
  jpeg("results/DEG_edgeR.jpeg")
  DEG_list <- filterDEGs(degDF=edgeDF, filter, plot = TRUE)
  dev.off()

  #edgeDF <- read.table("./results/edgeRcomp.xls")
  #DEG_list <- filterDEGs(degDF=edgeDF, filter, plot = FALSE)

  x <- 1:20
  if (interactive()) {
    for (val in x)
    {
      v1 <- readline("Number of starting for vennplot: ")
      v2 <- readline("Number of last point for vennplot: ")
      vennsetup <- overLapper(DEG_list$Up[v1:v2], type="vennsets")
      vennsetdown <- overLapper(DEG_list$Down[v1:v2], type="vennsets")
      jpeg("results/Venndiagm.jpeg")
      vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
      dev.off()

      n <- readline(prompt="Enter 0 TO STOP EXECUTION/ 1 TO CONTINUE")
      if (n == 0)
      {
        break
      }
    }

  } else {
    vennsetup <- overLapper(DEG_list$Up[1:4], type="vennsets")
    vennsetdown <- overLapper(DEG_list$Down[1:4], type="vennsets")
    jpeg("results/venndiagm_1_4.jpeg")
    vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
    dev.off()

    vennsetup <- overLapper(DEG_list$Up[3:6], type="vennsets")
    vennsetdown <- overLapper(DEG_list$Down[3:6], type="vennsets")
    jpeg("results/venndiagm_3_6.jpeg")
    vennPlot(list(vennsetup, vennsetdown), mymain="", mysub="", colmode=2, ccol=c("blue", "red"))
    dev.off()
  }


  rpkmDFeByg <- read.table("results/rpkmDFeByg.xls")
  d <- cor(rpkmDFeByg, method="spearman")
  hc <- hclust(dist(1-d))
  jpeg("results/Dendrogram_based_on_RPKM_values.jpeg")
  plot.phylo(as.phylo(hc), type="p", edge.color=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
  dev.off()
  
  
  args <- systemArgs(sysma, mytargets, type)
  countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
  colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
  condition <- con
 #dds <- DESeqDataSetFromMatrix(countDF, DataFrame(condition), ~ condition)
 #dds <- DESeq(dds)
  dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 5, ]
cts = counts(dds)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds = estimateSizeFactors(dds, geoMeans=geoMeans)
  
# dds <-exp(sum(log(dds[dds != 0]))/length(dds))
  d <- cor(assay(rlog(dds, blind=FALSE)), method="spearman")
  hc <- hclust(dist(1-d))
  jpeg("results/Dendrogram_based_on_rlog_values.jpeg")
  plot.phylo(as.phylo(hc), type="p", edge.color=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
  dev.off()
}
