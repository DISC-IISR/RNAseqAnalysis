#'@title Annotation based on GO terms
#'
#'@description Provides user with GO annotation and graphs for all the three GO terms i.e., Molecular Function, Biological Process and Cellular Components
#'@import biomaRt
#'@import grDevices
#'@importFrom systemPipeR filterDEGs
#'@importFrom utils write.table
#'@importFrom systemPipeR GOHyperGAll
#'@importFrom systemPipeR makeCATdb
#'@importFrom systemPipeR GOCluster_Report
#'@importFrom systemPipeR goBarplot
#'@importFrom systemPipeR run_edgeR
#'@param filter Named vector with filter cutoffs of format c(Fold=2, FDR=10) where Fold refers to the fold change cutoff (unlogged) and FDR to the p-value cutoff, default value is Fold=2, FDR=10.
#'@param id_type specifies type of IDs in input, can be assigned gene or affy, default is gene.
#'@param CLSZ minimum gene set (cluster) size to consider. Gene sets below this cutoff will be ignored, default value is 2.
#'@param cutoff p-value cutoff for GO terms to show in result data.frame, default is 0.01.
#'@export
#'@seealso \code{\link{useMart}}
#'@examples
#'
#'\dontrun{
#' annotate()
#' annotate(filter=c(Fold=2, FDR=70), id_type="gene", CLSZ=2, cutoff=0.01)
#'}
#'

annotate <- function (filter=c(Fold=2, FDR=10), id_type="gene", CLSZ=2, cutoff=0.01)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  dir.create("dataset/GO")
  dir.create("results/GO")

  if (interactive()) {
    ensembl=useMart(readline(prompt="enter a Mart name "), host = readline(prompt="enter a host id for ur mart"))
  } else {
    ensembl=useMart("metazoa_mart", host = "metazoa.ensembl.org")
  }
  
  print(listDatasets(ensembl))
  if (interactive()) {
    m <- useMart(readline(prompt="enter a Mart name as above"), host = readline(prompt="enter a host id as above"),  dataset=(readline(prompt="enter a datasets from the one listed above : ")))
  } else {
    m <- useMart(host = "metazoa.ensembl.org", biomart = "metazoa_mart", dataset = "celegans_eg_gene")
  }
  print(listAttributes(m))
  if (interactive()) {
    go <- getBM(attributes = c(readline(prompt="enter an atrribute for go id - "), readline(prompt="enter an attribute for gene id corresponding to that of your data - "), readline(prompt="enter an attribute for namespace - ")), mart=m)
  } else {
    go <- getBM(attributes = c("go_id", "ensembl_transcript_id", "namespace_1003"), mart = m)
  }
  go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
  write.table(go, "dataset/GO/GOannotationsBiomart.txt", append = FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  catdb <- makeCATdb(myfile="./dataset/GO/GOannotationsBiomart.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
  save(catdb, file="dataset/GO/catdb.RData")
  load("./dataset/GO/catdb.RData")
  edgeDF <- read.table("./results/edgeRcomp.xls")
  DEG_list <- filterDEGs(degDF=edgeDF, filter, plot=FALSE)
  up_down <- DEG_list$UporDown
  names(up_down) <- paste(names(up_down), "_up_down", sep="")
  up <- DEG_list$Up
  names(up) <- paste(names(up), "_up", sep="")
  down <- DEG_list$Down
  names(down) <- paste(names(down), "_down", sep="")
  DEGlist <- c(up_down, up, down)
  DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
  BatchResult <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="all", id_type, CLSZ, cutoff, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
  write.table(BatchResult, "results/GO/GO_BatchResult.xls", quote=FALSE, sep="\t", row.names = FALSE)
  goslimvec <- as.character(getBM(attributes=c("goslim_goa_accession", "goslim_goa_description"), mart=m)[,1])
  BatchResultslim <- GOCluster_Report(catdb=catdb, setlist=DEGlist, method="slim", id_type, myslimv=goslimvec, CLSZ, cutoff, gocats=c("MF", "BP", "CC"), recordSpecGO=NULL)
  write.table(BatchResultslim, "results/GO/GO_BatchResultslim.xls", quote=FALSE, sep="\t", row.names = FALSE)
  gos <- BatchResultslim
  jpeg("results/GO/GOslimbarplotMF.jpeg")
  goBarplot(gos, gocat="MF")
  dev.off()
  jpeg("results/GO/GOslimbarplotBP.jpeg")
  goBarplot(gos, gocat="BP")
  dev.off()
  jpeg("results/GO/GOslimbarplotCC.jpeg")
  goBarplot(gos, gocat="CC")
  dev.off()
}
