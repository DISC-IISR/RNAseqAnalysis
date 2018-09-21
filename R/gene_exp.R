#'@title Differential Gene Expression plots
#'
#'@description write a desc
#'@import edgeR
#'@import lattice
#'@import gplots
#'@import RColorBrewer
#'@import DESeq2
#'@import pheatmap
#'@import ggplot2
#'@import vsn
#'@import grDevices
#'@import graphics
#'@import SummarizedExperiment
#'@importFrom utils write.table
#'@importFrom utils read.table
#'@importFrom utils read.delim
#'@importFrom stats density
#'@importFrom systemPipeR systemArgs
#'@importFrom systemPipeR targetsin
#'@importFrom systemPipeR filterDEGs
#'@importFrom systemPipeR run_edgeR
#'@importFrom systemPipeR readComp
#'@importFrom easyRNASeq assay
#'@param sysma path to 'param' file; file structure follows a simple name/value syntax that converted into JSON format; for details about the file structure see sample files provided by package. Assign NULL to run the pipeline without 'param' file. This can be useful for running partial workflows, e.g. with pregenerated BAM files.
#'@param mytargets path to targets file
#'@param type type="SYSargs" returns SYSargs, type="json" returns param file content in JSON format (requires rjson library)
#'@param x used for ploting the data to specify weather the data is paired-end or single end
#'@param y used for plotting the data to indicate the count of data
#'@param con a formula which expresses how the counts for each gene depend on the variables in colData
#'@export
#'@seealso \code{\link{edgeR}}
#'@examples
#'
#'\dontrun{
#'pairs <- c( "paired-end", "paired-end", "paired-end","paired-end","paired-end","paired-end")
#'count <- c(1,1,2,2,3,3)
#'condition <- factor(c("N1","N1","N2","N2","M1","M1","M2","M2","M3","M3","M4","M4","U1","U1"))
#'gene_exp(sysma="./param/tophat.param", mytargets="targetsPE.txt", x = pairs, y = count, con = condition)
#'}
#'

gene_exp <- function (sysma, mytargets, type = "SYSargs", x, y, con)
{

  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  dir.create("results/expression")

  args <- systemArgs(sysma, mytargets, type)
  countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
  colData <- data.frame(row.names=targetsin(args)$SampleName, condition=targetsin(args)$Factor)
  condition <- con
  dds <- DESeqDataSetFromMatrix(countDF, DataFrame(condition), ~ condition)
  dds <- dds[ rowSums(counts(dds)) > 5, ]
  cts = counts(dds)
  geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds = estimateSizeFactors(dds, geoMeans=geoMeans)
  
  dds$condition <- droplevels(dds$condition)
  
  dds <- DESeq(dds)
  res <- results(dds)
  write.table(res, "results/expression/result.xls", col.names=NA, quote=FALSE, sep="\t")
  print(res)
  print(summary(res))
  resLFC <- lfcShrink(dds, coef=2, res=res)
  write.table(resLFC, "results/expression/shrunken_result.xls", col.names=NA, quote=FALSE, sep="\t")

  jpeg("results/expression/log2_fold_changes.jpeg")
  plotMA(res, ylim=c(-2,2))
  dev.off()
  jpeg("results/expression/log2_fold_changes_AfterNoiseRemoval.jpeg")
  plotMA(resLFC, ylim=c(-2,2))
  dev.off()



  d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
  ggplot(d, aes(condition,count)) + geom_point(position=position_jitter(width=0.1,height=0)) +scale_y_log10(breaks=c(25,100,400))
  ggsave("results/expression/counts_of_reads_across_groups.jpeg")


  rld <- rlog(dds, blind=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  plotPCA(vsd)
  ggsave("results/expression/PCA_plot.jpeg")


  # this gives log2(n + 1)
  ntd <- normTransform(dds)
  notAllZero <- (rowSums(counts(dds))>0)
  jpeg("results/expression/StandardDeviation_basedon_SLransformation.jpeg")
  meanSdPlot(assay(ntd)[notAllZero,])
  dev.off()

  jpeg("results/expression/StandardDeviation_basedon_RLtransformation.jpeg")
  meanSdPlot(assay(rld[notAllZero,]))
  dev.off()

  jpeg("results/expression/StandardDeviation_basedon_VSTransformation.jpeg")
  meanSdPlot(assay(vsd[notAllZero,]))
  dev.off()


  select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
  # user input must be taken for the below commands
  colData(dds)$type <- x
  df <- as.data.frame(colData(dds)[,c("condition","type")])
  jpeg("results/expression/count_matrix_based_on_shifted_logarithm_transformation.jpeg")
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  dev.off()

  jpeg("results/expression/count_matrix_based_on_regularized_log_transformation.jpeg")
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  dev.off()

  jpeg("results/expression/count_matrix_based_on_variance_stabilizing_transformation.jpeg")
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df)
  dev.off()


  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  jpeg("results/expression/heatmap_of_distance_matrix.jpeg")
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, color=colors)
  dev.off()

  countDF <- read.delim("./results/countDFeByg.xls", row.names=1, check.names=FALSE)
  print(dim(countDF))
  cpm <- cpm(countDF)
  lcpm <- cpm(countDF, log=TRUE)

  samplenames.txt <- substring(colnames(countDF), 1)
  nsamples <- ncol(countDF)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
       main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=c(samplenames.txt), text.col=col, bty="n")
  table(rowSums(countDF==0)==3)
  keep.exprs <- rowSums(cpm>1)>=2
  countDF <- countDF[keep.exprs,,]
  print(dim(countDF))
  lcpm <- cpm(countDF, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
       main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=c(samplenames.txt), text.col=col, bty="n")
  dev.copy(png,"results/expression/read_density_based_on_log-cpm.png")
  dev.off()


  countDF <- read.delim("./results/countDFeByg.xls", row.names=1, check.names=FALSE)
  x <- rowSums(countDF==0)!=ncol(countDF)
  newCountDF <- countDF[x,]
  conds <- y
  y <- DGEList(counts=newCountDF, group=conds)
  keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ]
  y <- calcNormFactors(y)
  x2 <- y
  x2$samples$norm.factors <- 1
  x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
  x2$counts[,2] <- x2$counts[,2]*5
  par(mfrow=c(1,2))
  lcpm <- cpm(y, log=TRUE)
  boxplot(lcpm, las=2, col=col, main="")
  title(main="A. Example: Unnormalised data",ylab="Log-cpm")
  x2 <- calcNormFactors(x2)
  lcpm <- cpm(x2, log=TRUE)
  boxplot(lcpm, las=2, col=col, main="")
  title(main="B. Example: Normalised data",ylab="Log-cpm")
  dev.copy(png,"results/expression/box_plot.png")
  dev.off()

  
}
