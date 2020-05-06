
WGCNA_6_1 <- function(expr=expr, pheno=pheno, cutoff=8000, cutHeight=200, project.name="WGCNA_net"){

  pdffile   <- paste0(project.name, "_01_input_inf.pdf");
  pdf(pdffile, 8,6);

  # Loading gene expression data
  Data       <- expr;
  gene.num   <- cutoff;
  dat.E      <- as.matrix(Data);
  var1       <- function(x) var(x, na.rm=T);
  vargenes   <- apply(dat.E, 1,var1)
  rankvar      <- rank(-vargenes)
  restVariance <- rankvar <= gene.num
  Data         <- dat.E[restVariance,]

  datExpr0   <- as.data.frame(t(Data))

  gsg        <- WGCNA::goodSamplesGenes(datExpr0, verbose = 3);

  sampleTree <- hclust(dist(datExpr0), method = "average");

  par(cex = 0.6);
  par(mar = c(0,12,6,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  # Plot a line to show the cut, remove outlier sample
  cutHeight <- cutHeight;
  abline(h = cutHeight, col = "red");
  # Determine cluster under the line
  clust     <- WGCNA::cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 6)
  table(clust)
  
  # check if sample(s) is(are) going to be dropped
  # clust 1 contains the samples we want to keep.
  # we are going to keep all samples?
  keepSamples <- (clust==1)
  datExpr     <- datExpr0[keepSamples, ]

  # remove the sample from pheno
  Samples   <- rownames(datExpr);
  phenoRows <- match(Samples, rownames(pheno));
  datpheno  <- pheno[phenoRows,];
  
  # Re-cluster samples
  sampleTree2 <- hclust(dist(datExpr), method = "average")
  traitColors <- numbers2colors(datpheno, signed=F);
  WGCNA::plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = colnames(pheno), marAll = c(1, 12, 3, 1), 
                      main        = "Sample dendrogram and pheno heatmap", cex.dendroLabels = 0.6);
  dev.off();
  return(list(expr=datExpr, pheno=pheno));
}