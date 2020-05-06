
hclustplotn=function(hier1, Color, RowLabels=NULL, cex.RowLabels = 0.9, ...) 
{
  options(stringsAsFactors=FALSE);
  if (length(hier1$order) != dim(Color)[[1]] ) 
  { 
    stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
  } else {
    No.Sets = dim(Color)[[2]];
    C = Color[hier1$order, ]; 
    step = 1/dim(Color)[[1]];
    ystep = 1/No.Sets;
    barplot(height=1, col = "white", border=F,space=0, axes=F, ...)
    for (j in 1:No.Sets)
    {
      ind = (1:(dim(C)[1]));
      xl = (ind-1) * step; xr = ind * step; 
      yb = rep(ystep*(j-1), dim(C)[1]); yt = rep(ystep*j, dim(C)[1]);
      rect(xl, yb, xr, yt, col = as.character(C[,j]), border = as.character(C[,j]));
      if (is.null(RowLabels))
      {
        text(as.character(j), pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      } else {
        text(RowLabels[j], pos=2, x=0, y=ystep*(j-0.5), cex=cex.RowLabels, xpd = TRUE);
      }
    }
    for (j in 1:No.Sets) lines(x=c(0,1), y=c(ystep*j,ystep*j));
  }
}

WGCNA_6_3 <- function(expr=expr, softPower=8, minM=minM, cutH=cutH, project.name="WGCNA_net"){
  softPower <- softPower; # from 02_TOM
  
  # the pdf figure size should to be done by pdf("pdffile", sizeX, sizeY)
  pdffile   <- paste(project.name, "_03_TOM_TreeCut.pdf", sep="");
  pdf(pdffile, 12, 8);
  dissTOM.f <- paste(project.name, "_dissTOM.xls", sep="");

  # Co-expression similarity and adjacency
  # calculate the adjacencies, using the soft thresholding power 8:
  datExpr    <- expr;
  adjacency  <- WGCNA::adjacency(datExpr, power = softPower);
  
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  TOM <- WGCNA::TOMsimilarity(adjacency);
  dissTOM <- 1-TOM
  # save dissTOM 
  outtab  <- dissTOM
  colnames(outtab)  = colnames(datExpr)
  row.names(outtab) = colnames(datExpr)
  write.table(outtab, dissTOM.f, sep = "\t", quote=F)
  rm(outtab)
  
  # Clustering using TOM
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(dissTOM), method = "average");
  
  # We like large modules, so we set the minimum module size relatively high:
  # Module identification using dynamic tree cut:
  deeS <- c(0, 1, 2, 3);
  
  DetectColors <- c();
  Row.L        <- c();
  for (i in 1:length(cutH)){
    cutHeight          <- cutH[i];
    for (j in 1:length(minM)){
      minModuleSize  <- minM[j];
      for (k in 1:length(deeS)){
        deepSplit <- deeS[k];
        dynamicMods <- cutreeDynamic(
          dendro = geneTree, 
          distM = dissTOM,
          deepSplit = deepSplit, 
          pamRespectsDendro = FALSE, 
          cutHeight=cutHeight,
          minClusterSize = minModuleSize
        );
        #table(dynamicMods)
        DetectColors <- cbind(DetectColors, labels2colors(dynamicMods));
        Row.L <- cbind(Row.L, paste("cH=", cutHeight, " mM=", minModuleSize, " dS=", deepSplit, sep=""))
      }
    }
  }
  
  layout(matrix(c(1,2)), heights=c(1,3))
  par(mar=c(0.2,6.0,2,1.2))
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = F, hang = 0.04);
  
  par(mar=c(0.8,6.0,0,1.2))
  hclustplotn(geneTree, DetectColors, RowLabels=Row.L, main="", cex.RowLabels=0.5)
  
  minM = c(25);
  rm(DetectColors)
  rm(Row.L)
  DetectColors = c();
  Row.L        = c();
  for (i in 1:length(cutH)){
    cutHeight          = cutH[i];
    for (j in 1:length(minM)){
      minModuleSize      = minM[j];
      for (k in 1:length(deeS)){
        deepSplit = deeS[k];
        dynamicMods = cutreeDynamic(
          dendro = geneTree, 
          distM = dissTOM,
          deepSplit = deepSplit, 
          pamRespectsDendro = FALSE, 
          cutHeight=cutHeight,
          minClusterSize = minModuleSize
        );
        #table(dynamicMods)
        DetectColors = cbind(DetectColors, labels2colors(dynamicMods));
        Row.L = cbind(Row.L, paste("cH=", cutHeight, " mM=", minModuleSize, " dS=", deepSplit, sep=""))
      }
    }
  }
  
  layout(matrix(c(1,2)), heights=c(1,3))
  par(mar=c(0.2,6.0,2,1.2))
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = F, hang = 0.04);
  
  par(mar=c(0.8,6.0,0,1.2))
  hclustplotn(geneTree, DetectColors, RowLabels=Row.L, main="", cex.RowLabels=0.5)

  dev.off()

}