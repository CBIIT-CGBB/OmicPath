
WGCNA_6_2 <- function(expr=expr, project.name="WGCNA_net"){
  pdffile   <- paste(project.name , "_02_power.pdf", sep="");
  outfile   <- paste(project.name , "_02_power.xls", sep="");
  datExpr   <- expr;
  pdf(pdffile, 8,6)
  # network construction and module detection
  # Choosing the soft-thresholding power: analysis of network topology
  # Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to=40, by=2))
  
  # Call the network topology analysis function
  sft   <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # write out
  write.table(sft$fitIndices, outfile, sep="\t", row.names=F, quote=F);
  # Plot the results:
  par(mfrow = c(1,2));
  cex1 = 0.90;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  r2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2];
  r2.max <- max(r2);
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h=r2.max, col="blue", lwd=2);
  r2.i <- which(r2==r2.max);
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(v=r2.i, col="blue", lwd=2);
  dev.off()
}