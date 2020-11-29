
WGCNA_6_4 <- function(expr=expr, pheno=pheno, mdir="module", cutHeight=200, softPower=8, minM=25, 
                      cutH=0.996, deepSplit=1,  project.name="WGCNA_net"){
  
  ###
  softPower     <- softPower; # from 02_TOM
  # We like large modules, so we set the minimum module size relatively high:
  # from 03_TOM cH=0.996 mM=15 dS=1
  minModuleSize <- minM;
  cutHeight     <- cutH;
  deepSplit     <- deepSplit;
  datExpr       <- expr;
  datpheno      <- pheno;
  subDir        <- mdir;
  if (dir.exists(subDir)){
  } else {
    dir.create(subDir);
  }
  # the pdf figure size should to be done by pdf("pdffile", sizeX, sizeY)
  pdffile  = paste(project.name, "_04_TOM_Association.pdf", sep="");

  
  # outfile
  nc.file   <- paste(project.name , "-networkConstruction-stepByStep.RData", sep="")
  saveB     <- paste(project.name , "-TOM-blockwise", sep="")
  mtp.file  <- paste(project.name, "_moduleTraitPvalue.xls", sep="")
  mtc.file  <- paste(project.name, "_moduleTraitCor.xls", sep="")
  mes.file  <- paste(project.name, "_moduleEigengenes.xls", sep="")
  pm.file   <- paste(project.name, "_ProbeModule.xls", sep="")
  MEDiss.file <- paste(project.name, "_ModuleCor.xls", sep="")
  
  dissTOM.f <- paste(project.name, "_dissTOM.xls", sep="");
  print("reading dissTOM ...")
  dissTOM   <- read.table(dissTOM.f, header=T, row.names=1);
  print("read dissTOM is done.")
  # Clustering using TOM
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(dissTOM), method = "average");
  
  # Module identification using dynamic tree cut:
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit   = deepSplit, 
                               pamRespectsDendro = FALSE, 
                               cutHeight         = cutHeight,
                               minClusterSize    = minModuleSize);
  print(table(dynamicMods));
  
  # Convert numeric lables into colors
  dynamicColors <- WGCNA::labels2colors(dynamicMods)
  print(table(dynamicColors));
  
  # Plot the dendrogram and colors underneath
  #sizeGrWindow(8,6)
  WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide     = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  # Calculate eigengenes ### re-run R if Error Error in colMeans(x, na.rm = TRUE) : 'x' must be numeric
  MEList <- WGCNA::moduleEigengenes(datExpr, colors = dynamicColors)
  MEs    <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1-cor(MEs);
  write.table(MEDiss, MEDiss.file, quote=F, sep="\t");
  # Cluster module eigengenes
  METree <- hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "", cex=1.5)
  
  MEDissThres = 0.02
  # Plot the cut line into the dendrogram
  # abline(h=MEDissThres, col = "red")
  
  # Call an automatic merging function
  merge <- WGCNA::mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors <- merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs    <- merge$newMEs;
  
  #
  WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05);
  
  # Rename to moduleColors
  moduleColors <- mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder   <- c("grey", standardColors(50));
  moduleLabels <- match(moduleColors, colorOrder)-1;
  MEs <- mergedMEs;
  # Save module colors and labels for use in subsequent parts
  # save(MEs, moduleLabels, moduleColors, geneTree, file = nc.file)
  
  # 
  bwnet <- WGCNA::blockwiseModules(datExpr, 
                            maxBlockSize = 2000,
                            power = softPower, 
                            minModuleSize = minModuleSize,
                            reassignThreshold = 0, 
                            mergeCutHeight = 0.1,
                            numericLabels = TRUE,
                            saveTOMs = FALSE,
                            saveTOMFileBase = saveB,
                            verbose = 3)
  
  # Relabel blockwise modules
  bwLabels <- WGCNA::matchLabels(bwnet$colors, moduleLabels);
  # Convert labels to colors for plotting
  bwModuleColors <- WGCNA::labels2colors(bwLabels)
  
  singleBlockMEs <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes;
  blockwiseMEs   <- WGCNA::moduleEigengenes(datExpr, bwModuleColors)$eigengenes;
  
  single2blockwise <- match(names(singleBlockMEs), names(blockwiseMEs))
  #signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)
  
  # Quantifying module{trait associations
  # Define numbers of genes and samples
  nGenes   <- ncol(datExpr);
  nSamples <- nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs  <- WGCNA::orderMEs(MEs0)
  moduleTraitCor    <- cor(MEs, datpheno, use = "p");
  moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples);
  write.table(moduleTraitPvalue, mtp.file, sep="\t", quote=F)
  write.table(moduleTraitCor, mtc.file, sep="\t", quote=F)
  write.table(MEs, mes.file, sep="\t", quote=F)
  
  #sizeGrWindow(16,10)
  ##############################################
  # Will display correlations and their p-values
  ##############################################
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  #
  par(mar = c(6, 10, 3, 3));
  # Display the correlation values within a heatmap plot
  pdf(pdffile, 8, 8)
  par(mar=c(8,10,4,2))
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(datpheno),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off();

  high     <- as.data.frame(datpheno);
  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  names(MMPvalue) <- paste("p.MM", modNames, sep="");
  geneTraitSignificance <- as.data.frame(cor(datExpr, high, use = "p"));
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  # output
  names(geneTraitSignificance) <- paste("GS.", names(high), sep="");
  names(GSPvalue) <- paste("p.GS.", names(high), sep="");
  probeColor <- cbind(colnames(datExpr), moduleColors);
  probeColor.colname <- c("Probe", "Module");
  write.table(probeColor, pm.file, col.names=probeColor.colname, row.names=F,
              sep="\t", quote=F)
  
  #names(datExpr)
  # Gene significance for Mets
  for (i in 1:length(modNames)){
    outfile = paste(subDir, "/", "module_", modNames[i], ".txt", sep="")
    write.table (names(datExpr)[moduleColors==modNames[i]], outfile, quote=F, row.names=F)
    #print (names(datExpr)[moduleColors=="blue"]);
  }
}
