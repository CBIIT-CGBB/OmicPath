
ES <- function(db=db, gene=gene, pdffile=pdffile, fdr=T){
  GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
    tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))
    no.tag.indicator <- 1 - tag.indicator 
    N <- length(gene.list) 
    Nh <- length(gene.set) 
    Nm <-  N - Nh 
    if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
    norm.tag    <- 1.0/sum.correl.tag
    norm.no.tag <- 1.0/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > - min.ES) {
      #      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      #      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
  }
  if(fdr){
    my.col <- c(2, 7);
  }
  if (is.list(db)==T) {
    # based on the input GSEA.db
    g.list  <- db$pathway.list;
    g.table <- db$pathway.table;
    g.ver   <- db$pathway.version;
  } else if (nchar(db) >= 1){
    # by known genome
    pofile.s <- paste(db, ".RData", sep="");
    pofile   <- system.file("data", pofile.s, package="OmicPath");
    if (file.exists(pofile)){
    } else {
      stop ("db name?");
    }
    GSEA.db <- local(get(load(pofile)));
    g.list  <- GSEA.db$pathway.list;
    g.table <- GSEA.db$pathway.table;
    g.ver   <- GSEA.db$pathway.version;
  } else {
    stop("db name ??")
  }
  
  gene.db.sym   <- unique(unlist(g.list));
  gene.db.num   <- length(gene.db.sym);
  gene          <- gene[!is.na(gene)];
  gene.num      <- length(gene);
  gene.list     <- rank(gene);
  
  pdf(pdffile, 6, 4);
  par(mar=c(4,5,4,2));
  
  out.s <- c();
  for (i in 1:length(g.list)){
    g.sub     <- unlist(g.list[[i]]);
    gene.set  <- which(gene %in% g.sub);
    gene.n    <- gene[gene.set];
    if (length(gene.set)==0){
      next;
    }
    the.i <- intersect(gene.list, gene.set);
    if (length(the.i) == 0){
      next;
    }
    out       <- GSEA.EnrichmentScore(gene.list=gene.list, gene.set=gene.set, weighted.score.type = 0, correl.vector = NULL);
    RES.m     <- max(abs(out$RES));
    ranking   <- rank(gene.list);
    ind       <-  gene.set;
    geneset   <-  ranking[ind]  
    background <- ranking[-ind]  
    out2  <- ks.test(geneset, background)  
    p.val <- sprintf("%.3e", out2$p.value);
    ## output
    g.n    <- paste(as.character(gene.n), collapse=" ");
    p.sym  <- as.character(g.table[i,1]);
    p.name <- as.character(g.table[i,3]);
    out.s <- rbind(out.s, c(p.val, g.n, p.sym, p.name));
    
    if (out2$p.value < 0.05){
      ## plot
      y.i <- which(out$indicator==1);
      y1  <- y2 <- out$indicator;
      RES.m1 <- RES.m/8;
      
      y1[y.i] <-  RES.m1;
      y2[y.i] <- -RES.m1;
      
      xlab  <- paste0("P value: ", p.val);
      ylab  <- "Running Enrichment Score\n(RES)";
      cols  <- rainbow(10, alpha=0.8);
      col1  <- rgb(0.5, 0.5, 0.5, 0.5);
      
      plot(1:length(gene.list), out$RES, type="l", ylim=c(-RES.m, RES.m), axes=F, 
           col=cols[7], ylab="", xlab="", lwd=3, main=p.name);
      points(1:length(gene.list), y1, type="h", col=col1, lwd=2);
      points(1:length(gene.list), y2, type="h", col=col1, lwd=2);
      mtext(xlab, 1, 0, cex=1.2);
      mtext(ylab, 2, 2, cex=1.2);
      axis(2);
    }
  }
  dev.off();
  
  colnames(out.s) <- c("Pvalue",	"Gene", "db.type",	"db.name")
  if (fdr){
    p.adjust.M <- p.adjust.methods[c(4,7)];
    p.adj      <- sapply(p.adjust.M, function(meth) p.adjust(out.s[,1], meth));
    out.p      <- cbind(out.s, p.adj);
    m.res.out  <- out.p[order(as.numeric(out.p[,1])),];
  } else {
    m.res.out <- out.s
  }
  return(m.res.out);
}





