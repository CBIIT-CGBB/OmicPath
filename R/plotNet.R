
plotNet <- function (db=db, dat=dat, Taxonomy="mm", cluster=TRUE, pdffile=pdffile){
  if (Taxonomy=="mm"){
    if (db == "biosystems"){
      infile.s <- paste0("biosystems_mm.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else if (db == "PPI"){
      infile.s <- paste0("PPI_mm.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else if (db == "KEGG") {
      infile.s <- paste0("relation_mmu.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else {
      infile.s <- paste0("all_interaction_mm.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    }
  } else if (Taxonomy=="hs"){
    if (db == "biosystems"){
      infile.s <- paste0("biosystems_hs.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else if (db == "PPI"){
      infile.s <- paste0("PPI_hs.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else if (db == "KEGG") {
      infile.s <- paste0("relation_hsa.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    } else {
      infile.s <- paste0("all_interaction_hs.txt")
      infile   <- system.file("extdata",  infile.s, package="OmicPath");
    }
  } else {
    stop("Taxonomy", Taxonomy);
  }
  
  gene  <- as.character(dat[,1]);
  dat.l <- read.table(infile);
  r1.i  <- which(dat.l[,1] %in% gene);
  r2.i  <- which(dat.l[,2] %in% gene);
  r.i   <- intersect(r1.i, r2.i);

  dat.s <- dat.l[r.i,];
  dat.s <- na.omit(dat.s);
  g   <- igraph::graph.data.frame(dat.s, directed=T);
  V(g)$label <- V(g)$name;
  
  col1       <- rainbow(10, alpha=1);
  col2       <- rainbow(10, alpha=0.4);
  l          <- igraph::layout.kamada.kawai(g);
  
  ## 
  if (cluster){
    wc   <- igraph::walktrap.community(g);
    col3 <- rainbow(length(wc), alpha=0.2);
  } 
  
  dat.i   <- which(dat[,1] %in% V(g)$label);
  if (length(dat.i)>0){
    vertex.size <- -log10(dat[dat.i,3])/5;
  } else {
    vertex.size <- rep(1, length(V(g)$label));
  }

  vertex.size[vertex.size > 4]  <- 4;
  vertex.color       <- rep(col2[7], length(V(g)$label));
  vertex.color.i     <- which(dat[V(g)$label,2] > 0);
  vertex.color[vertex.color.i] <- col2[1];
  vertex.frame.color <- vertex.color;
  
  vertex.label.cex   <- log2(c(degree(g) + 1));
  vertex.label.cex[vertex.label.cex < 0.4] <- 0.4;
  vertex.label.cex[vertex.label.cex > 0.6] <- 0.6;
  vertex.label.color <- rep(rgb(0,0,0,1), length(V(g)));
  
  pdf(pdffile, 6, 6);
  if (cluster){
    plot(g,
         mark.groups=wc, mark.col=col3, mark.border=col3, mark.shape=1,
         layout=l,
         ##
         edge.color         = col2[7], 
         edge.width         = 1,
         edge.arrow.size    = 0.2,
         ##
         vertex.size        = vertex.size, 
         vertex.shape       = "circle",
         vertex.color       = vertex.color,
         vertex.frame.color = vertex.frame.color,
         ##
         vertex.label.dist = 0.1, 
         vertex.label.cex   = vertex.label.cex, 
         vertex.label.color = vertex.label.color,
         ## 
         main               = "",
         xlab               = ""
    )
  } else {
    plot(g,
         layout=l,
         ##
         edge.color         = col2[7], 
         edge.width         = 1,
         edge.arrow.size    = 0.2,
         ##
         vertex.size        = vertex.size, 
         vertex.shape       = "circle",
         vertex.color       = vertex.color,
         vertex.frame.color = vertex.frame.color,
         ##
         vertex.label.dist = 0.1, 
         vertex.label.cex   = vertex.label.cex, 
         vertex.label.color = vertex.label.color,
         ## 
         main               = "",
         xlab               = ""
    )
  }

  dev.off();
}
