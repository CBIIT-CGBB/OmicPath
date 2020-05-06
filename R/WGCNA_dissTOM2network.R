
WGCNA_dissTOM2network <- function(m.files=m.files, indir=indir, outdir=outdir, project.name="ALL", softPower=6, 
                              cutoff=cutoff, label=FALSE, label.cex=0.4){
  if (cutoff=="power"){
    infile    <- paste0(project.name, "_02_power.xls");
    inf       <- read.table(infile, header=T);
    inf.i     <- which(inf[,1]==softPower);
    # get Mean connectivity
    mean.k    <- as.numeric(inf[inf.i, 5]);
    dissf     <- paste0(project.name, "_dissTOM.xls");
    ## 
    print("reading dissTOM ...")
    dissTOM   <- read.table(dissf, header=T, row.names=1);
    dissTOM   <- as.matrix(dissTOM);
    print("read dissTOM is done.")
    ## 
    dissTOM.tri   <- dissTOM[lower.tri(dissTOM)];
    dissTOM.v     <- as.vector(dissTOM.tri);
    dissTOM.w     <- sort(dissTOM.v);
    gene.n    <- dim(dissTOM)[1];
    exp.k     <- gene.n * mean.k;
    cutoff    <- dissTOM.w[exp.k];
  }
  if (length(m.files) < 1){
    warning("you did not give me a file!")
  }
  for (f in 1:length(m.files)){
    root    <- gsub(".matrix", "", m.files[f]);
    infile  <- paste0(indir, "/", m.files[f]);
    dat     <- read.table(infile, header=T, row.names=1);
    colnames(dat) <- row.names(dat);
    if (cutoff=="auto"){
      tmp    <- summary(as.vector(dat));
      cutoff <- tmp[2];
    } 

    dat[dat > cutoff] <- 0;
    dat     <- as.matrix(dat);
    g       <- igraph::graph_from_adjacency_matrix(dat, weighted=TRUE, mode="upper");
    out.s   <- igraph::get.edgelist(g, names=TRUE);
    pdff    <- paste0(outdir, "/", root, ".pdf");
    outf    <- paste0(outdir, "/", root, ".xls");
    dff     <- paste0(outdir, "/", root, "_df.txt");
    out.d   <- data.frame(gene=V(g)$name, df=degree(g, mode ="all", loops=F));
    out.i   <- order(out.d[,2], decreasing = T);
    out.d   <- out.d[out.i,];
    
    #out.d   <- out.d[!duplicated(out.d[,1]),];
    pdf(pdff, 8,8);
    cols <- rainbow(10, alpha=0.9);
    col2 <- rainbow(10, alpha=0.2);
    
    write.table(out.s, outf, quote=F, sep="\t", row.names=F);
    write.table(out.d, dff,  quote=F, sep="\t", row.names=F);
    
    g2   <- igraph::minimum.spanning.tree(g, weights=E(g)$weight);
    l    <- igraph::layout.fruchterman.reingold(g2);
    l    <- igraph::layout.norm(l, -1, 1, -1, 1);
    df   <- as.data.frame(get.edgelist(g));
    plot(l, type="n", axes=F, xlab="", ylab="", main="");
    if (is.null(E(g)$weight)){
      maxw <- 1;
    } else {
      maxw <- max(E(g)$weight);
    }
    V(g2)$label <- V(g2)$name;
    w           <- 4*E(g)$weight/maxw;
    i1 <- which(V(g)$name %in% df[,1]);  
    i2 <- which(V(g)$name %in% df[,2]);
    segments(l[i1,1], l[i1,2], l[i2,1], l[i2,2], col=col2[7], lwd=w);

    if (label==TRUE){
      plot.igraph(g2, vertex.size=6, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[1],
                  vertex.label.cex = label.cex, layout=l, add=T, edge.width=4*E(g2)$weight/maxw);
    } else {
      plot.igraph(g2, vertex.size=6, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[1],
                  vertex.label = NULL, layout=l, add=T, edge.width=4*E(g2)$weight/maxw);
    }
    dev.off();

  }
}