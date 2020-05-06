
WGCNA_6_5 <- function(mdir="module", project.name="WGCNA_net", cutoff=20000){

  files  <- dir(mdir, pattern = "^module", full.names = TRUE);
  file2  <- dir(mdir, pattern = "^module");
  file.n <- length(files);

  ## 
  print("reading dissTOM ...")
  dissf     <- paste0(project.name, "_dissTOM.xls");
  dissTOM   <- read.table(dissf, header=T, row.names=1);
  col.names <- colnames(dissTOM); 
  print("read dissTOM is done.")
  for (i in 1:file.n){
    infile     <- files[i];
    m.name     <- sub("module_", "", file2[i]);
    m.name     <- sub(".txt", "", m.name);
    m.col      <- read.table(infile, header=T);
    node.name  <- m.col$x;
    #node.name  <- paste("X", node.name, sep="");
    node.num   <- length(node.name);
    out.s      <- paste("file = ", infile, " node number = ", node.num, sep=""); 
    print(out.s);
    if (node.num > cutoff){
      next;
    }
    m.n.index  <- col.names %in% node.name;
    diss.v     <- dissTOM[m.n.index, m.n.index];
    diss.v     <- as.matrix(diss.v);
    g0         <- igraph::graph_from_adjacency_matrix(diss.v, weighted=TRUE, mode ="upper");
    m.mst      <- igraph::mst(g0, weights=E(g0)$weight);
    g          <- m.mst;
    outfile  <- paste(mdir, "/MST_", m.name, "_edges.txt", sep="");
    outfile2 <- paste(mdir, "/MST_", m.name, "_degrees.txt", sep="");
    outfile3 <- paste(mdir, "/ALL_", m.name, "_edges.txt", sep="");
    out2     <- data.frame(gene=V(g)$name, degree=degree(g));
    out2.i   <- order(out2[,2], decreasing = T);
    out2     <- out2[out2.i,];
    write.table(out2, outfile2, quote=F, sep="\t", row.names=F);
    out.s   <- igraph::get.edgelist(g, names=TRUE);
    colnames(out.s) <- c("node1", "node2")
    write.table(out.s, outfile, quote=F, row.names=F, sep="\t");
    out.s0   <- igraph::get.data.frame(g0);
    colnames(out.s0) <- c("node1", "node2", "weight")
    write.table(out.s0, outfile3, quote=F, row.names=F, sep="\t");
    # dino.mst of fossil
    # mst_restart of nnclust
    # spantree {vegan}
    # m.mst    <- minimum.spanning.tree(graph.adjacency(as.matrix(diss.v))) #igraph
    ###########################################################
    # Visualization
    ###########################################################
    V(g)$label <- V(g)$name;
    #l          <- layout.kamada.kawai(g);  
    l          <- igraph::layout.fruchterman.reingold(g);   
    pdf.f      <- paste(mdir, "/MST_", m.name, ".pdf", sep="");
    
    edge.color   <- c();
    vertex.color <- c();
    edge.color[1:length(E(g))]   <- rgb(0,0,0,0.4);
    vertex.color[1:length(V(g))] <- rgb(0,0,1,0.8);
    m.name2     <- gsub("network/MST_", "", m.name);
    m.name2     <- gsub("_edges_gene", "", m.name2);
    main        <- paste("Model = ", m.name2, sep="");
    sub         <- paste("infor = ", m.name2, sep="");
    vertex.label.cex   <- .3*log((degree(g)+.1)) + 0.0;
    vertex.label.cex[vertex.label.cex < 0.4] <- 0.4;
    vertex.label.cex[vertex.label.cex > 0.8] <- 0.8;
    vertex.size        = .8*degree(g);
    vertex.size[vertex.size < 0.2] <- 0.2;
    vertex.size[vertex.size > 0.8] <- 0.8;
    main        <- paste("Model = ", m.name, sep="");
    sub         <- paste("infor = ", m.name, sep="");
    pdf(pdf.f, 6,6);
    plot(g,
         layout=l,
         vertex.size        = vertex.size, 
         edge.color         = edge.color, 
         #vertex.shape       = "circle",
         vertex.shape       = "none",
         edge.width         = 0.2,
         #asp                = FALSE,
         #vertex.label.dist = 0.1, 
         vertex.label.cex   = vertex.label.cex, 
         vertex.label.color = vertex.color,
         main               = main,
         sub                = sub
    )
    dev.off();
  }
}