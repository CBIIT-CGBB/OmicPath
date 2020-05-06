
steiner_net <- function(g=g, steiner.node=steiner.node, plot.net=F){

  mst.g    <- minimum.spanning.tree(g);
  steiner.points <- which(V(g)$name %in% steiner.node);

  Gi             <- graph.full(length(steiner.points));
  V(Gi)$name     <- V(g)$name[steiner.points];
  mst <- minimum.spanning.tree(Gi);
  
  ##  For each edge in mst, replace with shortest path:
  edge_list <- get.edgelist(mst)

  Gs <- mst
  message ("Doing Steiner edges ...");

  pb <- txtProgressBar(min = 0, max = nrow(edge_list), style = 3);
  for (n in 1:nrow(edge_list)) {
    i <- edge_list[n,2]
    j <- edge_list[n,1]
    setTxtProgressBar(pb, n)
    ##  If the edge of T' mst is shared by Gi, then remove the edge from T'
    ##    and replace with the shortest path between the nodes of g: 
    if (length(E(Gi)[which(V(mst)$name==i) %--% which(V(mst)$name==j)]) == 1) {
      ##  If edge is present then remove existing edge from the 
      ##    minimum spanning tree:
      Gs <- Gs - E(Gs)[which(V(mst)$name==i) %--% which(V(mst)$name==j)]
      
      ##  Next extract the sub-graph from g corresponding to the 
      ##    shortest path and union it with the mst graph:
      tmp <- shortest.paths(g, v=V(g)[i], to=V(g)[j]);
      if (is.infinite(tmp)){
        next;
      }
      ## check get.shortest.paths
      net.s <- get.shortest.paths(g, from=V(g)[i], to=V(g)[j]);
      if (length(net.s$vpath[[1]])==0){
        next;
      }
      g_sub <- induced.subgraph(g, (net.s$vpath[[1]]))
      Gs    <- graph.union(Gs, g_sub, byname=T)
    }
  }
  close(pb);
  
  ## 
  df    <- as.data.frame(get.edgelist(g));
  df.G  <- as.data.frame(get.edgelist(Gs));
  out.g <- Gs;

  if (plot.net){
    l    <- layout.fruchterman.reingold(g);
    l2   <- layout.fruchterman.reingold(Gs);
    l2   <- layout.norm(l2, -1, 1, -1, 1);
    
    cols <- rainbow(10, alpha=0.9);
    col2 <- rainbow(10, alpha=0.2);

    l    <- layout.fruchterman.reingold(mst.g);
    V(mst.g)$label.cex <- 0.1;
    plot(mst.g, vertex.size=0.1, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[7],
         vertex.label=V(g)$name, edge.width=0.8, edge.arrow.size=0.1)
    

    message ("Adding more edges on Steiner network ...");
    
    plot(l2, type="n", axes=F, xlab="", ylab="", main="");
    pb <- txtProgressBar(min = 0, max = length(V(Gs)), style = 3);
    for (k1 in 1:length(V(Gs))){
      n1 <- V(Gs)$name[k1];
      for (k2 in k1:length(V(Gs))){
        if (k1==k2){
          next;
        }
        n2 <- V(Gs)$name[k2];
        if (are.connected(Gs, n1, n2)){
          next;
        }
        if (are.connected(g, n1, n2)){
          i1 <- which(V(Gs)$name==n1);
          i2 <- which(V(Gs)$name==n2);
          out.g <- add_edges(out.g, c(n1, n2));
          segments(l2[i1,1], l2[i1,2], l2[i2,1], l2[i2,2], col=col2[7], lwd=4);  
        }
      }
      setTxtProgressBar(pb, k1)
    }
    close(pb);
    
    V(Gs)$label.cex <- 0.6;
    plot(Gs, vertex.size=1, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[1],
         vertex.label=V(Gs)$name, add=T, layout=l2, edge.width=1.2);
    
    V(out.g)$label.cex <- 0.6;
    l3 <- layout.fruchterman.reingold(out.g);
    plot(out.g, vertex.size=2, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[7],
         vertex.label=V(out.g)$name, layout=l3, edge.width=1.2);
    
    plot(out.g, vertex.size=2, vertex.color=cols[2], vertex.frame.color=cols[2], edge.color=cols[7],
         vertex.label=V(out.g)$name, edge.width=1.2);
    
    out.df  <- as.data.frame(get.edgelist(out.g));
    return(list(df=out.df,g=out.g));

  } else {
    message ("Adding more edges on Steiner network ...");
    pb <- txtProgressBar(min = 0, max = length(V(Gs)), style = 3);
    for (k1 in 1:length(V(Gs))){
      n1 <- V(Gs)$name[k1];
      for (k2 in k1:length(V(Gs))){
        if (k1==k2){
          next;
        }
        n2 <- V(Gs)$name[k2];

        if (are.connected(Gs, n1, n2)){
          next;
        }
        if (are.connected(g, n1, n2)){
          i1 <- which(V(Gs)$name==n1);
          i2 <- which(V(Gs)$name==n2);
          out.g <- add_edges(out.g, c(n1, n2));
        }
      }
      setTxtProgressBar(pb, k1)
    }
    close(pb);
    
    out.df  <- as.data.frame(get.edgelist(out.g));
    return(list(df=out.df,g=out.g));
  }

}