rm(list=ls());

library(igraph);

## set node number
node.n <- 60;
## generate a scale free network 
g     <- barabasi.game(node.n);

## get cluster algorithms
clusters <- grep("^cluster\\_", ls("package:igraph"), value=TRUE)
## get layout
l        <- layout.davidson.harel(g)

pdf("04do_igraph_cluster.pdf", 8, 8);
par(mfrow=c(3, 3), mar=c(2,2,2,2));
for (cl in 1:length(clusters)){
  ## cluster
  if (cl==3){
    wc   <- do.call(clusters[cl], list(as.undirected(g), nb.trials = 20));
  } else if (cl==4){
    wc   <- do.call(clusters[cl], list(as.undirected(g), initial = sample(E(g)$name)));
  } else if (cl==5){
    wc   <- do.call(clusters[cl], list(as.undirected(g), start=membership(wc)));
  } else if (cl==8){
    wc   <- do.call(clusters[cl], list(as.undirected(g), spins = sample(10:30, 1)));
  } else if (cl==9){
    wc   <- do.call(clusters[cl], list(as.undirected(g), steps = sample(2:8, 1)));
  } else {
    wc   <- do.call(clusters[cl], list(as.undirected(g)));
  }
  main <- gsub("cluster_", "", clusters[cl]);
  ## set colors
  col1  <- rainbow(length(wc), alpha=0.2);
  col2  <- sample(rainbow(length(wc), alpha=0.8));
  plot(g, layout=l,  main=main, edge.arrow.size=0.1, 
       vertex.label=NA,
       mark.groups=wc, mark.col=col1, mark.border=col1, mark.shape=1,
       vertex.size=degree(g)+6, vertex.color=col2[membership(wc)], vertex.frame.color=col2[membership(wc)], 
       edge.color="blue", edge.lty=1);
}
dev.off();
