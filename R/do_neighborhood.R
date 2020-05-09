
do_neighborhood <- function(g=g, genes=genes, order=1){
  g.i <- which(V(g)$name %in% genes);
  out <- make_ego_graph(graph=g, nodes = g.i, order=order);
  out.s <- NULL;
  for (i in 1:length(out)){
    if (i==1){
      out.s <- out[[i]]
    } else {
      out.s <- out.s %u% out[[i]];
    }
  }
  return(out.s);
}