

pathway_creator <- function(gene.table=gene.table, creator=creator, 
                           version=version, mydate=NULL){
  dat.s  <- gene.table;
  out.s   <- list();
  for (i in 1:nrow(dat.s)){
    n     <- as.character(dat.s[i,1]);
    g.l   <- strsplit(dat.s[i,2], "\\s+");
    out.s <- c(out.s, list(i=g.l));
  }
  names(out.s) <- as.character(dat.s[,1]);
  out.b   <- dat.s;
  gene.u.num <- length(unique(unlist(out.s))); 
  if (is.null(mydate)){
    GSEA.db <- list(pathway.list=out.s, pathway.table=out.b, 
                pathway.version=list("version:" = version, "Date:" = date(), 
                                     "Gene number:" = gene.u.num, "Creator:" = "Ying Hu"));
  } else {
    GSEA.db <- list(pathway.list=out.s, pathway.table=out.b, 
                pathway.version=list("version:" = version, "Date:" = mydate, 
                                     "Gene number:" = gene.u.num, "Creator:" = "Ying Hu"));
  }
  return(GSEA.db);
}
