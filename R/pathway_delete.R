
pathway_delete <- function(db=db, pathway.ID=pathway.ID,
                           creator=creator, 
                           version=version, mydate=NULL){
  out.s <- db$pathway.list;
  out.s[[pathway.ID]] <- NULL;
  out.b   <- db$pathway.table;
  out.i   <- which(out.b[,1]==pathway.ID);
  if (length(out.i)==0){
    stop("The pathway.ID ", pathway.ID, " was not found.");
  }
  out.b   <- out.b[-out.i,];
  creator <- creator;
  if (is.null(mydate)){
    mydate <- date();
  } else {
    mydate <- mydate;
  }
  gene.u.num   <- length(unique(unlist(out.s)))
  GSEA.db <- list(pathway.list=out.s, pathway.table=out.b, 
                     pathway.version=list("version:" = version, "Date:" = mydate, 
                                          "Gene number:" = gene.u.num, "Creator:" = creator));
  return(GSEA.db);
}

