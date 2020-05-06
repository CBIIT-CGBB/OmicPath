
pathway_update <- function(db=db, gene.list=gene.list, pathway.ID=pathway.ID,
                           pathway.name=pathway.name, creator=creator, 
                           version=version, mydate=NULL){
  
  out.s               <- db$pathway.list;
  out.s[[pathway.ID]] <- gene.list;
  out.b   <- db$pathway.table;

  out.i   <- which(out.b[,1]==pathway.ID);

  if (length(out.i)==0){
    stop("The pathway.ID ", pathway.ID, " was not found.");
  }
  out.b <- out.b[-out.i,];
  out.b   <- rbind(out.b, 
                   data.frame(pathway.ID = toString(pathway.ID), 
                              pathway.name = toString(pathway.name)));

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
  return(GSEA.db)
}

