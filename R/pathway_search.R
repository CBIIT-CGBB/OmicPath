
pathway_search <- function(db=db, gene=gene, pathway.ID=NULL){
  if (is.list(db)==T) {
    # based on the input GSEA.db
    out.s  <- db$pathway.list;
    out.b  <- db$pathway.table;
  } else if (nchar(db) >= 1){
    # by known genome
    pofile.s <- paste(db, ".RData", sep="");
    pofile   <- system.file("data", pofile.s, package="OmicPath");
    if (file.exists(pofile)){
    } else {
      stop ("db name?");
    }
    GSEA.db <- local(get(load(pofile)));
    out.s   <- GSEA.db$pathway.list;
    out.b   <- GSEA.db$pathway.table;
  } else {
    stop("db name ??")
  }

  if (!is.null(gene)){
    gene <- as.character(gene);
    p.i  <- rep(seq_along(out.s), sapply(out.s, length))
    p.j  <- p.i[match(gene, unlist(out.s))];
    out  <- unique(names(out.s)[p.j]);
  }
  if (!is.null(pathway.ID)){
    pathway.ID <- as.character(pathway.ID);
    p.i  <- which(names(out.s) %in% pathway.ID);
    out1 <- c();
    for (j in p.i){
      out1 <- c(out1, out.s[[j]]);
    }
    out1  <- unique(out1);
    j     <- which(out.b[,1] %in% pathway.ID);
    out2  <- out.b[j,];
    out   <- list(gene=out1, info=out2);
  }
  return(out);
}

