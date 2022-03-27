
doGSEA_all <- function(db=db, gene=gene, filter.num=0, fdr=FALSE){
  if(fdr){
    my.col <- c(2, 7);
  }
  db.i <- grep("MM", db);
  if (length(db.i)==0){
    mydb <- c("GO_hs", "KEGG", "PID_biocarta", "PID_KEGG", "PID_NCI_Nature_Curated", "ReactomePathways")
  } else {
    mydb <- c("GO_mm", "KEGG_mm", "PID_biocarta_mm", "PID_KEGG_mm", "PID_NCI_Nature_Curated_mm");
  }
  
  out.s <- c();
  for (db in mydb){
    # by known genome
    pofile.s <- paste(db, ".RData", sep="");
    pofile   <- system.file("data", pofile.s, package="OmicPath");
    GSEA.db <- local(get(load(pofile)));
    g.list  <- GSEA.db$pathway.list;
    g.table <- GSEA.db$pathway.table;
    g.ver   <- GSEA.db$pathway.version;
    gene.db.sym   <- unique(unlist(g.list));
    gene.db.num   <- length(gene.db.sym);
    gene          <- gene[!is.na(gene)];
    ## the gene list in GSEA.db 
    all.gene      <- unique(unlist(GSEA.db$pathway.list));
    gene.i        <- which(all.gene %in% gene);
    gene.num      <- length(gene.i);
    for (i in 1:length(g.list)){
      g.sub  <- unlist(g.list[[i]]);
      g.i    <- which(g.sub %in% gene);
      if (length(g.i)==0){
        next;
      }
      g.n    <- paste(g.sub[g.i], collapse=" ");
      p.sym  <- as.character(g.table[i,1]);
      p.name <- as.character(g.table[i,3]);
      db.mn  <- gene.db.num;
      db.group      <- length(g.sub);
      seleced.group <- length(g.i);
      seleced.k     <- gene.num;
      p.v    <- phyper (seleced.group, db.group, db.mn - db.group, seleced.k, lower.tail=F);
      out.s <- rbind(out.s, c(seleced.group, seleced.k, db.group, db.mn, p.v, g.n, p.sym, p.name, db));
    }
  }
  colname.s <- c("SigGeneNumInGroup",	"SigGeneNum",	"GeneNumInGroup",	"GeneNumInDB",	"Pvalue",	"Gene", "item.ID",	"item.name", "db");
  colnames(out.s) <- colname.s;
  if (filter.num > 0){
    drop.i <- which(out.s[,1] <= filter.num);
    if (length(drop.i)==nrow(out.s)){
      print("No result if the filter is done.")
    } else {
      if (length(drop.i) > 0){
        out.s <- out.s[-drop.i,];
      }
    }
  }
  tmp.d <- dim(out.s);
  if(length(tmp.d) < 1){
    print("FDR p value can not be calculated due to one row only.");
    out.s <- t(out.s);
    colnames(out.s) <- colname.s;
    return(out.s);
  }
  if (fdr){
    p.adjust.M <- p.adjust.methods[c(4,7)];
    p.adj      <- sapply(p.adjust.M, function(meth) p.adjust(out.s[,5], meth));
    out.p      <- cbind(out.s, p.adj);
    m.res.out  <- out.p[order(as.numeric(out.p[,5])),];
  } else {
    m.res.out <- out.s
  }
  return(m.res.out);
}

doGSEA <- function(db=db, gene=gene, filter.num=0, fdr=FALSE){
  if(fdr){
    my.col <- c(2, 7);
  }
  if (is.list(db)==T) {
    # based on the input GSEA.db
    g.list  <- db$pathway.list;
    g.table <- db$pathway.table;
    g.ver   <- db$pathway.version;
  } else if (nchar(db) >= 1){
    db1 <- toupper(db);
    if (db1 == "ALL" | db1 == "ALL.MM"){
      m.res.out <- doGSEA_all(db=db1, gene=gene, filter.num=filter.num, fdr=fdr);
      return(m.res.out);
    }
    # by known genome
    pofile.s <- paste(db, ".RData", sep="");
    pofile   <- system.file("data", pofile.s, package="OmicPath");
    if (file.exists(pofile)){
    } else {
      stop ("db name?");
    }
    GSEA.db <- local(get(load(pofile)));
    g.list  <- GSEA.db$pathway.list;
    g.table <- GSEA.db$pathway.table;
    g.ver   <- GSEA.db$pathway.version;
  } else {
    stop("db name ??")
  }
  
  gene.db.sym   <- unique(unlist(g.list));
  gene.db.num   <- length(gene.db.sym);
  gene          <- gene[!is.na(gene)];
  ## the gene list in db 
  all.gene      <- unique(unlist(g.list));
  gene.i        <- which(all.gene %in% gene);
  gene.num      <- length(gene.i);
  out.s <- c();
  for (i in 1:length(g.list)){
    g.sub  <- unlist(g.list[[i]]);
    g.i    <- which(g.sub %in% gene);
    if (length(g.i)==0){
      next;
    }
    g.n    <- paste(g.sub[g.i], collapse=" ");
    p.sym  <- as.character(g.table[i,1]);
    p.name <- as.character(g.table[i,3]);
    db.mn  <- gene.db.num;
    db.group      <- length(g.sub);
    seleced.group <- length(g.i);
    seleced.k     <- gene.num;
    p.v    <- phyper (seleced.group, db.group, db.mn - db.group, seleced.k, lower.tail=F);
    out.s  <- rbind(out.s, c(seleced.group, seleced.k, db.group, db.mn, p.v, g.n, p.sym, p.name));
  }
  colname.s <- c("SigGeneNumInGroup",	"SigGeneNum",	"GeneNumInGroup",	"GeneNumInDB",	"Pvalue",	"Gene", "item.ID",	"item.name");
  if (class(out.s)=="NULL"){
    print.s <- paste("No result: no gene in the item");
    return()
  } else {
    colnames(out.s) <- colname.s;
  }
  if (filter.num > 0){
    drop.i <- which(out.s[,1] <= filter.num);
    if (length(drop.i)==nrow(out.s)){
      print("No result if the filter is done.")
    } else {
      if (length(drop.i) > 0){
        out.s <- out.s[-drop.i,];
      }
    }
  }
  tmp.d <- dim(out.s);
  if(length(tmp.d) < 1){
    print("FDR p value can not be calculated due to one row only.");
    out.s <- t(out.s);
    colnames(out.s) <- colname.s;
    return(out.s);
  }
  if (fdr){
    p.adjust.M <- p.adjust.methods[c(4,7)];
    p.adj      <- sapply(p.adjust.M, function(meth) p.adjust(out.s[,5], meth));
    out.p      <- cbind(out.s, p.adj);
    m.res.out  <- out.p[order(as.numeric(out.p[,5])),];
  } else {
    m.res.out <- out.s
  }
  return(m.res.out);
}


