
WGCNA_6_6 <- function(gdir="GSEA", mdir="module", db=db){
  subDir        <- gdir;
  if (dir.exists(subDir)){
  } else {
    dir.create(subDir);
  }
  files  <- dir(mdir, pattern = "^module", full.names = TRUE);
  file2  <- dir(mdir, pattern = "^module");
  file.n <- length(files);
  for (k in 1:length(db)){
    the.db <- db[k];
    for (i in 1:file.n){
      infile     <- files[i];
      g.name     <- sub("module_", "", file2[i]);
      g.name     <- sub(".txt", "", g.name);
      dat        <- read.table(infile, header=T);
      outfile    <- paste0(gdir, "/", g.name, "_", db[k], ".xls");
      out        <- doGSEA(db=the.db, gene=dat[,1], filter.num=1, fdr=T);
      write.table(out, outfile, quote=F, sep="\t", row.names=F);
    }
  }
}


