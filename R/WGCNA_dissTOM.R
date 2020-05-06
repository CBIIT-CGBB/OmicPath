
WGCNA_dissTOM <- function(tom=tom, modules=modules, indir=indir, outdir=outdir){
  packageStartupMessage("reading TOM, be patient please ...", appendLF = FALSE);
  dat.t <- read.table(tom, row.names=1, header=T);
  colnames(dat.t) <- row.names(dat.t)
  packageStartupMessage(" done")
  pb <- txtProgressBar(min = 0, max = length(modules), style = 3);
  for (f in 1:length(modules)){
    infile  <- paste0(indir, "/", modules[f]);
    outfile <- gsub(".txt", ".matrix", modules[f]);
    outfile <- paste0(outdir, "/", outfile);
    gene    <- read.table(infile, header=T);
    gene    <- unlist(gene);
    col.i   <- which(colnames(dat.t) %in% gene);
    row.i   <- which(row.names(dat.t) %in% gene);
    write.table(dat.t[row.i, col.i], outfile, quote=F, sep="\t");
    setTxtProgressBar(pb, f)
  }
  close(pb);
}