
mm2hs <- function(gene=gene){
  infile.s <- paste0("homologene.data_f.hs_mm")
  infile   <- system.file("extdata",  infile.s, package="OmicPath");
  dat      <- read.table(infile);
  gene.i   <- which(dat[,3] %in% gene);
  out      <- dat[gene.i,];
  return(out);
}