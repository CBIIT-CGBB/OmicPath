KEGGplot <- function(pid=pid, Taxonomy="hs", plot=TRUE){
  root    <- pid;
  if (Taxonomy == "mm"){
    file.c.in <- paste0("mmu/coord/", root, ".txt"); 
    file.p.in <- paste0("mmu/image/", root, ".png"); 
  } else {
    file.c.in <- paste0("hsa/coord/", root, ".txt"); 
    file.p.in <- paste0("hsa/image/", root, ".png"); 
  }

  file.c    <- system.file("extdata",  file.c.in, package="OmicPath");
  if(!file.exists(file.c)){
    stop("There is not ", file.c.in, " at extdata!")
  }
  file.p  <- system.file("extdata",  file.p.in, package="OmicPath");
  dat.c   <- read.table(file.c, header=T, sep="\t");
  dat2    <- dat.c;
  
  Image   <- EBImage::readImage(file.p);
  max.x  <- dim(Image)[1];
  max.y  <- dim(Image)[2];
  ima    <- png::readPNG(file.p);

  if (plot){
    plot(0:max.x, seq(0, max.y, length.out = (max.x+1)), type='n', main="", xlab="", ylab="", axes=F);
    lim <- par()
    rasterImage(ima, 0, 0, max.x, max.y);
  }

  out.s <- NULL;
  for (g in 1:nrow(dat2)){
    if (is.na(dat2[g,3])){
      next;
    }
    xy.co <- unlist(strsplit(as.character(dat2[g,2]), " "));
    xy.co <- as.numeric(xy.co);
    x <- (xy.co[1]+xy.co[3])/2;
    y <- max.y - (xy.co[2]+xy.co[4])/2;
    x <- as.numeric(x);
    y <- as.numeric(y);
    out.s <- rbind(out.s, c(as.character(dat2[g, 4]), x, y));
  }
  colnames(out.s) <- c("gene", "x", "y");
  return(out.s);
}

