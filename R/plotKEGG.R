plotKEGG <- function(pid=pid, pname=pname, Taxonomy="hs", gene=gene, glab=F, fc.value=fc.value, 
                     p.value=p.value, pdffile=pdffile){
  cols    <- myColorRamp(c("blue", "white", "red"), fc.value); 
  #cols    <- paste0(cols, "80");
  names(cols)     <- gene;
  names(p.value)  <- gene;
  names(fc.value) <- gene;

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
  
  Image   <- EBImage::readImage(file.p);
  max.x  <- dim(Image)[1];
  max.y  <- dim(Image)[2];
  ima    <- png::readPNG(file.p);
  
  col1  <- rainbow(20, alpha=0.5);
  dat2  <- dat.c;
  pdf(pdffile, 12, 8);
  par(mar=c(4,4,4,5));
  
  plot(0:max.x, seq(0, max.y, length.out = (max.x+1)), type='n', main=pname, xlab="", ylab="", axes=F);
  lim <- par()
  rasterImage(ima, 0, 0, max.x, max.y);
  
  for (g in 1:nrow(dat2)){
    if (is.na(dat2[g,3])){
      next;
    }
    xy.co <- unlist(strsplit(dat2[g,2], " "));
    xy.co <- as.numeric(xy.co);
    x <- (xy.co[1]+xy.co[3])/2;
    y <- max.y - (xy.co[2]+xy.co[2])/2;  

    cex <- 0.4 + rnorm(1)/10;
    rect(xy.co[1], max.y-xy.co[2], xy.co[3], max.y-xy.co[4], 
         col= NA, border = cols[as.character(dat2[g,4])], lwd=-log10(p.value[as.character(dat2[g,4])]));
    if (glab){
      text(x, y, dat2[g,4]);
    }
  }
  par(xpd=NA);
  coord <- c(50+max.x, max.y, 100+max.x, max.y*0.8);
  color_bar(cols=cols, val=fc.value, coord=coord, lab.s=20, horizontal=F);
  legend(max.x*0.3, -1, c(2,3,4,5,6),  lwd=c(2,3,4,5,6), horiz=T, bty="n", title="-log10 Pvalue", col="blue");
  legend(max.x*0.3, -60, c(2,3,4,5,6), lwd=c(2,3,4,5,6), horiz=T, bty="n", col="red");
  dev.off();
  
}

