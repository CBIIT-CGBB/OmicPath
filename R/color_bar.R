
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], names = NULL, maxColorValue = 255)
}

#dat  <- matrix(rnorm(200), ncol=2);
#cols <-  myColorRamp(c("red", "white", "blue"), dat[,2]) 
#plot(dat[,1], dat[,2], pch=19, col=cols, main="colors by Y");

color_bar <- function(cols=cols, val=val, coord=coord, lab.s=lab.s, cex=0.4, horizontal=T){
  x1   <- coord[1];
  y1   <- coord[2];
  x2   <- coord[3];
  y2   <- coord[4];
  s.i  <- order(val);
  col1 <- cols[s.i];
  co.i <- seq(1, length(val), length.out=60);
  col2 <- col1[co.i];
  #col2 <- rev(col2);
  if (horizontal){
    step <- (x2-x1)/60;
    min  <- round(min(val), 2);
    max  <- round(max(val), 2);
    text(x1, y1+lab.s, min, cex=cex);
    text(x2, y1+lab.s, max, cex=cex);
    x.1 <- x1;
    x.2 <- x1 + step;
    for (i in 1:60){
      rect(x.1, y1, x.2, y2, col=col2[i], border=NA);
      x.1 <- x.2;
      x.2 <- x.1 + step;
    }
  } else {
    step <- (y2-y1)/60;
    min  <- round(min(val), 2);
    max  <- round(max(val), 2);
    text(x1-lab.s, y1, max, cex=cex);
    text(x1-lab.s, y2, min, cex=cex);
    y.2 <- y2;
    y.1 <- y.2 - step;
    for (i in 1:60){
      rect(x1, y.1, x2, y.2, col=col2[i], border=NA);
      y.1 <- y.2;
      y.2 <- y.1 - step;
    }
  }
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

