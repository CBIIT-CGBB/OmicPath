rm(list=ls());

library(OmicPath);

cols  <- rainbow(10, alpha=0.3)[c(1,2,7)]

pid <- "hsa04014";
pdf("02do_KEGGplot.pdf", 12, 8);

## it will return gene coordinates.
out <- KEGGplot(pid=pid, Taxonomy="hs", plot=T);

set.seed(1234);
## random select some genes and plot
out.i <- sample(1:nrow(out), 30);
x     <- as.numeric(out[out.i,2]);
y     <- as.numeric(out[out.i,3]);
points(x, y, cex=5, col=cols, pch=19);

dev.off();
