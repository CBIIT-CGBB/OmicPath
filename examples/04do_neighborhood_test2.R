rm(list=ls());

library(OmicPath);
library(wordcloud); 
library(qgraph);

options(stringsAsFactors = FALSE);

pid      <- "hsa04014";
infile.s <- paste0("relation_pathway_hsa.txt")
infile   <- system.file("extdata",  infile.s, package="OmicPath");
dat      <- read.table(infile, header=T);
dat.i    <- which(dat[,3]==pid);
dat.s    <- dat[dat.i,];
## remove some relationahips
dat.s    <- dat.s[-c(1:1000),];
dat.s    <- dat.s[-c(dat.s[,1]==dat.s[,2]),];
set.seed(123);
gene.i   <- sample(nrow(dat.s), 100);
## select 10 random genes 
gene.d  <- unique(dat.s[gene.i, 2])[c(1:10)];

order.i <- 2;
pdffile <- paste0("04do_neighborhood_test", order.i, ".pdf");
g0    <- igraph::graph.data.frame(dat.s, directed=T);
g     <- do_neighborhood(g0, genes = gene.d, order=order.i);
V(g)$label <- V(g)$name;

wc   <- igraph::walktrap.community(g);
col3 <- rainbow(length(wc), alpha=0.2);
col2 <- rainbow(10, alpha=0.4);

g4    <- delete_vertex_attr(g, "name");
e     <- get.edgelist(g4);
l     <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(g4), weights=E(g)$weight,
                                           area=8*(vcount(g4)^2),repulse.rad=(vcount(g4)^3.1))

vertex.size        <- 4;
vertex.color       <- col2[7];
vertex.frame.color <- vertex.color;

vertex.label.cex   <- log2(c(degree(g)+.1));
vertex.label.cex[vertex.label.cex < 0.4] <- 0.4;
vertex.label.cex[vertex.label.cex > 0.6] <- 0.6;
vertex.label.color <- rgb(0,0,0,1);

pdf(pdffile, 6, 6);
plot(g,
     mark.groups=wc, mark.col=col3, mark.border=col3, mark.shape=1,
     layout=l,
     ##
     edge.color         = col2[7], 
     edge.width         = 1,
     edge.arrow.size    = 0.2,
     ##
     vertex.size        = vertex.size, 
     vertex.color       = vertex.color,
     vertex.frame.color = vertex.frame.color,
     ##
     vertex.label.dist = 0.1, 
     vertex.label.cex   = vertex.label.cex, 
     vertex.label.color = vertex.label.color,
     ## 
     main               = "",
     xlab               = ""
)
dev.off();
