# OmicPath: an R package for gene set enrichment analysis and pathway network analysis 
#### Installation
```r
library(devtools)
install_github("CBIIT-CGBB/OmicPath")
``` 
 
 
# Examples:
### Perform gene set analysis (GSA)
```r
set.seed(1234);
## get some genes from biosystems database of NCBI
dat   <- read.table(system.file("extdata/biosystems_hs.txt", package = "OmicPath"));
genes <- sample(dat[,1], 200);

## select databases (GO) from the package
db  <- db_names()[1];
## do GSA
out <- doGSEA(db=db, gene=genes, filter.num=2, fdr=T);
``` 
 
### Map gene information on the pathway
<img src="examples/02do_KEGGplot.png" width="650" height="360">
  
[code](examples/02do_KEGGplot.R)

### Extract gene relationships in the pathway and plot the pathway as network clusters
  
KEGG: Ras signaling pathway
  
<img src="examples/03data_network.png" width="200" height="400">  <img src="examples/03plot_network.png" width="400" height="400">
  
[R codes](examples/03plot_network.R)

### Construct networks by gene neighborhoods in the pathway
KEGG: Ras signaling pathway: construct networks by the genes ("PAK2", "BAD", "RASGRF1", "RAP1A", "NRAS", "HRAS", "TIAM1", "RRAS2", "KRAS", "SHC3") and their neighborhoods. The networks contain genes with the neighborhood order = 1 and 2. The seed genes were colored in red.

<img src="examples/04do_neighborhood_test1.png" width="300" height="300">  <img src="examples/04do_neighborhood_test2.png" width="300" height="300">
  
[R codes (order=1)](examples/04do_neighborhood_test1.R)
  
[R codes (order=2)](examples/04do_neighborhood_test2.R)


### Additional example: cluster algorithms in igraph
<img src="examples/04do_igraph_cluster.png" width="500" height="500">
  
[R codes](examples/04do_igraph_cluster.R)
