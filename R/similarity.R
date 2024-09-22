
similarity <- function(s1, s2){
  if (length(s1) < length(s2)){
    set1 <- s2;
    set2 <- s1;
  } else {
    set1 <- s1;
    set2 <- s2;
  }
  ov.s  <- intersect(set1, set2);
  ov    <- length(ov.s);
  s1.l  <- length(set1);
  s2.l  <- length(set2);
  jac   <- jaccard_sim(set1, set2);
  jac.d <- 1 - jac;
  ove   <- overlap_coef(set1, set2);
  cos   <- cosine_sim(set1, set2);
  fte   <- ftest2sets(set1, set2)$p.value;
  hte   <- htest2sets(set1, set2);
  out.v <- c(as.integer(ov), as.integer(s1.l), as.integer(s2.l), jac, jac.d, ove, cos, fte, hte);
  out.s <- c("overlapping.number", "longer.set.number", "shorter.set.number", 
             "jaccard.similarity.index", "jaccard.distance", "overlap.coefficient",
             "cosine.similarity.coefficient", "fishers.test.pvalue", 
             "hypergeometric.test.pvalue")
  return(list(results=out.v, names=out.s));
}

## Jaccard similarity index:
jaccard_sim <- function(x, y){
  intersect <- length(intersect(x, y))
  union <- length(union(x, y))
  return(intersect/union)
}

## Jaccard distance == 1 - jaccard_sim
jaccard_distance <- function(set1, set2) {
  1 - length(intersect(set1, set2)) / length(union(set1, set2))
}

## Overlap coefficient:
overlap_coef <- function(x, y){
  intersect <- length(intersect(x, y))
  min_size <- min(length(x), length(y))
  return(intersect/min_size)
}

## Cosine similarity:
cosine_sim <- function(x, y){
  dot_product <- length(intersect(x, y))
  mag_x <- sqrt(length(x))
  mag_y <- sqrt(length(y))
  return(dot_product/(mag_x*mag_y))
}

## fisher's test
ftest2sets <- function(gene_set_1, gene_set_2, universe_gene_number=20000){
  # Calculate the overlap between the gene sets
  overlap_genes <- intersect(gene_set_1, gene_set_2)
  
  # Calculate the number of genes unique to each set
  unique_genes_set_1 <- length(setdiff(gene_set_1, gene_set_2))
  unique_genes_set_2 <- length(setdiff(gene_set_2, gene_set_1))
  
  # Calculate the number of genes in the universe not in either set
  non_set_genes <- universe_gene_number - length(unique(c(gene_set_1, gene_set_2)))
  
  # Create a contingency table
  contingency_table <- matrix(c(length(overlap_genes), unique_genes_set_1,
                                unique_genes_set_2, non_set_genes),
                              nrow = 2,
                              byrow = TRUE,
                              dimnames = list(c("Set 1", "Set 2"),
                                              c("Overlap", "No Overlap")))
  # Perform Fisher's exact test
  test_result <- fisher.test(contingency_table)
  return(test_result)
}

## hypergeometric test
htest2sets <- function(gene_set_1, gene_set_2, universe_gene_number=20000){
  
  # Calculate the overlap between the gene sets
  overlap_genes <- length(intersect(gene_set_1, gene_set_2))
  
  # Calculate the number of genes unique to each set
  unique_genes_set_1 <- length(setdiff(gene_set_1, gene_set_2))
  unique_genes_set_2 <- length(setdiff(gene_set_2, gene_set_1))
  
  # Calculate the number of genes in the universe not in gene set 1
  non_set1_genes <- universe_gene_number - unique_genes_set_1
  
  # Perform the hypergeometric test using the phyper function
  p_value <- phyper(overlap_genes, unique_genes_set_1, non_set1_genes, unique_genes_set_2, lower.tail = FALSE)
  
  return(p_value)
}

