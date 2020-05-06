

db_names <- function(Taxonomy="hs"){
  db.names <- data(package="OmicPath")$results[,3];
  db.names <- gsub("GSEA.db \\(", "", db.names);
  db.names <- gsub(")", "", db.names);
  db.i     <- grep("Mouse", db.names);
  db.j     <- grep("mm", db.names);
  db.k     <- c(db.i, db.j);
  if (Taxonomy=="mm"){
    db.names <- db.names[db.k];
  } else {
    db.names <- db.names[-db.k];
  }
  return(db.names);
}