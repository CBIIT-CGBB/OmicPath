\name{WGCNA_6_4}
\alias{WGCNA_6_4}

\title{
 Step 4: generating modules for WGCNA 
}

\description{
 Step 4: generating modules for WGCNA 
}

\usage{
  WGCNA_6_4(expr=expr, pheno=pheno, mdir="module", cutHeight=200, 
  softPower=8, minM=25, cutH=0.996, deepSplit=1,  project.name="WGCNA_net");
}

\arguments{
   \item{expr}{expr: gene expression data}
   \item{pheno}{pheno: phenotype data}
   \item{mdir}{mdir: output directory for module}
   \item{cutHeight}{cutHeight: outlier threshold}
   \item{softPower}{softPower: power value from step 2}
   \item{minM}{minM: minimum node number per module from step 3}
   \item{cutH}{cutH: threshold for module selecting from step 3}
   \item{deepSplit}{deepSplit: deep split value from step 3}
   \item{project.name}{project.name: project name}
}

\author{
 Ying Hu <yhu@mail.nih.gov>
 Chunhua Yan <yanch@mail.nih.gov>
}

\references{
##
}

\examples{
##
}

