################################################################################
######## This file includes R code for obtaining the gene sets data
################################################################################

## Install packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!require("breastCancerMAINZ", quietly = TRUE)) 
  BiocManager::install("breastCancerMAINZ")
if (!require("limma", quietly = TRUE)) BiocManager::install("limma")
if (!require("hgu95av2.db", quietly = TRUE)) BiocManager::install("hgu95av2.db")

## load packages
library(breastCancerMAINZ)
library(Biobase)
library(limma)
library(hgu95av2.db)

################################################################################
######## search for gene sets corresponding to different GO terms
######## and contains at least 25 genes
################################################################################

## load the mainz expression data from package "breastCancerMAINZ"
data(mainz)
dat <- Biobase::exprs(mainz)
## get the group id of each observation 
id <- Biobase::pData(mainz)[, 9]

## normalize the data using 'normalizeQuantiles' in the 'limma' package
dat <- limma::normalizeQuantiles(dat, ties=TRUE)
dat0 <- t(dat)

#### functions to identify genes sets corresponding to different GO terms.
#### Use 'hgu95av2GO2ALLPROBES' in 'hgu95av2.db' package to maps between 
#### manufacturer IDs and Gene Ontology (GO) IDs.
GOTFfun <- function(GOID) {
  x <- hgu95av2.db::hgu95av2GO2ALLPROBES[[GOID]]
  unique(x)
}

GOTFDat <- function(GOID, dat0, id) {
  TFs <- unique(unlist(lapply(GOID, GOTFfun)))
  
  inSel <- match(TFs, colnames(dat0), nomatch=0)
  
  if (sum(inSel!=0) > 1) {
    dat_G1 <- dat0[which(id==1), inSel]
    dat_G2 <- dat0[which(id==2), inSel]
    dat_G3 <- dat0[which(id==3), inSel]
  } else if (sum(inSel!=0) == 1) {
    dat_G1 <- matrix(dat0[which(id==1), inSel], sum(id==1), 1)
    dat_G2 <- matrix(dat0[which(id==2), inSel], sum(id==2), 1)
    dat_G3 <- matrix(dat0[which(id==3), inSel], sum(id==3), 1)
  }
  return(list("dat_G1" = dat_G1, "dat_G2" = dat_G2, "dat_G3" = dat_G3))
}

#### obtain all GO IDs that link to manufacturer IDs in 'dat0'
## The manufacturer IDs are given by the column names of 'dat0'
geneIDs <- colnames(dat0)
tmp <- AnnotationDbi::select(hgu95av2.db, keys=geneIDs, column="GO")
GOIDs <- unique(tmp$GO)
GOIDs <- GOIDs[-which(is.na(GOIDs))]
# In total, we obtain 2576 candidate GO IDs that link to manufacturer IDs

#### obtain gene sets corresponding to each GO ID with at least 25 genes
geneSetDat <- vector("list", length(GOIDs))
i <- 1;  j <- 1
while (j <= length(GOIDs)) {
  GOID <- GOIDs[j]
  dat_tmp <- GOTFDat(GOID, dat0, id)
  # print(ncol(dat_tmp$dat_G1))
  if (ncol(dat_tmp$dat_G1) >= 25) {
    dat_G1 <- dat_tmp$dat_G1
    dat_G2 <- dat_tmp$dat_G2
    dat_G3 <- dat_tmp$dat_G3
    
    dat <- list(dat_G1, dat_G2, dat_G3)
    N <- c(nrow(dat_G1), nrow(dat_G2), nrow(dat_G3))
    p <- ncol(dat_G1)
    
    geneSetDat[[i]] = list("GOID"=GOID, "dat"=dat, "N"=N, "p"=p)
    i <- i + 1
  }
  j <- j + 1
}
geneSetDat1 <- vector("list", i-1)
for (ss in 1:(i-1)) {
  geneSetDat1[[ss]] <- geneSetDat[[ss]]
}
geneSetDat <- geneSetDat1
# In total, we have 149 gene sets with at least 25 genes.

saveRDS(geneSetDat, "geneSetDat.rds")

## get the number of genes in each gene set
P0 <- unlist(lapply(geneSetDat, function(x) ncol(x$dat[[1]])))
summary(P0)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 25.00   29.00   40.00   57.44   63.00  271.00