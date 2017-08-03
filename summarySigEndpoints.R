# Set R terminal widht to actual terminal size
#options(width = Sys.getenv("COLUMNS"))
# Make all repositories available for package download:
setRepositories(ind = c(1,2,3,4,5,6,7,8,9))

library(limma)
library(genefilter)
library(xlsx)
library(dplyr)
library(reshape2)


# define Platform and initial path ----------------------------
if ( grepl(x = R.Version()["os"], pattern = "linux" ) ) {
  pathPrefix <- "/data/SKDngs01/pg/"
} else {
  pathPrefix <- "X:/pg/"
}
# Set your working directory where output files will be saved ---------
#
workingDir <- paste(pathPrefix , "Projects/Pima/expression/DEG", sep = "");
rawDataDir <- paste(pathPrefix , "Projects/Pima/expression/array", sep = "");
#
pvalue = 0.05


setwd(dir = workingDir)
# Import PIMA dge ---------------------------------------------------------
# Pima
pimaComps <- readRDS(paste(pathPrefix , "Projects/Pima/expression/DEG/DEGTableWithAllComps.rds", sep = ""))

# Process dge -------------------------------------------------------------
# define comps of interestc
# s1 <- pimaComps[,intersect(grep(x = colnames(pimaComps), pattern = "_GFR_"),
#                            grep(x = colnames(pimaComps), pattern = "P.Val"))]
s1 <- pimaComps[,grep(x = colnames(pimaComps), pattern = "P.Val")]

summAll <- data.frame(apply(s1,2, function(x) {length(which(x < 0.05))}))
colnames(summAll) <- "SigGenes"
summAll$Comp <- rownames(summAll)
summAll$Biopsy <- substr(x = summAll$Comp,start = 3,stop = 3)
summAll$Tussue <- substr(x = summAll$Comp,start = 5,stop = 5)
summAll$Trait <- gsub(pattern = "_adj.P.Val",replacement = "",ignore.case = F,fixed = T,x = substr(summAll$Comp,7,1000))
WriteXLS::WriteXLS(x = c("summAll"),
                   ExcelFileName = "summarySigEndpoints.xlsx",
                   SheetNames = "summAll",
                   verbose = T,
                   col.names = T,
                   row.names = F,
                   BoldHeaderRow = T)

resColl <- list()

for (col in colnames(s1)) {
  message("Processing:",col)
  sigGenes <- which(s1[,col] < pvalue)
  if (length(sigGenes) > 0) {
    message("\tsig gene:", length(sigGenes))
    name <- gsub(x = col, pattern = "_adj.P.Val", replacement = "")
    resColl[[name]] <- as.character(pimaComps[sigGenes, "Symbol"])
  }
}


