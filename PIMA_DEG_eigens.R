setRepositories(ind = c(1,2,3,4,5,6,7,8,9))

library(limma)
library(genefilter)
library(readxl)
library(dplyr)
library(reshape2)
library(splitstackshape)

##set work and data directory
workingDir<-"//eu.boehringer.com/users/rdg/users4/jwang6/Documents/My Received Files/PIMA/scripts/expression/DEG"
rawDataDir <- "//eu.boehringer.com/users/rdg/users4/jwang6/Documents/My Received Files/PIMA/scripts/expression/array"
setwd(dir = rawDataDir)

IQR.threshold <- 0.5
pvalue = 0.1

# read in expression data
b1 <- readRDS(file = "biopsy1.rds")
b2 <- readRDS(file = "biopsy2.rds")

# read in phenotype data
phenoData <- read_excel(path = "Calculations_EigenTraits_rmSampl.xlsx",
                        sheet = 1,
                        col_names = T)
phenoData <- as.data.frame(phenoData)
# formatting sample ids and add T and G annotation to match the expression data
phenoData$index = 1:nrow(phenoData)
phenoData = rbind(phenoData,phenoData)
phenoData[c(1:135),9] <- "T" 
phenoData[c(136:270),9] <- "G" 
colnames(phenoData)[9] <- "Tissue"
for (i in c(1:dim(phenoData)[1])){
    phenoData[i,1] <- paste("PIMA",phenoData[i,1],phenoData[i,9],sep="_")
}

# divide pheno data by biopsy number
phenoData$Target <- paste("biopsy",phenoData$BiopsyNumber,phenoData$Tissue, sep = "_")
ph1 <- subset(phenoData, BiopsyNumber == 1)
rownames(ph1) <- ph1$NIH
ph2 <- subset(phenoData, BiopsyNumber == 2)
rownames(ph2) <- ph2$NIH

# create collection of biopsies
collBiopsies <- list(b1 = list("biopsy_1",b1,ph1),
                     b2 = list("biopsy_2",b2,ph2))
#trait variables
CovOI <- colnames(phenoData)[c(3:(ncol(phenoData)-2))]

#loop through each tissue&biopsy group
for (currBT in unique(phenoData$Target) ) {
  dmIndex <- paste("b",substr(x = currBT, 8, 8 ), sep = "")
  ph <- subset(collBiopsies[[dmIndex]][[3]], Target == currBT)
  dataMatrix <- collBiopsies[[dmIndex]][[2]]
  dataMatrix <- as.matrix(dataMatrix[,intersect(colnames(dataMatrix),as.character(ph$NIH))])
  message("Samples available for ",currBT,":",nrow(ph))
  message("\tdata matrix dimension:", paste(dim(dataMatrix), collapse = "x"))
  message("\n")
  
  #loop throuth each covariate
  for (myCurrCov in CovOI) {
    message("Processing effect of:", myCurrCov)
    
    myRegMod <- as.formula(paste(" ~", paste(myCurrCov), sep = ""))
    message("Using following model:",myRegMod)
    
    myDesign <- try(model.matrix( myRegMod, data = ph[intersect(rownames(ph),colnames(dataMatrix)),])) 
    message("Running limma...")
    fit <- lmFit(dataMatrix[,intersect(rownames(myDesign),colnames(dataMatrix))],myDesign)
    fit <- eBayes(fit)
    message("done")
      
      message("Running regression...")
      fullFmla <- as.formula(paste("value ~", paste(myCurrCov), sep = ""))
      #redFmla <- as.formula(paste("value ~", paste( CovFixedVec, collapse = "+"), sep = ""))
      dmLongAll <- dataMatrix[,rownames(myDesign)] %>% melt()
      colnames(dmLongAll) <- c("geneName", "SampleID", "value")
      dmLongAll[,CovOI] <- ph[as.character(dmLongAll$SampleID),CovOI];
      fullModel = dmLongAll %>% group_by(geneName) %>% do(fitGene = glm(fullFmla, data = ., family = "gaussian"))
      #redModel = dmLongAll %>% group_by(geneName) %>% do(fitGene = glm(redFmla, data = ., family = "gaussian"))
      
      #stopifnot(unique(fullModel$geneName == redModel$geneName) == TRUE)
      
      message("Regression completed")
      
      fitFileName <- paste("fit_",myCurrCov,"_",currBT,".RData", sep = "")
      # if is ordinal, then extract the linear coefficient
      myCoeff <- character()
      if (is.ordered(ph[,myCurrCov])) {
        if (length(levels(ph[,myCurrCov])) > 1 ) {
          #linear
          myCoeff[1] <- paste(myCurrCov,"L",sep = ".")
          #quadratic
          myCoeff[2] <- paste(myCurrCov,"Q",sep = ".")      
        } else {
          myCoeff[1] <- paste(myCurrCov,levels(ph[,myCurrCov])[2], sep = "")
        }
      } else {
        myCoeff[1] <- myCurrCov
      }
      
      for (myCurrCoeff in myCoeff) {
        message("Extracting following coefficient:", myCurrCoeff)
        # RStat <- sapply(1:length(fullModel$geneName), function(kk) { 
        #   pr2 <- sqrt((sum(resid(redModel$fitGene[[kk]]) ^ 2) - sum(resid(fullModel$fitGene[[kk]]) ^ 2))/sum(resid(redModel$fitGene[[kk]]) ^ 2)) *
        #     sign(fullModel$fitGene[[kk]]$coefficients[myCurrCoeff])
        #   return(setNames(nm = fullModel$geneName[[kk]], object = pr2))
        # })
        
        output <- topTable(fit, coef = myCurrCoeff, number = Inf, adjust.method = "BH", sort.by = "p")
        output <- as.data.frame(output)
        output <- cbind(geneName = rownames(output),
                        output,
                        # rstat = RStat[rownames(output)],
                        stringsAsFactors = F)
        compName <- paste(currBT,"_",myCurrCoeff, sep = "")
        outfileName <-  paste("deg_",compName,".txt", sep = "")
        
        write.table(x = output,
                    file = outfileName,
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = T
        )
        output <- subset(output, adj.P.Val < pvalue)
        message("Significant genes:",nrow(output))
      }
  }

}






PIMAIDtoNIH <- function(x) {
  for (i in c(1:length(colnames(x)))){
    pimaID <- colnames(x)[i]
    k <- strsplit(pimaID,"_")
    colnames(x)[i] <- k[[1]][[2]]
  }
}
