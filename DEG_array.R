# Set R terminal widht to actual terminal size
#options(width = Sys.getenv("COLUMNS"))
# Make all repositories available for package download:
setRepositories(ind = c(1,2,3,4,5,6,7,8,9))

library(limma)
library(genefilter)
library(readxl)
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
workingDir <- paste(pathPrefix , "collaborators/junke/eigenTrait_DEG/Raw", sep = "");
rawDataDir <- paste(pathPrefix , "collaborators/junke/eigenTrait_DEG/Raw", sep = "");
#
IQR.threshold <- 0.5
pvalue = 0.1

#set number of Cores available
if ( grepl(x = R.Version()[["os"]], pattern = "linux" )  ) {
  threads <- as.numeric(system(command = "grep -c ^processor /proc/cpuinfo",intern = T))
} else {
  threads <- sum(suppressWarnings(as.numeric(system(command = "wmic cpu get NumberOfCores",intern = T))),na.rm = T)
}
cat(paste("detected ",threads," cores\n",sep = ""))

# Read in data -------------------------------------------------------
setwd(dir = rawDataDir)
# read in expression for pima ONLY!
b1 <- readRDS(file = "biopsy1.rds")
b2 <- readRDS(file = "biopsy2.rds")


phenoData <- read_excel(path = "sampleInfo_Array.xlsx",
                        sheet = 1,
                        col_names = T)
phenoData <- as.data.frame(phenoData)
# phenoData$uACR_Class <- ordered(phenoData$uACR_Class, levels = c("Norm","MicroAlbuminuria","MacroAlbuminuria"))
phenoData$uACR_Class <- NULL
phenoData$Subclass_5RG <- NULL
phenoData$Target <- paste("biopsy",phenoData$BiopsyNumber,phenoData$Tissue, sep = "_")

ph1 <- subset(phenoData, BiopsyNumber == 1)
rownames(ph1) <- ph1$SampleID
ph2 <- subset(phenoData, BiopsyNumber == 2)
rownames(ph2) <- ph2$SampleID

setwd(dir = workingDir)

# create collection of biopsies
collBiopsies <- list(b1 = list("biopsy_1",b1,ph1),
                     b2 = list("biopsy_2",b2,ph2))

CovOI <- colnames(phenoData)[c(15:(ncol(phenoData) - 2))]
# CovFixed <- c("+ AgeAtVisit + SEX")
# CovFixedVec <- c("AgeAtVisit","SEX")
CovFixed <- c("+ DurationOfDibetes + SEX")
CovFixedVec <- c("DurationOfDibetes","SEX")

# loop tissues within biopsy -----------------------------------------------

currBT <- "biopsy_2_G"
#for (currBT in unique(phenoData$Target) ) {
  dmIndex <- paste("b",substr(x = currBT, 8, 8 ), sep = "")
  ph <- subset(collBiopsies[[dmIndex]][[3]], Target == currBT)
  dataMatrix <- collBiopsies[[dmIndex]][[2]]
  dataMatrix <- as.matrix(dataMatrix[,as.character(ph$SampleID)])
  message("Samples available for ",currBT,":",nrow(ph))
  message("\tdata matrix dimension:", paste(dim(dataMatrix), collapse = "x"))
  message("\n")
  for (myCurrCov in CovOI) {
    message("Processing effect of:", myCurrCov)
    
    # hack to skip
    if (file.exists(paste("deg_",paste(currBT,"_",myCurrCov, sep = ""),".txt", sep = ""))) {
      message("output file already present, skipping")
      next;
    }
    
    
    myRegMod <- as.formula(paste(" ~", paste(myCurrCov, CovFixed, sep = ""), sep = ""))
    message("Using following model:",myRegMod)

    myDesign <- try(model.matrix( myRegMod, data = ph)) 
    if (inherits(myDesign, "try-error")) { 
      message("Skipping due to lack of enough patients in one of the covariates to correct for")
      next
    } else {
      myDesign
    }
        message("Patients available with this model:",nrow(myDesign))
    if (nrow(myDesign) > 0 ) {
      message("Patients discarded because missing model values:\n",
              paste(setdiff(rownames(ph),rownames(myDesign)), collapse = ","))
      chkLev <- apply(myDesign[,2:ncol(myDesign)],2,FUN = function(x) {
        # chech if has less than 2 levels 
        length(levels(as.factor(x))) < 2
      })
      if (any(chkLev)) {
        message("Skipping due to lack of enough patients in one of the covariates to correct for")
        next;
      }
      message("Running limma...")
      fit <- lmFit(dataMatrix[,rownames(myDesign)],myDesign)
      fit <- eBayes(fit)
      message("done")
      
      message("Running regression...")
      fullFmla <- as.formula(paste("value ~", paste(myCurrCov, CovFixed, sep = ""), sep = ""))
      redFmla <- as.formula(paste("value ~", paste( CovFixedVec, collapse = "+"), sep = ""))
      dmLongAll <- dataMatrix[,rownames(myDesign)] %>% melt()
      colnames(dmLongAll) <- c("geneName", "SampleID", "value")
      dmLongAll[,c(CovOI,CovFixedVec)] <- ph[as.character(dmLongAll$SampleID),c(CovOI,CovFixedVec)];
      fullModel = dmLongAll %>% group_by(geneName) %>% do(fitGene = glm(fullFmla, data = ., family = "gaussian"))
      redModel = dmLongAll %>% group_by(geneName) %>% do(fitGene = glm(redFmla, data = ., family = "gaussian"))
      
      stopifnot(unique(fullModel$geneName == redModel$geneName) == TRUE)
      
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
        RStat <- sapply(1:length(fullModel$geneName), function(kk) { 
          pr2 <- sqrt((sum(resid(redModel$fitGene[[kk]]) ^ 2) - sum(resid(fullModel$fitGene[[kk]]) ^ 2))/sum(resid(redModel$fitGene[[kk]]) ^ 2)) *
            sign(fullModel$fitGene[[kk]]$coefficients[myCurrCoeff])
          return(setNames(nm = fullModel$geneName[[kk]], object = pr2))
        })
  
        output <- topTable(fit, coef = myCurrCoeff, number = Inf, adjust.method = "BH", sort.by = "p")
        output <- as.data.frame(output)
        output <- cbind(geneName = rownames(output),
                        output,
                        rstat = RStat[rownames(output)],
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
    } else {
      message("Nothing to do here...")
    }
  }
#}
system("mv pheno_Plot.rds pheno_Plot_Old.rds", wait = T, show.output.on.console = T)
saveRDS(phenoData[,c(1,3,7:9,13:18,20:68,72,119:127)],file = "pheno_Plot.rds")














