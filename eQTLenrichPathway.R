library(readxl)
options(stringsAsFactors = F)
setwd("/data/SKDngs01/pg/collaborators/junke/PIMA/GeneList/data")

## cleaning WGCNA data
sheets <- c("B1T_magenta","B1G_greenyellow","B2T_brown","B2G_magenta","B2G_salmon","B2G_turquoise","B1T_green",
            "B1T_turquoise","B1T_blue","B1T_magenta","B1G_greenyellow","B1G_grey60")

genes <- lapply(sheets, function(X) readxl::read_excel("WGCNA_pathAnalysis.xlsx", sheet = X)[1])
pathways <- lapply(sheets, function(X) readxl::read_excel("WGCNA_pathAnalysis.xlsx", sheet = X)[-1,c(2,3,5)])

names(genes)<- sheets
names(pathways)<-sheets

pathways <- lapply(pathways, na.exclude)

dataDir <-"/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/data"
workDir <- "/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/"

setwd(dataDir)

table <- data.frame(Target=character(),
                    Pathway = character(),
                    a=integer(),
                    b=integer(),
                    c=integer(),
                    d=integer(),
                    PVal=double(),
                    OR=double())
df <- read.table("/data/SKDngs01/pg/Projects/Genetics/eQTL_DB.nodir.txt",sep = "\t",header = T)
df$Source <- gsub("B1_G","B1G",df$Source)
df$Source <- gsub("B2_G","B2G",df$Source)
df$Source <- gsub("B1_T","B1T",df$Source)
df$Source <- gsub("B2_T","B2T",df$Source)
somegenes<- character()

for(module in sheets){
  target <- unlist(strsplit(module,"_"))[1]
  color <- unlist(strsplit(module,"_"))[2]
  qtl <- df[grepl(paste(target),df$Source),]
  load(paste0(target,".RData"))
  eqtls <- qtl[qtl$eGENE %in% colnames(datExpr[,moduleColors==color]),]
  nsiggene <- length(unique(eqtls$eGENE))
  ngene <- sum(moduleColors==color)
  result <- list()
  
  for(i in 1:length(as.data.frame(pathways[paste(module)])[,3])){
    molegrp = as.data.frame(pathways[paste(module)])[i,3]
    moles <- strsplit(molegrp,",")
    moles <- intersect(unlist(moles), unlist(genes[paste(module)]))
    if(i>20){next}else{
      pathway <- as.data.frame(pathways[paste(module)])[i,1]
      a <- length(intersect(moles,unique(eqtls$eGENE)))
      b <- nsiggene - a
      c <- length(moles) - a 
      d <- ngene -a - b -c
      result[[paste(target,color,sep="_")]] <- fisher.test(matrix(c(a,b,c,d),nrow = 2,ncol = 2),alternative = "greater")
      table <- rbind(table,c(module,pathway,a,b,c,d,result[[paste(target,color,sep="_")]]$p.value,result[[paste(target,color,sep="_")]]$estimate))
    }
  }
} 
colnames(table)<-c("Target","Pathway","a","b","c","d","P.value","OR")
setwd(workDir)
WriteXLS::WriteXLS(table,"eQTLenriched_modulePathway20.xls",row.names = F,col.names = T)
