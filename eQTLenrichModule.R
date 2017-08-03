library(readxl)
options(stringsAsFactors = F)
dataDir <-"/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/data"
workDir <- "/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/"

setwd(dataDir)


table <- data.frame(Target=character(),
                    Color = character(),
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

for(target in c("B1T","B2G","B1G","B2T")){
##prepare snp counts fro raw eqtl table
eqtls <- df[grepl(paste(target),df$Source),]

load(paste0(target,".RData"))
result <- list()
nsiggene <- length(unique(eqtls$eGENE))
ngene <- ncol(datExpr)

for(color in unique(moduleColors)){
  a <- length(unique(eqtls[eqtls$eGENE %in% colnames(datExpr[,moduleColors==color]),2]))
  b <- nsiggene - a
  c <- sum(moduleColors==color) - a 
  d <- ngene -a - b -c
  result[[paste(target,color,sep="_")]] <- fisher.test(matrix(c(a,b,c,d),nrow = 2,ncol = 2),alternative = "greater")
  table <- rbind(table,c(target,color,a,b,c,d,result[[paste(target,color,sep="_")]]$p.value,result[[paste(target,color,sep="_")]]$estimate))
}
}
colnames(table)<-c("Tissue","ModuleColor","a","b","c","d","P.value","OR")
setwd(workDir)
WriteXLS::WriteXLS(table,"eQTLenriched_WGCNAmodules.xls",row.names = F,col.names = T)



