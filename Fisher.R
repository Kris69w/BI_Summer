library(readxl)
options(stringsAsFactors = F)
dataDir <-"/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/data"
workDir <- "/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/"

setwd(dataDir)
target <- read_excel("SigTraits&targetWGCNA.xlsx")
##prepare snp counts fro raw eqtl table
eqtls <- read.table("rs_b2glom.broad.emmax.cis_extract.txt",header = T,sep = "\t",stringsAsFactors = F)
#count=lapply(unique(eqtls$MarkerID), function(x){length(eqtls[eqtls$MarkerID==x,1])})
#rawMarkercount <- data.frame(MarkerID = unique(eqtls$MarkerID),count=unlist(count))
#write.csv(rawMarkercount,"B1T.count.csv",row.names = F,col.names = T,quote = F)
rawMarkercount <- read.csv("B2G.count.csv")

load("PI_G_Bx2-networkConstruction-signpwr7.RData")
result <- list()
mastertab <- eqtls[,c(1,5,7,15)]

for(color in unlist(unique(target[target$Target=="B2G",3]))){
  eqtlsub <- eqtls[eqtls$GENE %in% colnames(datExpr[,moduleColors==color]),]
  table <- data.frame(MarkerID = character(),
                      a=integer(),
                      b=integer(),
                      c=integer(),
                      d=integer(),
                      PVal=double(),
                      OR=double())
  for(marker in unique(eqtlsub$MarkerID)){
    a=length(eqtlsub[eqtlsub$MarkerID==marker,1])
    b=nrow(eqtlsub)-a
    c=rawMarkercount[rawMarkercount$MarkerID==marker,2]-a
    d=nrow(eqtls)-nrow(eqtlsub)-c
    pval <- fisher.test(matrix(c(a,b,c,d),nrow = 2,ncol = 2),alternative = "greater")
    table <- rbind(table,c(marker,a,b,c,d,pval$p.value,pval$estimate))
  }
  colnames(table)<-c("MarkerID","a","b","c","d","P.value","OR")
  result[[paste(color)]] <- table
  colnames(table)[6]<- paste(color,"Pval",sep="_")
  mastertab <- merge(mastertab,table[,c(1,6)],by="MarkerID",all.x = T)
}
setwd(workDir)
write.csv(mastertab,"B2G_SigLocus.csv",quote = F,row.names = F)
save(rawMarkercount,result,mastertab,file = "B2G_Fisher.RData")


Genelist <- list()
for(i in 5:ncol(mastertab)){
  mastertab[,i] <- as.numeric(mastertab[,i])
  sub <- mastertab[mastertab[,i]<0.05,c(3,i)]
  sub <- sub[complete.cases(sub),]
  Genelist[[paste(colnames(mastertab)[i])]] <- unique(sub$GENE)
}
