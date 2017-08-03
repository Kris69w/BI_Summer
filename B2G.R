library(readxl)
options(stringsAsFactors = F)
dataDir <-"/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/data"
workDir <- "/data/SKDngs01/pg/collaborators/junke/PIMA/WGCNA_eQTL/run2"

setwd(dataDir)
target <- read_excel("SigTraits&targetWGCNA.xlsx")
##prepare snp counts fro raw eqtl table
eqtls <- read.table("b2glom.broad.emmax.cis_extract",header = T,sep = "\t",stringsAsFactors = F)
count=lapply(unique(eqtls$MARKER), function(x){length(eqtls[eqtls$MARKER==x,1])})
rawMarkercount <- data.frame(MARKER = unique(eqtls$MARKER),count=unlist(count))
#rawMarkercount <- read.csv("B2G.count.csv")

load("PI_G_Bx2-networkConstruction-signpwr7.RData")
result <- list()
mastertab <- eqtls[,c(4,6)]

for(color in unlist(unique(target[target$Target=="B2G",3]))){
  eqtlsub <- eqtls[eqtls$GENE %in% colnames(datExpr[,moduleColors==color]),]
  table <- data.frame(MARKER = character(),
                      a=integer(),
                      b=integer(),
                      c=integer(),
                      d=integer(),
                      PVal=double(),
                      OR=double())
  for(marker in unique(eqtlsub$MARKER)){
    a=length(eqtlsub[eqtlsub$MARKER==marker,1])
    b=nrow(eqtlsub)-a
    c=rawMarkercount[rawMarkercount$MARKER==marker,2]-a
    d=nrow(eqtls)-nrow(eqtlsub)-c
    pval <- fisher.test(matrix(c(a,b,c,d),nrow = 2,ncol = 2),alternative = "greater")
    table <- rbind(table,c(marker,a,b,c,d,pval$p.value,pval$estimate))
  }
  colnames(table)<-c("MARKER","a","b","c","d","P.value","OR")
  result[[paste(color)]] <- table
  colnames(table)[6]<- paste(color,"Pval",sep="_")
  mastertab <- merge(mastertab,table[,c(1,6)],by="MARKER",all.x = T)
}
setwd(workDir)
write.csv(mastertab,"B2G_SigLocus.csv",quote = F,row.names = F)
write.csv(rawMarkercount,"B2G.count.csv",row.names = F,col.names = T,quote = F)
save(rawMarkercount,result,mastertab,file = "B2G_Fisher.RData")

