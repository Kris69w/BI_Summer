library(R2HTML)
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
options(stringsAsFactors = F)

setwd("/data/SKDngs01/pg/collaborators/junke/PIMA/GeneList/data")

## cleaning WGCNA data
sheets <- c("B1T_magenta","B1G_greenyellow","B2T_brown","B2G_magenta","B2G_salmon","B2G_turquoise","B1T_green",
            "B1T_turquoise","B1T_blue","B1T_magenta","B1G_greenyellow","B1G_grey60")

genes <- lapply(sheets, function(X) readxl::read_excel("WGCNA_pathAnalysis.xlsx", sheet = X)[1])
pathways <- lapply(sheets, function(X) readxl::read_excel("WGCNA_pathAnalysis.xlsx", sheet = X)[-1,c(5,6,8)])

names(genes)<- sheets
names(pathways)<-sheets

pathways <- lapply(pathways, na.exclude)

df <- data.frame(Gene=character(),
                 Tissue=character(),
                 Pathway=character(),
                 Pval=double(),
                 stringsAsFactors=FALSE)
for(module in sheets){
  for(i in 1:length(as.data.frame(pathways[paste(module)])[,3])){
    molegrp = as.data.frame(pathways[paste(module)])[i,3]
    moles <- strsplit(molegrp,",")
    moles <- intersect(unlist(moles), unlist(genes[paste(module)]))
    for(m in moles){
      tmp <- data.frame(Gene=paste(m),Tissue=substr(module,3,3),Pathway=as.data.frame(pathways[paste(module)])[i,1],
                        Pval=as.data.frame(pathways[paste(module)])[i,2],stringsAsFactors=FALSE)
      df <- rbind(df,tmp) 
    }
  }
}
GeneToPath <- unique(df)
GeneToPath$PValue <- as.numeric(GeneToPath$Pval)
GeneToPath <- GeneToPath[with(GeneToPath,order(Gene,Tissue,-PValue)),]
GeneToPath <- GeneToPath[,-4]

##read in all other files
eqtlraw <- read.table("eQTL_DB.txt",header = T,sep="\t",stringsAsFactors = F)
mastertab <- readxl::read_excel("SigGenemastertable.xlsx",sheet=1)
NASHraw <- read.table("LOF_Nash_Ftest.txt",header = T,sep="\t",stringsAsFactors = F)
NASHtab <- NASHraw[NASHraw$F_pval<0.1,][,c(1,4)]
DNraw <- readxl::read_excel("BI_DN_results.xlsx",sheet=1)
DNtab <- as.data.frame(DNraw[DNraw$snp_pval<0.1,][,c(2,10)])
GWAS_CKD <- read.table("GWAS_Significant.txt",header=T,sep="\t",stringsAsFactors = F)

##clean GWAS_NASH data
GWAS_NASH_raw <- readxl::read_excel("GWAS_NASH.xlsx",sheet="nash")[,c(13,53)]
GWAS_NASH_raw <- GWAS_NASH[complete.cases(GWAS_NASH),]
GWAS_NASH <- data.frame(Gene=character(),
                 Phenotype=character(),
                 stringsAsFactors=FALSE)
for(i in 1:length(GWAS_NASH_raw$InGene)){
  if(grepl(")(",GWAS_NASH_raw$InGene[i],fixed = T)){
    g <- unlist(strsplit(GWAS_NASH_raw$InGene[i],split = ")(",fixed = T))
    g[1] <- substr(g[1],2,nchar(g[1]))
    g[length(g)] <- substr(g[length(g)],1,nchar(g[length(g)])-1)
    for(j in 1:length(g)){
      tmp <- data.frame(Gene=g[j],Phenotype=GWAS_NASH_raw$Phenotype[i])
      GWAS_NASH <- rbind(GWAS_NASH,tmp)
    }
  }else{
    g <- GWAS_NASH_raw$InGene[i]
    g <- substr(g,2,nchar(g[length(g)])-1)
    tmp <- data.frame(Gene=g,Phenotype=GWAS_NASH_raw$Phenotype[i])
    GWAS_NASH <- rbind(GWAS_NASH,tmp)
  }
} 
GWAS_NASH<-GWAS_NASH[complete.cases(GWAS_NASH),]

for(i in  1:nrow(mastertab)){
  symbol <- mastertab[i,1]
  ## appending pathway
  pathall <- GeneToPath[GeneToPath$Gene==paste(symbol),][,c(2,3,4)]
  pathTwp <- pathall[grepl(x = pathall$Tissue,pattern = "T"),]
  pathGwp <- pathall[grepl(x = pathall$Tissue,pattern = "G"),]
  path <- pathTwp[with(pathTwp,order(-PValue)),][,2]
  mastertab$Pathway_T[i] <- ifelse(length(unique(path)) == 0, "",paste(unique(path),collapse = "<BR>"))
  path <- pathGwp[with(pathGwp,order(-PValue)),][,2]
  mastertab$Pathway_G[i] <- ifelse(length(unique(path)) == 0, "",paste(unique(path),collapse = "<BR>"))
  
  ##appending eqtl
  eqtl <- eqtlraw[eqtlraw$eGENE==paste(symbol),][,5]
  mastertab$eQTL[i] <- ifelse(length(unique(eqtl)) == 0, "",paste(unique(eqtl),collapse = "<BR>"))
  
  ##appending NASH
  nash <- NASHtab[NASHtab$gene==paste(symbol),][,2]
  mastertab$NASH_OR[i] <- ifelse(length(unique(nash)) == 0, "",paste(unique(nash),collapse = "<BR>"))
  
  ##appending DN
  dn <- DNtab[DNtab$snp==paste(symbol),][,2]
  mastertab$DN_OR[i] <- ifelse(length(unique(dn)) == 0, "",paste(unique(dn),collapse = "<BR>"))
  
  ##appending DE_in
  dgeall <- mastertab[i,c(2:9)]
  dgein <- colnames(dgeall[,!unlist(as.list(is.na(dgeall[1,])))])
  mastertab$DE_in[i] <- ifelse(length(unique(dgein)) == 0, "",paste(unique(dgein),collapse = "<BR>"))
  
  ##appending GWAS_CKD
  gwas_ckd <- GWAS_CKD[GWAS_CKD$SymbolHg19==paste(symbol),][,5]
  mastertab$GWAS_CKD[i] <- ifelse(length(unique(gwas_ckd)) == 0, "",paste(unique(gwas_ckd),collapse = "<BR>"))
  
  ##appending GWAS_NASH
  gwas_nash <- GWAS_NASH[GWAS_NASH$Gene==paste(symbol),][,2]
  mastertab$GWAS_NASH[i] <- ifelse(length(unique(gwas_nash)) == 0, "",paste(unique(gwas_nash),collapse = "<BR>"))

}


for(i in  1:nrow(mastertab)){
  mastertab$eigenTrait2[i] <- ifelse(grepl("eT_2",mastertab$DE_in[i]), "Yes","")
  mastertab$eigenTrait6[i] <- ifelse(grepl("eT_6",mastertab$DE_in[i]), "Yes","")
  mastertab$Tubule[i] <- ifelse(grepl("T.",mastertab$DE_in[i],fixed = T), "Yes","")
  mastertab$Glom[i] <- ifelse(grepl("G.",mastertab$DE_in[i],fixed = T), "Yes","")
}

##appending query
mastertab$QueryGene = paste("<a href=https://www.ncbi.nlm.nih.gov/gene?term=",as.character(mastertab$Symbol),
                            "%5BSymbol%5D%20AND%209606%5Btaxid%5D&cmd=DetailsSearch>",as.character(mastertab$Symbol),
                            "</a>", sep = "")

## make result readable
mastertab$DE_in = gsub(x = mastertab$DE_in,pattern = "eT",replacement = "eigentrait")
mastertab$DE_in = gsub(x = mastertab$DE_in,pattern = "T",replacement = "Tubule")
mastertab$DE_in = gsub(x = mastertab$DE_in,pattern = "G",replacement = "Glom")


##get final table summary
summary <- cbind(mastertab$Symbol,mastertab$QueryGene,mastertab$nSig_eT,mastertab$DE_in,mastertab$eigenTrait2,
                 mastertab$eigenTrait6,mastertab$Tubule,mastertab$Glom, mastertab$Pathway_T,mastertab$Pathway_G,
                 mastertab$eQTL,mastertab$NASH_OR,mastertab$DN_OR,mastertab$GWAS_CKD,mastertab$GWAS_NASH)
colnames(summary) <- c("Symbol","QueryGene","nSig_eT","DE_in","eigenTrait2",
                       "eigenTrait6","Tubule","Glom","Pathway_T","Pathway_G",
                       "eQTL","LoF_NASH","LoF_DN","GWAS_CKD","GWAS_NASH")
summary <- as.data.frame(summary)

outFileName <- "eigenTraitDEG_summary_GWASeQTLWGCNA.html"
cat("<html><style>br {mso-data-placement:same-cell;}
    table.topalign td { vertical-align: top } </style>", file = outFileName)

HTML(x = summary,
     file = outFileName,
     Border = 0,
     innerBorder = 1,
     align = "left",
     classtable = "topalign",
     append = T
)
HTMLEndFile()






