# load libraries ----------------------------------------------------------

library(readxl)
library(xlsx)
library(R2HTML)
library(readr)
library(dplyr)
library(RSQLite)
library(VennDiagram)

# define Platform and initial path ----------------------------------------
if ( grepl(x = R.Version()["os"], pattern = "linux" )  ) {
  pathPrefix <- "/data/SKDngs01/pg/"
} else {
  pathPrefix <- "X:/pg/"
}

setwd(paste(pathPrefix , "Projects/Genetics/GWAS_ol_eQTL_PimaRadar", sep = ""))

# load reference databases ------------------------------------------------

GWAS <- read_delim(file = paste(pathPrefix , "database/GWAS_DB/GWAS_Significant.txt", sep = ""), 
                   delim = "\t", 
                   escape_double = FALSE, 
                   trim_ws = TRUE, 
                   progress = T)

# eQTLDB <- read_excel(paste(pathPrefix , "Projects/Genetics/eQTLDB.xlsx", sep = ""))
eQTLDB <- read_delim(file = paste(pathPrefix , "Projects/Genetics/eQTL_DB.txt", sep = ""),
                     delim = "\t",
                     quote = "",
                     trim_ws = T,
                     col_names = T,
                     progress = T)

MSSM_LoF_Full <- read_excel(path = "MSSM Biobank LoF Oct 11 2016.xlsx",
                            sheet = 1,
                            col_names = T)

MSSM_LoF = unique(MSSM_LoF_Full[,c(1,8:10)])

# # Create cpnn with grasp
# driver = dbDriver("SQLite")
# graspdb = dbConnect(driver, dbname = paste(pathPrefix , "database/grasp/grasp.db", sep = ""))
# test <- dbGetQuery(conn = graspdb,
#                    statement = "select * from grasp limit 10;")
# rm(test)

# SNP ---------------------------------------------------------------------

SNPsIsectGWAS <- intersect(unique(GWAS$SNP_ID),unique(eQTLDB$rsID))

# Gene Symbol -------------------------------------------------------------

MatchingSymbols <- intersect(unique(GWAS$Symbol),unique(eQTLDB$eGENE))


# Master list SNP section -------------------------------------------------

subsetSNP <- unique(eQTLDB[which(eQTLDB$rsID %in% SNPsIsectGWAS),])
subsetSNP <- unique(subsetSNP[, c("rsID","eGENE", "Source")]) # few records coming from radar as well

subsetGWAS <- unique(GWAS[which(GWAS$SNP_ID %in% SNPsIsectGWAS),c("SNP_ID","Symbol")])
# rename cols
colnames(subsetGWAS) <- c("rsID","GWAS_Symbol")
# write.table(x = subsetPima,file = "clipboard", row.names = F, col.names = T)
# write.table(x = subsetGWAS,file = "clipboard", row.names = F, col.names = T)
subsetCombined <- merge(x = subsetSNP,
                        y = subsetGWAS,
                        by = "rsID",
                        sort = F)


masterListSNP = subsetCombined %>%
  group_by(eGENE,GWAS_Symbol) %>%
  summarise_each(funs(toString(unique(.)))) 


# combine Gene names, they will be splitted later
masterListSNP$ProjectName <- apply(X = masterListSNP, MARGIN = 1, FUN = function(x) {
  res <- na.omit(unique(as.character(x[c("eGENE","GWAS_Symbol")])))
  paste(res, collapse = ";")
})
rm(subsetSNP,subsetGWAS,subsetCombined)

# Master list by gene symbol section ---------------------------------------------

subsetGENE <- unique(eQTLDB[which(eQTLDB$eGENE %in% MatchingSymbols),])
subsetGENE <- unique(subsetGENE[, c("eGENE", "Source")])

subsetCombined <- subsetGENE
subsetCombined$GWAS_Symbol = subsetGENE$eGENE
subsetCombined$rsID = NULL

masterListSymbol = subsetCombined %>%
  group_by(eGENE,GWAS_Symbol) %>%
  summarise_each(funs(toString(unique(.)))) 

masterListSymbol$ProjectName <- masterListSymbol$GWAS_Symbol

rm(subsetGENE,subsetCombined)

# Master list join --------------------------------------------------------
masterList = full_join(masterListSNP, masterListSymbol, by = colnames(masterListSymbol))
# simplify redundancies
masterList = masterList %>%
  group_by(ProjectName,eGENE,GWAS_Symbol) %>%
  summarise_each(funs(toString(unique(na.omit(.))))) 

masterList <- as.data.frame(masterList)


# Import RADAR dge --------------------------------------------------------

radarComps <- readRDS(paste(pathPrefix , "Projects/RADAR/DEG/DEGTableWithAllComps.rds", sep = ""))

colnames(radarComps) <- gsub(x = colnames(radarComps), pattern = "DEG_deg_", replacement = "")
colnames(radarComps) <- gsub(x = colnames(radarComps), pattern = "_DEG_", replacement = "_")
colnames(radarComps) <- gsub(x = colnames(radarComps), pattern = ".L_", replacement = "_")
colnames(radarComps) <- gsub(x = colnames(radarComps), pattern = "__", replacement = "_")
colnames(radarComps) <- gsub(x = colnames(radarComps), pattern = "_Interaction_", replacement = "_")

# Import PIMA dge ---------------------------------------------------------
# Pima
pimaComps <- readRDS(paste(pathPrefix , "Projects/Pima/expression/DEG/DEGTableWithAllComps.rds", sep = ""))

# Process dge -------------------------------------------------------------

# common <- intersect(rownames(pimaComps), as.character(radarComps$Symbol))
outList <- list()
outPVal <- list()
outRsq <- list()
#split multiple genes
for (i in 1:nrow(masterList)) {
  message("Processing project:", masterList$ProjectName[i])
  splittedGenes <- strsplit(x = masterList$ProjectName[i],split = ";", fixed = T)
  splittedGenes <- unlist(splittedGenes)
  for (k in splittedGenes) {
    message("\tsplitting gene:", k)
    # match comps in RADAR
    stats1 = radarComps[which(as.character(radarComps$Symbol) == k),3:ncol(radarComps)]
    stats2 = pimaComps[which(rownames(pimaComps) == k),2:ncol(pimaComps)]
    if (nrow(stats1) > 1 | nrow(stats2) > 1) {
      stop(paste(k, "has multiple matches", sep = ""))
    } 
    if (nrow(stats1) == 0) {
      stats1[2,] <- NA
    }
    if (nrow(stats2) == 0) {
      stats2[1,] <- NA
    }
    # append stats
    stats <- cbind(stats1,stats2)
    # message(length(as.character(c(masterList[i,],k, stats))))
    outList[[length(outList) + 1]] <- as.character(c(masterList[i,],k, stats))
    # filterOnly p-value
    outPVal[[length(outPVal) + 1]] <- as.character(c(masterList[i,],k, stats[,grep(x = colnames(stats), pattern = "adj.P.Val", ignore.case = T)]))
    outRsq[[length(outRsq) + 1]] <- as.character(c(masterList[i,],k, stats[,grep(x = colnames(stats), pattern = "rstat", ignore.case = T)]))
  }
}

# res <- matrix(data = unlist(outList),
#               nrow = length(outList),
#               ncol = (length(outList[[1]])), 
#               byrow = T)
# res = unique(res)


compCols <- c(colnames(radarComps)[3:ncol(radarComps)],
              colnames(pimaComps)[2:ncol(pimaComps)])

# colnames(res) <- c(colnames(masterList),"Symbol",compCols)

# make a table with pvalue only
resP <- matrix(data = unlist(outPVal),
               nrow = length(outPVal),
               ncol = (length(outPVal[[1]])), 
               byrow = T)
resP = unique(resP)

C1 <- compCols

colnames(resP) <- c(colnames(masterList),"Symbol",C1[grep(x = C1, pattern = "adj.P.Val", ignore.case = T)])

# make a table with R squared only
resRsq <- matrix(data = unlist(outRsq),
                 nrow = length(outRsq),
                 ncol = (length(outRsq[[1]])), 
                 byrow = T)
resRsq = unique(resRsq)

C2 <- compCols
colnames(resRsq) <- c(colnames(masterList),"Symbol",C2[grep(x = C2, pattern = "rstat", ignore.case = T)])


colnames(resP) <- gsub(colnames(resP), pattern = "_T_", replacement = "_Tubules_")
colnames(resP) <- gsub(colnames(resP), pattern = "_G_", replacement = "_Gloms_")

colnames(resRsq) <- gsub(colnames(resRsq), pattern = "_T_", replacement = "_Tubules_")
colnames(resRsq) <- gsub(colnames(resRsq), pattern = "_G_", replacement = "_Gloms_")


SigTubColl <- list()
SigGlomColl <- list()
assTubColl <- list()
assGlomColl <- list()

for (i in 1:nrow(resP)) {
  message("Processing gene:", resP[i,"Symbol"])
  vectorToTest <- resP[i, 7:ncol(resP)]
  vectorOfAss <- resRsq[i, 7:ncol(resRsq)]
  
  Tubules <- vectorToTest[grep(names(vectorToTest),pattern = "Tubules_",ignore.case = F)]
  Gloms <- vectorToTest[grep(names(vectorToTest),pattern = "Gloms_",ignore.case = F)]
  
  tubVecAss <- vectorOfAss[grep(names(vectorOfAss),pattern = "Tubules_",ignore.case = F)]
  glomVecAss <- vectorOfAss[grep(names(vectorOfAss),pattern = "Tubules_",ignore.case = F)]
  
  tidx <- which(as.numeric(Tubules) < 0.05)
  SigTubules <- paste(gsub(names(Tubules)[tidx],pattern = "_adj.P.Val", replacement = ""),collapse = "<BR>")
  assTubules <- paste(ifelse(as.numeric(tubVecAss[tidx]) > 0,"Positive","Negative"), collapse = "<BR>")
  
  gidx <- which(as.numeric(Gloms) < 0.05)
  SigGloms <- paste(gsub(names(Gloms)[gidx],pattern = "_adj.P.Val", replacement = ""),collapse = "<BR>")
  assGloms <- paste(ifelse(as.numeric(glomVecAss[gidx]) > 0,"Positive","Negative"), collapse = "<BR>")
  
  SigTubColl[[i]] <- SigTubules
  SigGlomColl[[i]] <- SigGloms
  
  assTubColl[[i]] <- assTubules
  assGlomColl[[i]] <- assGloms
}

Summary <- as.data.frame(resP[,1:6])

Summary$QueryGene = paste("<a href=https://www.ncbi.nlm.nih.gov/gene?term=",as.character(Summary$Symbol),"%5BSymbol%5D%20AND%209606%5Btaxid%5D&cmd=DetailsSearch>",as.character(Summary$Symbol),"</a>", sep = "")


Summary[,"SigTubules"] <- as.character(unlist(SigTubColl))
Summary[,"assTubules"] <- as.character(unlist(assTubColl))
Summary[,"SigGloms"] <- as.character(unlist(SigGlomColl))
Summary[,"assGloms"] <- as.character(unlist(assGlomColl))

# format eQTL results -----------------------------------------------------
Summary$Has_eQTLs <- sapply(as.character(Summary$Source),function(x) {
    vec <- unlist(strsplit(x = x, split = ", ", fixed = T))
    paste(unique(vec),collapse = "<BR>")
},simplify = T)

Summary$Source <- NULL
Summary = unique(Summary)
# # check additional GWAS ---------------------------------------------------
# sqlStat = 
#   paste("SELECT * FROM GRASP WHERE", 
#         paste(
#           paste("InGene = \'",Summary$Symbol,"\'", sep = ""
#           ), collapse = " OR "
#         ),
#         sep = " ")
# inGeneGrasp <- dbGetQuery(conn = graspdb,
#                           statement = sqlStat)
# sqlStat2 = 
#   paste("SELECT * FROM GRASP WHERE", 
#         paste(
#           paste("NearestGene = \'",Summary$Symbol,"\'", sep = ""
#           ), collapse = " OR "
#         ),
#         sep = " ")
# nearestGeneGrasp <- dbGetQuery(conn = graspdb,
#                                statement = sqlStat2)
# 
# 
# Summary$HasGeneOtherGWAS 


# Traits ------------------------------------------------------------------

# loop is the simplest
for (i in 1:nrow(Summary)) {
  currSym <- as.character(Summary$Symbol[i])
  traits <- as.character()
  PMIDs <- as.character()
  Has1GWASPaperFromUPenn <- character()
  if (currSym %in% GWAS$Symbol) {
    # get traits for the gene
    traits <- GWAS[which(GWAS$Symbol == currSym),]$TRAIT
    PMIDs <- GWAS[which(GWAS$Symbol == currSym),]$PMID
    Has1GWASPaperFromUPenn <- GWAS[which(GWAS$Symbol == currSym),]$InUPenn
  } else {
    for (rs in as.character(Summary[i,"rsID"])) {
      traits <- c(traits,GWAS[which(GWAS$SNP_ID == rs),]$TRAIT)
      PMIDs <- c(PMIDs,GWAS[which(GWAS$SNP_ID == rs),]$PMID)
      Has1GWASPaperFromUPenn <- c(Has1GWASPaperFromUPenn,GWAS[which(GWAS$SNP_ID == rs),]$InUPenn)
    }
  }
  # combine all the traits
  Summary$GWAS_Traits[i] <- ifelse(length(unique(traits)) == 0, "",paste(unique(traits),collapse = "<BR>"))
  Summary$GWAS_Traits_PMIDs[i] <- ifelse(test = length(unique(PMIDs)) == 0, 
                                         yes = "",
                                         no = paste(paste("<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/",
                                                          unique(PMIDs),"\">",unique(PMIDs) ,"</a>",
                                                          sep = "" ),
                                                    collapse = ";")
  )
  Summary$Has1GWASPaperFromUPenn[i] <- ifelse(any(Has1GWASPaperFromUPenn == "True"), "Yes","No")
}


# Sinai LoF data ----------------------------------------------------------

SummaryAnn <- merge(x = Summary,
                    y = MSSM_LoF,
                    by.x = "Symbol",
                    by.y = "gene",
                    all.x = T,
                    sort = F)
SummaryAnn = unique(SummaryAnn)
outFileName <- "eQTL_GWAS_Overlap_SummaryStats.html"
cat("<html><style>br {mso-data-placement:same-cell;}
    table.topalign td { vertical-align: top } </style>", file = outFileName)

HTML(x = SummaryAnn,
     file = outFileName,
     Border = 0,
     innerBorder = 1,
     align = "left",
     classtable = "topalign",
     append = T
)

# dbDisconnect(graspdb)
# dbUnloadDriver(driver)



# # Venn by SNP -------------------------------------------------------------
# 
# 
# futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# VennDiagram::venn.diagram(
#   list(
#     GWAS = unique(GWAS$SNP_ID),
#     PIMA = unique(c(eQTLDB$rsID[which(eQTLDB$Source == "Pima_B1_G")], 
#                     eQTLDB$rsID[which(eQTLDB$Source == "Pima_B2_G")],
#                     eQTLDB$rsID[which(eQTLDB$Source == "Pima_B1_T")],
#                     eQTLDB$rsID[which(eQTLDB$Source == "Pima_B2_T")])),
#     RADAR = unique(c(eQTLDB$rsID[which(eQTLDB$Source == "RADAR_T")], 
#                      eQTLDB$rsID[which(eQTLDB$Source == "RADAR_G")]))
#     
#   ),
#   fill = c("red", "green", "yellow"),
#   rotation.degree = 45,
#   main = "Mapping by SNP rsID",
#   na = "remove",
#   alpha = 0.5, 
#   force.unique = T,
#   filename = "GWAS_eQTL_combined_overlap_SNP.jpg")
# 
# 
# 
# 
# # Venn by gene symbol ------------------------------------------------
# 
# 
# 
# # Genes in at least 1 overlap
# futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# VennDiagram::venn.diagram(
#   list(
#     GWAS = unique(GWAS$Symbol),
#     PIMA = unique(c(eQTLDB$eGENE[which(eQTLDB$Source == "Pima_B1_G")], 
#                     eQTLDB$eGENE[which(eQTLDB$Source == "Pima_B2_G")],
#                     eQTLDB$eGENE[which(eQTLDB$Source == "Pima_B1_T")],
#                     eQTLDB$eGENE[which(eQTLDB$Source == "Pima_B2_T")])),
#     RADAR = unique(c(eQTLDB$eGENE[which(eQTLDB$Source == "RADAR_T")], 
#                      eQTLDB$eGENE[which(eQTLDB$Source == "RADAR_G")]))
#   ),
#   fill = c("red", "green", "yellow"),
#   rotation.degree = 45,
#   main = "Mapping by gene name",
#   na = "remove",
#   alpha = 0.5, 
#   force.unique = T,
#   filename = "GWAS_eQTL_combined_overlap.jpg")
# 
# 
