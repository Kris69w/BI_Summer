library(WGCNA)
options(stringsAsFactors = F)

##load in previous session
lnames = load(file = "PI_G_Bx2-networkConstruction-signpwr7.RData");

##load in trait data
df <- read.table("/data/SKDngs01/pg/collaborators/junke/PIMA/data/eigenTrait_JSD_formatted.txt",header = T,sep = "\t")
df$Target <- paste("b",df$BiopsyNumber,df$Tissue,sep="_")

numT <- grep("eigen*",colnames(df))
datTrait0 <- df[,c(1,numT)]
pimaSamples = rownames(datExpr);
traitRows = match(pimaSamples, datTrait0$NIH);
datTraits = datTrait0[traitRows, -1];
rownames(datTraits) = datTrait0[traitRows, 1];
collectGarbage();

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
#MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#MEs = orderMEs(MEs0)
# names(MEs)
# names(MEs0)
# write.csv(MEs, file = "PI_TA_Bx1_samplesMEs_signpwr12.csv")

moduleTraitCor = round(cor(MEs, datTraits, use = "p"),2);
moduleTraitPvalue = round(corPvalueStudent(moduleTraitCor, nSamples),2);

# write.csv(moduleTraitCor , file = "moduleTraitCorPI_TA_Bx1samples_signpwr12.csv")
# write.csv(moduleTraitPvalue , file = "moduleTraitPvaluePI_TA_Bx1samples_signpwr12.csv")

allgeneshub= chooseTopHubInEachModule(datExpr, moduleColors,omitColors = FALSE)
moduleallgeneshub<-data.frame(allgeneshub)
write.table(moduleallgeneshub,file = "PI_TA_Bx1_hubgenes_signpwr12.txt",quote=F,sep="\t",row.names=T)

get_stars = function(p) {
  stars = findInterval(p, c(0, 0.001, 0.01, 0.05,0.1))
  codes = c("****" , "***","**", "*", "")
  codes[stars]
}

# sizeGrWindow(15,8)
# # Will display correlations and their p-values
# ##textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
# signif(moduleTraitPvalue, 1), ")", sep = "");
# 
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n", get_stars(moduleTraitPvalue));	
# yColorWidth =0.03,
# xColorWidth =1.5,
# 

# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.8, 3, 2.2)); 
# 
# #        ySymbols = names(MEs[c(1,6,10,11)]),
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.58,
#                cex.axis = 0.75,
#                cex.lab=0.75,
#                cex.lab.x = 0.6,
#                zlim = c(-1,1),
#                main = paste("Module Association"))		 
# 
# 

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# names (colors) of the modules
modNames = substring(names(MEs), 3)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# Define variable  containing the  column of datTrait
eigen6 = as.data.frame(datTraits$eigenTrait6);
names(eigen6) = "EigenTrait6"
eigen2 = as.data.frame(datTraits$eigenTrait2);
names(eigen2) = "EigenTrait2"



##find genes with high connectivity and high gene significance
ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)

##for trait 6
GS1 = as.data.frame(cor(eigen6,datExpr, use = "p"));
GeneSignificance=abs(GS1)

##for trait 6
colorlevels=c("green","turquoise","blue","magenta")
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
FilterGenes= abs(GS1)> .6 & abs(datKME$MM.magenta)>.8
colnames(FilterGenes)[which(FilterGenes==TRUE)]

# [1] "C12orf29" "CAPN7"    "DLAT"     "DLD"      "EBAG9"    "FBXO3"    "LARP4"    "MRS2"     "RABGGTB"  "RNF14"    "RPRD1A"  
# [12] "TNFRSF25" "TPD52"  

FilterGenes= abs(GS1)> .6 & abs(datKME$MM.blue)>.8
colnames(FilterGenes)[which(FilterGenes==TRUE)]

# [1] "ACSF2"    "ACY1"     "AKR1A1"   "ALDH1L1"  "ALDH4A1"  "ALDH6A1"  "ALOX5"    "AZGP1"    "C2orf47"  "CCL20"    "COQ9"    
# [12] "CRADD"    "CRYL1"    "CXCL1"    "DCXR"     "ECHDC2"   "ECI2"     "ENPP2"    "EPHX2"    "FAH"      "FBP1"     "FHL1"    
# [23] "FOLR1"    "FSTL1"    "GHITM"    "GLYAT"    "GPD1"     "GPR137B"  "HAGH"     "HMGCL"    "HPD"      "HS3ST1"   "MAOB"    
# [34] "MIPEP"    "PC"       "PCK2"     "PEPD"     "PINK1"    "POLDIP2"  "SLC13A3"  "SLC22A6"  "SLC22A8"  "SLC25A44" "SLC7A9"  
# [45] "SUGCT"    "TMBIM6"   "TNFAIP6"  "UPB1"     "VCAN"   

FilterGenes= abs(GS1)> .6 & abs(datKME$MM.green)>.8
colnames(FilterGenes)[which(FilterGenes==TRUE)]

# [1] "FSTL1"  "GPD1"   "HS3ST1" "MAP7"   "PMEPA1" "TLE1"  

FilterGenes= abs(GS1)> .6 & abs(datKME$MM.turquoise)>.8
colnames(FilterGenes)[which(FilterGenes==TRUE)]

# [1] "ALOX5"   "CALHM2"  "CD1C"    "CRADD"   "CRYL1"   "EFHD1"   "FLI1"    "GHITM"   "GPR137B" "IGSF6"   "PCCA"    "PIK3CG" 

##for trait 2
GS1 = as.data.frame(cor(eigen2,datExpr, use = "p"));
GeneSignificance=abs(GS1)

##for trait 2
colorlevels=c("magenta")
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

FilterGenes= abs(GS1)> .2 & abs(datKME$MM.magenta)>.8
colnames(FilterGenes)[which(FilterGenes==TRUE)]
#############################################################################################
##    "DLAT"     "EBAG9"    "FBXO3"    "LARP4"    "RABGGTB"  "RNF14"    "TNFRSF25" "TPD52" ##
#############################################################################################

