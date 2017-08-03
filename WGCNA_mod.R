library(WGCNA);
library(flashClust);
library(readxl)
enableWGCNAThreads()

##############################DATA Input#############################

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the expression data set
pimab1Data = readRDS("/data/SKDngs01/pg/Projects/Pima/expression/array/biopsy2.rds");

##prepare the required data format
Tlist<-c();Glist <- c();
for(i in names(pimab1Data)){
  if(grepl("*T$",i)){
    Tlist<-c(Tlist,i)
  }else{
    Glist <- c(Glist,i)
  }
}

pimab1Data <- as.data.frame(t(pimab1Data))
datExpr0 <- pimab1Data[Tlist,] ##b1T expression data

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(15,9)
# pdf(file = "expr_b1T_avg_sampleClustering.pdf", width = 12, height = 9);
# par(cex = 0.6);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "expr_b1T_avg_Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#      cex.axis = 1.5, cex.main = 2)
# 
# dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 125, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr= datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##########################Network#######################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sftsign = pickSoftThreshold(datExpr, networkType = "signed",powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftsign$fitIndices[,1], -sign(sftsign$fitIndices[,3])*sftsign$fitIndices[,2],
     xlab="Soft Threshold (power)sign",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sftsign$fitIndices[,1], -sign(sftsign$fitIndices[,3])*sftsign$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="blue")
abline(h=0.95,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sftsign$fitIndices[,1],sftsign$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sftsign$fitIndices[,1], sftsign$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 7;
adjacency = adjacency(datExpr,type = "signed",  power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene Averageclustering on TOM-based dissimilarity_signpwr7",
     labels = FALSE, hang = 0.04);

### Cut the tree	 
minModuleSize = 80;
mColorh=NULL
# Module identification using dynamic tree cut:
dynamicMods4 = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit =4 , pamRespectsDendro = FALSE,cutHeight = 0.98,
                minClusterSize = minModuleSize);
table(dynamicMods4)

dynamicMods=dynamicMods4
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

### Merge modules
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "geneDendro_merge_PI_B2T_signpwr7.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

module<-data.frame(gene=colnames(datExpr),module=mergedColors)
write.table(module,file = "PI_B2T-genemodule_signpwr7.txt",quote=F,sep="\t",row.names=F)

save(MEs, datExpr,moduleLabels, moduleColors, geneTree, file = "PI_T_Bx2-networkConstruction-signpwr7.RData")
