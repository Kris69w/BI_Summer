library(scater, quietly = TRUE)
library(knitr)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
require(scales)
options(stringsAsFactors = FALSE)
studyName <- "PBMC_treat+ctrl"
setwd(paste("/data/SKDngs01/pg/collaborators/junke/SingleCell/",studyName,"/", sep = ""))

#reads <- readr::read_delim(file = "RNASeqData.Count.txt", delim = "\t",progress = T)
reads <- readxl::read_excel(path = "RNASeqData.Count.xlsx", sheet = 1, col_names = T)
ph <- readr::read_tsv("sampleInfo.txt")
colnames(ph)[1] <- "SampleName"
ann <- reads[,c(1,(nrow(ph) + 2):ncol(reads))]
rmSample <- lapply(c("00"),function(x){paste0("MS00011_",x)})
ph <- ph[-which(ph$SampleName %in% rmSample),]

matReads <- as.matrix(reads[,ph$SampleName])
rownames(matReads) <- reads$GeneID

# read in aligment stats
al <- readr::read_delim(file = "AlignmentReport.txt",delim = "\t",
                        quote = "",col_names = T,
                        trim_ws = T)

ph$GROUPLABEL[1:48] <- "Control"
ph$GROUPLABEL[49:96] <- "Treated"
pd <- new("AnnotatedDataFrame", data = as.data.frame(ph))
rownames(pd) <- pd$SampleName

fd <- new("AnnotatedDataFrame", data = as.data.frame(ann))
rownames(fd) <- fd$GeneID


reads <- scater::newSCESet(
  countData = matReads,
  phenoData = pd,
  featureData = fd
)

keep_feature <- rowSums(counts(reads) > 0) > 1


#remove empty samples
keep_samples <- colSums(counts(reads) > 0 ) > 0
reads <- reads[keep_feature, ]

qc <- calculateQCMetrics(reads)
rmSample <- pData(qc)$total_features > 500
reads <- reads[,rmSample]
qc <- calculateQCMetrics(reads)
# scater_gui(qc)
sum(pData(reads)$GROUPLABEL=="Treated")



# This method plots the cumulative proportion of each cells library that is accounted 
# for by the top highest-expressed features (by default showing the cumulative proportion 
# across the top 500 features).This type of plot gives an overall idea of differences 
# in expression distributions for different cells. It is used in the same way as per-sample 
# boxplots are for microarray or bulk RNA-seq data. Due to the large numbers of zeroes in expression 
# values for single-cell RNA-seq data, boxplots are not as useful, so instead we focus 
# on the contributions from the most expressed features for each cell.
pdf("qc_cumulativePropOfLib.pdf", width = 10, height = 4)
p <- plot(qc, block1 = "GROUPLABEL",
     colour_by = "EffectiveAlignmentCount", 
     nfeatures = 1000, exprs_values = "counts")
print(p)
dev.off()

# Plot expression of markers
pdf("qc_markerGenes.pdf", width = 10, height = 4)
plotExpression(qc, c("GPR65","CD4", "CD3D"),
               x = "GROUPLABEL", 
               exprs_values = "counts")
dev.off()

# highest expression genes
pdf("qc_highestExpression.pdf", width = 4,height = 4)
plotQC(qc, type = "highest-expression", exprs_values = "counts")
dev.off()

# plot the frequency of expression (that is, number of cells with expression for the gene 
# above the defined threshold (default is zero)) against mean expression expression level .
pdf("qc_expressionVsFreq.pdf")
plotQC(qc, type = "exprs-freq-vs-mean")
dev.off()

pdf("qc_phenoData.pdf")
plotPhenoData(qc, aes(x = GROUPLABEL, y = total_features,
                                  colour = log10_total_counts))
dev.off()

pdf(file = "qc_tSNE.pdf")
plotTSNE(qc, colour_by = "GROUPLABEL", size_by="total_features")
dev.off()

p1 <- plotQC(qc, type = "find-pcs", variable = "total_features",
             plot_type = "pcs-vs-vars")
p2 <- plotQC(qc, type = "find-pcs", variable = "GROUPLABEL",
             plot_type = "pcs-vs-vars")
pdf("qc_PCvsVar.pdf")
multiplot(p1, p2, cols = 2)
dev.off()



al_t <- t(al[,-1])
colnames(al_t) <- al$Pair
# al_t <- type_convert(tibble::as.tibble(al_t))
al_t <- type_convert(tibble::as_tibble(al_t))
al_t$SampleName <- colnames(al)[-1]

qc_reads <- pData(qc)
qc_reads <- qc_reads %>% inner_join(al_t, by = c("SampleName" = "SampleName"))
qc_reads_m1 <- reshape2::melt(qc_reads[,c("SampleName","Total read#","total_counts")])
qc_reads_m2 <- reshape2::melt(qc_reads[,c("SampleName","Uniquely mapped paired read#","Non-uniquely mapped paired read#")])
colnames(qc_reads_m2)[2] <- "MapType"
qc_reads_m2$variable <- "Mapped reads"
qc_reads_m1$MapType <- ""
qc_reads_m <- rbind(qc_reads_m1[,colnames(qc_reads_m1)], qc_reads_m2[,colnames(qc_reads_m1)])
qc_reads_m <- inner_join(x = qc_reads_m, y = qc_reads[,c("SampleName","total_features")], by = c("SampleName" = "SampleName"))
qc_reads_m$variable <- plyr::revalue(qc_reads_m$variable,c("Total read#" = "Total Reads",
                                                           "Mapped reads" = "Mapped Reads",
                                                           "total_counts" = "Exonic Reads"), warn_missing = T)
qc_reads_m$variable <- factor(qc_reads_m$variable, levels = levels(qc_reads_m$variable)[c(1,3,2)])
qc_reads_m[qc_reads_m$variable == "Total Reads",]$MapType = "Total"
qc_reads_m[qc_reads_m$variable == "Exonic Reads",]$MapType = "Exonic"
reads_sums <- qc_reads_m %>% group_by(variable) %>% summarise(., sum = sum(value))


pdf(file = "qc_readsPerSample.pdf", height = 8, width = 15)

p1 <- ggplot(qc_reads_m, aes(SampleName, value, fill = MapType, group = variable)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 1,color = "black") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = c(0,0), labels = comma) +
  ylab("Number of reads") +
  xlab("Name of the sample") +
  labs(title = studyName, fill = "Type of reads:") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")

p2 <- ggplot(qc_reads_m, aes(SampleName, "Expressed genes")) +
  geom_raster(aes(fill = total_features)) + 
  scale_fill_gradient(low = "black", high = "yellow") +
  theme_minimal() + 
  xlab("Name of the sample") +
  labs(
    fill = "Expressed genes",
    caption = paste(
         paste("Total Reads:",round(reads_sums[1,2], digits = 0), sep = ""),
         paste("Total Mapped:",round(reads_sums[2,2], digits = 0), sep = ""),
         paste("Total Exonic:",round(reads_sums[3,2], digits = 0), sep = ""),
         sep = "; ")) +
  # theme(axis.text.x = element_text(angle = 90, vjust = .5 )) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5))  

plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(7,3) )

dev.off()


pdf(file = "qc_genesExpressedPerCell.pdf")
ggplot(pData(qc), aes(SampleName, total_features,  fill = total_features)) +
  geom_col(colour = "black") +
  labs(title = studyName, subtitle = "Genes expressed per cell", 
       caption = paste("Dataset total counts:",round(sum(pData(qc)$total_count), digits = 0))) +
  scale_fill_gradient(low = "black", high = "yellow") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        legend.position = "none")
dev.off()


qc_reads_sort <- qc_reads_m[qc_reads_m$MapType=="Total",]
pdf(file = "qc_genesExpressedPerCell_sorted.pdf",width = 12,height = 10)
ggplot(qc_reads_sort, aes(reorder(SampleName,value), total_features, fill =value)) +
  geom_col(colour = "black") +
  geom_text(aes(label=value),size=3,angle=90,hjust=0.001)+
  labs(title = studyName, subtitle = "Genes expressed per cell", 
       caption = paste("Dataset total counts:",round(sum(pData(qc)$total_count), digits = 0)),
       x="Sample Name",y="Number of detected Genes") +
  scale_fill_gradient(name="Total read#",low = "black", high = "yellow") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.key.width = unit(2,"cm"),
        axis.text.x = element_text(angle = 90, vjust = .5))
dev.off()

