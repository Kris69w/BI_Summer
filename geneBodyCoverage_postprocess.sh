dir="/data/SKDngs01/pg/collaborators/junke/SingleCell/PBMC_treat+ctrl" ;
file=${dir}/o_geneBodyCoverageAll.R
plot=${dir}/o_geneBodyCoverage.heatMap.pdf
> ${file}
ls -t ${dir}/plot/*.r |\

while read i;
do
	head -1 $i >>$file
done

/bin/cat <<EOM >> ${file}

cols <- ls()
data_matrix <- matrix(unlist(
  lapply(cols, get)
  ), byrow = F, nrow = 100)
colnames(data_matrix) <- cols
library(reshape2)
dm <- reshape2::melt(data_matrix)
data_matrix <- t(data_matrix)
rowLabel <- cols

library(grDevices)
myColors <- colorRampPalette(c("black", "yellow"))(256)

pdf(file = "${dir}/o_geneBodyCoverage.heatMap.pdf")
rc <- colorRampPalette(c("black", "yellow"))(ncol(data_matrix))

heatmap(data_matrix, scale = c("none"),keep.dendro = F, 
labRow  =  rowLabel ,Colv  =  NA,Rowv  =  NA,labCol = NA,
col = myColors,margins  =  c(6, 8),
ColSideColors  =  rc,
cexRow = 1,cexCol = 1,
xlab = "Gene body percentile (5'->3')", 
add.expr = x_axis_expr <- axis(side = 1,at = c(1,10,20,30,40,50,60,70,80,90,100),
labels = c("1","10","20","30","40","50","60","70","80","90","100")))
dev.off()

library(ggplot2)
library(ggplus)
pl <- ggplot(dm, aes(Var1, value)) + 
  geom_line(color = "navyblue") +
  xlab("Gene body percentile (5'->3')") +
  ylab("Coverage") +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "navyblue"),
        strip.text = element_text(color = "white"),
        legend.position = "none")
pdf(file =  "${dir}/o_geneBodyCoverage.profiles.pdf")
facet_multiple(plot = pl, 
               facets = 'Var2', 
               ncol = 2)
dev.off()
EOM

echo "Running script.."
Rscript ${file}
echo -e "\tdone"