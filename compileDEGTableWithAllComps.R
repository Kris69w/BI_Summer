# Set R terminal widht to actual terminal size
#options(width=Sys.getenv("COLUMNS"))
# Make all repositories available for package download:
setRepositories(ind = c(1,2,3,4,5,6,7,8,9))

# load libraries

# define Platform and initial path ----------------------------
if ( grepl(x = R.Version()["os"], pattern = "linux" )  ) {
  pathPrefix <- "/data/SKDngs01/pg/"
} else {
  pathPrefix <- "X:/pg/"
}


# Set your working directory where output files will be saved ---------
workingDir <- paste(pathPrefix , "collaborators/junke/eigenTrait_DEG/Raw", sep = "");
#
setwd(dir = workingDir)

# read in file names 
# b1_tub_comps <- dir(pattern = "_T_*", all.files = F,full.names = T,include.dirs = F, recursive = F)
# b2_tub_comps <- dir(pattern = "_T_*", all.files = F,full.names = T,include.dirs = F, recursive = F)
# glom_comps <- dir(pattern = "*.txt", all.files = F,full.names = T,include.dirs = F, recursive = T)
#combined_comps <- c(tub_comps,glom_comps)
combined_comps <- dir(pattern = "deg_biopsy*", all.files = F,full.names = T,include.dirs = F, recursive = F)
# collection of resulting outputs from analysis
res <- list()
colsToExtract <- c("adj.P.Val","rstat")
# isFirst <- TRUE
genesOrder <- character()

message("Preprocessing comparisons...")
for (i in combined_comps) {
  message("Reading file: ",i)
  if (grepl(pattern = "dge_Summary.txt", x = i, fixed = T)) {
    message("Skipping")
    next;
  }
  output <- read.table(file = i,header = T,sep = "\t",quote = "",as.is = T)
  rownames(output) <- output[,"geneName"]
  # parse file name
  i <- gsub(pattern = ".txt",replacement = "", x = i, fixed = T)
  i <- gsub(pattern = "./deg_",replacement = "", x = i, fixed = T)
  i <- gsub(pattern = "biopsy_",replacement = "B_", x = i, fixed = T)
  i <- gsub(pattern = "_X.",replacement = "_", x = i, fixed = T)
  message("Name parsed as: ",i)
  # if (isFirst) {
  #   genesOrder <- setNames(object = as.character(output[,"geneName"]), nm = as.character(row.names(output)))
  #   isFirst <- FALSE 
  # }
  # res[[i]] <- output[names(genesOrder), colsToExtract]
  genesOrder <- unique(c(genesOrder,as.character(output[,"geneName"])))
  res[[i]] <- output[, colsToExtract]
}

for (i in names(res)) {
  message("Resorting: ", i)
  res[[i]] <- res[[i]][genesOrder,]
}

FirstHeader <- rep(names(res), each = 2)
SecondHeader <- rep(colsToExtract, length(res))
FinalHeader <- paste(FirstHeader,SecondHeader, sep = "_")
resDF <- as.data.frame(matrix(unlist(res), nrow = length(genesOrder), ncol = length(res)*length(colsToExtract)))
colnames(resDF) <- FinalHeader
FinalDF <- cbind(Symbol = genesOrder,
                 resDF)
write.table(x = FinalDF,
            file = "DEGTableWithAllComps.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
saveRDS(object = FinalDF, compress = T, file = "DEGTableWithAllComps.rds")







