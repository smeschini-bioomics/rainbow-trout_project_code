library(QFeatures)
library(msqrob2)
library(tidyverse)
library(gt)
library(plotly)
library(ggplot2)
library(stageR)
library(poolr)
library(RColorBrewer)
library(seqinr)
library(stringr)
library(ExploreModelMatrix)

# import data
df <- df <- read.csv("01_input_files/20260420_TroutGut_peptides_ions_raw_ab_rank3_msqrob-input.csv",
                     header=TRUE,
                     skip=0)
ecols = grep("240816", colnames(df))
#order dataframe by protein, for the normalisation step
df = df[order(df$Protein),]
df$peptidoform <- paste(df$Sequence,
                        df$Variable.modifications...position..description.,
                        sep = "_")
pe <- readQFeatures(df,quantCols= ecols, name="peptidoformRaw", fnames = "peptidoform")

## Preprocessing
print(nrow(pe[["peptidoformRaw"]]))
#Get number of columns values per row that are larger than zero 
rowData(pe[["peptidoformRaw"]])$nNonZero <- rowSums(assay(pe[["peptidoformRaw"]]) > 0)
#Change the zeroes into NA values, because they represent missing values, not values of zero intensity
pe <- zeroIsNA(pe, i = "peptidoformRaw")

# inspect missigness
mat <- assay(pe[["peptidoformRaw"]])
miss_feature <- rowMeans(is.na(mat)) * 100  # percentage
# total number of features
length(miss_feature)
plot(
  seq_along(miss_feature),
  miss_feature,
  pch = 16, cex = 0.4,
  xlab = "Feature index",
  ylab = "Missingness (%)",
  main = "Feature-level missingness"
)

abline(h = 50, col = "red", lty = 2)
length(miss_feature)


#Filter out peptides with too many missing values
n_samples <- ncol(assay(pe[["peptidoformRaw"]]))
min_prop <- 0.5
min_nonzero <- ceiling(min_prop * n_samples)
pe <- pe[rowData(pe[["peptidoformRaw"]])$nNonZero >= min_nonzero, , ]
#Log transform the intensity values
pe <- logTransform(pe, i = "peptidoformRaw", base = 2, name = "peptidoformLog")
pe <- normalize(pe, i = "peptidoformLog", method = "center.median", name = "peptidoform")
print(nrow(pe[["peptidoform"]]))

master_proteins <- c("H3", "H4")
rowData(pe[["peptidoform"]])$Master_Protein <- case_when(
  grepl(pattern = "H3", rowData(pe[["peptidoform"]])$Protein) ~ "H3",
  grepl(pattern = "H4", rowData(pe[["peptidoform"]])$Protein) ~ "H4"
)

## Normalisation via robust summarisation
#list to store summarised values per protein
pe_robustval <- list()

pe <- aggregateFeatures(pe,
                        i = "peptidoform",
                        fcol = "Master_Protein",
                        na.rm = TRUE,
                        maxit = 100,
                        name = "proteinRobust",
                        fun = MsCoreUtils::robustSummary)

#Fill out list with the summarised values
for (prot in unique(rowData(pe)[["peptidoform"]]$Master_Protein)){
  print(prot)
  pe_robustval[[prot]] = assay(pe[["proteinRobust"]])[prot,]
}


y <- as_tibble(assay(pe[["peptidoform"]]))
y_new <- tibble()
rowdata <- tibble()
for (prot in unique(rowData(pe)[["peptidoform"]]$Master_Protein)){
  pe_sub <- pe[["peptidoform"]][grepl(prot, rowData(pe[["peptidoform"]])$Master_Protein, fixed = T),]
  y_ <- as_tibble(assay(pe_sub))
  #center assay based on the corresponding protein value
  y_scale <- base::scale(y_, center = pe_robustval[[prot]], scale = FALSE) 
  rownames(y_scale) <- rownames(assay(pe_sub))
  y_new <- rbind(y_new, y_scale)
  rowdata <- rbind(rowdata, as.data.frame(rowData(pe_sub)))
}


#Add the normalized assay as a new assay to the existing pe
y_assay <- SummarizedExperiment(assays=as.matrix(y_new), rowData=rowdata, colData=colData(pe[["peptidoform"]]))
#Filter out peptidoforms that now have 0 intensities everywhere
rowData(y_assay)$nNonZero2 <- rowSums(assay(y_assay)!=0, na.rm = T)
y_assay <- y_assay[rowData(y_assay)$nNonZero2>0,]
pe <- addAssay(pe, y_assay, name = "peptidoformNorm2")
pe <- QFeatures::normalize(pe, method = "center.median", i = "peptidoformNorm2",
                           name = "peptidoformNorm")
rowData(pe[["peptidoformNorm"]])$peptidoform <- rownames(pe[["peptidoformNorm"]])


pe <- zeroIsNA(pe, "peptidoformNorm")
limma::plotDensities(assay(pe[["peptidoformNorm"]]),
                     legend = FALSE)
                     #col = as.numeric(colData(pe)$diet))
boxplot(assay(pe[["peptidoformNorm"]]),
        main = "Peptide distributions after Normalisation", ylab = "intensity")



## ptm summarisation
#In this step we will summarise the peptidoform level data to ptm level data.
### Get location of modification in protein
# For this we need the location information of each modification in the protein.
### Get ptm level intensity matrix
# In our case, a ptm is a unique protein - modification(+location) combination.

#Add ptm variable = protein + modification
rowData(pe[["peptidoformNorm"]])$ptm <- ifelse(rowData(pe[["peptidoformNorm"]])$"Variable.modifications...position..description." != "",
                                               paste(rowData(pe[["peptidoformNorm"]])$Master_Protein,
                                                     rowData(pe[["peptidoformNorm"]])$"Variable.modifications...position..description.",
                                                     sep="_"),
                                               "")
prots <- unique(rowData(pe[["peptidoformNorm"]])$Master_Protein)
#Do for each protein
ptms <- sapply(prots, function(i) {
  pe_sub <- pe[["peptidoformNorm"]][grepl(i, rowData(pe[["peptidoformNorm"]])$Master_Protein, fixed = T),]
  #Get all unique modifications present on that protein
  mods <- unique(unlist(strsplit(rowData(pe_sub)$"Variable.modifications...position..description.", split = "|", fixed = TRUE)))
  #Add protein info to mods
  ptm <- paste(rep(i, length(mods)), mods)
  #return all the protein-mods combinations
  ptm
})
ptms <- as.vector(unlist(ptms))



# For every unique ptm, we take all its associated peptidoforms and aggregate all its intensity values per sample into one intensity value per sample for that ptm.
#For each ptm do
ptm_x_assay <- sapply(seq(1:length(ptms)), function(i){ 
  x <- ptms[i]
  #Get current protein and mod from ptm
  prot <- str_split(x, " ", 2)[[1]][1]
  current_ptm <- str_split(x, " ", 2)[[1]][2]
  #filter on that protein and on that mod to obtain all peptidoforms that correspond to the ptm
  pe_sub <- pe[["peptidoformNorm"]][grepl(prot, rowData(pe[["peptidoformNorm"]])$Master_Protein, fixed = T),]
  ptm_sub <- pe_sub[grepl(current_ptm, rowData(pe_sub)$"Variable.modifications...position..description.", fixed = T),]
  #Get intensity values of those peptidoforms
  y <- assay(ptm_sub)
  #And summarise them to 1 row of intensity values: 1 value per sample for that ptm
  if (any(is.finite(y))){
    ptm_y <- MsCoreUtils::robustSummary(y)
  }
  else {ptm_y <- rep(NA, ncol(y))}
  ptm_y
})
#Then we get the intensity assay on ptm level
ptm_x_assay <- t(ptm_x_assay)
rownames(ptm_x_assay) <- ptms
colnames(ptm_x_assay) <- colnames(assay(pe[["peptidoformNorm"]]))



### Make new QFeatures object
# From this assay, we make a new QFeatures object and fill out its rowData and colData
ptm <- readQFeatures(as.data.frame(ptm_x_assay),
                     ecol= 1:ncol(ptm_x_assay),
                     name="ptm")
ptm <- renamePrimary(ptm, rownames(colData(pe)))
colData(ptm) <- colData(pe)
ptm <- renamePrimary(ptm, rownames(colData(pe)))
rownames(ptm[["ptm"]]) <- rownames(ptm_x_assay)


rowData(ptm[["ptm"]])$Protein <- sapply(str_split(rownames(ptm[["ptm"]]),
                                                  pattern=" ", n=2),
                                        function(x) x[1])
rowData(ptm[["ptm"]])$modification <- sapply(str_split(rownames(ptm[["ptm"]]),
                                                       pattern=" ", n=2),
                                             function(x) x[2])

# Output directory
OUTDIR <- "01_output_files"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Output tables
write.csv(
  assay(pe[["peptidoformNorm"]]),
  file = file.path(OUTDIR, "TroutGut_peptidoform_H3H4.csv")
)

write.csv(
  ptm_x_assay,
  file = file.path(OUTDIR, "TroutGut_hPTM_H3H4.csv")
)

# Reproducibility record
writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTDIR, "sessionInfo.txt")
)
