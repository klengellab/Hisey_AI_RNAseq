# WGCNA analysis on the juvenile CSDS for Erin Hisey

library(WGCNA)
library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(VennDiagram)
library(EnhancedVolcano)
library(topGO)
library(knitr)
library(kableExtra)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(limma)
library(edgeR)
library(EDASeq)
library(QCEWAS)
library(ggrepel)
library(msigdbr)
library(scater)
library(openxlsx)
library(parallel)
options(stringsAsFactors = FALSE)
set.seed(1234)

# Functions
source("/PHShome/je637/gitlab/general_functions/WGCNA_functions.R")
ensembl2gene <- function(ensembl_name){
  if(length(ensembl_name > 1)){
    gene <- list()
    for(i in ensembl_name){
      ens_id <- genes %>% dplyr::filter(ensembl_gene_id == i) %>%
        dplyr::select(external_gene_name) %>% pull()
      if(length(ens_id) == 0){
        gene <- append(gene, i)
      } else if(ens_id == ""){
        gene <- append(gene, i)
      } else if(length(ens_id != 0)){
        gene <- append(gene, ens_id)
      } else {
        gene <- append(gene, i)
      }
    }
  } else {
    ens_id <- genes %>% dplyr::filter(ensembl_gene_id == ensembl_name) %>%
      dplyr::select(external_gene_name) %>% pull()
    if(length(ens_id) == 0){
      gene <- append(gene, i)
    } else if(ens_id == ""){
      gene <- append(gene, i)
    } else if(length(ens_id != 0)){
      gene <- append(gene, ens_id)
    } else {
      gene <- append(gene, i)
    }
  }
  return(gene)
}
genes <- read.csv("/PHShome/je637/general/tables/ensembl_w_description.mouse.csv")

# Set Image path, output this is due to the different functions in WGCNA
image_path <- "/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/output/"
output <- "/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/output/"

# Load Data and Metadata 
norm.data <- readRDS("/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/workenvironment/norm_data.RData")
metadata <- read.csv("/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/data/metadata_hisey.csv", header = TRUE)
rownames(metadata) <- paste0("X", metadata$description)

# data needs to be in a format wherein the colnames are the genes and the rows the samples
data <- as.data.frame(t(norm.data))

# Signed or unsigned network
direction = "signed"
data <- data_preprocessing(data)

# clustering samples to detect outliers
data <- clust_outliers(data, filter == FALSE)

metadata$num_group <- as.numeric(as.factor(metadata$Group))
metadata$num_condition <- as.numeric(as.factor(metadata$Condition))
metadata$SI <- as.numeric(as.factor(metadata$SI))
traitData = metadata[,c(5:7)]
x <- match(rownames(data), rownames(traitData))
traitData <- traitData[x,]

# recluster samples with traitdata, make sure that traitData is numeric
cluster_samples(data, traitData)

# check the power_calc.pdf image to be sure for power calculation
power <- power_calc(data)

# Automatic module detection via dynamic tree cutting
modules <- WGCNA_module(Expr = data, power = power, direction = direction, split = 2)

table(modules$colors)
moduleLabels <- modules$colors
moduleColors = labels2colors(moduleLabels)
MEs <- modules$MEs;
geneTree = modules$dendrograms[[1]]

# To identify which genes are in which modules located
num <- c(0:21)
modulecolours_labels <- as.data.frame(moduleColors)
modulecolours_labels$label <- moduleLabels
rownames(modulecolours_labels) <- names(moduleLabels)
# check for each experiment if the order of columns to numbers is the same so for example for grey the number is 0 and for turquoise it's 1 etc.
col <- as.list(c("grey", "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink",
                       "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", 
                       "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred",
                       "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white"))
col <- col[1:22]
names(col) <- num

# for loop to create a list of genes present in each module. Notice that when you transfer the ensembl names to external gene names
# that some genes will be lost. If you want to overcome this also deliver the ensembl names next to the external gene names.

for(i in num){
  message(i)
  colors_module <- as.data.frame(modules[["unmergedColors"]])
  colnames(colors_module) <- "colors"
  colors_module$genes <- rownames(colors_module)
  x <- which(colors_module$colors == i)
  colors_module <- colors_module[x,]
  # transfer ensembl names to gene names
  gene_names <- as.character(ensembl2gene(colors_module$genes))
  y <- as.character(i)
  module <- as.character(col[[y]])
  if(i == 0){
    n <- module
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, gene_names)
    saveWorkbook(wb, file = paste0(output, "WGCNA_genes_per_module.xlsx"), overwrite = TRUE)
  } else{
    n <- module
    wb <- loadWorkbook(paste0(output, "WGCNA_genes_per_module.xlsx"))
    addWorksheet(wb, n)
    writeData(wb, n, gene_names, startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = paste0(output, "WGCNA_genes_per_module.xlsx"), overwrite = TRUE)
    }
}

# Plotting the dendrogram and module colors

# Metadata variables
names <- colnames(traitData)
colors <- list()
groups <- list()
for(i in names){
  x <- which(colnames(traitData) == i)
  x <- as.data.frame(traitData[,x])
  names(x) = i
  # next use this trait to define a gene significant variable
  x <- as.numeric(stats::cor(data, x, use = "p", method = c("pearson")))
  groups <- append(groups, list(x))
  color <- numbers2colors(x, signed = T)
  colors <- append(colors, list(color))
}
names(colors) <- names
names(groups) <- names

# plot the dendrogram and the module colors underneath 
# change each colors variable for every experiment and grouplabels!
blocknumber = 1
datColors = data.frame(moduleColors, colors[["num_group"]], colors[["num_condition"]], colors[["SI"]])[modules$blockGenes[[1]],]

pdf(paste0(image_path, "cluster_dendrogram.pdf"))
plotDendroAndColors(modules$dendrograms[[1]], colors = datColors,
                    groupLabels = c("Module colors", "ctrl vs defeated", "S vs NS", "SI"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Calculation eigengenes
ME1 <- eigengene_network(Expr = data, module = modules, metadata = traitData)

# Relate modules to phenotype data
MEs0 <- module_trait(data, traitData)

# Gene relationship to traint and important modules: Gene significance and module membership
## calculate module membership values (aka. module eigengene based connectivity kME)
datKME = signedKME(data, MEs0, corOptions = "use = 'p', method = 'pearson'")

# groups of interest/metadata
group = colnames(traitData)
intramodular_analyses(kme = datKME, group = group, variables = groups)
# The pvalue from the intramodular plots and the module trait relationships are not the same. 
# This is due to a difference since the intramodular connectivity is the correlation of the genes within a module
# whereas the module trait relationship is correlating the module with the metadata

## The yellow and pink module are interesting to investigate further

# Hub genes
hubgenes(Expr = data, moduleColors = moduleColors, power = power, direction = direction, MEs1 = ME1, output = output)

# Enrichment analysis
Label <- as.data.frame(moduleLabels)

results_enrichment <- enrichment_analysis(moduleLabels = moduleLabels, moduleColors = moduleColors, Expr = data, Labels1 = Label, output = output,
                                          universe = rownames(datKME))

col <- c("yellow", "pink")
for(i in col){
  message(i)
  for(database in names(results_enrichment[[i]])){
    data <- results_enrichment[[i]]
    if(database == "GO"){
      n <- data[[database]]
      wb <- createWorkbook()
      addWorksheet(wb, database)
      writeData(wb, database, n)
      saveWorkbook(wb, file = paste0(output, "enrichment_results_module", i, ".xlsx"), overwrite = TRUE)
    }else {
      n <- data[[database]]
      wb <- loadWorkbook(paste0(output, "enrichment_results_module", i, ".xlsx"))
      addWorksheet(wb, database)
      writeData(wb, database, n, startRow = 1, startCol = 1, colNames = TRUE)
      saveWorkbook(wb, file = paste0(output, "enrichment_results_module", i, ".xlsx"), overwrite = TRUE)
    }
  }
}
 

save.image("/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/workenvironment/WGCNA_analysis.RData")
