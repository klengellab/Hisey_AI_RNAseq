# RNAseq dataset from male juvenile CSDS mice at the insula from Erin Hisey

#libraries

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
library(patchwork)
options(stringsAsFactors = FALSE)
set.seed(1234)

# functions:
ggPCA <- function(pca, metadata, variables){
  PCA_out <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(colnames(PCA_out), " (", percentage, "%", ")")
  theme <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 strip.background = element_blank(), 
                 axis.text.x = element_text(colour = "black"), 
                 axis.text.y = element_text(colour = "black"),
                 axis.ticks = element_line(colour = "black"), 
                 plot.margin = unit(c(1, 1, 1, 1), "line"))
  for(i in variables){
    p <- ggplot(PCA_out, aes(x = PC1, y = PC2, color = metadata[, i])) + 
      geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2]) + 
      labs(color = i)
    print(p)
  }
}
source("/PHShome/je637/gitlab/general_functions/plot/customPCA.R")
biomart_ensembl_to_gene <- function(gene){
  ensembl=useMart("ensembl")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = gene, 
                 mart = ensembl)
  return(names)
}
gene_boxplot_V2 <- function(count, gene, group){
  gene_count <- as.data.frame(t(count[gene,]))
  gene_count <- cbind(gene_count, group)
  colnames(gene_count) <- c("c","g")
  gene_count <- as.data.frame(gene_count)
  #gene_name <- biomart_ensembl_to_gene(gene)
  ggplot(gene_count, aes(x = factor(g), y = c)) + 
    geom_boxplot()+
    geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.5) + 
    geom_text(label = rownames(gene_count), nudge_x = 0.25) +
    labs(x = "Group", y = "Log2(normalized count)", title = gene_name) + 
    scale_x_discrete(labels= c("control", "defeated")) + theme_classic()
}
ensembl2gene <- function(genes){
  ensembl=useMart("ensembl", host = "https://useast.ensembl.org/")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = genes, 
                 mart = ensembl)
  return(names)
}

# output
output <- "/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/output/"

# data
data <- read.csv("/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/data/tximport-counts.csv",
                 header = TRUE)
metadata <- read.csv("/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/data/metadata_hisey.csv", header = TRUE)

# before filtering we have 54532 genes
rownames(data) <- data$gene
data <- data[,-1]

# PCA analysis before normalisation to check for unknown co-variates
seq_depth <- colSums(data)
metadata$seq_depth <- seq_depth

PCA1 <- as.matrix(t(data))
PCA1 <- log2(PCA1 + 1)
PCA1 <- prcomp(PCA1)

ggPCA(PCA1, metadata, colnames(metadata)[3:6])

# custom PCA plot before normalisation
pd <- metadata[,c(3:6)]
log_data <- log2(data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA before normalisation")
# There are signifiant covariables in condition and Group but the majority is explained by sequencing depth
# this should be gone after normalisation

# Filtering counts that a gene needs a read with 10 or more counts
keep <- apply(data,1,function(x){any(x >=10)})
filtered_data <- data[keep, ]
filtered_data <- as.matrix(filtered_data)
dim(filtered_data) # 19443 28
log_data <- log2(filtered_data + 1)
## there are 32649 genes filtered out
barplot(colSums(log_data), main = "Total reads", xlab = "sample", 
        ylab = "log2 read counts")

# Differential expression analysis with Deseq2
#-----
coldata <- metadata[,c(3:6)]
rownames(coldata) <- metadata$description
rownames(coldata) <- paste0("X", rownames(coldata))
all(rownames(coldata) == colnames(filtered_data)) # FALSE
# need to match to ensure that the rows and columns have the same order
x <- match(rownames(coldata), colnames(filtered_data))
filtered_data <- filtered_data[,x]
all(rownames(coldata) == colnames(filtered_data)) # TRUE
design <- as.data.frame(metadata[,3])
rownames(design) <- rownames(metadata)
coldata$Condition <- as.factor(coldata$Condition)
coldata$Group <- as.factor(coldata$Group)
dds <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata, design = ~ Group)
dds$condition <- factor(dds$Group, levels = c("control", "defeat"))                       
dds <- DESeq(dds)
resultsNames(dds)


# Check if data is normal distributes
# See if every sample has comparable counts
norm.data <-counts(dds, normalized=TRUE)
log_norm <- as.data.frame(log2(norm.data + 1))
barplot(colSums(log_norm), main = "Normalized Total reads", xlab = "sample", ylab = "log read counts")
saveRDS(norm.data, "/PHShome/je637/RNAseq/RNAseq_Hisey_CSDS/workenvironment/norm_data.RData")


# PCA analysis after normalisation to check for unknown co-variates
seq_depth_norm <- colSums(norm.data)
PCA1 <- as.matrix(t(norm.data))
PCA1 <- prcomp(PCA1)

ggPCA(PCA1, metadata, colnames(metadata)[3:6])

# custom PCA plot after normalisation
pd <- metadata[,c(3:6)]
custom.PCA(beta = norm.data, pd = pd, plot.title = "PCA after normalisation")
# Correcting the data with sequencing depth it seems that Deseq2 is able to normalise the data correctly

# Differential expression analysis with Deseq2
# Differential expression analysis
control_vs_defeat <- results(dds, pAdjustMethod = "BH", contrast=c("Group", "defeated", "control"))
control_vs_defeat <- as.data.frame(control_vs_defeat) #37 FDR significant DEGs

sig_DEGs <- dplyr::filter(control_vs_defeat, padj <= 0.05)
characters <- ensembl2gene(rownames(sig_DEGs))

sig_DEGs$gene_names <- as.character(characters$external_gene_name)
write.csv(sig_DEGs, paste0(output, "sig_DEGS_CSDS.csv"))

#Volcanoplots
library(ggrepel)
control_vs_defeat$dif <- "NS"
# Make sure that the up and downregulated are underneath the nominal ones otherwise 
# the FDR significant genes are taken up by the nominal ones
control_vs_defeat$dif[control_vs_defeat[,2] > 0 & control_vs_defeat$pvalue < 0.05] <- "UP-nom"
control_vs_defeat$dif[control_vs_defeat[,2] < 0 & control_vs_defeat$pvalue < 0.05] <- "DOWN-nom"
control_vs_defeat$dif[control_vs_defeat[,2] > 0 & control_vs_defeat$padj < 0.05] <- "UP"
control_vs_defeat$dif[control_vs_defeat[,2] < 0 & control_vs_defeat$padj < 0.05] <- "DOWN"
mycolors <- c("#4169E1","#7d9bf5", "#AFADB3", "#c97171","#B22222")
names(mycolors) <- c("DOWN","DOWN-nom", "NS", "UP-nom", "UP")
control_vs_defeat$symbol <- rownames(control_vs_defeat)
#control_vs_defeat$delabel <- biomart_ensembl_to_gene(control_vs_defeat$symbol)
control_vs_defeat$delabel <- NA
control_vs_defeat$delabel[control_vs_defeat$dif %in% c("UP", "DOWN")] <- control_vs_defeat$symbol[control_vs_defeat$dif %in% c("UP", "DOWN")]
#test <- biomart_ensembl_to_gene(control_vs_defeat$delabel)
test <- ensembl2gene(control_vs_defeat$delabel)
x <- which(control_vs_defeat$delabel %in% test$ensembl_gene_id)
control_vs_defeat[x,9] <- test$external_gene_name
ggrepel.max.overlaps = Inf
p <- ggplot(data = control_vs_defeat, aes(x=control_vs_defeat[,2], y=-log10(pvalue), col=dif, label = delabel)) + geom_point() + theme_classic()
p2 <- p + scale_color_manual(values = mycolors) + geom_text_repel(aes(label=delabel), max.overlaps = Inf)
p3 <- p2 + labs(title = paste0("control vs. defeat")) + xlab("Log2FC") + ylab("-log10 p.value")
pdf(paste0(output, "volcanoplots.pdf"))
plot(p3)
dev.off()
# ------

# Boxplots of the top 3 p.FDR significant genes
rownames(metadata) <- paste0("X", metadata$description)
gene_boxplot_V2(log_norm, "ENSMUSG00000059751", metadata$Group)
gene_boxplot_V2(log_norm, "ENSMUSG00000083166", metadata$Group)
gene_boxplot_V2(log_norm, "ENSMUSG00000104802", metadata$Group)

# Boxplots of the top3 log2FC significant genes
gene_boxplot_V2(log_norm, "ENSMUSG00000112023", metadata$Group)
gene_boxplot_V2(log_norm, "ENSMUSG00000051457", metadata$Group)
gene_boxplot_V2(log_norm, "ENSMUSG00000022483", metadata$Group)


# Enrichment analysis based on the FDR significant genes from DESeq2
# -----------
# Pathway analysis
enrichment_analysis <- function(genes, universe){
  enriched <- list()
  # Enrichment analysis with msigdbr/clusterProfiler
  # genes is a character vector of gene names
  # univere is a character vector of gene names
  # Databases are from GSEA
  message("GO analysis")
  GO_datasets <- msigdbr(species = "Mus musculus", category = "C5")
  # Filter GO_datasets for CC, BP, and MF categories
  filtered_GO_datasets <- GO_datasets[GO_datasets$gs_subcat %in% c("GO:CC", "GO:BP", "GO:MF"), ]
  # Perform enrichment analysis using filtered GO_datasets
  enriched_GO <- as.data.frame(enricher(
    genes,
    universe = universe,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE = filtered_GO_datasets[, c("gs_name", "gene_symbol")],
    TERM2NAME = filtered_GO_datasets[, c("gs_name", "gs_exact_source")]
  ))
  enriched <- append(enriched, list(enriched_GO))
  # KEGG analysis
  message("KEGG analysis")
  KEGG_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
  enriched_KEGG <- as.data.frame(enricher(genes, universe = universe, pvalueCutoff = 0.05, 
                                          qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE =
                                            KEGG_datasets[,c("gs_name", "gene_symbol")]))
  enriched_KEGG$Description <- str_sub(enriched_KEGG$ID,6) 
  enriched_KEGG$Description <- str_replace_all(enriched_KEGG$Description, "_", " ")
  
  enriched <- append(enriched, list(enriched_KEGG))
  # Reactome analysis
  message("Reactome analysis")
  reactome_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
  enriched_reactome <- as.data.frame(enricher(genes, universe = universe, pvalueCutoff = 0.05, 
                                              qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE =
                                                reactome_datasets[,c("gs_name", "gene_symbol")]))
  enriched_reactome$Description <- str_sub(enriched_reactome$Description, 10)
  enriched_reactome$ID <- str_sub(enriched_reactome$ID, start = 1, end = 8)
  
  enriched <- append(enriched, list(enriched_reactome))
  
  
  names(enriched) <- c("GO", "KEGG", "Reactome")
  
  return(enriched)
}
genes_down <- which(control_vs_defeat$padj < 0.05 & control_vs_defeat$log2FoldChange < 0)
genes_down <- control_vs_defeat[genes_down,]
genes_up <- which(control_vs_defeat$padj < 0.05 & control_vs_defeat$log2FoldChange > 0)
genes_up <- control_vs_defeat[genes_up,]

universe <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(control_vs_defeat),
  mart = ensembl
)
genes_down <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(genes_down),
  mart = ensembl
)

genes_up <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(genes_up),
  mart = ensembl
)

enrichment_down <- enrichment_analysis(genes = genes_down$external_gene_name, universe = universe$external_gene_name)
enrichment_up <- enrichment_analysis(genes = genes_up$external_gene_name, universe = universe$external_gene_name)

all_results <- list()
all_results <- append(all_results, list(enrichment_down))
all_results <- append(all_results, list(enrichment_up))
names(all_results) <- c("down", "up")

# Visualization enrichment analysis
results <- data.frame()
for(i in names(all_results)){
  if(i == "down"){
    tmp <- all_results[[i]]
    result <- data.frame()
    for(database in names(tmp)){
      data <- data.frame()
      t <- tmp[[database]]
      if(dim(t)[1] != 0){
        if(dim(t)[1] > 5){
          d <- t[1:5,]
          data <- rbind(data, d)
          data$regulation <- "down"    
        }
        else {
          l <- dim(t)[1]
          d <- t[c(1:l),]
          data <- rbind(data, d)
          data$regulation <- "down"
        }
      }
      result <- rbind(result, data)
    }
  }
  results <- rbind(results, result)
  if(i == "up"){
    tmp <- all_results[[i]]
    result <- data.frame()
    for(database in names(tmp)){
      data <- data.frame()
      t <- tmp[[database]]
      if(dim(t)[1] != 0){
        if(dim(t)[1] > 5){
          d <- t[1:5,]
          data <- rbind(data, d)
          data$regulation <- "up"    
        }
        else {
          l <- dim(t)[1]
          d <- t[c(1:l),]
          data <- rbind(data, d)
          data$regulation <- "up"
        }
      }
      result <- rbind(result, data)
    }
    results <- rbind(results, result)
  }
}

pdf(paste0(output, "FDR_enriched_terms.pdf"))
if(dim(results)[1] != 0){
  results$log10 <- -log10(results$p.adjust)
  plot <- results %>% mutate(log10 = ifelse(regulation == "down", log10(p.adjust), log10))
  dot_plot <- ggplot(plot, aes(x = log10, y = reorder(Description, log10), size = Count)) +
    geom_point() + 
    scale_size_continuous(range = c(1, 10)) +
    labs(x = "-log10 FDR", y = "Enriched Terms") +
    ggtitle(paste0("Top 5 Enriched Terms")) +
    theme_classic() + theme(axis.text.y = element_text(size = 4))
  plot(dot_plot)
}
dev.off()

library(openxlsx)
# Writing the FDR enriched terms to excel sheets
for (c in names(all_results)){
  data <- all_results[[c]]
  tmp_data <- list()
  tmp_name <- list()
  
  for (i in names(data)) {
    if (nrow(data[[i]]) == 0) {
      message("empty")
    } else {
      d <- data[[i]]
      tmp_name <- append(tmp_name, i)
      tmp_data <- append(tmp_data, list(d))
    }
  }
  
  tmp_name <- as.character(tmp_name)
  names(tmp_data) <- tmp_name
  
  if (length(tmp_data) == 0) {
    print("No enrichment")
  } else {
    for (databases in 1:length(tmp_data)) {
      name <- names(tmp_data)
      message(databases)
      
      if (databases == 1) {
        n <- name[[databases]]
        wb <- createWorkbook()
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]])
        saveWorkbook(wb, file = paste0(output, "Enrichment_analysis_FDR_", c, "_regulated.xlsx"), overwrite = TRUE)
      } else {
        n <- name[[databases]]
        wb <- loadWorkbook(paste0(output, "Enrichment_analysis_FDR_", c, "_regulated.xlsx"))
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]], startRow = 1, startCol = 1, colNames = TRUE)
        saveWorkbook(wb, file = paste0(output, "Enrichment_analysis_FDR_", c, "_regulated.xlsx"), overwrite = TRUE)
      }
    }
  }
}



#-----------

# Enrichment analysis based on the nominal significant genes from DESeq2
nominal_DEGs <- filter(control_vs_defeat, pvalue <= 0.05)
control_vs_defeat$dif <- "NO"
control_vs_defeat$dif[control_vs_defeat[,2] > 0 & control_vs_defeat$pvalue < 0.05] <- "UP"
control_vs_defeat$dif[control_vs_defeat[,2] < 0 & control_vs_defeat$pvalue < 0.05] <- "DOWN"

genes_down <- which(control_vs_defeat$pvalue < 0.05 & control_vs_defeat$log2FoldChange < 0)
genes_down <- control_vs_defeat[genes_down,]
genes_up <- which(control_vs_defeat$pvalue < 0.05 & control_vs_defeat$log2FoldChange > 0)
genes_up <- control_vs_defeat[genes_up,]

genes_down <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(genes_down),
  mart = ensembl
)

genes_up <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(genes_up),
  mart = ensembl
)

enrichment_down <- enrichment_analysis(genes = genes_down$external_gene_name, universe = universe$external_gene_name)
enrichment_up <- enrichment_analysis(genes = genes_up$external_gene_name, universe = universe$external_gene_name)

all_results <- list()
all_results <- append(all_results, list(enrichment_down))
all_results <- append(all_results, list(enrichment_up))
names(all_results) <- c("down", "up")
# Writing the nominal enriched terms to excel sheets
for (c in names(all_results)){
  data <- all_results[[c]]
  tmp_data <- list()
  tmp_name <- list()
  
  for (i in names(data)) {
    if (nrow(data[[i]]) == 0) {
      message("empty")
    } else {
      d <- data[[i]]
      tmp_name <- append(tmp_name, i)
      tmp_data <- append(tmp_data, list(d))
    }
  }
  
  tmp_name <- as.character(tmp_name)
  names(tmp_data) <- tmp_name
  
  if (length(tmp_data) == 0) {
    print("No enrichment")
  } else {
    for (databases in 1:length(tmp_data)) {
      name <- names(tmp_data)
      message(databases)
      
      if (databases == 1) {
        n <- name[[databases]]
        wb <- createWorkbook()
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]])
        saveWorkbook(wb, file = paste0(output, "Enrichment_analysis_nominal_", c, "_regulated.xlsx"), overwrite = TRUE)
      } else {
        n <- name[[databases]]
        wb <- loadWorkbook(paste0(output, "Enrichment_analysis_nominal_", c, "_regulated.xlsx"))
        addWorksheet(wb, n)
        writeData(wb, n, tmp_data[[n]], startRow = 1, startCol = 1, colNames = TRUE)
        saveWorkbook(wb, file = paste0(output, "Enrichment_analysis_nominal_", c, "_regulated.xlsx"), overwrite = TRUE)
      }
    }
  }
}

# Visualization
results <- data.frame()
for(i in names(all_results)){
  if(i == "down"){
    tmp <- all_results[[i]]
    result <- data.frame()
    for(database in names(tmp)){
      data <- data.frame()
      t <- tmp[[database]]
      if(dim(t)[1] != 0){
        if(dim(t)[1] > 5){
          d <- t[1:5,]
          data <- rbind(data, d)
          data$regulation <- "down"    
        }
        else {
          l <- dim(t)[1]
          d <- t[c(1:l),]
          data <- rbind(data, d)
          data$regulation <- "down"
        }
      }
      result <- rbind(result, data)
    }
  }
  results <- rbind(results, result)
  if(i == "up"){
    tmp <- all_results[[i]]
    result <- data.frame()
    for(database in names(tmp)){
      data <- data.frame()
      t <- tmp[[database]]
      if(dim(t)[1] != 0){
        if(dim(t)[1] > 5){
          d <- t[1:5,]
          data <- rbind(data, d)
          data$regulation <- "up"    
        }
        else {
          l <- dim(t)[1]
          d <- t[c(1:l),]
          data <- rbind(data, d)
          data$regulation <- "up"
        }
      }
      result <- rbind(result, data)
    }
    results <- rbind(results, result)
  }
}
pdf(paste0(output, "nominal_enriched_terms.pdf"))
if(dim(results)[1] != 0){
  x <- grep("GO:", results$Description)
  results[x,2] <- results[x,1]
  results$log10 <- -log10(results$p.adjust)
  plot <- results %>% mutate(log10 = ifelse(regulation == "down", log10(p.adjust), log10))
  dot_plot <- ggplot(plot, aes(x = log10, y = reorder(Description, log10), size = Count)) +
    geom_point() + 
    scale_size_continuous(range = c(1, 10)) +
    labs(x = "-log10 FDR", y = "Enriched Terms") +
    ggtitle(paste0("Top 5 Enriched Terms")) +
    theme_classic() + theme(axis.text.y = element_text(size = 4))
  plot(dot_plot)
}
dev.off()

# Correlation analysis of the FDR significant genes and social interaction score
gene2ensembl <- function(x){
  ensembl=useMart("ensembl", host = "https://useast.ensembl.org/")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "external_gene_name",
    values = x,
    mart = ensembl)
  return(as.character(genes$ensembl_gene_id))
}

biomart_ensembl_to_gene <- function(gene){
  biomartCacheClear()
  ensembl=useMart("ensembl", host = "https://useast.ensembl.org/")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = gene, 
                 mart = ensembl)
  return(names$external_gene_name)
}
library(patchwork)

genes_for_correlation <- as.character(rownames(sig_DEGs))

pdf(paste0(output, "correlation_plots.pdf"))
for(gene_name in genes_for_correlation){
  message(gene_name)
  gene <- gene2ensembl(gene_name)
  x <- which(rownames(log_norm) == gene_name)
  gene_count <- as.data.frame(log_norm[x,])
  x <- which(metadata$Group == "defeated")
  metadata_defeated <- metadata[x,]
  x <- which(rownames(metadata_defeated) == colnames(gene_count))
  gene_count_df <- gene_count[,x]
  gene_count_df <- as.data.frame(cbind(t(gene_count_df), metadata_defeated$SI))
  colnames(gene_count_df) <- c("counts", "SI")
  gene_count_df <- as.data.frame(gene_count_df)
  gene_count_df <- gene_count_df %>% arrange(SI)
  test.lm <- lm(gene_count_df$counts ~ gene_count_df$SI)
  sum.lm <- summary(test.lm)
  cor_value <- cor.test(gene_count_df$SI, gene_count_df$counts)
  cor_coefficient <- cor_value$estimate
  p.val <- cor_value$p.value
  p1 <- ggplot(gene_count_df, aes(x=SI, y = counts)) +
    geom_point() + geom_smooth(method = 'lm') +
    geom_text(label = rownames(gene_count_df), nudge_x = 0.1) +
    labs(x = "Social interaction score", y = "Log2 normalized counts", title = gene,
         subtitle = paste("cor value", round(cor_coefficient,5), "", "p.value ", round(p.val, 3))) +
    theme_classic()
   p2 <- gene_boxplot_V2(log_norm, gene, metadata$Group)
   plot(wrap_plots(p1, p2))
}
dev.off()

