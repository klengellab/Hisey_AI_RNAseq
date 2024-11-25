# WGCNA analysis functions

# Joy Otten
# 01/25/2021; last updated 02/02/2021

##==========================================================================###
# Version of champ.SVD that uses PCA and has more options for the plot         #
# Modified by Roy Lardenoije                                                   #
# Version 1, last updated 2019-12-10                                           #
###==========================================================================###

# This function not just replaces SVD with PCA, but also performs an anova(lm())
# test for both numeric and factor variables, instead of an lm and kruskal test,
# as done by the original function. It also adds percentage of variance 
# explained for each PC to the axis label, and the F value (a measure of effect 
# size) to each box. Besides the centering of the beta values as done by the 
# original function, this function also scales and centers the numeric variables
# of pd

custom.PCA <- function (beta, rgSet = NULL, pd, RGEffect = FALSE, 
                        PDFplot = TRUE, Rplot = TRUE, 
                        resultsDir = "./PCAimages/",
                        cex.axis = 1,
                        return.p = FALSE) 
{
  message("[===========================]")
  message("[<<<<< CUSTOM.PCA START >>>>>]")
  message("-----------------------------")
  GenPlot <- function(thdens.o, estdens.o, evalues.v) {
    minx <- min(min(thdens.o$lambda), min(evalues.v))
    maxx <- max(max(thdens.o$lambda), max(evalues.v))
    miny <- min(min(thdens.o$dens), min(estdens.o$y))
    maxy <- max(max(thdens.o$dens), max(estdens.o$y))
  }
  EstDimRMTv2 <- function(data.m) {
    M <- data.m
    for (c in 1:ncol(M)) M[, c] <- (data.m[, c] - mean(data.m[, c]))/sqrt(var(data.m[, c]))
    sigma2 <- var(as.vector(M))
    Q <- nrow(data.m)/ncol(data.m)
    thdens.o <- thdens(Q, sigma2, ncol(data.m))
    C <- 1/nrow(M) * t(M) %*% M
    eigen.o <- eigen(C, symmetric = TRUE)
    estdens.o <- density(eigen.o$values, from = min(eigen.o$values), 
                         to = max(eigen.o$values), cut = 0)
    GenPlot(thdens.o, estdens.o, eigen.o$values)
    intdim <- length(which(eigen.o$values > thdens.o$max))
    return(list(cor = C, dim = intdim, estdens = estdens.o, 
                thdens = thdens.o))
  }
  thdens <- function(Q, sigma2, ns) {
    lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
    delta <- lambdaMAX - lambdaMIN
    roundN <- 3
    step <- round(delta/ns, roundN)
    while (step == 0) {
      roundN <- roundN + 1
      step <- round(delta/ns, roundN)
    }
    lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
    dens.v <- vector()
    ii <- 1
    for (i in lambda.v) {
      dens.v[ii] <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - i) * (i - lambdaMIN))/i
      ii <- ii + 1
    }
    return(list(min = lambdaMIN, max = lambdaMAX, step = step, 
                lambda = lambda.v, dens = dens.v))
  }
  runpca <- function(featuresbysamples, npca) { # modified function from DEGreport
    pca.res <- prcomp(t(featuresbysamples), center = FALSE, retx = TRUE)
    pc.var <- pca.res$sdev^2L
    pve <- 100L * (pc.var/sum(pc.var))
    #npca <- max(1L, length(which(pve >= min_pve_pct_pc)))
    samplepcvals <- pca.res$x[, 1L:npca, drop = FALSE]
    list(samplepcvals = samplepcvals, pve = pve)
  }
  if (!file.exists(resultsDir)) 
    dir.create(resultsDir)
  message("custom.PCA Results will be saved in ", resultsDir, " .\n")
  if (length(which(is.na(beta))) > 0) 
    message(length(which(is.na(beta))), 
            " NAs are detected in your beta value data set, which may cause errors in the PCA. You may want to impute NAs first.")
  message("[PCA analysis will proceed with ", dim(beta)[1], " probes and ", 
          dim(beta)[2], " samples.]\n")
  message("\n[ custom.PCA() will only check the dimensions between data and pd, instead of checking if sample names are correctly matched, thus please make sure your pd file is in accord with your data sets (beta values and rgSet).]\n")
  if (is.null(pd) | class(pd) == "list") 
    stop("pd parameter in data frame or matrix format is necessary and must contain at least two factors. If your pd is a list, please change its format.")
  if (class(pd) == "matrix") 
    pd <- as.data.frame(pd)
  PhenoTypes.lv_tmp <- pd[, !colnames(pd) %in% c("Sample_Name", "Project", 
                                                 "filenames", "Basename") & apply(pd, 2, function(x) length(unique(x))) != 1]
  PhenoTypes.lv <- PhenoTypes.lv_tmp
  if (!is.null(rownames(pd))) 
    rownames(PhenoTypes.lv) <- rownames(pd)
  if (ncol(PhenoTypes.lv) >= 2) {
    message("<< The following factors in your pd will be analyzed: >>")
    sapply(colnames(PhenoTypes.lv_tmp), 
           function(x) message("<", x, ">(", class(PhenoTypes.lv[[x]]), "):", paste(unique(PhenoTypes.lv_tmp[, x]), collapse = ", ")))
    message("[custom.PCA has automatically selected ALL factors that contain at least two different values from pd. If you don't want to analyze some of them, please remove them from pd.]")
  }
  else {
    stop("There are no factors with at least two values to be analyzed. Maybe your factors contain only one value?")
  }
  if (ncol(pd) > ncol(PhenoTypes.lv)) {
    message("\n<< The following factors in pd will not be analyzed: >>")
    sapply(setdiff(colnames(pd), colnames(PhenoTypes.lv)), 
           function(x) message("<", x, ">"))
    message("[Factors are ignored because they only indicate name or project, or they contain only one value across all samples.]")
  }
  if (RGEffect == TRUE & is.null(rgSet)) 
    message("If you want to check the effects of control probes, you must provide an rgSet. Now custom.SVD can only analyze factors in pd.")
  if (!is.null(rgSet) & RGEffect) {
    if (rgSet@annotation[1] == "IlluminaHumanMethylation450k") 
      data(ControlProbes450K)
    else data(ControlProbesEPIC)
    dataC2.m <- as.data.frame(log2(apply(ControlProbes, 1, 
                                         function(x) if (x[3] == "Grn") getGreen(rgSet)[x[2],] else getRed(rgSet)[x[2], ])))
    PhenoTypes.lv <- cbind(PhenoTypes.lv, dataC2.m)
    message("\n<< The following rgSet information has been included: >>")
    sapply(colnames(dataC2.m), function(x) message("<", x, ">"))
  }
  if (nrow(PhenoTypes.lv) == ncol(beta)) 
    message("\n<< PhenoTypes.lv generated successfully. >>")
  else stop("Dimension of your pd file (and rgSet information) is not equal to your beta value matrix.")
  beta <- as.matrix(beta)
  tmp.m <- beta - rowMeans(beta)
  rmt.o <- EstDimRMTv2(tmp.m)
  pca.o <- runpca(tmp.m, rmt.o$dim)
  if (ncol(pca.o$samplepcvals) > 20) 
    topPCA <- 20
  else topPCA <- ncol(pca.o$samplepcvals)
  pcaPV.m <- matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  pcaFV.m <- matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  colnames(pcaPV.m) <- colnames(pcaFV.m) <- colnames(PhenoTypes.lv)
  # for (c in 1:topPCA) for (f in 1:ncol(PhenoTypes.lv)) if (class(PhenoTypes.lv[, f]) != "numeric") 
  #     pcaPV.m[c, f] <- kruskal.test(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
  # else pcaPV.m[c, f] <- summary(lm(pca.o$samplepcvals[, c] ~ PhenoTypes.lv[[f]]))$coeff[2, 4]
  for (c in 1:topPCA) for (f in 1:ncol(PhenoTypes.lv)) if (class(PhenoTypes.lv[, f]) != "numeric") {
    pcaPV.m[c, f] <- anova(lm(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]])))$`Pr(>F)`[1]
    pcaFV.m[c, f] <- round(anova(lm(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]])))$`F value`[1], 0)
  } else {
    pcaPV.m[c, f] <- anova(lm(pca.o$samplepcvals[, c] ~ scale(PhenoTypes.lv[[f]])))$`Pr(>F)`[1]
    pcaFV.m[c, f] <- round(anova(lm(pca.o$samplepcvals[, c] ~ scale(PhenoTypes.lv[[f]])))$`F value`[1], 0)
  }
  message("<< Calculated PC matrix successfully. >>")
  myPalette <- c("darkred", "red", "orange", "pink", "white")
  breaks.v <- c(-200, -10, -5, -2, log10(0.05), 0)
  if (Rplot) {
    par(mar = c(5, 15, 2, 1))
    image(x = 1:nrow(pcaPV.m), y = 1:ncol(pcaPV.m), z = log10(pcaPV.m), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = "Principal Components Analysis (PCA)")
    text(x = row(pcaFV.m), y = col(pcaFV.m), pcaFV.m, cex = 0.5)
    axis(1, at = 1:nrow(pcaPV.m), labels = paste0("PC-", 1:nrow(pcaPV.m), 
                                                  " (", round(pca.o$pve[1:topPCA], 2), 
                                                  "%)"), 
         las = 2, cex.axis = 0.6)
    suppressWarnings(axis(2, at = 1:ncol(pcaPV.m), 
                          labels = colnames(pcaPV.m), las = 2, 
                          cex.axis = cex.axis))
    legend(x = -(topPCA/2.2), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-5}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = c("darkred", "red", "orange", "pink", "white"), 
           par("usr")[2], 
           par("usr")[4], xpd = NA,
           bty = "n",
           cex = 0.7)
  }
  if (PDFplot) {
    pdf(paste(resultsDir, "PCAsummary.pdf", sep = ""), width = 10, 
        height = 8)
    par(mar = c(5, 15, 2, 1), xpd = TRUE)
    image(x = 1:nrow(pcaPV.m), y = 1:ncol(pcaPV.m), z = log10(pcaPV.m), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = "Principal Components Analysis (PCA)")
    text(x = row(pcaFV.m), y = col(pcaFV.m), pcaFV.m, cex = 0.5)
    axis(1, at = 1:nrow(pcaPV.m), labels = paste0("PC-", 1:nrow(pcaPV.m), 
                                                  " (", round(pca.o$pve[1:topPCA], 2), 
                                                  "%)"), 
         las = 2, cex.axis = 0.6)
    suppressWarnings(axis(2, at = 1:ncol(pcaPV.m), 
                          labels = colnames(pcaPV.m), las = 2, 
                          cex.axis = cex.axis))
    legend(x = -(topPCA/2.2), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-5}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = c("darkred", "red", "orange", "pink", "white"), 
           par("usr")[2], 
           par("usr")[4], xpd = NA,
           bty = "n",
           cex = 0.7)
    dev.off()
  }
  message("<< Plotted PCA matrix successfully. >>")
  message("[<<<<<< custom.PCA END >>>>>>]")
  message("[===========================]")
  if(return.p){
    pcaPV.m <- apply(pcaPV.m, 2, function(x) sprintf("%.2e", x))
    rownames(pcaPV.m) <- paste0("PC", 1:topPCA)
    return(t(pcaPV.m)[ncol(pcaPV.m):1,])
  } 
}


data_preprocessing <- function(Expr, verbose = FALSE){
  # Expr data is dataframe columns are the samples and rows the genes
  # This function filters on that the median variance is higher than 0,
  # if not this is filtered out. 
  #keep <- apply(Expr,2,function(x){any(mad(x) > 0)})
  #Expr <- Expr[, keep]
  
  # checks for genes and samples with too many missing values. 
  print(dim(Expr))
  gsg = goodSamplesGenes(Expr, verbose = 3)
  print(gsg$allOK)
  
  # removes genes with too many missing values
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(Expr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(Expr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    Expr = Expr[gsg$goodSamples, gsg$goodGenes]
    dim(Expr)
  }
  return(Expr)
}

clust_outliers <- function(Expr, filter, verbose = FALSE){
  # filter is an object wherein you choose if you want to filter out the outliers
  # from the rest of the data or that you want to take them along.
  # plots the outliers in the samples
  sampleTree <- hclust(dist(Expr), method = "average")
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  pdf(paste0(image_path, "sample_outliers.pdf"))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="",
       xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  abline(h = 15, col = "red")
  dev.off()
  clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
  print(table(clust))
  if(filter == TRUE){
    keepSamples = (clust==1)
    Expr <- Expr[keepSamples, ]
  }
  nGenes = ncol(Expr)
  nSamples = nrow(Expr)
  
  return(Expr)
}

cluster_samples <- function(Expr, datTraits, verbose = FALSE){
  sampleTree2 = hclust(dist(Expr), method = "average")
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  pdf(paste0(image_path, "sampledendro.pdf"))
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}

power_calc <- function(Expr){
  # Calculates the power needed for WGCNA check always the plots
  # it could be that the power estimate gives another value back than 
  # that you think according to the plots.
  powers = c(c(1:10), seq(from = 12, to = 20, by =2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5)
  #Plot the results
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9
  pdf(paste0(image_path, "power_calc.pdf"))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2", type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  dev.off()
  pdf(paste0(image_path, "mean_connectivity.pdf"))
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab ="Soft Threshold (power)", ylab = "Mean connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  print(sft$fitIndices)
  power <- sft$powerEstimate
  
  return(power)
}
  
WGCNA_module <- function(Expr, power, direction, split){
  cor <- WGCNA::cor
  enableWGCNAThreads(8)
  modules = blockwiseModules(Expr, power = power, cortype = "pearson", 
                             networkType = direction,
                             TOMType = direction, minModuleSize = 30,
                             maxBlockSize = 30000, deepSplit = split,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             verbose = 3)
  disableWGCNAThreads()
  cor <- stats::cor
  
  # plotting the figures
  plot_modules <- function(module){
    mergedColors = labels2colors(module$colors)
    pdf(paste0(image_path, "modules.pdf"))
    plotDendroAndColors(module$dendrogram[[1]],
                        mergedColors[module$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
  }
  
  return(modules)
}

TOMplot <- function(Expr, module, power){
  # Plots a TOM adjacency plot, although I don't find it particulary handy
  blocknumber = 1
  m <- module$colors
  moduleColors <- labels2colors(m)
  dynamicColors <- data.frame(moduleColors)[module$blockGenes[[1]],]
  restGenes <- (dynamicColors != "grey")
  diss1 <- 1-TOMsimilarityFromExpr(Expr[,restGenes], power = power)
  Subgenenames = colnames(Expr)
  colnames(diss1) = rownmaes(diss1) = Subgenenames[restGenes]
  hier1 = flashClust(as.dist(diss1), method = "average")
  diag(diss1) = NA
  
  pdf(paste0(image_path, "Tomplot.pdf"))
  TOMplot(diss1, hier1, as.character(moduleColors[restGenes]))
  dev.off()
  
  module_colors = setdiff(unique(moduleColors), "grey")
  for (color in module_colors){
    module = Subgenenames[which(dynamicColors == color)]
    write.table(module, paste("module_", color, ".txt", sep = ""), sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  module.order <- unlist(tapply(1:ncol(Expr),as.factor(dynamicColors),I))
  m<-t(t(Expr[,module.order])/apply(Expr[,module.order],2,max))
  
  # plots expression patterns of genes in modules although not very handy if
  # you use dynamicTreecut.
  pdf(paste0(image_path, "expression_patterns.pdf"))
  heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
  dev.off()
}

eigengene_network <- function(Expr, module, metadata){
  moduleLabels <- module$colors
  moduleColors = labels2colors(moduleLabels)
  # plot eigeen gene network
  ME1 <- moduleEigengenes(Expr, moduleColors)$eigengenes
  MET = orderMEs(cbind(ME1, metadata))
  pdf(paste0(image_path, "eigengene_network.pdf"))
  plotEigengeneNetworks(MET, "", marDendro = c(4,4,4,4), marHeatmap = c(4,6,2,6),
                        cex.lab = 0.8, xLabelsAngle=90)
  dev.off()
  
  return(ME1)
}

module_trait <- function(Expr, metadata){
  # Visualisation of modules with the phenotype data
  nGenes = ncol(Expr)
  nSamples = nrow(Expr)
  MEs0 <- moduleEigengenes(Expr, moduleColors)$eigengenes
  MEs0 <- orderMEs(MEs0)
  modTraitCor = stats::cor(MEs0, metadata, use = "p", method = "pearson")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
  textMatrix = paste(signif(modTraitCor, 2), "\n(",
                     signif(modTraitP, 1), ")", se = "")
  dim(textMatrix) = dim(modTraitCor)
  par(mar = c(4,6,2,2))
  pdf(paste0(image_path, "module_trait_relationship.pdf"))
  labeledHeatmap(Matrix = modTraitCor, xLabels = names(metadata), 
                 yLabels = names(MEs0), ySymbols = names(MEs0),
                 colorLabels = F, colors = blueWhiteRed(50),
                 textMatrix = textMatrix, setStdMargins = F, cex.text = 0.5,
                 zlim = c(-1,1), main = paste("Module-trait relationships"))
  dev.off()
  
  return(MEs0)
}

intramodular_analyses <- function(kme, group, variables){
  # Intramodular analysis: identifying genes with high GS and MM
  colorOfColumn = substring(names(kme), 4)
  par(mfrow = c(2,2))
  color_column <- as.character(unique(colorOfColumn))
  selectModules = as.character(unique(colorOfColumn))
  par(mfrow=c(2,length(selectModules)/2))
  
  for(i in group){
    pdf(paste0(image_path, "intramodular_", i, ".pdf"))
    for (module in selectModules){
      column = match(module, colorOfColumn)
      restModule = moduleColors==module
      p1 <- verboseScatterplot(kme[restModule,column], variables[[i]][restModule],
                               xlab=paste("Module Membership", module, "module"),
                               ylab=i, main = paste("kME", module, "GS"),
                               col = module)
    print(p1)
  }
  dev.off()
  }
}

hubgenes <- function(Expr, moduleColors, power, direction, MEs1, ensembl = FALSE, output){
  # output is a path were the results will be saved
  # ensembl: If there are ensembl gene id's
  hub_genes <- as.data.frame(chooseTopHubInEachModule(Expr, moduleColors, omitColors = "grey", 
                                                      power = power, type = direction))
  names(hub_genes) <- "genes"
  
  list_hubgenes <- t(as.data.frame(mclapply(colnames(Expr), function(x){
    y <- numeric()
    for(i in colnames(MEs1)){
      y <- append(y, unlist(cor.test(Expr[, x], MEs1[, i], method = "pearson")[4:3]))
    }
    return(y)
  })))
  rownames(list_hubgenes) <- colnames(Expr)
  colnames(list_hubgenes) <- paste0(rep(colnames(MEs1), each = 2), rep(c(".r", ".p"), ncol(MEs1)))
  
  if(ensembl == TRUE){
    # convert Ensembl ID to NCBI gene ID
    ensembl <- useMart("ensembl")
    ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
    attributes = listAttributes(ensembl)
    
    x <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = hub_genes,
               mart = ensembl)
    y <- match(x$ensembl_gene_id, hub_genes$genes)
    x <- x[y,]
    hub_genes$gene_names <- x$external_gene_name
    
    print(hub_genes)
    write.csv(hub_genes, paste0(output, "hubgenes.csv"))
    
    x <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = rownames(list_hubgenes),
               mart = ensembl)
    
    y <- match(x$ensembl_gene_id, rownames(list_hubgenes))
    message("It could be that in your list of hubgenes are less genes compared to
          the expr data. This is due to that those ensembl names don't have a gene
          name")
    list_hubgenes <- list_hubgenes[y,]
    rownames(list_hubgenes) <- x$external_gene_name
  }
  
  write.csv(list_hubgenes, paste0(output, "list_hubgenes.csv"))
}

# GO analysismodule_genes, all_genes
GO_analysis <- function(genes, universe){
  require(dplyr)
  require(clusterProfiler)
  require(msigdbr)
  require(stringr)
  y <- enrichGO(
    gene = genes,
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universe,
    qvalueCutoff = 0.05,
    readable = TRUE)
  if(length(y) == 1){
    q <- as.data.frame(y@result)
    return(q)
  } else {
    message("empty data frame no data")
    q <- DataFrame()
    return(q)
    }
}

# Kegg pathway analysis
Kegg_analysis <- function(x, universe){
  enrichKEGG(gene = x, organism = "mmu", keyType = "kegg", 
             pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = universe)
}

enrichment_analysis <- function(moduleLabels, moduleColors, Expr, Labels1, output, universe){
  enrichment <- list()
  # keytype: character to identify the gene id's such as ENSEMBL, external_gene_name etc
  Labels1$moduleColors <- moduleColors
  Labels1$moduleColors = moduleColors
  Labels1$genes <- colnames(Expr)
  select_genes <- function(x){
    require(dplyr)
    y <- dplyr::filter(Labels1, moduleColors == x)
    genes <- y$genes
    return(as.data.frame(genes))
  }
  
  # Getting the genes per module
  Colours <- as.list(unique(Labels1$moduleColors))
  names(Colours) <- unique(Labels1$moduleColors)
  list_modules <- lapply(Colours, select_genes)
  
  for(i in Colours){
    enriched <- list()
    GO_datasets <- msigdbr(species = "Mus musculus", category = "C5")
    # Select for only the BP, CC and MF category in GO analysis
    x <- which(GO_datasets$gs_subcat %in% c("GO:BP", "GO:CC","GO:MF"))
    GO_datasets <- GO_datasets[x,]
    
    x <- Colours[[i]]
    genes <- list_modules[[x]][["genes"]]
    enriched_GO <- as.data.frame(enricher(genes, pvalueCutoff = 0.05, universe = universe,
                                          qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = GO_datasets[,c("gs_name", "ensembl_gene")], TERM2NAME = GO_datasets[,c("gs_name", "gs_exact_source", "ensembl_gene", "gene_symbol")]))
    enriched <- append(enriched, list(enriched_GO))
    n <- enriched_GO[c(1:10),]
    n <- na.omit(n)
    p1 <- ggplot(data=n,
                aes(x = ID, y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + theme_classic() + coord_flip() + 
      labs(x = "", y = "Count") + guides(fill=guide_legend(title = "p.value"))
    
    # KEGG analysis
    KEGG_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
    enriched_KEGG <- as.data.frame(enricher(genes, universe = universe, pvalueCutoff = 0.05, 
                                            qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = KEGG_datasets[,c("gs_name", "ensembl_gene")]))
    enriched_KEGG$Description <- str_sub(enriched_KEGG$ID,6) 
    enriched_KEGG$Description <- str_replace_all(enriched_KEGG$Description, "_", " ")
    n <- enriched_KEGG[c(1:10),]
    n <- na.omit(n)
    p2 <- ggplot(data=n,
                 aes(x = Description, y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + theme_classic() + coord_flip() + 
      labs(x = "", y = "Count") + guides(fill=guide_legend(title = "p.value"))
    
    enriched <- append(enriched, list(enriched_KEGG))
    
    # Reactome analysis
    reactome_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
    enriched_reactome <- as.data.frame(enricher(genes, universe = universe, pvalueCutoff = 0.05, 
                                                qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = reactome_datasets[,c("gs_name", "ensembl_gene")]))
    enriched_reactome$Description <- str_sub(enriched_reactome$Description, 10)
    enriched_reactome$ID <- str_sub(enriched_reactome$ID, start = 1, end = 8)
    n <- enriched_reactome[c(1:10),]
    n <- na.omit(n)
    p3 <- ggplot(data=n,
                 aes(x = Description, y = Count, fill = p.adjust)) +
      geom_bar(stat = "identity") + theme_classic() + coord_flip() + 
      labs(x = "", y = "Count") + guides(fill=guide_legend(title = "p.value"))
    
    
    enriched <- append(enriched, list(enriched_reactome))
    pdf(paste0(output, "module_", i, "enrichment_analysis.pdf"))
    plot(p1)
    plot(p2)
    plot(p3)
    dev.off()
    names(enriched) <- c("GO", "KEGG", "Reactome")
    enrichment <- append(enrichment, list(enriched))
  }
  names(enrichment) <- as.character(Colours)
  return(enrichment)
  
}




