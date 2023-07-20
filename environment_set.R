libraries_set <- function(){
  library(limma)
  library(NbClust)
  library(ggplot2)
  library(ggrepel)
  library(statmod)
  library(Hmisc)
  library(viper)
  library(caret)
  library(ComplexHeatmap)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggbeeswarm)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(genefilter)
  library(EnhancedVolcano)
  library(circlize)
  library(enrichplot)
  library(DOSE)
  library(GSVA)
  library(msigdbr)
  library(tidyverse)
  library(reactome.db)
  library(igraph)
  library(fgsea)
  library(plyr)
  library(viper)
  library(dorothea)
  library(dplyr)
  library(tidyverse)
  library(progeny)
  library(clusterProfiler)
  library(ggupset)
  library(ggfortify)
  library(pathview)
  library(ggridges)
  library(WGCNA)
  library(factoextra)
  library(corrplot)
  library(gprofiler2)
  library(RCy3)
  library(cluster)
  library(ReactomePA)
  library(ggvenn)
  library(forcats)
  library(ggtext)
  library(mctest)
  library(GGally)
  library(readxl)
  library(PCAtools)
  library(cowplot)
  library(maftools)
  library(ggplotify)
  library(pROC)
  library(purrr)
  library(tidyverse)
  library(dendextend)
  library(ggthemr)
  library(reshape2)
  library(randomForest)
  library(caret)
  library(MASS)
  library(Seurat)
  library(SeuratDisk)
  library(decoupleR)
  library(tibble)
  library(patchwork)
  library(OmnipathR)
  library(ExperimentHub)
  library(SingleR)
}

minMax <- function(x) {
  #columns: features
  x = data.matrix(x)
  for(i in 1:ncol(x)){
    x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))    
  }
  
  return(x)

}

#' @param data matrix; TPM values (genes as rows and samples as columns)
TPM_normalization <- function(data, log = FALSE, pseudo = 1) {
  
  # TPM normalization
  if(log){
    
    if (pseudo == 0)
      
      warning("Using 0 pseudo: Inf may be generated.\n")
    
    data <- log2(data + pseudo)
  }
  
  TPM_data <- t(t(data)*1e6/apply(data,2,sum))
  
  return(TPM_data)
}


# Convert FPKM values to transcripts per million (TPM).
fpkm2tpm <- function(fpkm){
  tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
  tpm[which(is.na(tpm))] <- 0
  return(tpm)
}

##Pathway analysis

compute_pathways_scores <- function(RNA.counts.normalized, file_name){
  gene_expr = RNA.counts.normalized
  
  ########################################################################################################################################
  
  # Pathways activity (Progeny package)
  Pathway_scores <- progeny::progeny(gene_expr, scale = FALSE, organism = "Human", verbose = TRUE)
  
  ########################################################################################################################################
  
  # GSVA pathways (Reactome database)
  reactome_gene_sets <- AnnotationDbi::select(reactome.db, keys=keys(reactome.db), columns=c("PATHNAME","REACTOMEID"))
  reactome_gene_sets_hs <- reactome_gene_sets[grep("Homo sapiens: ", iconv(reactome_gene_sets$PATHNAME)),]
  gene_set_list <- split(reactome_gene_sets_hs$ENTREZID, reactome_gene_sets_hs$PATHNAME)
  
  gene_expr = data.frame(gene_expr)
  gene_expr$Symbol = rownames(gene_expr)
  mapped_df <- data.frame("entrez_id" = mapIds(org.Hs.eg.db,
                                               keys = gene_expr$Symbol,
                                               keytype = "SYMBOL",
                                               column = "ENTREZID",
                                               multiVals = "first")
  ) %>% 
    dplyr::filter(!is.na(entrez_id)) %>%  
    tibble::rownames_to_column("Symbol") %>% 
    dplyr::inner_join(gene_expr, by = "Symbol") 
  
  mapped_df$Symbol <- NULL
  rownames(mapped_df) <- mapped_df$entrez_id
  mapped_df$entrez_id <- NULL
  mapped_df <- as.matrix(mapped_df)
  
  gsva_scores <- gsva(
    mapped_df,
    gene_set_list,
    method = "gsva", 
    kcdf = "Gaussian", 
    min.sz = 15, 
    max.sz = 500,
    mx.diff = TRUE,
    verbose = FALSE
  )
  
  ########################################################################################################################################
  
  #ORA enrichment
  mapped_df = data.frame(mapped_df)
  ora_results <- enrichKEGG(gene       = rownames(mapped_df),
                            organism     = "hsa",
                            pvalueCutoff = 0.05, 
                            use_internal_data = T)
  
  ########################################################################################################################################
  
  # Output list:
  progeny_results <- as.data.frame(Pathway_scores)
  reactome_results <- as.data.frame(gsva_scores)
  kegg_results <- as.data.frame(ora_results)
  
  message("\n Pathway scores computed \n")
  
  write.csv(progeny_results,paste0(getwd(),"/Output/Pathways/Pathways_Progeny_", file_name, ".csv"))
  write.csv(reactome_results,paste0(getwd(),"/Output/Pathways/Pathways_Reactome_", file_name, ".csv"))
  write.csv(kegg_results,paste0(getwd(),"/Output/Pathways/Pathways_KEGG_", file_name, ".csv"))
  
  return(list(progeny_results, reactome_results, kegg_results))
  
}

## TF analysis

compute_TFs_scores = function(RNA.counts.normalized,  file_name){

  dorothea2viper_regulons <- function(df) {
    regulon_list <- split(df, df$tf)
    viper_regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    
    return(viper_regulons)
  }
  
  data("dorothea_hs", package = "dorothea")
  regulons <- dorothea_hs %>%
    filter(confidence %in% c("A", "B", "C", "D"))
  
  regu <- dorothea2viper_regulons(regulons)
  RNA.counts.normalized = as.data.frame(RNA.counts.normalized)
  vpres<- viper(RNA.counts.normalized, regu, verbose = FALSE, minsize = 4)
  
  write.csv(vpres,paste0(getwd(),"/Output/TFs/TFs_scores_all_", file_name, ".csv"))
  
  return(vpres)
  
}

compute_msVIPER_scores = function(RNA.counts.normalized, test_name, ref_name, measure, file_name){
  
  dorothea2viper_regulons <- function(df) {
    regulon_list <- split(df, df$tf)
    viper_regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    
    return(viper_regulons)
  }
  
  data("dorothea_hs", package = "dorothea")
  regulons <- dorothea_hs %>%
    filter(confidence %in% c("A", "B", "C", "D"))
  
  regu <- dorothea2viper_regulons(regulons)
  RNA.counts.normalized = as.data.frame(RNA.counts.normalized)
  vpres<- viper(RNA.counts.normalized, regu, verbose = FALSE, minsize = 4)
  
  # Generating test and ref data
  test_i <- which(measure == test_name)
  ref_i <- which(measure == ref_name)
  
  mat_test <- as.matrix(RNA.counts.normalized[,test_i])
  mat_ref <- as.matrix(RNA.counts.normalized[,ref_i])
  
  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000)
  
  # Generating signature
  signature <- rowTtest(mat_test, mat_ref)
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))
  signature <- na.omit(signature)
  signature <- signature[,1]
  
  # Running msVIPER
  mra <- msviper(signature, regu, dnull, verbose = FALSE)
  
  # Plot DiffActive TFs
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  
}


##Correlations

correlation <- function(df) {
  
  M <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
  
  corr_df = Mdf %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
  
  corr_df = na.omit(corr_df)  #remove the ones that are the same TFs (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)
  
  return(corr_df)
  
}

##Regression plot 

ggplotRegression <- function (x, y, measure) {
  x = data.frame(x)
  for(i in 1:ncol(x)){
    dlm = cbind(x[,i], y)
    colnames(dlm)[1] = colnames(x)[i]
    colnames(dlm)[2] = measure
    dlm = as.data.frame(dlm)
    fit <- lm(dlm[,1] ~ dlm[,2], dlm)
    png(paste0(getwd(),"/Figures/Regression_plots/", colnames(x)[i], "vs", measure, ".png"))
    print(ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
            geom_point() +
            stat_smooth(method = "lm", col = "#1d4e89ff") +
            xlab(colnames(x)[i]) + 
            ylab(measure) +
            labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                               "Intercept =",signif(fit$coef[[1]],5 ),
                               " Slope =",signif(fit$coef[[2]], 5),
                               " P =",signif(summary(fit)$coef[2,4], 5))))
    dev.off()
  }
}

##Violin plot

compute_violin = function(data, feature, file_name){
  
  matrix = data.frame(data[1:nrow(data)-1,])
  matrix$measure = feature
  pval = data[nrow(data),]
  
  for (i in (1:(ncol(matrix)-1))) {
    violin = ggplot(matrix, aes(x=factor(measure), y=matrix[,i], fill=measure)) +
      geom_violin(width=0.6) +
      geom_boxplot(width=0.07, color="black", alpha=0.2) +
      scale_fill_brewer() +
      geom_smooth(aes(x=as.factor(measure), y=matrix[,i]), method = "loess") +
      ylab(colnames(matrix[i])) +
      xlab(file_name) +
      labs(title=colnames(matrix[i]),
           subtitle=paste0("ANOVA test\npvalue: ", pval[i])) +
      theme(axis.text.x = element_text(angle = 0),
            axis.title.y = element_text(size = 8, angle = 90))
    
    
    ggsave(violin, file=paste0(getwd(),"/Figures/",file_name, "/", colnames(matrix[i]),".png"))
  }
}


##Correlation plot

plot_correlation <- function(data, title){
  #row: features and columns: samples 
  d <- dist(data, method = "euclidean")
  hc1 <- hclust(d, method = "complete" )
  vec = hc1[["order"]]
  order = rownames(data)[vec]
  
  M <- Hmisc::rcorr(data.matrix(t(data)), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
  
  corr_df = Mdf %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
  
  g <- corr_df %>%
    mutate(measure1 = fct_relevel(measure1, order)) %>%
    mutate(measure2 = fct_relevel(measure2, order)) %>%
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title=title,
         subtitle="Only significant Pearson's correlation coefficients shown") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    rotate_x_text(angle = 45) + theme(axis.text.x=element_markdown()) + theme(axis.text.y=element_markdown())
  
  print(g)
  
}

#Compute PCA analysis
compute_pca_analysis <- function(data, coldata, n_components,  file_name){
  
  p <- pca(data, metadata = coldata, removeVar = 0.1)
  
  peigencor  <- eigencorplot(p,
                             components = getComponents(p, 1:n_components),
                             metavars = colnames(coldata),
                             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
                             cexCorval = 1.2,
                             fontCorval = 2,
                             posLab = 'all',
                             rotLabX = 45,
                             scale = TRUE,
                             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
                             plotRsquared = TRUE,
                             corFUN = 'pearson',
                             corUSE = 'pairwise.complete.obs',
                             corMultipleTestCorrection = 'BH',
                             signifSymbols = c('****', '***', '**', '*', ''),
                             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
  
  print(peigencor)
}

random_forest <- function(data, target, k_features) {
  # Step 1: Random Forest Feature Selection
  target = as.factor(target)
  # Step 2: Random Forest Classification
  rf_classifier <- randomForest(target ~ ., data = data)
  rf_importance <- importance(rf_classifier)
  rf_importance <- rf_importance[order(rf_importance, decreasing = TRUE), ]
  rf_selected_features <- data[, names(rf_importance[order(rf_importance, decreasing = TRUE)][1:k_features])]
  rf_importance = data.frame(rf_importance)
  colnames(rf_importance)[1] = "Importance"
  return(list(rf_selected_features, rf_importance))
}

heatmap_annotation = function(annotations){
  
  greyscale <- grey.colors(10, rev = T)
  col_fun = colorRamp2(c(0, 0.5, 1), c(brewer.pal(9, "YlOrBr")[1], brewer.pal(9, "YlOrBr")[4], brewer.pal(9, "YlOrBr")[9]))
  
  ann_colors <- HeatmapAnnotation(
    Stage = annotations$Stage, 
    CYT = annotations$CYT,
    Roh_IS = annotations$Roh_IS,
    Davoli_IS = annotations$Davoli_IS,
    Chemokines = annotations$chemokines,
    IFNy = annotations$IFNy,
    Ayers_expIS = annotations$Ayers_expIS,
    Tcell_inflamed = annotations$Tcell_inflamed,
    
    col = list(
      Stage = c("I" = greyscale[1], "II" = greyscale[4], "III" = greyscale[7], "IV" = greyscale[10]),
      CYT = col_fun,
      Roh_IS = col_fun,
      Davoli_IS = col_fun,
      Chemokines = col_fun,
      IFNy = col_fun,
      Ayers_expIS = col_fun,
      Tcell_inflamed = col_fun
    ),
    
    annotation_legend_param = list(labels_gp = gpar(fontsize = 10), legend_width = unit(12, "cm"), 
                                   legend_heigh = unit(12, "cm"), title_gp = gpar(fontsize = 12))
    
  )
  
  return(ann_colors)
  
}


