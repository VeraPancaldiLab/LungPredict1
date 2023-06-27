libraries_set <- function(){
  library(limma)
  library(NbClust)
  library(ggplot2)
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
  library("dendextend")
  library(ggthemr)
  library(reshape2)
  library(randomForest)
  library(caret)
  library(MASS)
  #library(gProfileR)
  #library(rsample)
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

##Function differential expression

diff_expression = function(ddsObj, measure, num, den, file_name){

  ddsObj = DESeq(ddsObj)
  resultsALL = results(ddsObj, contrast = c(measure, num, den), pAdjustMethod = "BH")

  png(paste0(getwd(),"/Figures/Volcanos/Volcano_", file_name, ".png"), width = 3000, height = 2000, res=100)
  print(EnhancedVolcano(resultsALL,
                        title = paste(measure, num, "vs", den),
                        subtitle = NULL,
                        caption = NULL,
                        lab = rownames(resultsALL),
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-20, 20),
                        ylim = c(0, 15), 
                        pCutoff = 10e-2,
                        FCcutoff = 1))
  dev.off()
  
  write.csv(resultsALL,paste0(getwd(),"/Output/DEG_", file_name, ".csv")) 
  
}

##Data vsd normalization

vsd_normalization = function(RNA.counts, file_name, coldata, remove.batch = T){

  #RNA.counts = Counts
  RNA.counts = data.frame(RNA.counts)
  raw_counts <- RNA.counts
  genes <- rownames(raw_counts)
  
  if (is.integer(raw_counts) == FALSE) {
    raw_counts.integer <- apply(raw_counts, 2, as.integer)
    rownames(raw_counts.integer) <- rownames(raw_counts)
  } else{
    raw_counts.integer <- raw_counts
  }
  
  message("DESeq2 -->  \n")
  
  dset <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts.integer,
                                         colData = coldata,
                                         design = ~ 1) #design =~1 to consider all samples as part of the same group   

  
  dset <- DESeq2::estimateSizeFactors(dset)
  dset <- DESeq2::estimateDispersions(dset)
  gene_expr <- DESeq2::getVarianceStabilizedData(dset)
  rownames(gene_expr) <- rownames(raw_counts.integer)
  
  write.csv(gene_expr,paste0(getwd(),"/Output/Counts_normalized", file_name, ".csv"))
  
  return(gene_expr)
  
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
  
  png(paste0(getwd(),"/Figures/Dotplots/ORA_PathwaysKEGG_", file_name, ".png"), width = 1500, height = 1000, res=100)
  print(dotplot(ora_results))
  dev.off()
  
  ########################################################################################################################################
  
  # Output list:
  progeny_results <- as.data.frame(Pathway_scores)
  reactome_results <- as.data.frame(gsva_scores)
  kegg_results <- as.data.frame(ora_results)
  
  message("\n Pathway scores computed \n")
  
  write.csv(progeny_results,paste0(getwd(),"/Output/Pathways/Pathways_Progeny_", file_name, ".csv"))
  write.csv(reactome_results,paste0(getwd(),"/Output/Pathways/Pathways_Reactome_", file_name, ".csv"))
  write.csv(kegg_results,paste0(getwd(),"/Output/Pathways/Pathways_KEGG_", file_name, ".csv"))
  
  return(list(progeny_results, reactome_results, kegg_results))#, kegg_results))
  
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
  
  write.csv(vpres,paste0(getwd(),"/Output/TFs/VIPER/TFs_scores_all_", file_name, ".csv"))
  
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
  
  # Extract Differential Active TFs
  TFs = data.frame(mra[["es"]][["p.value"]])
  colnames(TFs) = "pval"
  vec = TFs$pval < 0.05
  TFs_msViper_adj = TFs[vec,]
  TFs_msViper_adj = as.data.frame(TFs_msViper_adj)
  rownames(TFs_msViper_adj) = rownames(TFs)[vec]
  TFs_msviper = vpres[rownames(vpres)%in%rownames(TFs_msViper_adj),]
  
  write.csv(TFs_msviper,paste0(getwd(),"/Output/TFs/msVIPER/TFs_scores_DifferentialActive_", file_name, ".csv"))
  
  # Plot DiffActive TFs
  png(paste0(getwd(),"/Figures/msVIPER/DifferentialActive_TFs_", file_name, ".png"), width = 1000, height = 800, res=100)
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  dev.off()
  
  return(TFs_msviper)
  
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
##Wilcox analysis

compute_wilcox = function(data, test, ref, feature, file_name){
  
  data = data.frame(data)
  data$measure <- feature

  test_df = data[grep(test, data$measure), ]
  ref_df = data[grep(ref, data$measure), ]
  
  testVSref = list()
  
  for (i in (1:(ncol(data)-1))) {
    HL = (wilcox.test(test_df[,1], ref_df[,1], exact = FALSE))
    testVSref[[length(testVSref)+1]] = HL$p.value
    keepWilcox <- which(testVSref <0.05)
    wilcoxPass_test = test_df[keepWilcox]
    wilcoxPass_ref = ref_df[keepWilcox]
    wilcoxPass = rbind(wilcoxPass_test, wilcoxPass_ref)
    rank = c(test_df$measure, ref_df$measure)
    wilcoxPass$measure = rank
  }
  
  matrix = wilcoxPass
  
  for (i in (1:(ncol(matrix)-1))) {
    WILCOX = ggplot(matrix, aes(x=factor(measure), y=matrix[,i])) +
      geom_violin() +
      geom_jitter(aes(colour=measure), width = 0.15) +
      geom_smooth(aes(x=as.factor(measure), y=matrix[,i]), method = "loess") +
      ylab(colnames(matrix[i])) +
      ggtitle(paste(colnames(matrix[i]))) +
      theme(axis.text.x = element_text(angle = 0),
            axis.title.y = element_text(size = 8, angle = 90))
    
    ggsave(WILCOX, file=paste0(getwd(),"/Figures/Wilcox/", colnames(matrix[i]), "_", file_name, ".png"))
    
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
    
    
    ggsave(violin, file=paste0(getwd(),"/Figures/Heatmaps/Maha/", colnames(matrix[i]), "_", file_name, ".png"))
  }
}

##Unsupervised clustering

unsupervised_clustering = function(data, clinical.data, n_clusters, file_name){
  
  #Perform k-means algorithm: rows->features and columns->samples (kmeans affect rows)
  set.seed(240) 
  kmeans.samples <- kmeans(data, centers = n_clusters, nstart = 20)
  #Join cluster information as annotations  
  clinical.data$cluster <- kmeans.samples$cluster
  
  clusters = data.frame(cbind(rownames(clinical.data), clinical.data$cluster))
  rownames(clusters) = clusters$X1
  clusters$X1 = NULL
  colnames(clusters)[1] = "Cluster"
  
  colnames(clinical.data)[ncol(clinical.data)] = paste0(file_name, "_Cluster")
  
  return(clinical.data)
  
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

## Remove highly correlate features

remove_high_corr = function(df){
  #rows = features and columns = samples
  corr_pearson <- as.data.frame(rcorr((as.matrix(t(df))), type = "pearson")$r)
  for(i in 1:ncol(corr_pearson)){
    idx <- which((corr_pearson[,i]>0.9)&(corr_pearson[,i]!=1))
    if(is_empty(idx)==TRUE){
      corr_pearson <- corr_pearson
    }else{
      corr_pearson <- corr_pearson[-idx,-idx]
    }
    N <- ncol(corr_pearson)
    if (i == N){
      break
    }
  }
  idx <- which(colnames(corr_pearson)%in%rownames(df))
  df <- as.data.frame(df[idx,])
  return(df)
}

compute_subgroups = function(data, corr_thres, corr_if, iteration, file_name){
  #data: features columns and samples rows
  for (iter in 1:iteration) {
    if(iter == 1){
      corr_df <- correlation(data.matrix(data))
      vec = colnames(data)
      indice = 1
      subgroup_init = list()
      vector_init = c() 
      data_sub = c()
      
      while(length(vec)>0){
        if(vec[1] %in% corr_df$measure1){
          tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature
          tab = tab[tab$r>corr_thres,] #Select corr above the threshold
          if(nrow(tab)!=0){
            subgroup_init[[indice]] = c(vec[1], tab$measure2) #Create subgroup
            corr_df = corr_df[-which(corr_df$measure1 %in% subgroup_init[[indice]]),] #Remove varible already subgroupped from corr
            vec = vec[-which(vec%in%subgroup_init[[indice]])] #Remove feature already subgroupped 
            indice = indice + 1
          }else{ #Condition when there is no correlation above the first threshold
            tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature
            tab = tab[tab$r>corr_if,]  #Reduce corr threshold to corr_if
            if(nrow(tab)!=0){
              subgroup_init[[indice]] = c(vec[1], tab$measure2) #Create subgroup
              corr_df = corr_df[-which(corr_df$measure1 %in% subgroup_init[[indice]]),] #Remove varible already subgroupped from corr
              vec = vec[-which(vec%in%subgroup_init[[indice]])] #Remove feature already subgroupped
              indice = indice + 1
            }else{ #Condition when there is no correlation above the second threshold 
              corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable (outlier cause it doesnt subgroup with anything else)
              vector_init = c(vector_init, vec[1]) #Save variable for further analysis
              vec = vec[-1]
              indice = indice 
            }
          }
        }else{ #If feature is not in corr_df it means that there is no any significant correlation against it and other features 
          vector_init = c(vector_init, vec[1]) #Save variable for further analysis
          vec = vec[-1] 
          indice = indice 
        }
      }
      
      if(length(subgroup_init)==0){
        return("No subgroups above 0.5")
      }
      
      for (i in 1:length(subgroup_init)) { #Name Subgroups
        names(subgroup_init)[i] = paste0("Subgroup", i)
      }
      
      for(i in 1:length(subgroup_init)){ #Create data frame with features subgroupped
        sub = data.frame(data[,colnames(data)%in%subgroup_init[[i]]]) #Map features that are inside each group from input
        sub$average = rowMeans(sub) #Compute score of subgroup (average of features from the subgroup) 
        data_sub = data.frame(cbind(data_sub, sub$average)) #Save scores in the new data frame
        colnames(data_sub)[i] = paste0(file_name, "_Subgroup_",i)
      }
      
      rownames(data_sub) = rownames(data) #List of patients
      data = data_sub
      
      if(ncol(data)==1){ #If data has only one subgroup
        lis <- vector(mode = "list", length = 1)
        lis[[1]] = subgroup_init[[1]] #Save names of features inside the subgroup
        lis[[1]] = unique(unlist(lis[[1]])) #Take only unique features inside the subgroup
        names(lis) = colnames(data)
        return(list(data, lis)) #End function
      }
    }
    
    #Next iterations between subgroups (Same process but not reducing threshold)
    
    else{ 
      corr_df <- correlation(data.matrix(data))
      indice = 1
      subgroup = list()
      vector = c() 
      x = c()
      data_sub = c()
      vec = colnames(data)
      
      if(nrow(corr_df) == 0){ #Enter this point when subgroup is only one feature (impossible to make correlation --> end function)
        iter == iteration
      }
      
      while(length(vec)>0){
        if(vec[1] %in% corr_df$measure1){
          tab = corr_df[corr_df$measure1 == vec[1],]
          tab = tab[tab$r>corr_if,]
          if(nrow(tab)!=0){
            subgroup[[indice]] = c(vec[1], tab$measure2)
            corr_df = corr_df[-which(corr_df$measure1 %in% subgroup[[indice]]),]
            vec = vec[-which(vec%in%subgroup[[indice]])]
            indice = indice + 1
          }else{
            vector = c(vector, vec[1])
            corr_df = corr_df[-which(corr_df$measure1 == vec[1]),]
            vec = vec[-1]
            indice = indice 
          }
        }
        else{
          vector = c(vector, vec[1])
          vec = vec[-1]
          indice = indice 
        }
      }
      
      if(length(subgroup)>0){
        for(i in 1:length(subgroup)){
          sub = data.frame(data[,colnames(data)%in%subgroup[[i]]])
          sub$average = rowMeans(sub)
          data_sub = data.frame(cbind(data_sub, sub$average))
          colnames(data_sub)[i] = paste0(file_name, "_Subgroup_",str_c(sub(paste0(file_name, "_Subgroup_"), "", unlist(subgroup[[i]])), collapse = "_"))
        }
        rownames(data_sub) = rownames(data)
        if(length(vector) > 0){
          data_sub = cbind(data[,colnames(data)%in%vector], data_sub)
          if(length(vector)==1){
            colnames(data_sub)[1] = vector
          }
        }
        data = data_sub
      }else{ #Reach this point when there is no more to subgroup
        iter = iteration
      }
    }
    
    #Return data when iterations are done 
    if(iter == iteration){
      lis <- vector(mode = "list", length = ncol(data))
      for (i in 1:ncol(data)) {
        idx = unlist(strsplit(colnames(data)[i], split=paste0(file_name,'_Subgroup_')))[2]
        num = unlist(strsplit(idx, split='_'))
        for (j in 1:length(num)){
          lis[[i]] = c(lis[[i]], subgroup_init[[as.numeric(num[j])]]) 
        }
        lis[[i]] = unique(unlist(lis[[i]]))
      }
      names(lis) = colnames(data)
      saveRDS(lis, file=paste0(getwd(), "/Output/Subgroups/", file_name, "_subgroups.rds"))
      return(list(data, lis, vector_init))
    }
    
  }
}

compute_cell.types = function(data, dataset){
  #columns:features and rows:samples
  ################################# B cells
  Bcells <- data[,grep("B", colnames(data))]
  if(length(grep("M2", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("M2", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("ancer", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("ancer", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("asophils", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("asophils", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("CD4", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("CD4", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("Mono", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("Mono", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("Treg", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("Treg", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("Neu", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("Neu", colnames(Bcells))]}else{Bcells <- Bcells}
  if(length(grep("M1", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("M1", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("CD8", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("CD8", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("M0", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("M0", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("NK", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("NK", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("ndothelial", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("ndothelial", colnames(Bcells))]}else{Bcells <- Bcells}      
  if(length(grep("CAF", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("CAF", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("acrophage", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("acrophage", colnames(Bcells))]}else{Bcells <- Bcells}        
  if(length(grep("endritic", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("endritic", colnames(Bcells))]}else{Bcells <- Bcells}      
  if(length(grep("alignant", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("alignant", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("Mast", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("Mast", colnames(Bcells))]}else{Bcells <- Bcells}    
  if(length(grep("yocyte", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("yocyte", colnames(Bcells))]}else{Bcells <- Bcells}  
  if(length(grep("ibroblast", colnames(Bcells)))>0){Bcells <- Bcells[,-grep("ibroblast", colnames(Bcells))]}else{Bcells <- Bcells} 
  x = correlation(Bcells)
  if(length(which(x$r == 1))>0){
    Bcells = Bcells[,-which(colnames(Bcells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(Bcells, paste0(getwd(),"/Output/CellTypes/B_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/B_", dataset, ".png"), width = 1600, height = 1000, res=100)
  # plot_correlation(t(Bcells), paste0("B cells ", dataset))
  # dev.off()
  
  ################################# CD4 cells
  CD4cells <- data[,grep("CD4", colnames(data))]
  x = correlation(CD4cells)
  if(length(which(x$r == 1))>0){
    CD4cells = CD4cells[,-which(colnames(CD4cells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(CD4cells, paste0(getwd(),"/Output/CellTypes/CD4_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/CD4_", dataset, ".png"), width = 1800, height = 1000, res=100)
  # plot_correlation(t(CD4cells), paste0("CD4 cells ", dataset))
  # dev.off()
  
  ################################# CD8 cells
  CD8cells <- data[,grep("CD8", colnames(data))]
  x = correlation(CD8cells)
  if(length(which(x$r == 1))>0){
    CD8cells = CD8cells[,-which(colnames(CD8cells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(CD8cells, paste0(getwd(),"/Output/CellTypes/CD8_", dataset, ".txt"))
  # 
  # png(paste0(getwd(),"/Figures/Corrplots/CD8_", dataset, ".png"), width = 1800, height = 1000, res=100)
  # plot_correlation(t(CD8cells), paste0("CD8 cells ", dataset))
  # dev.off()
  # 
  ################################# M1 cells
  M1cells <- data[,grep("M1", colnames(data))]
  x = correlation(M1cells)
  if(length(which(x$r == 1))>0){
    M1cells = M1cells[,-which(colnames(M1cells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(M1cells, paste0(getwd(),"/Output/CellTypes/M1_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/M1_", dataset, ".png"), width = 1200, height = 600, res=100)
  # plot_correlation(t(M1cells), paste0("M1 cells ", dataset))
  # dev.off()
  
  ################################# M2 cells
  M2cells <- data[,grep("M2", colnames(data))]
  if(length(grep("LM22", colnames(M2cells)))>0){M2cells <- M2cells[,-grep("LM22", colnames(M2cells))]}else{M2cells <- M2cells} 
  x = correlation(M2cells)
  if(length(which(x$r == 1))>0){
    M2cells = M2cells[,-which(colnames(M2cells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(M2cells, paste0(getwd(),"/Output/CellTypes/M2_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/M2_", dataset, ".png"), width = 800, height = 600, res=100)
  # plot_correlation(t(M2cells), paste0("M2 cells ", dataset))
  # dev.off()
  
  ################################# NK cells
  NKcells <- data[,grep("NK", colnames(data))]
  x = correlation(NKcells)
  if(length(which(x$r == 1))>0){
    NKcells = NKcells[,-which(colnames(NKcells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(NKcells, paste0(getwd(),"/Output/CellTypes/NK_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/NK_", dataset, ".png"), width = 1500, height = 1000, res=100)
  # plot_correlation(t(NKcells), paste0("NK cells ", dataset))
  # dev.off()
  
  ################################# cancer cells
  cancer <- data[,grep("ancer", colnames(data))]
  x = correlation(cancer)
  if(length(which(x$r == 1))>0){
    cancer = cancer[,-which(colnames(cancer)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(cancer, paste0(getwd(),"/Output/CellTypes/Cancer_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/Cancer_", dataset, ".png"), width = 800, height = 500, res=100)
  # plot_correlation(t(cancer), paste0("Cancer cells ", dataset))
  # dev.off()
  
  ################################# Neutrophils cells
  Neutrophils <- data[,grep("Neutro", colnames(data))]
  x = correlation(Neutrophils)
  if(length(which(x$r == 1))>0){
    Neutrophils = Neutrophils[,-which(colnames(Neutrophils)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(Neutrophils, paste0(getwd(),"/Output/CellTypes/Neutrophils_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/Neutrophils_", dataset, ".png"), width = 1000, height = 800, res=100)
  # plot_correlation(t(Neutrophils), paste0("Neutrophils cells ", dataset))
  # dev.off()
  
  ################################# M0cells cells
  M0cells <- data[,grep("M0", colnames(data))]
  x = correlation(M0cells)
  if(length(which(x$r == 1))>0){
    M0cells = M0cells[,-which(colnames(M0cells)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(M0cells, paste0(getwd(),"/Output/CellTypes/M0_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/M0_", dataset, ".png"), width = 800, height = 500, res=100)
  # plot_correlation(t(M0cells), paste0("M0 cells ", dataset))
  # dev.off()
  
  ################################# Monocytes cells 
  Mono <- data[,grep("Mono", colnames(data))]
  x = correlation(Mono)
  if(length(which(x$r == 1))>0){
    Mono = Mono[,-which(colnames(Mono)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(Mono, paste0(getwd(),"/Output/CellTypes/Mono_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/Mono_", dataset, ".png"), width = 1500, height = 800, res=100)
  # plot_correlation(t(Mono), paste0("Monocytes cells ", dataset))
  # dev.off()
  
  ################################# CAF cells
  CAF <- data[,grep("CAF", colnames(data))]
  x = correlation(CAF)
  if(length(which(x$r == 1))>0){
    CAF = CAF[,-which(colnames(CAF)%in%x$measure1[which(x$r == 1)])]
  }
  
  write.csv(CAF, paste0(getwd(),"/Output/CellTypes/CAF_", dataset, ".txt"))

  # png(paste0(getwd(),"/Figures/Corrplots/CAF_", dataset, ".png"), width = 1000, height = 800, res=100)
  # plot_correlation(t(CAF), paste0("CAF cells ", dataset))
  # dev.off()
  
  ################################# Dendritic cells
  Dendritic <- data[,grep("endritic", colnames(data))]
  x = correlation(Dendritic)
  if(length(which(x$r == 1))>0){
    Dendritic = Dendritic[,-which(colnames(Dendritic)%in%x$measure1[which(x$r == 1)])]
  }
  
  write.csv(Dendritic, paste0(getwd(),"/Output/CellTypes/Dendritic_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/Dendritic_", dataset, ".png"), width = 1000, height = 800, res=100)
  # plot_correlation(t(Dendritic), paste0("Dendritic cells ", dataset))
  # dev.off()
  
  ################################# Endothelial cells
  Endothelial <- data[,grep("Endothelial", colnames(data))]
  x = correlation(Endothelial)
  if(length(which(x$r == 1))>0){
    Endothelial = Endothelial[,-which(colnames(Endothelial)%in%x$measure1[which(x$r == 1)])]
  }
  write.csv(Endothelial, paste0(getwd(),"/Output/CellTypes/Endothelial_", dataset, ".txt"))

  # png(paste0(getwd(),"/Figures/Corrplots/Endothelial_", dataset, ".png"), width = 1000, height = 800, res=100)
  # plot_correlation(t(Endothelial), paste0("Endothelial cells ", dataset))
  # dev.off()
  
  ################################# Macrophages cells
  Macrophages <- data[,grep("Macrophages", colnames(data))]
  if(length(grep("M2", colnames(Macrophages)))>0){Macrophages <- Macrophages[,-grep("M2", colnames(Macrophages))]}else{Macrophages <- Macrophages}
  if(length(grep("M1", colnames(Macrophages)))>0){Macrophages <- Macrophages[,-grep("M1", colnames(Macrophages))]}else{Macrophages <- Macrophages} 
  x = correlation(Macrophages)
  
  if(length(which(x$r == 1))>0){
    Macrophages = Macrophages[,-which(colnames(Macrophages)%in%x$measure1[which(x$r == 1)])]
  }
  
  write.csv(Macrophages, paste0(getwd(),"/Output/CellTypes/Macrophages_", dataset, ".txt"))
  
  # png(paste0(getwd(),"/Figures/Corrplots/Macrophages_", dataset, ".png"), width = 1000, height = 800, res=100)
  # plot_correlation(t(Macrophages), paste0("Macrophages cells ", dataset))
  # dev.off()
  
  cell_types = list(Bcells, CD4cells, CD8cells, M1cells, M2cells, NKcells, cancer, Neutrophils, M0cells, Mono, CAF, Dendritic, Endothelial, Macrophages)
  names(cell_types) = c("Bcells", "CD4cells", "CD8cells", "M1cells", "M2cells", "NKcells", "Cancer", "Neutrophils", "M0cells", "Monocytes", "CAF", "Dendritic", "Endothelial", "Macrophages")
  
  return(cell_types)
  
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
  
  ploadings <- plotloadings(p, components = getComponents(p, c(1:n_components)),
                            rangeRetain = 0.01,
                            labSize = 4.0,
                            title = 'Loadings plot',
                            caption = 'Top 1% variables',
                            shape = 24,
                            col = c('limegreen', 'black', 'red3'),
                            drawConnectors = TRUE)
  
  ############################################################################################
  
  png(paste0(getwd(),"/Figures/PCA/ploadings_", file_name, ".png"), width = 3500, height = 2500, res = 150)
  print(ploadings)
  dev.off()
  
  png(paste0(getwd(),"/Figures/PCA/peigencor_", file_name, ".png"), width = 2500, height = 1500, res = 150)
  print(peigencor)
  dev.off()
}

compute_wgcna = function(data, coldata, file_name){
  
  #data --> rows samples and genes as columns
  
  ################################################################################################################
  #Choose Soft Threshold (power) parameter
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
  
  png(paste0(getwd(),"/Figures/WGCNA/SoftThreshold_", file_name, ".png"), width = 800, height = 700, res=100)
  print(plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence")))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
  abline(h=0.90,col="red")
  dev.off()
  
  ################################################################################################################  
  #Compute WGCNA modules
  cor <- WGCNA::cor        
  picked_power = readline("Enter soft threshold: ")
  picked_power = as.integer(picked_power);
  netwk <- blockwiseModules(data,               
                            power = picked_power,              
                            networkType = "signed",
                            deepSplit = 2,
                            pamRespectsDendro = F,
                            minModuleSize = 10, 
                            reassignThreshold = 0,
                            mergeCutHeight = 0.15,
                            numericLabels = T,
                            verbose = 3)
  mergedColors = labels2colors(netwk$colors)
  module_df <- data.frame(gene_id = names(netwk$colors),
                          colors = labels2colors(netwk$colors))
  
  ################################################################################################################    
  # Get Module Eigengenes per cluster
  nGenes = ncol(data)
  nSamples = nrow(data)
  
  MEs0 <- moduleEigengenes(data, mergedColors)$eigengenes
  MEs0 <- orderMEs(MEs0)
  hubgenes = chooseTopHubInEachModule(data, mergedColors)
  module_order = names(MEs0) %>% gsub("ME","", .)
  MEs0$treatment = row.names(MEs0)
  mME = MEs0 %>%
    pivot_longer(-treatment) %>%
    mutate(
      name = gsub("ME", "", name),
      name = factor(name, levels = module_order)
    )
  p = mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0,
      limit = c(-1,1)) +
    theme(axis.text.x = element_text(angle=90)) +
    labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
  
  png(paste0(getwd(),"/Figures/WGCNA/Module_Eigengenes_", file_name, ".png"), width = 1200, height = 600, res=150)
  print(p)
  dev.off()
  
  ################################################################################################################  
  #Relate module eigengenes with features
  MEs0$treatment = NULL
  moduleTraitCor = cor(MEs0, coldata, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  png(paste0(getwd(),"/Figures/WGCNA/Module_trait_relationships_", file_name, ".png"), width = 3500, height = 1500, res=150)
  par(mar = c(15, 10, 4, 4)) #bottom, left, top, right
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(coldata),
                 yLabels = names(MEs0),
                 ySymbols = names(MEs0),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  ################################################################################################################
  #ORA enrichment analysis
  for (i in unique(module_df$colors)) {
    entrezID = gconvert(module_df$gene_id[module_df$colors==i], organism = "hsapiens", target = "ENTREZGENE_ACC", filter_na = F, mthreshold = 1)
    #GO analysis
    goBP <- as.data.frame(enrichGO(entrezID$target, OrgDb="org.Hs.eg.db", pvalueCutoff=0.05, ont = "BP"))
    write.csv(goBP,paste0(getwd(),"/Output/WGCNA/GO_BP_", i, "_", file_name, ".csv"))
    #KEGG analysis
    kegg_results = as.data.frame(enrichKEGG(gene = entrezID$target, organism = "hsa", pvalueCutoff = 0.05))
    write.csv(kegg_results,paste0(getwd(),"/Output/WGCNA/KEGG_", i, "_", file_name, ".csv"))
  }
  
  return(list(module_df, hubgenes))
  
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


##Functions scores from https://github.com/olapuentesantana/easier_manuscript

#' Compute Expanded Immune signature
#'
#' \code{compute_ayersEI} computes Expanded Immune signature score as the arithmetic mean of genes included
#' in the Expanded Immune signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Expanded Immune signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Ayers_expIS <- function(RNA.tpm){
  # Literature genes
  Ayers_expIS.read <- c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                        "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA","HLA-DRA", "LAG3", "IDO1", "TAGAP")
  match_Ayers_expIS.genes <- match(Ayers_expIS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Ayers_expIS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Ayers_expIS.read[!Ayers_expIS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Ayers_expIS.genes <- stats::na.omit(match_Ayers_expIS.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm  <- log2.RNA.tpm[match_Ayers_expIS.genes, ]
  
  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)
  
  message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}


#' Compute cytolytic activity score
#'
#' \code{compute_CYT} computes cytolytic activity score as the geometric mean of immune cytolytic genes
#' (Rooney et al., 2015).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.CYT <- function(RNA.tpm){
  
  # Literature genes
  CYT.read <- c("GZMA", "PRF1")
  match_CYT.genes <- match(CYT.read, rownames(RNA.tpm))
  
  if (anyNA(match_CYT.genes)){
    warning(paste0("differenty named or missing signature genes : \n", paste(CYT.read[!CYT.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_CYT.genes <- stats::na.omit(match_CYT.genes)
  }
  
  # Subset RNA.tpm
  subset_RNA.tpm <- RNA.tpm[match_CYT.genes, ]
  
  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))
  
  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}

#' Compute Davoli immune signature
#'
#' \code{compute_davoliIS} computes Davoli immune signature as the arithmetic mean of cytotoxic
#' immune infiltrate signature genes, after rank normalization (Davoli et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Davoli immune signature
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Davoli_IS <- function(RNA.tpm){
  
  # Literature genes
  Davoli_IS.read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  match_Davoli_IS.genes <- match(Davoli_IS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Davoli_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Davoli_IS.read[!Davoli_IS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Davoli_IS.genes <- stats::na.omit(match_Davoli_IS.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_Davoli_IS.genes, ]
  
  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)
  
  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)
  
  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)
  
  message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}

#' Compute IFNy signature score
#'
#' \code{compute_ayersIFNy} computes IFNy signature score as the arithmetic mean of genes included
#' in the IFN-γ signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IFNy signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.IFNy <- function(RNA.tpm){
  
  # Literature genes
  IFNy.read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA")
  match_IFNy.genes <- match(IFNy.read, rownames(RNA.tpm))
  
  if (anyNA(match_IFNy.genes)){
    warning(paste0("differenty named or missing signature genes : \n", IFNy.read[!IFNy.read %in% rownames(RNA.tpm)]))
    match_IFNy.genes <- stats::na.omit(match_IFNy.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)
  
  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_IFNy.genes, ]
  
  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)
  
  message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}


#' Compute MSI score
#'
#' \code{compute_MSI} computes MSI score by applying logical comparison of MSI-related gene pairs
#' (Fu et al., 2019).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=MSI score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.MSI <- function(RNA.tpm){
  
  # Literature genes: * (CCRN4L in tcga, NOCT approved symbol)
  MSI.basis <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE",
                                     "PRPF39","NOCT","AMFR"),
                          Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3",
                                     "ATF6","GRM8","DUSP18"))
  MSI.read <- unique(as.vector(as.matrix(MSI.basis))) # 20 genes
  
  # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
  if (any(rownames(RNA.tpm) %in% "CCRN4L")){
    cat("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "CCRN4L"] <- "NOCT"
  }
  
  # Subset RNA.tpm
  match_F_1 <- match(as.character(MSI.basis[,1]), rownames(RNA.tpm))
  match_F_2 <- match(as.character(MSI.basis[,2]), rownames(RNA.tpm))
  
  if (anyNA(c(match_F_1,match_F_2))){
    warning(c("differenty named or missing signature genes : \n", paste(MSI.read[!MSI.read %in% rownames(RNA.tpm)], collapse = "\n")))
  }
  
  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  F_pair_expr_B <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  MSI.matrix <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm)) ; colnames(MSI.matrix) <- colnames(RNA.tpm)
  remove_pairs <- vector("list", length = ncol(RNA.tpm)) ; names(remove_pairs) <- colnames(RNA.tpm)
  score <- vector("numeric", length = ncol(RNA.tpm)) ; names(score) <- colnames(RNA.tpm)
  
  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))
  
  # Calculation:
  F_pair_expr_A <- log2.RNA.tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA.tpm[match_F_2, ]
  
  if(anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }
  
  MSI.matrix <- F_pair_expr_A > F_pair_expr_B
  if(anyNA(MSI.matrix)){
    score <- colSums(MSI.matrix, na.rm = TRUE)
    score <- (score * nrow(MSI.matrix)) / (nrow(MSI.matrix) - length(remove_pairs))
  }else{
    score <- colSums(MSI.matrix)
  }
  
  message("MSI score computed")
  return(data.frame(MSI = score, check.names = FALSE))
}

#' Compute Roh immune score
#'
#' \code{compute_rohIS} computes Roh immune score as the geometric-mean of immune score genes
#' (Roh et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Roh immune score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Roh_IS <- function(RNA.tpm){
  
  # Literature genes
  Roh_IS.read <- c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                   "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
                   "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
                   "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
                   "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1")
  match_Roh_IS.genes <- match(Roh_IS.read, rownames(RNA.tpm))
  
  if (anyNA(match_Roh_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Roh_IS.read[!Roh_IS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Roh_IS.genes <- stats::na.omit(match_Roh_IS.genes)
  }
  
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[match_Roh_IS.genes, ]
  
  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01
  
  # Pseudocount of 1 for genes with 0 expr
  if(any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] + 1
  
  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))
  
  message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}

#' Compute tertiary lymphoid structures signature
#'
#' \code{compute_TLS} computes TLS signature as the geometric-mean of TLS signature genes
#' (Cabrita et al., 2020).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=TLS signature
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.TLS <- function(RNA.tpm){
  
  # Literature genes
  TLS.read <- c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS")
  match_TLS.read <- match(TLS.read, rownames(RNA.tpm))
  
  if (anyNA(match_TLS.read)){
    warning(c("differenty named or missing signature genes : \n", paste(TLS.read[!TLS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_TLS.read <- stats::na.omit(match_TLS.read)
  }
  
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[match_TLS.read, ]
  
  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log2(X + 1))))
  
  message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}

#' Compute T cell-inflamed signature score
#'
#' \code{compute_ayersTcellInfl} computes T cell-inflamed signature score by taking a weighted sum of
#'  the housekeeping normalized values of the T cell-inflamed signature genes
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=T cell-inflamed signature score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.Tcell_inflamed <- function(RNA.tpm){
  
  # Literature genes
  Tcell_inflamed.read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                           "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  Housekeeping.read <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2", "UBB", "TBP", "SDHA") # C14orf102 = NRDE2
  weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                        CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                        PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767, check.names = FALSE)
  
  # Some genes might have other name: case for "C14orf102", it's called "NRDE2", be careful
  if (any(rownames(RNA.tpm) %in% "C14orf102")){
    cat("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "C14orf102"] <- "NRDE2"
  }
  
  match_genes.housekeeping <- match(Housekeeping.read, rownames(RNA.tpm))
  match_genes.predictors <- match(Tcell_inflamed.read, rownames(RNA.tpm))
  
  if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
    tmp <- c(Tcell_inflamed.read, Housekeeping.read)
    warning(c("differenty named or missing signature genes : \n", paste(tmp[!tmp %in% rownames(RNA.tpm)], collapse = "\n")))
    match_genes.housekeeping <- stats::na.omit(match_genes.housekeeping)
    match_genes.predictors <- stats::na.omit(match_genes.predictors)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)
  
  # Subset log2.RNA.tpm
  ## housekeeping
  log2.RNA.tpm.housekeeping <- log2.RNA.tpm[match_genes.housekeeping, ]
  ## predictors
  log2.RNA.tpm.predictors <- log2.RNA.tpm[match_genes.predictors, ]
  weights <- weights[,rownames(log2.RNA.tpm.predictors)]
  
  # Housekeeping normalization
  average.log2.RNA.tpm.housekeeping <- apply(log2.RNA.tpm.housekeeping, 2, mean)
  log2.RNA.tpm.predictors.norm <- sweep(log2.RNA.tpm.predictors, 2, average.log2.RNA.tpm.housekeeping, FUN = "-")
  
  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.RNA.tpm.predictors.norm), colnames(weights))
  score <- t(log2.RNA.tpm.predictors.norm[tidy,]) %*% t(weights)
  
  message("Tcell_inflamed score computed")
  return(data.frame( Tcell_inflamed = score, check.names = FALSE))
}

#' Compute chemokine score
#'
#' \code{compute_chemokine} computes chemoine score as the PC1 score that results from applying PCA
#' to z-score expression of 12 chemokine genes (Messina et al., 2012).
#'
#' @importFrom stats na.omit prcomp
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
#-------------------------------------------------------------------------------------------------------------

compute.chemokines <- function(RNA.tpm){
  RNA.tpm = vsd_df
  # Literature genes
  chemokines.read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                       "CXCL9", "CXCL10", "CXCL11", "CXCL13")
  match_chemokines.genes <- match(chemokines.read, rownames(RNA.tpm))
  
  if (anyNA(match_chemokines.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(chemokines.read[!chemokines.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_chemokines.genes <- stats::na.omit(match_chemokines.genes)
  }
  
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)
  
  # Subset gene_expr
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_chemokines.genes, ]
  
  sub_log2.RNA.tpm = sub_log2.RNA.tpm[apply(sub_log2.RNA.tpm,1,sum)>0,]
  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine.pca <- stats::prcomp(t(sub_log2.RNA.tpm), center = TRUE, scale = TRUE)
  score <- chemokine.pca$x[, 1]
  
  message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}


#-------------------------------------------------------------------------------------------------------------


#' Compute the expression of the immune checkpoints genes
#'
#' \code{computation_ICB_genes} computes the scores for the immune checkpoint genes.
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#'
#' @return List with the expression of the immune checkpoint genes
#'
#-------------------------------------------------------------------------------------------------------

compute_ICB_genes <- function(RNA.tpm){
  
  # Extract position genes for GZMA and PRF1
  tmp <- match(c("CD274","CTLA4","PDCD1"), rownames(RNA.tpm))
  
  # PDL-1 calculation
  PDL1_expr = RNA.tpm[tmp[1],] ; rownames(PDL1_expr) <- "PDL1"
  
  # CTLA-4 calculation
  CTLA4_expr = RNA.tpm[tmp[2],] ; rownames(CTLA4_expr) <- "CTLA4"
  
  # PD-1 calculation
  PD1_expr = RNA.tpm[tmp[3],] ; rownames(PD1_expr) <- "PD1"
  
  ICB_genes_expr <- list(PDL1 = PDL1_expr , CTLA4 = CTLA4_expr , PD1 =PD1_expr)
  message("ICB genes expression computed")
  return(ICB_genes_expr)
}



#' \code{computation_gold_standards} computes the scores for the gold standards required by the user
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#' @param list_gold_standards string with gold standards names
#' @param cancertype string character
#'
#' @return List with the scores of all the gold standards specified.
#'
#-------------------------------------------------------------------------------------------------------
# Input: Transcriptomics data as tpm values
# A list with the names of the scores to be computed has to be provided
# Output: Gold standards scores
#-------------------------------------------------------------------------------------------------------
compute_gold_standards <- function(RNA.tpm,
                                   list_gold_standards, dataset){
  
  
  gold.standards <- sapply(list_gold_standards, function(X){
    
    if ("CYT" == X) {
      
      # calculate Cytolytic activity #
      CYT <- t(compute.CYT(RNA.tpm))
      return(list(CYT))
      
    }else if("Roh_IS" == X) {
      
      # calculate roh immune signature #
      Roh_IS <- t(compute.Roh_IS(RNA.tpm))
      return(list(Roh_IS))
      
    }else if("chemokines" == X) {
      
      # calculate chemokine signature #
      chemokines <- t(compute.chemokines(RNA.tpm))
      return(list(chemokines))
      
    }else if("Davoli_IS" == X) {
      
      # calculate davoli cytotoxic immune signature #
      Davoli_IS <- t(compute.Davoli_IS(RNA.tpm))
      return(list(Davoli_IS))
      
    }else if("IFNy" == X) {
      
      # calculate ayers IFNy #
      IFNy <- t(compute.IFNy(RNA.tpm))
      return(list(IFNy))
      
    }else if("Ayers_expIS" == X) {
      
      # calculate ayers expanded immune signature #
      Ayers_expIS <- t(compute.Ayers_expIS(RNA.tpm))
      return(list(Ayers_expIS))
      
    }else if("Tcell_inflamed" == X) {
      
      # calculate ayers T cell inflamed signature #
      Tcell_inflamed <- t(compute.Tcell_inflamed(RNA.tpm))
      return(list(Tcell_inflamed))
      
    }else if("MSI" == X) {
      
      # calculate MSI signature #
      MSI <- t(compute.MSI(RNA.tpm))
      return(list(MSI))
      
    }else if("RIR" == X) {
      
      # calculate MSI signature #
      RIR <- t(compute.RIR(RNA.tpm))
      return(list(RIR))
      
    }else if("TLS" == X) {
      
      # calculate MSI signature #
      TLS <- t(compute.TLS(RNA.tpm))
      return(list(TLS))
      
    }
    
  })
  
  gold.standards = lapply(gold.standards,function(x){t(x) %>% as.data.frame() %>% rownames_to_column("condition")})
  gold.standards = Reduce(function(x,y) {left_join(x,y ,by="condition")}, gold.standards)
  rownames(gold.standards) = gold.standards$condition 
  gold.standards$condition = NULL
  
  write.csv(gold.standards,paste0(getwd(),"/Output/Immunoscores/ImmuneScores_",dataset,".csv"), row.names = T)
  
  return(gold.standards)
  
}

gsva_analysis = function(Counts, category){
  
  ######################################################################Input Gene matrix
  EntrezID = gconvert(row.names(Counts), organism = "hsapiens", target = "ENTREZGENE_ACC", filter_na = F, mthreshold = 1)
  Counts$EntrezID = EntrezID$target
  Counts = na.omit(Counts)
  rownames(Counts) = NULL
  
  gene_means <- rowMeans(Counts %>% dplyr::select(-EntrezID))
  Counts <- Counts %>%
    dplyr::mutate(gene_means) %>%
    dplyr::select(EntrezID, gene_means, dplyr::everything()) %>%
    dplyr::arrange(dplyr::desc(gene_means)) %>%
    dplyr::distinct(EntrezID, .keep_all = TRUE) %>%
    dplyr::select(-gene_means) %>%
    tibble::column_to_rownames("EntrezID") %>%
    as.matrix()
  #####################################################################Input Gene sets
  gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens", 
    category = category 
  )
  
  gene_list <- split(gene_sets$entrez_gene, gene_sets$gs_name)
  
  #####################################################################Run GSVA
  gsva_results <- data.frame(gsva(Counts, gene_list, method = "gsva", kcdf = "Gaussian", min.sz = 15, max.sz = 500, mx.diff = TRUE, verbose = FALSE))
  return(gsva_results)
}
