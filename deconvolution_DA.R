
#Deconvolution function analysis
source(paste0(getwd(),"/environment_set.R"))
libraries_set()

#Input file
deconv <- read.csv("~/Downloads/all_deconvolutions_Counts_LPVan_TPM.txt_.csv", row.names=1)

#Look for outliers
SampleDist = dist(deconv)
SampleDistMatrix = as.matrix(SampleDist)
colors = colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(SampleDistMatrix, clustering_distance_rows = SampleDist, clustering_distance_cols = SampleDist, col = colors)
idx = c('R4163_YZ_8', 'R4163_YZ_26')
deconv = deconv[-which(rownames(deconv)%in%idx),]

#Split by method-signature
# idx = grep(paste0(methods[1], "_", signatures[1]), colnames(deconv))
# dt = View(deconv[, idx])

#Filtering equal method-signatures (corr=1)
deconv = data.matrix(deconv[,-which(duplicated(t(deconv)))])

#Remove columns with low variance
deconv = data.frame(deconv[,-which(colVars(deconv)<=summary(colVars(deconv))[2])])
#Remove columns with more than 70% of zeros
i <- colSums(deconv == 0, na.rm=TRUE) < round(0.7*nrow(deconv))
deconv <- data.frame(deconv[, i, drop=FALSE])

#Map methods, signatures and cells 
x = strsplit(colnames(deconv), "_")
methods = c()
signatures = c()
celltypes = c()
for(i in 1:length(x)){
  methods = c(methods, x[[i]][1])
  signatures = c(signatures, x[[i]][2])
  celltypes = c(celltypes, x[[i]][3])
}

methods = unique(methods)
signatures = unique(signatures)

#Split deconvolution matrix by cell type
compute_cell.types = function(data, celltypes,dataset){
  data = data.frame(data)
  #celltypes = celltypes
  B = grep("B", celltypes)
  B = data[, B, drop = FALSE]
  B.naive <- grep("NAIVE", colnames(B))
  B.naive <- B[, B.naive, drop = FALSE]
  B.memory = grep("MEMORY", colnames(B))
  B.memory <- B[, B.memory, drop = FALSE]
  B = B[,-which(colnames(B)%in%c(colnames(B.naive), colnames(B.memory)))]
  
  Macrophages = data[,grep("MACROPHAGES", celltypes)]
  M0 = data[,grep("M0", celltypes)]
  M1 = data[,grep("M1", celltypes)]
  M2 = data[,grep("M2", celltypes)]
  Monocytes = data[,grepl(paste(c("MONO", "MONOCYTES"), collapse="|"), celltypes)]
  Neutrophils = data[,grepl(paste(c("NEU", "NEUTROPHILS"), collapse="|"), celltypes)]
  
  NK = data[,grep("NK", celltypes)]
  NK.activated <- grep("ACTIVATED", colnames(NK), value = TRUE)
  NK.activated <- NK[, NK.activated, drop = FALSE]
  NK.resting = grep("RESTING", colnames(NK))
  NK.resting <- NK[, NK.resting, drop = FALSE]
  NKT = grep("NKT", colnames(NK))
  NKT <- NK[, NKT, drop = FALSE]
  NK = NK[,-which(colnames(NK)%in%c(colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]  
  
  CD4 = data[,grepl(paste(c("CD4", "FOLLICULAR.HELPER"), collapse="|"), celltypes)]
  CD4.activated = grep("ACTIVATED", colnames(CD4))
  CD4.activated = CD4[, CD4.activated, drop = FALSE]
  CD4.resting = grep("RESTING", colnames(CD4))
  CD4.resting = CD4[, CD4.resting, drop = FALSE]
  CD4.naive = grep("NAIVE", colnames(CD4))
  CD4.naive = CD4[, CD4.naive, drop = FALSE]
  CD4.non.regulatory = grep("NON.REGULATORY", colnames(CD4))
  CD4.non.regulatory = CD4[, CD4.non.regulatory, drop = FALSE]
  CD4.follicular.helper = grep("FOLLICULAR.HELPER", colnames(CD4))
  CD4.follicular.helper = CD4[, CD4.follicular.helper, drop = FALSE]
  CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(CD4.activated), colnames(CD4.resting), colnames(CD4.naive), colnames(CD4.non.regulatory), colnames(CD4.follicular.helper)))]  
  
  CD8 = data[,grep("CD8", celltypes)]
  Tregs = data[,grep("TREGS", celltypes)]
  Dendritic = data[,grep("DENDRITIC", celltypes)]
  Dendritic.activated = grep("ACTIVATED", colnames(Dendritic))
  Dendritic.activated = Dendritic[, Dendritic.activated, drop = FALSE]
  Dendritic.resting = grep("RESTING", colnames(Dendritic))
  Dendritic.resting = Dendritic[, Dendritic.resting, drop = FALSE]
  Dendritic = Dendritic[,-which(colnames(Dendritic)%in%c(colnames(Dendritic.activated), colnames(Dendritic.resting)))]
  
  Cancer = data[,grep("CANCER", celltypes)]

  Mast = data[,grep("MAST", celltypes)]
  Mast.activated = grep("ACTIVATED", colnames(Mast))
  Mast.activated = Mast[, Mast.activated, drop = FALSE]
  Mast.resting = grep("RESTING", colnames(Mast))
  Mast.resting = Mast[, Mast.resting, drop = FALSE]
  Mast = Mast[,-which(colnames(Mast)%in%c(colnames(Mast.activated), colnames(Mast.resting)))]
  
  Eosinophils = grep("EOSINOPHILS", celltypes)
  Eosinophils = data[, Eosinophils, drop = FALSE]
  
  Endothelial = grep("ENDOTHELIAL", celltypes)
  Endothelial = data[, Endothelial, drop = FALSE]
  
  Myocytes = grep("MYOCYTES", celltypes)
  Myocytes = data[, Myocytes, drop = FALSE]
  
  CAF = data[,grep("CAF", celltypes)]
  
  T_gamma_delta = grep("GAMMA.DELTA", celltypes)
  T_gamma_delta = data[, T_gamma_delta, drop = FALSE]
  
  Uncharacterized = grep("UNCHARACTERIZED", celltypes)
  Uncharacterized = data[, Uncharacterized, drop = FALSE]
  
  cell_types = list(B, B.naive, B.memory, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.activated, CD4.resting, CD4.naive, CD4.non.regulatory, CD4.follicular.helper,
                    CD8, Tregs, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Mast, Mast.activated, Mast.resting, Eosinophils, Endothelial, Myocytes, CAF, T_gamma_delta, Uncharacterized)
  names(cell_types) = c("B", "B.NAIVE", "B.MEMORY", "MACROPHAGES", "M0", "M1", "M2", "MONOCYTES", "NEUTROPHILS", "NK", "NK.ACTIVATED", "NK.RESTING", "NKT", "CD4", "CD4.ACTIVATED", "CD4.RESTING", "CD4.NAIVE", "CD4.NON.REGULATORY", "CD4.FOLLICULAR.HELPER",
                        "CD8", "TREGS", "DENDRITIC", "DENDRITIC.ACTIVATED", "DENDRITIC.RESTING", "CANCER", "MAST", "MAST.ACTIVATED", "MAST.RESTING", "EOSINOPHILS", "ENDOTHELIAL", "MYOCYTES", "CAF", "GAMMA.DELTA", "UNCHARACTERIZED")
  
  return(cell_types)
}
cells = compute_cell.types(data.frame(deconv),celltypes, "LP")

#Grouping cell features
remove_subgroups = function(groups){
  lis = c()
  for (pos in 1:length(groups)){
    x = c()
    for (i in 1:length(groups[[pos]])) {
      x =  c(x,str_split(groups[[pos]][[i]], "_")[[1]][1])
    }
    if(length(unique(x)) == 1){
      lis = c(lis, names(groups)[[pos]])
    }
  }
  
  if(length(lis)>0){
    groups = groups[-which(names(groups)%in%lis)]
  }
  
  return(groups)
}
compute_subgroups = function(data, cell, thres_similarity = 0.05, thres_corr = 0.75, thres_change = 0.01){
  # data  = data.frame(cells[[1]]) 
  # cell = names(cells)[1]
  data = data.frame(data)
  if (ncol(data) < 2) {
    warning("Data must have at least two columns for subgrouping.")
    cell_subgroups = list()
    subgroup = list()
    lis = c()
    return(list(data, cell_subgroups, subgroup, lis))
  }else{
    cell_subgroups = list()
    lis = list()
    ###################################################################################################################Similarity part  
    is_similar <- function(value1, value2, threshold) {return(abs(value1 - value2) <= threshold)}
    similarity_matrix <- matrix(FALSE, nrow = ncol(data), ncol = ncol(data), dimnames = list(names(data), names(data)))
    for (col1 in names(data)) {
      for (col2 in names(data)) {
        similarity <- all(mapply(is_similar, data[[col1]], data[[col2]], MoreArgs = list(thres_similarity)))
        similarity_matrix[col1, col2] <- similarity
      }
    }
    # Get upper triangle of the correlation matrix
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat, diag = T)]<- NA
      return(cormat)
    }
    
    upper_tri <- get_upper_tri(similarity_matrix)
    x <- melt(upper_tri) %>%
      na.omit() %>%
      mutate_all(as.character)
    indice = 1
    subgroup = list()
    vec = unique(x$Var1)
    while(length(vec)>0){
      sub = x[which(x$Var1%in%vec[1]),] 
      sub = sub[which(sub$value==T),]
      if(nrow(sub)!=0){
        subgroup[[indice]] = c(vec[1], sub$Var2)
        x = x[-which(x$Var1%in%subgroup[[indice]]),]
        x = x[-which(x$Var2%in%subgroup[[indice]]),]
        vec = vec[-which(vec%in%subgroup[[indice]])]
        indice = indice + 1
      }else{
        indice = indice
        vec = vec[-1]
      }
    }
    
    if(length(subgroup)!=0){
      #Name subgroups    
      for (i in 1:length(subgroup)){ 
        names(subgroup)[i] = paste0(cell, "_Subgroup.Similarity", i)
      }
      
      lis = remove_subgroups(subgroup)
      data_sub = c()
      #Take average expression of subgroups
      for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
        sub = data.frame(data[,colnames(data)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
        sub$average = rowMeans(sub) #Compute average of subgroup across patients 
        data_sub = data.frame(cbind(data_sub, sub$average)) #Save avg in a new data frame
        colnames(data_sub)[i] = names(subgroup)[i]
        name = colnames(data)[which(!(colnames(data)%in%subgroup[[i]]))]
        data = data[,-which(colnames(data)%in%subgroup[[i]])] #Remove from deconvolution features that are subgrouped
        if(ncol(data.frame(data))==1){
          data = as.data.frame(data)
          colnames(data)[1] = name
        }
      }
      
      rownames(data_sub) = rownames(data) #List of patients
      data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(lis)])
      colnames(data_sub) = names(lis)
      
      data = cbind(data, data_sub) #Join subgroups in deconvolution file
      k = 2
    }else{
      k = 3
    }
    
    if(ncol(data) == 1){
      cell_subgroups = list()
      return(list(data, cell_subgroups, subgroup, lis))  
    }
    ###################################################################################################################Correlation part
    
    if(k==2 | k==3){
      terminate = FALSE
      iteration = 1
      while (terminate == FALSE) {
        corr_df <- correlation(data.matrix(data))
        vec = colnames(data)
        indice = 1
        subgroup = list()
        data_sub = c() 
        while(length(vec)>0){ #Keep running until no features are left
          if(vec[1] %in% corr_df$measure1){ #Check if feature still no-grouped
            tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature against the others
            tab = tab[tab$r>thres_corr,] #Select features corr above the threshold = 0.75
            if(nrow(tab)!=0){ #If algorithm found features above corr
              subgroup[[indice]] = c(vec[1], tab$measure2) #Save features as subgroup
              corr_df = corr_df[-which(corr_df$measure1 %in% subgroup[[indice]]),] #Remove features already subgroupped
              corr_df = corr_df[-which(corr_df$measure2 %in% subgroup[[indice]]),] #Remove features already subgroupped
              vec = vec[-which(vec%in%subgroup[[indice]])] #Remove feature already subgroupped from vector  
              indice = indice + 1
            }else{ #Condition when there is no correlation above the threshold (features no subgroupped)
              corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable from corr matrix to keep subgrouping the others
              vec = vec[-1] #Remove variable from vector to keep analyzing the others 
              indice = indice #Not increase index cause no subgroup appeared
            }
          }else{ #If feature is not in corr matrix it means that there is no any significant correlation against it and other features 
            vec = vec[-1] #Remove variable from vector to keep analyzing the others
            indice = indice  #Not increase index cause no subgroup appeared
          }
        }
        
        if(length(subgroup)!=0){
          #Name subgroups
          if(iteration > 1){
            idx = c()
            for (i in 1:length(subgroup)) {
              for(j in 1:length(subgroup[[i]])){
                idx = as.numeric(c(idx, str_split(subgroup[[i]], "_")[[j]][3]))
              }
            }
            for (i in 1:length(subgroup)){ 
              names(subgroup)[i] = paste0(cell, "_Subgroup_", paste0(idx, collapse = "."))
            }
          }else{
            for (i in 1:length(subgroup)){ 
              names(subgroup)[i] = paste0(cell, "_Subgroup_", i)
            }
            cell_subgroups = subgroup 
          }
          
          #Take average expression of subgroups
          for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
            sub = data.frame(data[,colnames(data)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
            sub$average = rowMeans(sub) #Compute average of subgroup across patients 
            data_sub = data.frame(cbind(data_sub, sub$average)) #Save avg in a new data frame
            colnames(data_sub)[i] = names(subgroup)[i]
            name = colnames(data)[which(!(colnames(data)%in%subgroup[[i]]))]
            data = data.frame(data[,-which(colnames(data)%in%subgroup[[i]])]) #Remove from deconvolution features that are subgrouped
            if(ncol(data.frame(data))==1){
              data = as.data.frame(data)
              colnames(data)[1] = name
            }
          }
          
          rownames(data_sub) = rownames(data) #List of patients
          
          new_average_values = colMeans(data.matrix(data_sub))
          
          if(iteration == 1){
            data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_subgroups)])
            colnames(data_sub) = names(cell_subgroups)
          }
          
          #Compare averages and test if they are above certain threshold 
          if(iteration == 1){
            df = data
          }else{
            for (i in 1:length(subgroup)) {
              for (j in 1:length(idx)) {
                change = max(abs(average_values[idx[j]] - new_average_values[i]))
                if (change > thres_change) {
                  terminate <- TRUE
                }
              }
            }
          }
          average_values = new_average_values
          
          if(ncol(data_sub)>1){
            data = data_sub
          }else{
            terminate = TRUE
          }
          iteration = iteration + 1
        }else{
          terminate = TRUE
          if(!is.null(tryCatch(eval(parse(text = df)), error = function(e) NULL))==F){
            df = data
          }
        }
        
      }
      
      if(is.null(data_sub)==TRUE){
        data = cbind(df, data)
      }else{
        data = cbind(df, data, data_sub)
      }
      
      idx = which(duplicated(t(data)))
      if(length(idx)>0){
        data = data[,-idx]
      }
    }
    
    
    return(list(data, cell_subgroups, subgroup, lis))
  }
  
}
res = list()
groups = list()
groups_similarity = list()
for (i in 1:length(cells)) {
  x = compute_subgroups(cells[[i]], names(cells)[i])
  res = c(res, x[1])
  groups = c(groups, x[2])
  groups_similarity = c(groups_similarity, x[4])
}

names_cells = c("B", "B.NAIVE", "B.MEMORY", "MACROPHAGES", "M0", "M1", "M2", "MONOCYTES", "NEUTROPHILS", "NK", "NK.ACTIVATED", "NK.RESTING", "NKT", "CD4", "CD4.ACTIVATED", "CD4.RESTING", "CD4.NAIVE", "CD4.NON.REGULATORY", 
                "CD4.FOLLICULAR.HELPER", "CD8", "TREGS", "DENDRITIC", "DENDRITIC.ACTIVATED", "DENDRITIC.RESTING", "CANCER", "MAST", "MAST.ACTIVATED", "MAST.RESTING", "EOSINOPHILS", "ENDOTHELIAL", "MYOCYTES", "CAF", "T.GAMMA.DELTA", "UNCHARACTERIZED")
names(res) = names_cells
names(groups) = names_cells
names(groups_similarity) = names_cells

#Joining subgroups in a new deconvolution output
dt = c()
for (i in 1:length(res)) {
  dt = c(dt, res[[i]])
}
dt = data.frame(dt)
rownames(dt) = rownames(deconv)

write.csv(dt, "Deconvolution.csv")

ha <- HeatmapAnnotation(Cluster = target, col = list(Cluster = c("1" = "green", "2" = "purple", "3" = "blue")))
png("heatmap.png", width = 7000, height = 4000, res=250)
Heatmap(t(scale(dt)), clustering_method_columns = "ward.D2")
dev.off()

#Subgrouping between cell types
corr_df <- correlation(data.matrix(dt))
vec = colnames(dt)
indice = 1
subgroup = list()
data_sub = c() 
thres_corr = 0.75
while(length(vec)>0){ #Keep running until no features are left
  if(vec[1] %in% corr_df$measure1){ #Check if feature still no-grouped
    tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature against the others
    tab = tab[tab$r>thres_corr,] #Select features corr above the threshold = 0.75
    if(nrow(tab)!=0){ #If algorithm found features above corr
      subgroup[[indice]] = c(vec[1], tab$measure2) #Save features as subgroup
      corr_df = corr_df[-which(corr_df$measure1 %in% subgroup[[indice]]),] #Remove features already subgroupped
      corr_df = corr_df[-which(corr_df$measure2 %in% subgroup[[indice]]),] #Remove features already subgroupped
      vec = vec[-which(vec%in%subgroup[[indice]])] #Remove feature already subgroupped from vector  
      indice = indice + 1
    }else{ #Condition when there is no correlation above the threshold (features no subgroupped)
      corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable from corr matrix to keep subgrouping the others
      vec = vec[-1] #Remove variable from vector to keep analyzing the others 
      indice = indice #Not increase index cause no subgroup appeared
    }
  }else{ #If feature is not in corr matrix it means that there is no any significant correlation against it and other features 
    vec = vec[-1] #Remove variable from vector to keep analyzing the others
    indice = indice  #Not increase index cause no subgroup appeared
  }
}

#Name subgroups
for (i in 1:length(subgroup)){ 
  names(subgroup)[i] = paste0("Cells_group_", i)
}

data = dt
#Take average expression of subgroups
for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
  sub = data.frame(dt[,colnames(dt)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
  sub$average = rowMeans(sub) #Compute average of subgroup across patients 
  data_sub = data.frame(cbind(data_sub, sub$average)) #Save avg in a new data frame
  colnames(data_sub)[i] = names(subgroup)[i]
  name = colnames(dt)[which(!(colnames(dt)%in%subgroup[[i]]))]
  dt = data.frame(dt[,-which(colnames(dt)%in%subgroup[[i]])]) #Remove from deconvolution features that are subgrouped
  if(ncol(data.frame(dt))==1){
    dt = as.data.frame(dt)
    colnames(dt)[1] = dt
  }
}

rownames(data_sub) = rownames(dt) #List of patients

dt = cbind(data_sub, dt)

png("heatmap2.png", width = 6000, height = 6000, res=250)
Heatmap(t(scale(data)), clustering_method_columns = "ward.D2", width = unit(35, "cm"), height = unit(45, "cm"), column_dend_height = unit(10, "cm"))
dev.off()

####################################################################End of deconvolution analysis







#Patients clustering 
clinical.data = read.csv("~/LungPredict/ClinicalLPVan.csv", row.names=1)
hc1 = as.dendrogram(hclust(dist(scale(Deconvolution)), method = "ward.D2"))
clusters <- cutree(hc1, k=3)
plot(hc1)
rect.dendrogram(hc1, k=3,horiz=F)
clinical.data$cluster = clusters

#Feature selection
anova_rfe_random_forest <- function(data, target) {
  # Step 1: ANOVA Feature Selection
  target = as.factor(target)
  anova_result <- apply(data, 2, function(x) {
    fit <- glm(target ~ x, data = data, family = binomial())
    summary(fit)$coef[2, 4]  # Extract the p-value (change the index if needed)
  })
  anova_selected_features <- data[, names(anova_result)[anova_result < 0.05]]  # Select features with p-value less than 0.05 (adjust as needed)
  
  # Step 2: RFE Feature Selection
  rfe_result <- rfe(data, target, sizes = c(50), rfeControl = rfeControl(functions = rfFuncs))
  rfe_selected_features <- data[, rfe_result$optVariables]
  
  # Step 3: Random Forest Classification
  rf_classifier <- randomForest(target ~ ., data = data)
  rf_importance <- importance(rf_classifier)
  rf_selected_features <- data[, names(rf_importance[order(rf_importance, decreasing = TRUE),][1:20])]
  
  # Find the common features from ANOVA, RFE, and Random Forest
  common_features <- Reduce(intersect, list(colnames(anova_selected_features), colnames(rfe_selected_features), colnames(rf_selected_features)))
  pval = anova_result[which(names(anova_result)%in%common_features)]
  
  return(list(common_features, pval))
}
df = Deconvolution[which(clinical.data$cluster%in%c(1,2)),]
target = clinical.data$cluster[which(clinical.data$cluster%in%c(1,2))]
res = anova_rfe_random_forest(df, target)

#Plot FS results
deconv_sub = df[,colnames(df)%in%res[[1]]]
matrix = rbind(deconv_sub, res[[2]])

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
           subtitle=paste0("ANOVA_RFE_RANDOM.FOREST test\npvalue: ", pval[i])) +
      theme(axis.text.x = element_text(angle = 0),
            axis.title.y = element_text(size = 8, angle = 90))
    
    
    ggsave(violin, file=paste0(getwd(),"/Figures/Violin_plots/", colnames(matrix[i]), "_", file_name, ".png"))
  }
}
compute_violin(matrix, target, "Clusters")

ha <- HeatmapAnnotation(Cluster = target, col = list(Cluster = c("1" = "green", "2" = "purple")))
Heatmap(scale(t(matrix[1:nrow(matrix)-1,])), clustering_method_columns = "ward.D2", top_annotation = ha)

#TFs analysis
TFs <- read.csv(paste0(getwd(),"/Output/TFs/VIPER/TFs_scores_all_LPVan.csv"), row.names = 1)
TFs_top = TFs[order(rowMeans(TFs), decreasing=TRUE)[1:50],]
Heatmap(t(scale(t(TFs_top))), clustering_method_columns = "ward.D2", top_annotation = ha)

