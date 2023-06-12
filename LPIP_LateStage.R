#Deconvolution function analysis
source(paste0(getwd(),"/environment_set.R"))
libraries_set()

dataset = "LungPredict"
dataset = "ImmunoPredict"
dataset = "Response"
deconv <- as.matrix(read.csv(paste0(getwd(),"/Input/Deconvolution/Deconvolution_", dataset, ".csv"), row.names = 1))
clinical.data <- read.csv(paste0(getwd(),"/Input/ClinicalData/ClinicalData_", dataset, ".csv"), row.names = 1)
response <- read.csv(paste0(getwd(),"/Input/ClinicalData/ClinicalData_", dataset, ".csv"), row.names = 1)
response <- read.csv(paste0(getwd(),"/RawFiles/ColumnData_Response.csv"), row.names=1)
response = response[which(rownames(response)%in%rownames(clinical.data)),]
clinical.data = clinical.data[rownames(clinical.data)%in%rownames(response),]
clinical.data$response = response$Response
clinical.data = clinical.data[-which(clinical.data$Stages_simplified%in%c("III","IV")),]
clinical.data = clinical.data[which(clinical.data$Stages_simplified%in%"IV"),]
clinical.data = clinical.data[-which(clinical.data$response%in%"No data"),]
deconv = deconv[rownames(deconv)%in%rownames(clinical.data),]

#Filtering equal method-signatures (corr=1)
deconv = data.matrix(deconv[,-which(duplicated(t(deconv)))])

#Remove columns with low variance
deconv = data.frame(deconv[,-which(colVars(deconv)<=summary(colVars(deconv))[2])])
#Remove columns with more than 70% of zeros
i <- colSums(deconv == 0, na.rm=TRUE) < round(0.7*nrow(deconv))
deconv <- data.frame(deconv[, i, drop=FALSE])

#Split deconvolution matrix by cell type
compute_cell.types = function(data){
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
  B = Bcells
  
  Macrophages = data[,grep("acrophages", colnames(data))]
  M0 = data[,grep("M0", colnames(data))]
  M1 = data[,grep("M1", colnames(data))]
  M2 <- data[,grep("M2", colnames(data))]
  if(length(grep("LM22", colnames(M2)))>0){M2 <- M2[,-grep("LM22", colnames(M2))]}else{M2 <- M2} 
  test = data[,grep("LM22", colnames(data))]
  test = test[,grep("Macrophages.M2", colnames(test))]
  M2 = cbind(M2, test)
  Macrophages = Macrophages[,-which(colnames(Macrophages)%in%c(colnames(M0), colnames(M1), colnames(M2)))]
  
  Monocytes = data[,grep("Mono", colnames(data))]
  Neutrophils <- data[,grep("Neu", colnames(data))]
  
  NK = data[,grep("NK", colnames(data))]
  NK.activated <- grep("activated", colnames(NK), value = TRUE)
  NK.activated <- NK[, NK.activated, drop = FALSE]
  NK.resting = grep("resting", colnames(NK))
  NK.resting <- NK[, NK.resting, drop = FALSE]
  NKT = grep("NKT", colnames(NK))
  NKT <- NK[, NKT, drop = FALSE]
  NK = NK[,-which(colnames(NK)%in%c(colnames(NK.activated), colnames(NK.resting), colnames(NKT)))]  
  
  ################################# CD4 cells
  CD4 <- data[,grep("CD4", colnames(data))]
  CD4.activated = grep("activated", colnames(CD4))
  CD4.activated = CD4[, CD4.activated, drop = FALSE]
  CD4.resting = grep("resting", colnames(CD4))
  CD4.resting = CD4[, CD4.resting, drop = FALSE]
  CD4.naive = grep("naive", colnames(CD4))
  CD4.naive = CD4[, CD4.naive, drop = FALSE]
  CD4.non.regulatory = grep("regulatory", colnames(CD4))
  CD4.non.regulatory = CD4[, CD4.non.regulatory, drop = FALSE]
  #CD4.follicular.helper = grep("follicular", colnames(CD4))
  #CD4.follicular.helper = CD4[, CD4.follicular.helper, drop = FALSE]
  CD4 = CD4[,-which(colnames(CD4)%in%c(colnames(CD4.activated), colnames(CD4.resting), colnames(CD4.naive), colnames(CD4.non.regulatory)))]  
  
  ################################# CD8 cells
  CD8 <- data[,grep("CD8", colnames(data))]
  Tregs = data[,grep("regs", colnames(data))]
  Dendritic = data[,grep("endritic", colnames(data))]
  Dendritic.activated = grep("activated", colnames(Dendritic))
  Dendritic.activated = Dendritic[, Dendritic.activated, drop = FALSE]
  Dendritic.resting = grep("resting", colnames(Dendritic))
  Dendritic.resting = Dendritic[, Dendritic.resting, drop = FALSE]
  Dendritic = Dendritic[,-which(colnames(Dendritic)%in%c(colnames(Dendritic.activated), colnames(Dendritic.resting)))]
  
  Cancer = data[,grep("ancer", colnames(data))]
  Mast = data[,grep("ast", colnames(data))]
  Mast.activated = grep("activated", colnames(Mast))
  Mast.activated = Mast[, Mast.activated, drop = FALSE]
  Mast.resting = grep("resting", colnames(Mast))
  Mast.resting = Mast[, Mast.resting, drop = FALSE]
  Mast = Mast[,-which(colnames(Mast)%in%c(colnames(Mast.activated), colnames(Mast.resting)))]
  
  Eosinophils = grep("nophils", colnames(data))
  Eosinophils = data[, Eosinophils, drop = FALSE]
  
  Endothelial = grep("dothelial", colnames(data))
  Endothelial = data[, Endothelial, drop = FALSE]
  
  Myocytes = grep("Myocytes", colnames(data))
  Myocytes = data[, Myocytes, drop = FALSE]
  
  CAF = data[,grep("CAF", colnames(data))]
  
  T_gamma_delta = grep("delta", colnames(data))
  T_gamma_delta = data[, T_gamma_delta, drop = FALSE]
  
  Uncharacterized = grep("uncharacterized", colnames(data))
  Uncharacterized = data[, Uncharacterized, drop = FALSE]
  
  cell_types = list(B, Macrophages, M0, M1, M2, Monocytes, Neutrophils, NK, NK.activated, NK.resting, NKT, CD4, CD4.activated, CD4.resting, CD4.naive, CD4.non.regulatory,
                    CD8, Tregs, Dendritic, Dendritic.activated, Dendritic.resting, Cancer, Mast, Mast.activated, Mast.resting, Eosinophils, Endothelial, Myocytes, CAF, T_gamma_delta, Uncharacterized)
  names(cell_types) = c("B", "MACROPHAGES", "M0", "M1", "M2", "MONOCYTES", "NEUTROPHILS", "NK", "NK.ACTIVATED", "NK.RESTING", "NKT", "CD4", "CD4.ACTIVATED", "CD4.RESTING", "CD4.NAIVE", "CD4.NON.REGULATORY",
                        "CD8", "TREGS", "DENDRITIC", "DENDRITIC.ACTIVATED", "DENDRITIC.RESTING", "CANCER", "MAST", "MAST.ACTIVATED", "MAST.RESTING", "EOSINOPHILS", "ENDOTHELIAL", "MYOCYTES", "CAF", "GAMMA.DELTA", "UNCHARACTERIZED")
  
  return(cell_types)
  
}
cells = compute_cell.types(data.frame(deconv))
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
compute_subgroups = function(data, cell, thres_similarity = 0.05, thres_corr = 0.7, thres_change = 0.01){
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
names_cells = c("B", "MACROPHAGES", "M0", "M1", "M2", "MONOCYTES", "NEUTROPHILS", "NK", "NK.ACTIVATED", "NK.RESTING", "NKT", "CD4", "CD4.ACTIVATED", "CD4.RESTING", "CD4.NAIVE", "CD4.NON.REGULATORY", 
                "CD8", "TREGS", "DENDRITIC", "DENDRITIC.ACTIVATED", "DENDRITIC.RESTING", "CANCER", "MAST", "MAST.ACTIVATED", "MAST.RESTING", "EOSINOPHILS", "ENDOTHELIAL", "MYOCYTES", "CAF", "T.GAMMA.DELTA", "UNCHARACTERIZED")
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

png(paste0(getwd(),"/Figures/Heatmaps/New/DeconvolutionLPIP.png"), width = 6000, height = 4500, res=250)
#ha <- HeatmapAnnotation(Response = clinical.data$response, col = list(Response = c("Responders"="green", "Non responders" = "black")))
Heatmap(t(scale(dt)), clustering_method_columns = "ward.D2",  width = unit(30, "cm"), height = unit(32, "cm"), column_dend_height = unit(5, "cm"))
dev.off()

####################################################################End of deconvolution analysis

#Patients clustering 
hc1 = as.dendrogram(hclust(dist(scale(dt)), method = "ward.D2"))
clusters <- cutree(hc1, k=2)
plot(hc1)
rect.dendrogram(hc1, k=2,horiz=F)
clinical.data$cluster = NULL
clinical.data$cluster = clusters

#Feature selection
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
res = random_forest(dt, response$cluster, 20)

ggplot(res[[2]], aes(x = reorder(rownames(res[[2]]), -Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Features", y = "Importance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

response = clinical.data[-which(clinical.data$response=="No data"),]
dt = dt[rownames(dt)%in%rownames(response),]
png(paste0(getwd(),"/Figures/Heatmaps/New/DeconvolutionLPIP_FS.png"), width = 6000, height = 4500, res=250)
ha <- HeatmapAnnotation(Response = response$response, Cluster = response$cluster,col = list(Cluster = c("1" = "yellow", "2"= "purple"), Response = c("Responders"="green", "Non responders"="red")))
Heatmap(t(scale(res[[1]])), clustering_method_columns = "ward.D2",  top_annotation = ha, width = unit(30, "cm"), height = unit(20, "cm"), column_dend_height = unit(10, "cm"))
dev.off()

response <- read.csv(paste0(getwd(),"/Input/ClinicalData/ClinicalData_Response.csv"), row.names = 1)
response = response[rownames(response)%in%rownames(clinical.data),]
response_dt = dt[rownames(dt)%in%rownames(response),]
response_dt = response_dt[,colnames(response_dt)%in%colnames(res[[1]])]
response = response[-which(response$response=="No data"),]
response_dt = response_dt[rownames(response_dt)%in%rownames(response),]
png(paste0(getwd(),"/Figures/Heatmaps/New/DeconvolutionLP_LATE_STAGE_FS_RESPONSE_ONLY.png"), width = 6000, height = 4500, res=250)
ha <- HeatmapAnnotation(Response = response$response,col = list(Response = c("Responders" = "green", "Non responders"= "red")))
Heatmap(t(scale(response_dt)), clustering_method_columns = "ward.D2",  top_annotation = ha, width = unit(30, "cm"), height = unit(20, "cm"), column_dend_height = unit(10, "cm"))
dev.off()

compute_violin(res[[1]], clinical.data$response, "_RespondeIP")

compute_violin = function(data, feature, file_name){
  matrix = data
  feature = as.factor(feature)
  #matrix = data.frame(data[1:nrow(data)-1,])
  matrix$measure = feature
  #pval = data[nrow(data),]
  
  for (i in (1:(ncol(matrix)-1))) {
    violin = ggplot(matrix, aes(x=factor(measure), y=matrix[,i], fill=measure)) +
      geom_violin(width=0.6) +
      geom_boxplot(width=0.07, color="black", alpha=0.2) +
      scale_fill_brewer() +
      geom_smooth(aes(x=as.factor(measure), y=matrix[,i]), method = "loess") +
      ylab(colnames(matrix[i])) +
      xlab(file_name) +
      #labs(title=colnames(matrix[i]),
           #subtitle=paste0("ANOVA_RFE_RANDOM.FOREST test\npvalue: ", pval[i])) +
      theme(axis.text.x = element_text(angle = 0),
            axis.title.y = element_text(size = 8, angle = 90))
    
    
    ggsave(violin, file=paste0(getwd(),"/Figures/Violin_plots/", colnames(matrix[i]), "_", file_name, ".png"))
  }
}
compute_violin(matrix, clinical.data$cluster, "_Van")

##################################DEG between clusters
Counts <- read.csv(paste0(getwd(),"/Input/Counts/Counts_", dataset, ".csv"), row.names = 1)
Counts = Counts[,colnames(Counts)%in%rownames(clinical.data)]
m2 <- do.call(rbind, strsplit(rownames(Counts), split="_", fixed = TRUE))
Counts = as.matrix(Counts)
rownames(Counts) = m2[,2]

clinical.data$Batch <- as.factor(clinical.data$Batch)
clinical.data$Smoking_Status <- as.factor(clinical.data$Smoking_Status)
clinical.data$Gender <- as.factor(clinical.data$Gender)
clinical.data$Age_Bracket <- as.factor(clinical.data$Age_Bracket)
clinical.data$Stages_simplified <- as.factor(clinical.data$Stages_simplified)
clinical.data$cluster = as.factor(clinical.data$cluster)

dds <- DESeqDataSetFromMatrix(
  countData = round(Counts),
  colData = clinical.data,
  design = ~Batch + Gender + Smoking_Status + Age_Bracket + cluster)

dds <- estimateSizeFactors(dds)
ddsObj = DESeq(dds)
vsd <- vst(dds)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$Batch)

resultsALL = results(ddsObj, contrast = c("cluster", "1","2"), pAdjustMethod = "BH")
res = data.frame(resultsALL)
res = na.omit(res)
res = res[res$padj < 0.05,]

png(paste0(getwd(),"/Figures/Heatmaps/Volcano_DEG_Cluster_LP_LATE_STAGE.png"), width = 3000, height = 2000, res=200)
print(EnhancedVolcano(res,
                      title = "Cluster 1 vs Cluster 2",
                      subtitle = NULL,
                      caption = NULL,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlim = c(-15, 15),
                      ylim = c(0, 20), 
                      pCutoff = 10e-2,
                      FCcutoff = 1))
dev.off()

res = res[which((res$log2FoldChange>1.5)|(res$log2FoldChange<(-1.5))),]
res$Symbol = rownames(res)
res <- data.frame("entrez_id" = mapIds(org.Hs.eg.db,
                                       keys = rownames(res),
                                       keytype = "SYMBOL",
                                       column = "ENTREZID",
                                       multiVals = "first")
) %>% 
  dplyr::filter(!is.na(entrez_id)) %>%  
  tibble::rownames_to_column("Symbol") %>% 
  dplyr::inner_join(res, by = "Symbol") 

counts_norm = assay(vsd)
counts_DEG = counts_norm[rownames(counts_norm)%in%res$Symbol,]
counts_DEG = counts_DEG[order(rowVars(counts_DEG), decreasing=TRUE)[1:100],]
png(paste0(getwd(),"/Figures/Heatmaps/Heatmap_DEG_Van.png"), width = 6000, height = 5000, res=250)
ha <- HeatmapAnnotation(Cluster = clinical.data$cluster, col = list(Cluster = c("1" = "green", "2" = "purple")))
Heatmap(t(scale(t(counts_DEG))), clustering_method_columns = "ward.D2", top_annotation = ha,  width = unit(30, "cm"), height = unit(40, "cm"), column_dend_height = unit(5, "cm"))
dev.off()

#Enrichment of DEGs
# up = res[res$log2FoldChange>1.5,]
# down = res[res$log2FoldChange<1.5,]
# kegg_results = enrichKEGG(gene = down$entrez_id, organism = "hsa",  pvalueCutoff = 0.05, pAdjustMethod = "BH")
# print(dotplot(kegg_results))
# 
# x <- enrichPathway(gene=down$entrez_id, pvalueCutoff = 0.05, readable=TRUE)
# print(dotplot(x))

##################################TFs between clusters
compute_msVIPER_scores = function(RNA.counts.normalized, test_name, ref_name, measure){
  
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
  mra <- ledge(mra)
  # Extract Differential Active TFs
  TFs = data.frame(mra[["es"]][["p.value"]])
  colnames(TFs) = "pval"
  vec = TFs$pval < 0.05
  TFs_msViper_adj = TFs[vec,]
  TFs_msViper_adj = as.data.frame(TFs_msViper_adj)
  rownames(TFs_msViper_adj) = rownames(TFs)[vec]
  TFs_msviper = vpres[rownames(vpres)%in%rownames(TFs_msViper_adj),]
  
  # write.csv(TFs_msviper,paste0(getwd(),"/Output/TFs/msVIPER/TFs_scores_DifferentialActive_", file_name, ".csv"))
  # 
  # # Plot DiffActive TFs
  # png(paste0(getwd(),"/Figures/msVIPER/DifferentialActive_TFs_", file_name, ".png"), width = 1000, height = 800, res=100)
  # print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  # dev.off()
  
  return(list(TFs_msviper, mra))
  
}
res_TFs = compute_msVIPER_scores(assay(vsd), "1", "2", clinical.data$cluster)

png(paste0(getwd(),"/Figures/Heatmaps/Heatmap_TFs_LP_LATE_STAGE.png"), width = 6000, height = 4500, res=250)
ha <- HeatmapAnnotation(Cluster = clinical.data$cluster, col = list(Cluster = c("1" = "green", "2" = "purple")))
Heatmap(t(scale(t(res_TFs[[1]]))), clustering_method_columns = "ward.D2", top_annotation = ha,  width = unit(30, "cm"), height = unit(20, "cm"), column_dend_height = unit(10, "cm"))
dev.off()

TFs = data.frame(res_TFs[[1]])
TFs$Symbol = rownames(TFs) 
TFs <- data.frame("entrez_id" = mapIds(org.Hs.eg.db, 
                                       keys = rownames(TFs),
                                       keytype = "SYMBOL",
                                       column = "ENTREZID",
                                       multiVals = "first")
) %>% 
  dplyr::filter(!is.na(entrez_id)) %>%  
  tibble::rownames_to_column("Symbol") %>% 
  dplyr::inner_join(TFs, by = "Symbol") 

###TFs enrichment

#Split by UP-DOWN regulation

#Vanderbilt
upTFs = TFs[which(TFs$Symbol%in%c("ZNF207", "TBX21", "IRF1", "USF1", "FOXP2", "KMT2B", "NFYB")),]
downTFs = TFs[-which(TFs$Symbol%in%c("ZNF207", "TBX21", "IRF1", "USF1", "FOXP2", "KMT2B", "NFYB")),]

#LungPredict
upTFs = TFs[which(TFs$Symbol%in%c("NR2C2", "E2F4", "BCL6", "ETV6", "DEAF1", "ZEB2")),]
downTFs = TFs[-which(TFs$Symbol%in%c("NR2C2", "E2F4", "BCL6", "ETV6", "DEAF1", "ZEB2")),]

#LPVan
upTFs = TFs[which(TFs$Symbol%in%c("MYC", "FOXM1", "E2F4", "E2F1")),]
downTFs = TFs[-which(TFs$Symbol%in%c("MYC", "FOXM1", "E2F4", "E2F1")),]

#Add targets 
mra = res_TFs[[2]]
x = mra[["ledge"]][which(names(mra[["ledge"]])%in%upTFs$Symbol)]
x = mra[["ledge"]][which(names(mra[["ledge"]])%in%downTFs$Symbol)]
genes = c()
for (i in 1:length(x)) {
  genes = c(genes,x[[i]])
}
genesID <- data.frame("entrez_id" = mapIds(org.Hs.eg.db,
                                           keys = genes,
                                           keytype = "SYMBOL",
                                           column = "ENTREZID",
                                           multiVals = "first")
) %>% 
  dplyr::filter(!is.na(entrez_id))

#Filter only DEG targets
genesID = genesID[which(genesID$entrez_id%in%res$entrez_id),]

kegg_results = enrichKEGG(gene =c(downTFs$entrez_id, genesID), organism = "hsa",  pvalueCutoff = 0.05, pAdjustMethod = "BH")
print(dotplot(kegg_results))

x <- enrichPathway(gene=c(upTFs$entrez_id, genesID), pvalueCutoff = 0.05, readable=TRUE)
print(dotplot(x))

##########################################

p <- pca(t(scale(matrix[1:nrow(matrix)-1,])), metadata = clinical.data, removeVar = 0.1)

pbiplot_loadings <- biplot(p, x = "PC1", y = "PC2",
                           showLoadings = TRUE,
                           lengthLoadingsArrowsFactor = 1.5,
                           sizeLoadingsNames = 3,
                           colLoadingsNames = 'red4',
                           lab = NULL,
                           colby = "cluster",
                           hline = 0, vline = c(-25, 0, 25),
                           vlineType = c('dotdash', 'solid', 'dashed'),
                           gridlines.major = FALSE, gridlines.minor = FALSE,
                           pointSize = 3,
                           legendPosition = 'left', legendLabSize = 14, legendIconSize = 8.0,
                           drawConnectors = FALSE,
                           title = 'PCA loadings',
                           subtitle = "PC1 versus PC2")

png(paste0(getwd(),"/Figures/Heatmaps/Biplot_deconvolutionIP_LATE_STAGE.png"), width = 2000, height = 1500, res=200)
print(pbiplot_loadings)
dev.off()

x = dt[rownames(dt)%in%rownames(bootstrap_data),]
# Set the number of bootstrapping iterations
num_bootstraps <- 10

# Create empty lists to store the bootstrap results
bootstrap_results <- list()
df_responders = clinical.data[which(clinical.data$response=="Responders"),]
df_Non_responders = clinical.data[which(clinical.data$response=="Non responders"),]
# Perform random bootstrapping
for (i in 1:num_bootstraps) {
  # Create a bootstrap sample
  bootstrap_data <- rbind(df_responders, df_Non_responders[sample(nrow(df_Non_responders), size = 7, replace = FALSE), ])
  target = as.factor(bootstrap_data[,'response'])
  # Perform random forest classification
  rf_classifier <- randomForest(target ~ ., data = dt[rownames(dt)%in%rownames(bootstrap_data),])
  rf_importance <- importance(rf_classifier)
  rf_importance <- rf_importance[order(rf_importance, decreasing = TRUE), ]
  rf_selected_features <- dt[, names(rf_importance[order(rf_importance, decreasing = TRUE)][1:5])]
  # Store the random forest model in the bootstrap_results list
  bootstrap_results[[i]] <- rf_selected_features
}

x = c()
for (i in 1:length(bootstrap_results)) {
  x = c(x, names(bootstrap_results[[i]]))
}

intersection = unique(x)
dt_response = dt[,colnames(dt)%in%intersection]
compute_violin = function(data, feature, file_name){
  matrix = data 
  feature = as.factor(feature)
  #matrix = data.frame(data[1:nrow(data)-1,])
  matrix$measure = feature
  #pval = data[nrow(data),]
  
  for (i in (1:(ncol(matrix)-1))) {
    violin = ggplot(matrix, aes(x=factor(measure), y=matrix[,i], fill=measure)) +
      geom_violin(width=0.6) +
      geom_boxplot(width=0.07, color="black", alpha=0.2) +
      scale_fill_brewer() +
      geom_smooth(aes(x=as.factor(measure), y=matrix[,i]), method = "loess") +
      ylab(colnames(matrix[i])) +
      xlab(file_name) +
      #labs(title=colnames(matrix[i]),
      #subtitle=paste0("ANOVA_RFE_RANDOM.FOREST test\npvalue: ", pval[i])) +
      theme(axis.text.x = element_text(angle = 0),
            axis.title.y = element_text(size = 8, angle = 90))
    
    
    ggsave(violin, file=paste0(getwd(),"/Figures/Violin_plots/", colnames(matrix[i]), "_", file_name, ".png"))
  }
}
dt_response$measure = NULL
compute_violin(dt_response, clinical.data$response, "response")

png(paste0(getwd(),"/Figures/Heatmaps/DeconvolutionResponse_FS_TEST2.png"), width = 6000, height = 4500, res=250)
ha <- HeatmapAnnotation(Response = clinical.data$response, col = list(Response = c("Responders"="green", "Non responders" = "black")))
Heatmap(t(scale(dt_response)), clustering_method_columns = "ward.D2",  top_annotation = ha, width = unit(30, "cm"), height = unit(20, "cm"), column_dend_height = unit(10, "cm"))
dev.off()

correlation <- function(df) {
  df = cbind(dt, clinical.data$response)
  M <- Hmisc::rcorr(as.matrix(df), type = "spearman")
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

clinical.data$response[which(clinical.data$response=="Responders")] = 1
clinical.data$response[which(clinical.data$response=="Non responders")] = 2

moduleTraitCor = cor(dt, clinical.data, use = "p")
z = correlation(cbind(dt, clinical.data$response))
z[z$measure1=="clinical.data$response",]
