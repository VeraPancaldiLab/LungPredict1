
Counts1 = Counts1[rownames(Counts1)%in%rownames(Counts2),]
Counts2 = Counts2[rownames(Counts2)%in%rownames(Counts1),]

entrz = gconvert(row.names(Counts2), organism = "hsapiens", target = "ENTREZGENE_ACC", filter_na = F, mthreshold = 1)

entrz = entrz[,c(2,5)]
colnames(entrz)[c(1,2)] = c("ENSEMBL","SYMBOL") 

Counts2$ENSEMBL <- rownames(Counts2)
Counts2 <- Counts2 %>%
  inner_join(., entrz, by="ENSEMBL") %>%
  dplyr::filter(!is.na(SYMBOL))
gene_means <- rowMeans(Counts2 %>% dplyr::select(-ENSEMBL, -SYMBOL))
Counts2 <- Counts2 %>%
  dplyr::mutate(gene_means) %>% 
  dplyr::select(ENSEMBL, SYMBOL, gene_means, dplyr::everything()) 
Counts2 <- Counts2 %>%
  dplyr::arrange(dplyr::desc(gene_means)) %>% 
  dplyr::distinct(SYMBOL, .keep_all = TRUE)  
Counts2 <- Counts2 %>%
  dplyr::select(-ENSEMBL, -gene_means) %>%   
  tibble::column_to_rownames("SYMBOL") %>%   
  as.matrix() 

Counts1 = data.frame(Counts1)
Counts2 = data.frame(Counts2)

####TPM normalization
Counts1 = TPM_normalization(Counts1)
Counts2 = TPM_normalization(Counts2)

##Immunoscores
list_gold_standards=c("CYT","Roh_IS","Davoli_IS", "IFNy","Ayers_expIS","Tcell_inflamed")
gold.standards <- compute_gold_standards(RNA.tpm = data.frame(Counts2), list_gold_standards, "VANDERBILT")
gold.standards = minMax(gold.standards)
gold.standards = gold.standards[rownames(gold.standards)%in%rownames(clinical.data),]
immunoscores = rowMeans(gold.standards)
clinical.data$immuno = immunoscores
###Batch effect
clinical.data <- read.csv(paste0(getwd(),"/Input/ClinicalData/ClinicalData_Vanderbilt.csv"), row.names = 1)
Counts2 = Counts2[,colnames(Counts2)%in%rownames(clinical.data)]

Counts2 <- limma::removeBatchEffect(Counts2, batch=clinical.data$Batch) #Samples as columns

write.csv(Counts1, "TPM_LUNGPREDICT.csv", row.names = T)
write.csv(Counts2, "TPM_VANDERBILT.csv", row.names = T)
