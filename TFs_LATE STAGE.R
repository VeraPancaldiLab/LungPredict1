Counts1 <- read.csv(paste0(getwd(),"/Input/Counts/Counts_LungPredict.csv"), row.names = 1)
Counts2 <- read.csv(paste0(getwd(),"/Input/Counts/Counts_ImmunoPredict.csv"), row.names = 1)
Counts = cbind(Counts1, Counts2)
Counts = Counts[,colnames(Counts)%in%rownames(clinical.data)]

clinical.data$Batch <- as.factor(clinical.data$Batch)
clinical.data$Smoking_Status <- as.factor(clinical.data$Smoking_Status)
clinical.data$Gender <- as.factor(clinical.data$Gender)
clinical.data$Age_Bracket <- as.factor(clinical.data$Age_Bracket)
clinical.data$Stages_simplified <- as.factor(clinical.data$Stages_simplified)

dds <- DESeqDataSetFromMatrix(
  countData = round(Counts),
  colData = clinical.data,
  design = ~Batch  + Gender + Smoking_Status + Age_Bracket)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds)
