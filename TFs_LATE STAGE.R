
source(paste0(getwd(),"/environment_set.R"))
libraries_set()

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
Counts_normalized = assay(vsd)

m2 <- do.call(rbind, strsplit(rownames(Counts_normalized), split="_", fixed = TRUE))
Counts_normalized = as.matrix(Counts_normalized)
rownames(Counts_normalized) = m2[,2]

TFs = compute_TFs_scores(Counts_normalized, "LateStage")
