# LungPredict 

Profiling the Non-Small Cell Lung Cancer (NSCLC) microenvironment by integrating transcriptomics to uncover potential  phenotypic profiles associated to patterns in immune infiltration

## Summary
Here, we applied a computational immunology approach involving differential expression and pathway analyses, immune cell proportion quantification via deconvolution, transcription factor inference and immune score estimation in order to better characterize bulk transcriptomics samples of lung adenocarcinoma. This analysis allowed us to identify biomarkers of disease progression and potential immune infiltration patterns across disease stages and patient groups enabling us to compute and evaluate hallmarks of immune response across patients. Through our methodology and feature selection pipeline, we identified niches of immune cells and report a duality in behavior of NK cells in the TME corresponding to either immune response or dysfunctional/exhausted states. We validate our results using an independent cohort with bulk and scRNAseq samples, allowing us to further investigate the duality of the NK cell profiles identified.

## Data 
Bulk RNASeq data from NSCLC patients 
- LungPredict dataset: 62 patients with Lung Cancer (Adenocarcinoma) - Stage I, II, III, IV
- Vanderbilt dataset: 76 patients NSCLC (Adenocarcinoma) - Early Stage
- ImmunoPredict dataset: 19 patients NSCLC (Adenocarcinoma) - Late Stage
- Vanderbilt Single-cell RNAseq data: 15 patients NSCL (Adenocarcinoma) - Early Stage

![image](https://github.com/VeraPancaldiLab/LungPredict1/assets/37853385/2641fa06-91e4-46f5-bc6f-4f83baacb035)

## Project organization
- Scripts: Codes used to develop the project's methods 
- Figures: Images results from analysis in Early and Late Stage patients
- Output: Result files from analysis in Early and Late Stage patients
- 
## Methodology
Analysis was divided in two parts for Early (Stage I, II) and Late Stage (Stage IV) patients. Some of the methods used:
- Unsupervised and supervised analysis
- Differential Gene expression (DESeq2)
- Immune cell type Deconvolution (GEMDECan)
- TFs inference for bulk and single cell data 
- Feature selection (Random forest, mutual information)
- Clustering algorithms (PCA, UMAP)
- Data visualization
- Univariate Linear model 
- Pseudobulk analysis for SingleCell
  
## Contributing
If you are interested in contributing to this project, please contact the project leader for more information (vera.pancaldi@inserm.fr)

## Authors
- Marcelo Hurtado (https://github.com/Marcelo1308)
- Leila Khajavi (https://github.com/LeilaKhajavi)
- Mouneem Essabar (https://github.com/mouneem)
- Vera Pancaldi (https://github.com/VeraPancaldi)


