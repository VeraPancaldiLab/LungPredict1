# LungPredict 

Profiling the Non-Small Cell Lung Cancer (NSCLC) microenvironment by integrating transcriptomics to uncover potential  phenotypic profiles associated to patterns in immune infiltration

## Summary
Here, we applied a computational immunology approach involving differential expression and pathway analyses, along with the quantification of immune cell proportions using computational deconvolution methods, inference of transcription factor activities, and specification of immune score types to better characterize bulk transcriptomics samples of lung adenocarcinoma. This comprehensive analysis enabled us to identify biomarkers of disease progression and immune infiltration patterns across disease stages and patient groups, thus allowing us to compute and evaluate specific hallmarks of immune response across patients. Through our methodology and feature selection pipeline, which we briefly describe, we identified niches of immune cells and report on a duality in the behavior of NK cells in the Tumor Microenvironment (TME), corresponding to either immune response or characteristic dysfunctional/exhausted states. We validated our results using an independent cohort, detailed in terms of size and nature, with both bulk and single-cell RNA sequencing (scRNAseq) samples. This validation enabled us to further investigate the NK cell profiles that we had identified.


## Data 
Bulk RNASeq data from NSCLC patients 
- LungPredict dataset: 62 patients with Lung Cancer (Adenocarcinoma) - Stage I, II, III, IV
- Vanderbilt dataset: 76 patients NSCLC (Adenocarcinoma) - Early Stage
- Vanderbilt Single-cell RNAseq data: 15 patients NSCL (Adenocarcinoma) - Early Stage

![image](https://github.com/VeraPancaldiLab/LungPredict1/assets/37853385/2641fa06-91e4-46f5-bc6f-4f83baacb035)

## Project organization
- Scripts: Codes used to develop the project's methods. 
- Figures: Images results from analysis in Early and Late Stage patients
- Output: Result files from analysis in Early and Late Stage patients
For producing the figures of the paper, use file `Scripts/Analysis_EarlyStage.rmd`

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
- Marcelo Hurtado (https://github.com/mhurtado13)
- Leila Khajavi (https://github.com/LeilaKhajavi)
- Mouneem Essabar (https://github.com/mouneem)
- Vera Pancaldi (https://github.com/VeraPancaldi)


