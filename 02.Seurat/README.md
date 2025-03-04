# Spatial Transcriptomics data from Yadav et al., 2023. 

4 slices from the patient 19 [19]
4 samples from differnet patients (1 = 19, 2 = 43, 3 = 47, 4 = 45) [MIX]
The slice/sample 1 from patient 19 is always used as the refernce.

- saveSeurat_forAnnData.R & h5ad_create.py: save the csv files from a rds file (Seurat, R) to create an AnnData file (python) in AnnData19 & AnnDataMIX

- R: R scripts for flip functions, and for merging the four slices/samples into a unique rds file