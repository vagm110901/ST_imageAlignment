# ST_imageAlignment
Image alignment methods for Spatial Transcriptomics data of Spinal Cord

* 02.Seurat: Spatial Transcriptomics data from Yadav et al., 2023. 
    4 slices from the patient 19 [19]
    4 samples from differnet patients (1 = 19, 2 = 43, 3 = 47, 4 = 45) [MIX]
    The slice/sample 1 from patient 19 is always used as the refernce.
    - saveSeurat_forAnnData.R & h5ad_create.py: save the csv files from a rds file (Seurat, R) to create an AnnData file (python) in AnnData19 & AnnDataMIX
    - R: R scripts for flip functions, and for merging the four slices/samples into a unique rds file

* 04.ImageAlignment: image alignment by the approaches: 
    - Geometric Transformation Estimation Model [GTEM] 
        (TR = Translation+Rotation, 
        TR_S = Translation+Rotation and Scaling after them, 
        STR = Scaling+Translation+Rotation simultaneously, 
        RT = Rotation+Translation)
    - Procrustes Transformation [Proc] 
        (with and without scaling)
    - ImageJ plugin Register Virtual Stack Slices [ImageJ] 
        results saved in * 05.ImageJ

* 06.RCTD: single cell deconvolution for ST with RCTD
    - RCTD_deconvolution.R : 
        1. Reference using single-nucleus RNA-seq data from Yadav et al., 2023
        2. Spatial using the three-terms objects (Reference + NotAlign + YesAlign)
        3. Cell type deconvolution
            Complete tissue
            Bottom left butterfly region
            Central canal region
    - RCTD_visualization.R : cell type composition visualization 

* 07.PASTE2: ST alignment using PASTE2 (histological + gene expression alignment)
    - paste2_align.py/ipynb: ST alignment 
    - RGB_values.py/ipynb: estimate the RGB mean value for all the pixels in each spot 

