import scanpy as sc

import pandas as pd
import numpy as np
import anndata as ad

#carpetaData = "/home/vgaya/ST_imageAlignment/07.PASTE2/alignAnnData/"
carpetaData = "/home/vgaya/ST_imageAlignment/07.PASTE2/alignAnnData/09/"

patient = "MIX"
N = "2"

patient_1_name = (f"{carpetaData}Paciente{patient}_merge_1_align1{N}_h.h5ad")
patient_1 = sc.read_h5ad(patient_1_name)

patient_prob_name = (f"{carpetaData}Paciente{patient}_merge_{N}_align1{N}_h.h5ad")
patient_prob = sc.read_h5ad(patient_prob_name)

# 2 columns: imagerow, imagecol
np.savetxt(f"{carpetaData}Paciente{patient}_1_align1{N}_h_coord.csv",
           patient_1.obsm['spatial'], delimiter=",", fmt="%.2f", header="imagerow,imagecol", comments="")
np.savetxt(f"{carpetaData}Paciente{patient}_{N}_align1{N}_h_coord.csv",
           patient_prob.obsm['spatial'], delimiter=",", fmt="%.2f", header="imagerow,imagecol", comments="")



