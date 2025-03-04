# Packages

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squidpy as sq
import pandas as pd
import numpy as np
from paste2 import PASTE2, projection

# Data
carpetaData = "/home/vgaya/ST_imageAlignment/02.Seurat/"

options = ["19", "MIX"]
print("Select which patient do you want to align:\n(patient 19 refers to three different slices from patient 19 and\npatient MIX refers to three different samples from different patients,\nin both choices the alignment would be performed to the same slice from patient 19\nwhich is considered as the reference)\n")
for i, option in enumerate(options, start=1):
    print(f"{i}. {option}") 
choice = int(input("\nEnter the number of your choice: "))
while choice < 1 or choice > len(options):
    print("\nInvalid choice. Please try again.")
    choice = int(input("\nEnter the number of your choice: "))

patient = options[choice - 1]
        
adatas = []
for i in range(1,4):
    adata_name = f"Paciente{patient}_merge_{i}.h5ad"
    adata = sc.read_h5ad(adata_name)
    adatas.append(adata)

adatas[0]



sq.pl.spatial_scatter(
    adatas[0],
    size = 1.5,
    color='percent.butterfly',
    figsize= (5, 5)
)
