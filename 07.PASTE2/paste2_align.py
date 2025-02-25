# Packages

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squidpy as sq
import pandas as pd
import numpy as np
import anndata as ad 
from paste2 import PASTE2, projection

# Data
carpetaData = "/home/vgaya/ST_imageAlignment/02.Seurat/"

options = ["19", "MIX"]
print("Select which patient do you want to align:\n(patient 19 refers to three different slices from patient 19 and\npatient MIX refers to three different samples from different patients,\nin both choices the alignment would be performed to the same slice from patient 19\nwhich is considered as the reference)\n")
for i, option in enumerate(options, start=1):
    print(f"{i}. {option}") 
choice = 0
while choice not in (1:len(options)):
    print("\nInvalid choice. Please try again.")
    choice = int(input("\nEnter the number of your choice: "))

patient = options[choice - 1]
        
# Call the function
patient = select_patient()
