import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from paste2 import PASTE2, projection, model_selection


def plot_slices_overlap(slices, name):
    plt.figure(figsize=(7,7))

    all_values = np.concatenate([adata.obs['percent.butterfly'].values for adata in slices])
    vmin, vmax = all_values.min(), all_values.max()

    '''
    adata_img = slices[0]  
    image_name = str(adata_img.obs.name.iloc[0])
    img = adata_img.uns['spatial'][image_name]['images']['lowres']
    scale = adata_img.uns['spatial'][image_name]['scalefactors']['tissue_lowres_scalef']

    plt.imshow(img, extent=[
        0, img.shape[1] / scale,
        img.shape[0] / scale, 0
    ])
    '''

    for i in range(len(slices)):
        adata = slices[i]

        sc = plt.scatter(
            adata.obsm['spatial'][:, 0],
            adata.obsm['spatial'][:, 1],
            c=adata.obs['percent.butterfly'],
            cmap='viridis',  # Puedes cambiar el colormap si quieres (ej: 'plasma', 'magma', etc.)
            vmin=vmin,       # Mismo mínimo para todas las figuras
            vmax=vmax,       # Mismo máximo para todas las figuras
            s=60,
            marker=".",
            linewidth=0,
        )

    plt.colorbar(sc, label='Percent Butterfly')
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.savefig(f"{name}.png", dpi=300, bbox_inches='tight', transparent=True)
    plt.show()


carpetaData = "/home/vgaya/ST_imageAlignment/07.PASTE2/alignAnnData/"

patient = "19"
N = "2"

patient_1_name = (f"{carpetaData}Paciente{patient}_merge_1_align1{N}_h.h5ad")
patient_1_h = sc.read_h5ad(patient_1_name)

patient_prob_name = (f"{carpetaData}Paciente{patient}_merge_{N}_align1{N}_h.h5ad")
patient_prob_h = sc.read_h5ad(patient_prob_name)

plot_slices_overlap([patient_1_h], "grafico1")
#plot_slices_overlap([patient_prob_h], "graficoProb")