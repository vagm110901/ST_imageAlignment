import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import anndata as ad
import squidpy as sq
from skimage import io
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

adatas = []

for i in range(1, 5):
    # Load expression matrixes and metadata for each slice
    expression = pd.read_csv(f'expression_matrix_slice{i}.csv', index_col=0)
    cell_metadata = pd.read_csv(f'cell_metadata_slice{i}.csv', index_col=0)
    gene_metadata = pd.read_csv('gene_metadata.csv', index_col=0)  # Los metadatos de genes son comunes
    
    # Cargar las coordenadas espaciales (row, col, imagerow, imagecol)
    image_coords = pd.read_csv(f"image_coordinates_slice{i}.csv", index_col=0)
    spatial_coords = image_coords[['row', 'col']]  # Ajustar según los nombres de las columnas
    scale_factors = pd.read_csv(f'scale_factors_slice{i}.csv', index_col=0)

    # Crear objeto AnnData para el slice
    adata = sc.AnnData(X=expression.values, obs=cell_metadata, var=gene_metadata)

    # Incluir las coordenadas espaciales
    adata.obsm['spatial'] = spatial_coords.values

    # Cargar la imagen espacial para el slice
    spatial_image = io.imread(f"spatial_image_slice{i}.png")

    # Asociar la imagen al objeto AnnData
    image_name = cell_metadata.name.iloc[0]
    adata.uns['spatial'] = {image_name: {}}
    adata.uns['spatial'][image_name]['images'] = {}
    adata.uns['spatial'][image_name]["images"] = {"hires": spatial_image}
    adata.uns['spatial'][image_name]["scalefactors"] = {
        "tissue_hires_scalef": scale_factors.tissue_hires_scalef.iloc[0],
        "spot_diameter_fullres": scale_factors.spot_diameter_fullres.iloc[0]
    }
    

    # Agregar el objeto AnnData a la lista
    adatas.append(adata)

# Combined the AnnData objects into a unique one
adata_combined = adatas[0].concatenate(adatas[1:], join='outer')

adata_combined.write(f"Paciente19_merge.h5ad")