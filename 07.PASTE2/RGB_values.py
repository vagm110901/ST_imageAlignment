import numpy as np
import pandas as pd
from skimage.draw import disk
import scanpy as sc
import cv2

carpetaData = "../02.Seurat/"
carpetaSave = "../07.PASTE2/"
patient = "19"
adatas = []
ims = []
for i in range(1,5):
    adata_name = f"{carpetaData}Paciente{patient}_merge_{i}.h5ad"
    adata = sc.read_h5ad(adata_name)
    adatas.append(adata)
    im_name = f"{carpetaData}AnnData{patient}/spatial_image_slice{i}.png"
    im = cv2.imread(im_name)
    ims.append(im)

for i in range(0,len(adatas)):
    adata = adatas[i]
    im = ims[i]
    # I am using the lowres image, because is the one that was available 
    image_name = str(adata.obs.name.iloc[0])
    spot_diam_fullres = float(adata.uns['spatial'][image_name]['scalefactors']['spot_diameter_fullres'])
    scale_lowres = float(adata.uns['spatial'][image_name]['scalefactors']['tissue_lowres_scalef'])
    spot_diam_lowres = round(spot_diam_fullres * scale_lowres)
    spot_rad_lowres = round(spot_diam_lowres / 2)

    rgb = []
    for j in range(adata.n_obs):
        x, y = round(adata.obsm['spatial'][j][0] * scale_lowres), round(adata.obsm['spatial'][j][1] * scale_lowres)
        spot_mask = np.zeros((im.shape[0], im.shape[1]), np.uint8)
        cv2.circle(spot_mask, (x, y), radius=int(spot_rad_lowres), color=(255, 255, 255), thickness=-1)
        rgb.append(tuple(int(round(c)) for c in cv2.mean(im, spot_mask)[::-1][1:]))

        #rgb.append(cv2.mean(im, spot_mask)[::-1][1:])
        
    adata.obsm['rgb'] = np.array(rgb)

    adata.write(f'{carpetaSave}Paciente{patient}_merge_{i+1}.h5ad')

