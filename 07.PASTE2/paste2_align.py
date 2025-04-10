# Packages

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squidpy as sq
import pandas as pd
import numpy as np
from paste2 import PASTE2, projection, model_selection

###
from scipy.spatial import distance
from paste2.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix, generalized_kl_divergence, \
    high_umi_gene_distance, pca_distance, glmpca_distance

def partial_pairwise_align_histology(sliceA, sliceB, alpha=0.1, s=None, armijo=False, dissimilarity='glmpca', use_rep=None, G_init=None, a_distribution=None,
                   b_distribution=None, norm=True, return_obj=False, verbose=False, **kwargs):
    """
    Optimal partial alignment of two slices using both gene expression and histological image information.

    sliceA, sliceB must be AnnData objects that contain .obsm['rgb'], which stores the RGB value of each spot in the histology image.
    """
    m = s
    print("PASTE2 starts...")

    # subset for common genes
    common_genes = intersect(sliceA.var.index, sliceB.var.index)
    sliceA = sliceA[:, common_genes]
    sliceB = sliceB[:, common_genes]
    # print('Filtered all slices for common genes. There are ' + str(len(common_genes)) + ' common genes.')

    # Calculate spatial distances
    D_A = distance.cdist(sliceA.obsm['spatial'], sliceA.obsm['spatial'])
    D_B = distance.cdist(sliceB.obsm['spatial'], sliceB.obsm['spatial'])

    # Calculate expression dissimilarity
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, use_rep)), to_dense_array(extract_data_matrix(sliceB, use_rep))
    if dissimilarity.lower() == 'euclidean' or dissimilarity.lower() == 'euc':
        M_exp = distance.cdist(A_X, B_X)
    elif dissimilarity.lower() == 'kl':
        s_A = A_X + 0.01
        s_B = B_X + 0.01
        M_exp = kl_divergence(s_A, s_B)
    elif dissimilarity.lower() == 'glmpca':
        M_exp = glmpca_distance(A_X, B_X, latent_dim=50, filter=True, verbose=verbose)
    else:
        print("ERROR")
        exit(1)

    # Calculate RGB dissimilarity
    # sliceA_rgb = (sliceA.obsm['rgb'] - np.mean(sliceA.obsm['rgb'], axis=0)) / np.std(sliceA.obsm['rgb'], axis=0)
    # sliceB_rgb = (sliceB.obsm['rgb'] - np.mean(sliceB.obsm['rgb'], axis=0)) / np.std(sliceB.obsm['rgb'], axis=0)
    M_rgb = distance.cdist(sliceA.obsm['rgb'], sliceB.obsm['rgb'])
    # M_rgb = distance.cdist(sliceA_rgb, sliceB_rgb)

    # Scale M_exp and M_rgb, obtain M by taking half from each
    M_rgb /= M_rgb[M_rgb > 0].max()
    M_rgb *= M_exp.max()
    # M_exp /= M_exp[M_exp > 0].max()
    # M_rgb /= M_rgb[M_rgb > 0].max()
    
    M = 0 * M_exp + 1 * M_rgb               ### line 345
    # M = 0.5 * M_exp + 0.5 * M_rgb 

    # init distributions
    if a_distribution is None:
        a = np.ones((sliceA.shape[0],)) / sliceA.shape[0]
    else:
        a = a_distribution

    if b_distribution is None:
        b = np.ones((sliceB.shape[0],)) / sliceB.shape[0]
    else:
        b = b_distribution

    if norm:
        D_A /= D_A[D_A > 0].min().min()
        D_B /= D_B[D_B > 0].min().min()

        """
        Code for normalizing distance matrix
        """
        D_A /= D_A[D_A>0].max()
        D_A *= M.max()
        D_B /= D_B[D_B>0].max()
        D_B *= M.max()
        """
        Code for normalizing distance matrix ends
        """

    # Run OT
    pi, log = partial_fused_gromov_wasserstein(M, D_A, D_B, a, b, alpha=alpha, m=m, G0=G_init, loss_fun='square_loss', armijo=armijo, log=True, verbose=verbose)

    if return_obj:
        return pi, log['partial_fgw_cost']
    return pi

###

# Data
carpetaData = "/home/vgaya/ST_imageAlignment/07.PASTE2/"

patient = "19"

patient_1_name = f"{carpetaData}Paciente{patient}_merge_1.h5ad"
patient_1 = sc.read_h5ad(patient_1_name)

patient_2_name = f"{carpetaData}Paciente{patient}_merge_2.h5ad"
patient_2 = sc.read_h5ad(patient_2_name)

#patient_3_name = f"{carpetaData}Paciente{patient}_merge_3.h5ad"
#patient_3 = sc.read_h5ad(patient_3_name)

#patient_4_name = f"{carpetaData}Paciente{patient}_merge_4.h5ad"
#patient_4 = sc.read_h5ad(patient_4_name)

patient_1.X = patient_1.X.astype('int32')
patient_2.X = patient_2.X.astype('int32')

s_value = model_selection.select_overlap_fraction(patient_1, patient_2)

print("The estimated overlap percentage of the two slices is ", s_value)

pi = PASTE2.partial_pairwise_align(patient_1, patient_2, s = s_value)

np.savetxt("pi_12.csv", pi, delimiter=",")
print("Partial alignment of spots correctly saved")

"""
sq.pl.spatial_scatter(
    adatas[0],
    size = 1.5,
    color='percent.butterfly',
    figsize= (5, 5)
)
"""