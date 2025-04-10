

###############################################################################
#                     packages and reticulate.                                #        
###############################################################################

# Load a custom function for alignment evaluation
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_selectCoord.R")

library(ggplot2)
library(reticulate)
use_python("/opt/anaconda3/bin/python3")
py_config()
#knitr::knit_engines$set(python = reticulate::eng_python)

py_run_string("
import sys
sys.path.append('/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/08.STalign/STalign/STalign/')
import STalign
import matplotlib.pyplot as plt 
import numpy as np
import torch
")

###############################################################################
#           function transform_points_source_with_A(A, pointsI)                #        
#                                                                             #
# Transform points with an affine matrix.                                     #
# Note points are in row column order, not xy.                                # 
#                                                                             #
# Parameters                                                                  #
# ----------                                                                  #
#                                                                             #
# A  : torch tensor                                                           #
#      Affine transform matrix                                                #
# pointsI : torch tensor                                                      #
#      N x 2 set of corresponding points for matching in source image.        #
#      Default None (no points).                                              #
#                                                                             #
# Returns                                                                     #
# -------                                                                     #
# pointsI : torch tensor                                                      #
#      N x 2 set of points I after affine transformation A                    #
###############################################################################
py_run_string("
def transform_points_source_with_A(A, pointsI):
    if isinstance(pointsI,torch.Tensor):
        pointsIt = torch.clone(pointsI)
    else:
        pointsIt = torch.tensor(pointsI)
        
    pointsIt = (A[:2,:2]@pointsIt.T + A[:2,-1][...,None]).T
    return pointsIt
")

###############################################################################
#     Code for reading the datasets from Seurat objects:                      #        
#                                                                             #
# * tissue position list                                                      #
# * barcodes list                                                             #
# * scalefactors - tissue hires/lowres scalef                                 #
# * images                                                                    #
###############################################################################

options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/08.STalign/"
  
paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))

listSTalignResults <- list()
listSTalignResults[[1]] <- list()

ref_image_coords <- paciente.merge@images[["slice1.1"]]@coordinates
ref_pos <- as.data.frame(ref_image_coords)  # Convertir en un data.frame
ref_barcodes <- rownames(ref_pos)
ref_scalefactor <- paciente.merge@images[["slice1.1"]]@scale.factors$lowres
ref_image <- paciente.merge@images[["slice1.1"]]@image

posCalc_ref <- ref_pos[ref_barcodes,5:4] * ref_scalefactor
colnames(posCalc_ref) <- c("x", "y")
max_yref <- max(posCalc_ref$y, na.rm = TRUE)
min_yref <- min(posCalc_ref$y, na.rm = TRUE)

listSTalignResults[[1]][["image"]] <- ref_image
listSTalignResults[[1]][["coords"]] <- ref_pos

py_run_string("
## source (reference)
im_ref_norm = STalign.normalize(r.ref_image)
")

###############################################################################
#                 Reference landmarks selection                               #
#                                                                             #
#   To assist with the spatial alignment, I will place a few landmark points. #
#   This will help mitigate the likelihood of us falling into a local minimum #
#   in the gradient descent and arrive at a suboptimal solution. I will just  #
#   manually create them by eyeballing the image. They can be very approximate#
#   as STalign will integrate these landmarks with other imaging features in  #
#   its optimization.                                                         #
###############################################################################

coordenadas_im_ref <- selectoCoordSTalign(py$im_ref_norm) # Execute in command line
for (j in seq_along(coordenadas_im_ref)) {
  coordenadas_im_ref[[j]] <- round(coordenadas_im_ref[[j]])
}

###############################################################################
###############################################################################
####                     IMPORTANT INFORMATION                             ####
#### It must be applied each image manually, because the combination of    ####
#### python code using reticulate, R code and the landmark selection,      ####
#### does NOT work correctly                                               ####
####                                                                       ####
###############################################################################
###############################################################################

im <- 2

listSTalignResults[[im]] <- list()
  
prob_image_coords <- paciente.merge@images[[paste0("slice1.",im)]]@coordinates
prob_pos <- as.data.frame(prob_image_coords)  # Convertir en un data.frame
prob_barcodes <- rownames(prob_pos)
prob_scalefactor <- paciente.merge@images[[paste0("slice1.",im)]]@scale.factors$lowres
prob_image <- paciente.merge@images[[paste0("slice1.",im)]]@image
  
posCalc_prob <- prob_pos[prob_barcodes,5:4] * prob_scalefactor
colnames(posCalc_prob) <- c("x", "y")
max_yprob <- max(posCalc_prob$y, na.rm = TRUE)
min_yprob <- min(posCalc_prob$y, na.rm = TRUE)

# df <- data.frame(rbind(cbind(posCalc_ref, time='ref'), cbind(posCalc_prob, time='prob')))
# ggplot(df, aes(x=x, y=y, col=time)) + geom_point(alpha=0.5) + theme_minimal() + scale_y_reverse()
  
py_run_string("
## target (problem)
im_prob_norm = STalign.normalize(r.prob_image)
")
  
###############################################################################
#                 Visualization before alignment                              #
###############################################################################
  
par(mfrow=c(1,2))
  
## source (reference)
plot(c(0,0), xlim=c(0,ncol(py$im_ref_norm)), ylim=c(0,nrow(py$im_ref_norm)),
     xlab=NA, ylab=NA, xaxt="n", yaxt="n", bty="n", type="n")
rasterImage(py$im_ref_norm, xleft = 0, xright = ncol(py$im_ref_norm),
            ytop = 0, ybottom = nrow(py$im_ref_norm), interpolate = FALSE)
points(posCalc_ref$x, max_yref - posCalc_ref$y + min_yref, pch=16, cex=0.2, col='blue')
  
## target (problem)
plot(c(0,0), xlim=c(0,ncol(py$im_prob_norm)), ylim=c(0,nrow(py$im_prob_norm)),
     xlab=NA, ylab=NA, xaxt="n", yaxt="n", bty="n", type="n")
rasterImage(py$im_prob_norm, xleft = 0, xright = ncol(py$im_prob_norm),
            ytop = 0, ybottom = nrow(py$im_prob_norm), interpolate = FALSE)
points(posCalc_prob$x, max_yprob - posCalc_prob$y + min_yprob, pch=16, cex=0.2, col='red')
  
###############################################################################
#                 Problem landmarks selection                               #
#                                                                             #
#   To assist with the spatial alignment, I will place a few landmark points. #
#   This will help mitigate the likelihood of us falling into a local minimum #
#   in the gradient descent and arrive at a suboptimal solution. I will just  #
#   manually create them by eyeballing the image. They can be very approximate#
#   as STalign will integrate these landmarks with other imaging features in  #
#   its optimization.                                                         #
###############################################################################

coordenadas_im_prob <- selectoCoordSTalign(py$im_prob_norm) # Execute in command line
for (j in seq_along(coordenadas_im_prob)) {
  coordenadas_im_prob[[j]] <- round(coordenadas_im_prob[[j]])
}
  
print(coordenadas_im_ref); print(coordenadas_im_prob)
  
# Visualization
  
par(mfrow=c(1, 2))
## source (reference)
points_im_ref <- t(data.frame(a = c(coordenadas_im_ref$x[[1]], coordenadas_im_ref$y[[1]]), 
                              b = c(coordenadas_im_ref$x[[2]], coordenadas_im_ref$y[[2]]), 
                              c = c(coordenadas_im_ref$x[[3]], coordenadas_im_ref$y[[3]]),
                              d = c(coordenadas_im_ref$x[[4]], coordenadas_im_ref$y[[4]]),
                              e = c(coordenadas_im_ref$x[[5]], coordenadas_im_ref$y[[5]])))
colnames(points_im_ref) <- c('x', 'y')
  
plot(c(0,0), xlim=c(0,ncol(py$im_ref_norm)), ylim=c(0,nrow(py$im_ref_norm)),  xlab=NA, ylab=NA,
      xaxt="n", yaxt="n", bty="n", type="n")
rasterImage(py$im_ref_norm, xleft = 0, xright = ncol(py$im_ref_norm),
            ytop = 0, ybottom = nrow(py$im_ref_norm), interpolate = FALSE)
points(points_im_ref, col='red', pch=16)
  
## target (problem)
points_im_prob <- t(data.frame(a = c(coordenadas_im_prob$x[[1]], coordenadas_im_prob$y[[1]]), 
                               b = c(coordenadas_im_prob$x[[2]], coordenadas_im_prob$y[[2]]), 
                               c = c(coordenadas_im_prob$x[[3]], coordenadas_im_prob$y[[3]]),
                               d = c(coordenadas_im_prob$x[[4]], coordenadas_im_prob$y[[4]]),
                               e = c(coordenadas_im_prob$x[[5]], coordenadas_im_prob$y[[5]])))
colnames(points_im_prob) <- c('x', 'y')
  
plot(c(0,0), xlim=c(0,ncol(py$im_prob_norm)), ylim=c(0,nrow(py$im_prob_norm)),  xlab=NA, ylab=NA,
     xaxt="n", yaxt="n", bty="n", type="n")
rasterImage(py$im_prob_norm, xleft = 0, xright = ncol(py$im_prob_norm),
            ytop = 0, ybottom = nrow(py$im_prob_norm), interpolate = FALSE)
points(points_im_prob, col='red', pch=16)
  
###############################################################################
#                 R - Python correspondence                                   #
#                                                                             #
#   STalign with tensors in Python use row-column designation while position  #
#   coordinates in R use x-y.                                                 #
###############################################################################
  
# switch to row col order
points_im_prob_rc = points_im_prob[,2:1]
points_im_ref_rc = points_im_ref[,2:1]
  
# split so we can access them more easily in Python later
x_im_prob <- posCalc_prob[,1]
y_im_prob <- posCalc_prob[,2]
x_im_ref <- posCalc_ref[,1]
y_im_ref <- posCalc_ref[,2]
  
###############################################################################
#                 Alignment with STalign                                      #
#                                                                             #
#   We can access the landmark points points_im_prob_rc and points_im_ref_rc  #  
#   we created in R now in Python using reticulate. We will use these         # 
#   landmarks to compute an initial affine transformation.                    #
###############################################################################
  
py_run_string("
# compute initial affine transformation from points
L,T = STalign.L_T_from_points(r.points_im_prob_rc, r.points_im_ref_rc)
    
im_ref_trans = im_ref_norm.transpose(2,0,1)
Y_im_ref = np.array(range(im_ref_trans.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
X_im_ref = np.array(range(im_ref_trans.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extent_im_ref = STalign.extent_from_x((Y_im_ref,X_im_ref))

# transpose matrices into right dimensions, set up extent
im_prob_trans = im_prob_norm.transpose(2,0,1)
Y_im_prob = np.array(range(im_prob_trans.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
X_im_prob = np.array(range(im_prob_trans.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extent_im_prob = STalign.extent_from_x((Y_im_prob,X_im_prob))
")
  
  
# There are a few parameters that may need to be tuned depending on the size of 
# your tissue and the degree of local diffeomorphism you want to tolerate 
# (such as the smoothness scale of velocity field a), and how long you want to 
# run the algorithm for (such as the number of iterations in the gradient 
# descent niter and the step size epV).

# We are running the transformation with values from a tutorial: Aligning 10X 
# Visium spatial transcriptomics datasets using STalign with Reticulate in R
# {https://jef.works/blog/2023/11/05/aligning-visium-using-STalign-with-reticulate/}

# After aligning, we can apply the learned transform, which includes both an 
# affine and diffeomorphic component, to move the spatial positions of our 
# Visium spots and the histological images from the source dataset to be 
# aligned with the target


py_run_string("
torch.set_default_device('cpu')
device = 'cpu'

# keep all other parameters default
params = {'L':L,'T':T,
       'niter':400,
       'pointsI':r.points_im_prob_rc,
       'pointsJ':r.points_im_ref_rc,
       'device':device,
       'sigmaM':0.15,
       'sigmaB':0.2,
       'sigmaA':0.3,
       'a':300,
       'epV':10,
       'muB': torch.tensor([255,255,255]) # white is background
}

# run LDDMM
out = STalign.LDDMM([Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref],im_ref_trans,**params)

# get necessary output variables
A = out['A']
v = out['v']
xv = out['xv']
")

transformations <- c("affine","affine_points","affine+diffeo")
for (transf in transformations) {
  print(transf)
  if (transf == 'affine') {
    py_run_string("
# modify the source
align_points_im_prob = transform_points_source_with_A(A,np.stack([r.y_im_prob, r.x_im_prob], -1)) # transform the points
align_points_im_prob = align_points_im_prob.numpy() # convert from tensor to numpy to access in R later

im_prob_trans = im_prob_trans.astype(np.float64)
align_image_im_prob = STalign.transform_image_source_with_A(A,[Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref])
align_image_im_prob = align_image_im_prob.permute(1,2,0).numpy()
")
  }
  else if (transf == 'affine_points') {
    py_run_string("
# Generate the A matrix using L, T parameters affine transformation from points
device='cpu'
dtype = torch.float64
Lpoints = torch.tensor(L,device=device,dtype=dtype,requires_grad=False)
Tpoints = torch.tensor(T,device=device,dtype=dtype,requires_grad=False)
  
Apoints = STalign.to_A(Lpoints,Tpoints)
  
# modify the source
align_points_im_prob = transform_points_source_with_A(Apoints,np.stack([r.y_im_prob, r.x_im_prob], -1)) # transform the points
align_points_im_prob = align_points_im_prob.numpy() # convert from tensor to numpy to access in R later

#im_prob_trans = im_prob_trans.astype(np.float64)
align_image_im_prob = STalign.transform_image_source_with_A(Apoints,[Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref])
align_image_im_prob = align_image_im_prob.permute(1,2,0).numpy()
")
  }
  else if (transf == 'affine+diffeo') {
    py_run_string("
# modify the source
align_points_im_prob = STalign.transform_points_source_to_target(xv,v,A,np.stack([r.y_im_prob, r.x_im_prob], -1)) # transform the points
align_points_im_prob = align_points_im_prob.numpy() # convert from tensor to numpy to access in R later

im_prob_trans = im_prob_trans.astype(np.float64)
align_image_im_prob = STalign.transform_image_source_to_target(xv,v,A,[Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref])
align_image_im_prob = align_image_im_prob.permute(1,2,0).numpy()
")
  }
  
  listSTalignResults[[im]][[transf]][["image"]] <- py$align_image_im_prob
  
  posAligned <- py$align_points_im_prob ## currently in row-col order
  posAligned <- data.frame(posAligned[,2:1]) ## put into x-y order
  rownames(posAligned) <- prob_barcodes
  colnames(posAligned) <- c('x', 'y')
  
  posAligned_adapt <- posAligned[prob_barcodes,2:1] / prob_scalefactor
  colnames(posAligned_adapt) <- c("imagerow", "imagecol")
  
  listSTalignResults[[im]][[transf]][["coords"]] <- posAligned_adapt
  
  ###############################################################################
  #           Visualizing the results after alignment                           #
  ###############################################################################
  
  # Graph 1
  df <- data.frame(rbind(cbind(posCalc_prob, time='Prob'),
                         cbind(posAligned, time='Prob aligned'),
                         cbind(posCalc_ref, time='Ref')))
  ggplot(df, aes(x=x, y=y, col=time)) + geom_point(alpha=0.5) + scale_y_reverse() + theme_minimal()
  
  # Graph 2
  par(mfrow=c(1,1))
  plot(c(0,0), xlim=c(0,ncol(py$im_ref_norm)), ylim=c(0,nrow(py$im_ref_norm)),  xlab=NA, ylab=NA, 
       xaxt="n", yaxt="n", bty="n", type="n")
  rasterImage(py$align_image_im_prob, xleft = 0, xright = ncol(py$align_image_im_prob),
              ytop = 0, ybottom = nrow(py$align_image_im_prob), interpolate = FALSE)
  points(posAligned$x, max_yprob - posAligned$y + min_yprob, cex=0.4, col='green', pch=16)
  
  # Graph 3
  par(mfrow=c(1,2))
  
  ## source (reference)
  plot(c(0,0), xlim=c(0,ncol(py$im_ref_norm)), ylim=c(0,nrow(py$im_ref_norm)),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n", bty="n", type="n")
  rasterImage(py$im_ref_norm, xleft = 0, xright = ncol(py$im_ref_norm),
              ytop = 0, ybottom = nrow(py$im_ref_norm), interpolate = FALSE)
  points(posCalc_ref$x, max_yref - posCalc_ref$y + min_yref, pch=16, cex=0.2, col='blue')
  
  ## target (problem)
  plot(c(0,0), xlim=c(0,ncol(py$align_image_im_prob)), ylim=c(0,nrow(py$align_image_im_prob)),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n", bty="n", type="n")
  rasterImage(py$align_image_im_prob, xleft = 0, xright = ncol(py$align_image_im_prob),
              ytop = 0, ybottom = nrow(py$align_image_im_prob), interpolate = FALSE)
  points(posAligned$x, max_yprob - posAligned$y + min_yprob, pch=16, cex=0.2, col='red')
}

im <- im + 1

###############################################################################
###############################################################################
####                     IMPORTANT INFORMATION                             ####
####                                                                       ####
#### Repeat the code of image alignment with STalign until you complete    ####
#### the alignment with all your problem images {go to line 128}           ####
####                                                                       ####
###############################################################################
###############################################################################

###############################################################################
#                   Save the list of aligned results                          #
###############################################################################

saveRDS(listSTalignResults, paste0(saveDir,"patient",patient,"_list_STalign.rds"))
