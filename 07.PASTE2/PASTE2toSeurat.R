library(ggplot2)

# Load a custom function for alignment evaluation
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_evalAlign.R")
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_selectCoord.R")

escale <- function(new_coords, patient, N, objSeurat) {
  if (patient == "MIX") {
    # Select the max values changed (row = col & col = row)
    newCOLmax <- max(new_coords$imagerow)
    newROWmax <- max(new_coords$imagecol)
    
    origCOLmax <- max(objSeurat@images[[paste0("slice1.", N)]]@coordinates$imagecol)
    origROWmax <- max(objSeurat@images[[paste0("slice1.", N)]]@coordinates$imagerow)
    
    new_coords$imagerow <- new_coords$imagerow / newCOLmax * origCOLmax
    new_coords$imagecol <- new_coords$imagecol / newROWmax * origROWmax
  }
  return(new_coords)
}

options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
saveDir <- paste0("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/07.PASTE2/")

paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))
paciente.merge@images[["slice1"]] <- paciente.merge@images[["slice1.1"]]

Ns = c("2","3","4")
for (N in Ns) {
  #new_coords1 <- read.csv(file = paste0(saveDir, "Paciente", patient, "_1_align1", N, "_h_coord.csv"))
  new_coordsProb <- read.csv(file = paste0(saveDir, "Paciente", patient, "_", N, "_align1", N, "_h_coord.csv"))
  
  #paciente.merge@images[["slice1.1"]]@coordinates$imagerow <- unlist(new_coords1$imagecol)
  #paciente.merge@images[["slice1.1"]]@coordinates$imagecol <- unlist(new_coords1$imagerow)
  
  new_coordsProb <- escale(new_coordsProb, patient, N, paciente.merge)
  
  paciente.merge@images[[paste0("slice",N)]] <- paciente.merge@images[[paste0("slice1.", N)]]
  paciente.merge@images[[paste0("slice", N)]]@coordinates$imagerow <- unlist(new_coordsProb$imagecol)
  paciente.merge@images[[paste0("slice", N)]]@coordinates$imagecol <- unlist(new_coordsProb$imagerow)
}

saveRDS(paciente.merge, file = paste0(saveDir,"patient",patient,"_merge_PASTE2.rds"))

############################################################################### 
#                Create objects for deconvolution                             #
############################################################################### 

# Obtain the individual object of the slide
# Splits the Seurat object 'paciente19.merge' into a list of Seurat objects based on the 'name' metadata
pacientes.seurat.split.orig <- SplitObject(readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds")), split.by = "name")

# Create a copy of the split objects for transformed images
pacientes.seurat.split.trans <- pacientes.seurat.split.orig

# Assign transformed images to the corresponding objects in the new split list
pacientes.seurat.split.trans[[2]]@images$slice1.2 <- paciente.merge@images[[2]]
pacientes.seurat.split.trans[[3]]@images$slice1.3 <- paciente.merge@images[[3]]
pacientes.seurat.split.trans[[4]]@images$slice1.4 <- paciente.merge@images[[4]]

# Loop through each of the split original objects starting from the second one
for (i in 2:length(pacientes.seurat.split.orig)) {
  listaObjDeconv <- list()
  
  # Create a list to hold reference and problematic objects for deconvolution
  listaObjDeconv[["reference"]] <- pacientes.seurat.split.orig[[1]]
  listaObjDeconv[["problemNOTalign"]] <- pacientes.seurat.split.orig[[i]]
  listaObjDeconv[["problemYESalign"]] <- pacientes.seurat.split.trans[[i]]
  
  # Update the name in the metadata for each object to reflect its status
  listaObjDeconv[["reference"]]@meta.data$name <- gsub(listaObjDeconv[["reference"]]@meta.data$name[[1]], "reference", listaObjDeconv[["reference"]]@meta.data$name)
  listaObjDeconv[["problemNOTalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemNOTalign"]]@meta.data$name[[1]], "problemNOTalign", listaObjDeconv[["problemNOTalign"]]@meta.data$name)
  listaObjDeconv[["problemYESalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemYESalign"]]@meta.data$name[[1]], "problemYESalign", listaObjDeconv[["problemYESalign"]]@meta.data$name)
  
  savename <- paste0("patient",patient,"_merge_PASTE2_list_im", i, ".rds")
  saveRDS(listaObjDeconv, paste0(saveDir, savename))
}

####################### TRY

patient <- "MIX"
i <- 4
savename <- paste0("patient",patient,"_merge_PASTE2_list_im", i, ".rds")
listaObjDeconv <- readRDS(paste0(saveDir, savename))


###############################################################################
#           Visualization                                                     #
###############################################################################
library(Seurat)
SpatialDimPlot(paciente.merge, alpha = 0.5, crop = FALSE, ncol = 2)

# Plot spatial expression of the "percent.butterfly" feature
SpatialPlot(object = paciente.merge, 
            features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 2) & 
  theme(legend.position = "none")

pl1 <- SpatialDimPlot(listaObjDeconv[[3]], alpha = 0, crop = FALSE, ncol = 2) & 
  theme(legend.position = "none")
pl2 <- SpatialPlot(object = listaObjDeconv[[3]], 
            features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 1) & 
  theme(legend.position = "none")
pl1 + pl2

###############################################################################
#           Evaluation of the alignment                                       #
###############################################################################

listaCoordenadas <- list()
for (i in seq_along(1:4)) {
  coordenadas <- selectCoordSTalign(paciente.merge@images[[paste0("slice1.",i)]]@image)
  for (j in seq_along(coordenadas)) {
    coordenadas[[j]] <- round(coordenadas[[j]])
  }
  listaCoordenadas[[i]] <- coordenadas
}

listaCoordenadasNEW <- listaCoordenadas

listaRawImages <- list()
listaTransImages <- list()
for (i in 1:length(listaCoordenadasNEW)) {
  listaRawImages[[i]] <- paciente.merge@images[[paste0("slice1.",i)]]@image
  listaTransImages[[i]] <- paciente.merge@images[[paste0("slice",i)]]@image
}

if (patient == "19") {
  patientType <- "unique"
} else if (patient == "MIX") {
  patientType <- "multiple"
}

Evaluation <- evaluationComplete(listaCoordenadasNEW, listaCoordenadas,
                                 listaRawImages, listaTransImages, patient = patientType)

saveRDS(Evaluation, paste0(saveDir, "/patient",patient,"_EVALUATION_merge_PASTE2.rds"))
