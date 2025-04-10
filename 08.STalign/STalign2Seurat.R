library(Seurat)
library(ggplot2)

# Load a custom function for alignment evaluation
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_evalAlign.R")
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_selectCoord.R")

options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/08.STalign/"

paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))
listSTalignResults <- readRDS(file = paste0(saveDir, "patient",patient,"_list_STalign.rds"))

############################################################################### 
#           Save the aligned data (image & spots) to Seurat object            #
############################################################################### 

paciente.merge@images[["slice1"]] <- paciente.merge@images[["slice1.1"]]

for (i in c(2:4)) {
  paciente.merge@images[[paste0("slice",i)]] <- paciente.merge@images[[paste0("slice1.",i)]]
  paciente.merge@images[[paste0("slice",i)]]@image <- listSTalignResults[[i]][["affine+diffeo"]][["image"]]
  paciente.merge@images[[paste0("slice",i)]]@coordinates[,4:5] <- listSTalignResults[[i]][["affine+diffeo"]][["coords"]]
}

saveRDS(paciente.merge, file = paste0(saveDir,"patient",patient,"_merge_STalign.rds"))

############################################################################### 
#                Create objects for deconvolution                             #
############################################################################### 

# Obtain the individual object of the slide
# Splits the Seurat object 'paciente19.merge' into a list of Seurat objects based on the 'name' metadata
pacientes.split.orig <- Seurat::SplitObject(paciente.merge, split.by = "name")

# Create a copy of the split objects for transformed images
pacientes.split.trans <- pacientes.split.orig

# Assign transformed images and coordinates to the corresponding objects in the new split list
pacientes.split.trans[[2]]@images$slice1.2@image <- listSTalignResults[[2]][["affine+diffeo"]][["image"]]
pacientes.split.trans[[3]]@images$slice1.3@image <- listSTalignResults[[3]][["affine+diffeo"]][["image"]]
pacientes.split.trans[[4]]@images$slice1.4@image <- listSTalignResults[[4]][["affine+diffeo"]][["image"]]

pacientes.split.trans[[2]]@images$slice1.2@coordinates[,4:5] <- listSTalignResults[[2]][["affine+diffeo"]][["coords"]]
pacientes.split.trans[[3]]@images$slice1.3@coordinates[,4:5] <- listSTalignResults[[3]][["affine+diffeo"]][["coords"]]
pacientes.split.trans[[4]]@images$slice1.4@coordinates[,4:5] <- listSTalignResults[[4]][["affine+diffeo"]][["coords"]]

# Loop through each of the split original objects starting from the second one
for (i in 2:length(pacientes.split.trans)) {
  listaObjDeconv <- list()
  
  # Create a list to hold reference and problematic objects for deconvolution
  listaObjDeconv[["reference"]] <- pacientes.split.orig[[1]]
  listaObjDeconv[["problemNOTalign"]] <- pacientes.split.orig[[i]]
  listaObjDeconv[["problemYESalign"]] <- pacientes.split.trans[[i]]
  
  # Update the name in the metadata for each object to reflect its status
  listaObjDeconv[["reference"]]@meta.data$name <- gsub(listaObjDeconv[["reference"]]@meta.data$name[[1]], "reference", listaObjDeconv[["reference"]]@meta.data$name)
  listaObjDeconv[["problemNOTalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemNOTalign"]]@meta.data$name[[1]], "problemNOTalign", listaObjDeconv[["problemNOTalign"]]@meta.data$name)
  listaObjDeconv[["problemYESalign"]]@meta.data$name <- gsub(listaObjDeconv[["problemYESalign"]]@meta.data$name[[1]], "problemYESalign", listaObjDeconv[["problemYESalign"]]@meta.data$name)
  
  savename <- paste0("patient",patient,"_merge_STalign_list_im", i, ".rds")
  saveRDS(listaObjDeconv, paste0(saveDir, savename))
}

########################################## TRY

i <- 3
savename <- paste0("patient",patient,"_merge_STalign_list_im", i, ".rds")
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

pl1 <- SpatialDimPlot(listaObjDeconv[[2]], alpha = 0, crop = FALSE, ncol = 2) & 
  theme(legend.position = "none")
pl2 <- SpatialPlot(object = listaObjDeconv[[2]], 
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

listaCoordenadasNEW <- list()
for (i in seq_along(1:4)) {
  coordenadas <- selectCoordSTalign(paciente.merge@images[[paste0("slice",i)]]@image)
  for (j in seq_along(coordenadas)) {
    coordenadas[[j]] <- round(coordenadas[[j]])
  }
  listaCoordenadasNEW[[i]] <- coordenadas
}

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

saveRDS(Evaluation, paste0(saveDir, "/patient",patient,"_EVALUATION_merge_STalign.rds"))
