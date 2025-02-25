#load required libraries
#library(Seurat)
#library(SeuratObject)
library(SeuratDisk)
library(hdf5r)

# Define the files' folders 
carpetaData <- "/home/vgaya/ST_imageAlignment/02.Seurat/"

# Load the merged Seurat object 
options <- c("19", "MIX")
#choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
#patient <- options[choice]
patient <- options[1]
paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))

paciente.merge[["Spatial"]] <- as(object = paciente.merge[["Spatial"]], Class = "Assay")

# Convert Seurat to AnnData
SaveH5Seurat(paciente.merge, filename = paste0(carpetaData,"Paciente19_merge.h5Seurat"))
Convert("Paciente19_merge.h5Seurat", dest = "h5ad")

patient <- options[2]
paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))

paciente.merge[["Spatial"]] <- as(object = paciente.merge[["Spatial"]], Class = "Assay")

# Convert Seurat to AnnData
SaveH5Seurat(paciente.merge, filename = paste0(carpetaData,"PacienteMIX_merge.h5Seurat"), overwrite = TRUE)
Convert("PacienteMIX_merge.h5Seurat", dest = "h5ad")