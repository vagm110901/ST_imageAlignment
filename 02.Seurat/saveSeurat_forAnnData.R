options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
saveDir <-     paste0("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/AnnData",patient,"/")

paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))

samples <- c("19","43","47","45")
# Expression matrixes for each slice
for (i in 1:4) {
  
  # Expression data = counts
  counts_layer <- paciente.merge@assays$Spatial@layers[[paste0("counts.", i)]]
  counts_layer_T <- t(counts_layer)
  write.csv(as.matrix(counts_layer_T), paste0(saveDir,"expression_matrix_slice", i, ".csv"))
  
  # Cell metadata 
  if (patient == "MIX") {
    sample <- samples[i]
    cell_metadata <- paciente.merge@meta.data[paciente.merge@meta.data$name == paste0("slice",sample,"_1"), ]
  } else if (patient == "19") {
    cell_metadata <- paciente.merge@meta.data[paciente.merge@meta.data$name == paste0("slice19_", i), ]
  }
  write.csv(cell_metadata, paste0(saveDir, "cell_metadata_slice", i, ".csv"))
  
  # Spatial coordinates (row, col, imagerow, imagecol)
  image_coords <- paciente.merge@images[[paste0("slice1.", i)]]@coordinates
  image_coords_df <- as.data.frame(image_coords)  # Convertir en un data.frame
  write.csv(image_coords_df, paste0(saveDir, "image_coordinates_slice", i, ".csv"))
  
  # Extraer los factores de escala (scale.factors)
  scale_factors <- list()
  scale_factors$tissue_lowres_scalef <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$lowres
  scale_factors$tissue_hires_scalef <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$hires
  scale_factors$spot_diameter_fullres <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$spot
  scale_factors$fiducial_diameter_fullres <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$fiducial
  write.csv(scale_factors, paste0(saveDir, "scale_factors_slice", i, ".csv"))
  
  # Gene metadata (shared for all the slices)
  write.csv(rownames(paciente.merge@assays$Spatial@features[which(paciente.merge@assays$Spatial@features[,i] == TRUE),]), 
            paste0(saveDir, "gene_metadata", i, ".csv"))
}

# Images for each slice
# I manually copy them from the original directory.


