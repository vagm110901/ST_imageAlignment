library(Matrix)


options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

carpetaData <- "./"
saveDir <-     paste0("./AnnData",patient,"/")

paciente.merge <- readRDS(file = paste0(carpetaData, "Paciente", patient, "_merge.rds"))

# Expression matrixes for each slice
for (i in 1:4) {
  # Expression data = counts
  counts_layer <- paciente.merge@assays$Spatial@layers[[paste0("counts.", i)]]
  counts_layer_T <- t(counts_layer)
  write.csv(as.matrix(counts_layer_T), paste0(saveDir,"expression_matrix_slice", i, ".csv"))
  
  # Cell metadata 
  cell_metadata <- paciente.merge@meta.data[paciente.merge@meta.data$name == paste0("slice19_", i), ]
  write.csv(cell_metadata, paste0(saveDir, "cell_metadata_slice", i, ".csv"))
  
  # Spatial coordinates (row, col, imagerow, imagecol)
  image_coords <- paciente.merge@images[[paste0("slice1.", i)]]@coordinates
  image_coords_df <- as.data.frame(image_coords)  # Convertir en un data.frame
  write.csv(image_coords_df, paste0(saveDir, "image_coordinates_slice", i, ".csv"))
  
  # Extraer los factores de escala (scale.factors)
  scale_factors <- list()
  scale_factors$tissue_hires_scalef <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$hires
  scale_factors$spot_diameter_fullres <- paciente.merge@images[[paste0("slice1.", i)]]@scale.factors$spot
  write.csv(scale_factors, paste0(saveDir, "scale_factors_slice", i, ".csv"))
}

# Gene metadata (shared for all the slices)
write.csv(rownames(paciente.merge@assays$Spatial@features), paste0(saveDir, "gene_metadata.csv"))

# Images for each slice
# I manually copy them from the original directory.


