library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

# Set the maximum size for global objects
options(future.globals.maxSize = 100000 * 1024^2)

counts <- read.csv("./GSE190442_aggregated_counts_postqc.csv", row.names = 1)
counts_matrix <- as.matrix(counts)

metadata <- read.csv("./GSE190442_aggregated_metadata_postqc.csv", row.names = 1)

seurat_object <- CreateSeuratObject(counts = counts_matrix, project = "Yadav")

# Verificar coincidencias
all(rownames(metadata) %in% colnames(seurat_object))
all(rownames(metadata) == colnames(seurat_object))

# Agregar los metadatos
seurat_object <- AddMetaData(object = seurat_object, metadata = metadata)

# Revisar las cuentas
#head(seurat_object@assays$RNA@counts)

# Revisar los metadatos
#head(seurat_object@meta.data)

saveRDS(seurat_object, file.path("./",'snRNA_QC_Yadav.rds'))
