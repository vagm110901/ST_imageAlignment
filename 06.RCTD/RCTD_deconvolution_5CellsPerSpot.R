library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)

library(dplyr)
# remotes::install_version("Seurat", version = "5.0.2")
library(Seurat) 
library(patchwork)

library(ggplot2)
library(sctransform)

###############################################################################
as_AssayObject <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$sub_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table::data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = r[[i]]$cell_type_list,
          weight = r[[i]]$sub_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      swm <- as.matrix(swd[, -1])
      rownames(swm) <- swd$barcode
      swm <- t(swm)
      swm <- spacexr::normalize_weights(swm)
      #swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      swm <- as(swm, "sparseMatrix")
      return(CreateAssayObject(data = swm))
    }
  } else if (length(r) == 1) {
    m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
  }
}

as_AssayObject_complete <- function(object) {
  if (is(object, "RCTD")) {
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      stop("Install spacexr.")
    }
  } else {
    stop("Only RCTD objects supported")
  }
  r <- object@results
  if (length(r) > 1) {
    if (!is.null(r[[1]]$all_weights)) {
      sw <- data.table::rbindlist(lapply(seq_along(r), function(i)
        data.table::data.table(
          barcode = colnames(object@spatialRNA@counts)[i],
          cell_type = names(r[[i]]$all_weights),
          weight = r[[i]]$all_weights
        )), fill = TRUE)
      sw$cell_type[is.na(sw$cell_type)] <- "unassigned"
      swd <- data.table::dcast(sw, barcode ~ cell_type, value.var = "weight", fill = 0)
      swm <- as.matrix(swd[, -1])
      rownames(swm) <- swd$barcode
      swm <- t(swm)
      swm <- spacexr::normalize_weights(swm)
      #swm <- rbind(swm, max = apply(swm[!rownames(swm) %in% "unassigned", ], 2, max))
      swm <- as(swm, "sparseMatrix")
      return(CreateAssayObject(data = swm))
    }
  } else if (length(r) == 1) {
    m <- t(spacexr::normalize_weights(as.matrix(r$weights)))
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
  }
}

###############################################################################

saveDir <- "/home/vgaya/ST_imageAlignment/06.RCTD"

reference <- readRDS(paste0(saveDir,'/SCRef_YadavAnnot.rds'))

patients <- c("19","MIX")
#patients <- c("19")
#modes <- c("alignTR","procrustes","imageJ","alignTR_S","procrustes_S","imageJ_S")
modes <- c("alignTR")
#ims <- c("2","3","4")
ims <- c("3")

for (patient in patients) {
  for (mode in modes) {
    if (mode == "imageJ" | mode == "imageJ_S") {
      dataDir <- "/home/vgaya/ST_imageAlignment/05.ImageJ" } else ( 
        dataDir <- "/home/vgaya/ST_imageAlignment/04.ImageAlignment")
    for (im in ims) {
      # Spatial Transcriptomics data 
      listaObjST <- readRDS(paste0(dataDir, "/patient", patient, "_merge_", 
                                   mode, "_list_im", im, ".rds"))
      
      for (i in 1:length(listaObjST)) {
        objectST <- listaObjST[[i]]
        # we normalize the objects
        objectST <- SCTransform(objectST, assay = "Spatial", verbose = FALSE) %>%
          RunPCA(verbose = FALSE)
        
        objname <- names(listaObjST[i])
        listaObjST[[objname]] <- objectST
      }
      
      m1 <- SpatialDimPlot(listaObjST[[1]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      m2 <- SpatialDimPlot(listaObjST[[2]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      m3 <- SpatialDimPlot(listaObjST[[3]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      print(m1 + m2 + m3)
      
      listaObjST_RCTD <- list()
      
      for (i in 1:length(listaObjST)) {
        # 1. coords: A numeric data.frame (or matrix) representing the spatial pixel 
        # locations. rownames are barcodes/pixel names, and there should be two columns 
        # for ‘x’ and for ‘y’.
        coords <- listaObjST[[i]]@images[[1]]@coordinates[,c("imagecol", "imagerow")]
        colnames(coords) <- c("x", "y") 
        
        
        # 2. counts: A matrix (or dgCmatrix) representing Digital Gene Expression (DGE). 
        # Rownames should be genes and colnames represent barcodes/pixel names. 
        # Counts should be untransformed count-level data.
        counts <- listaObjST[[i]]@assays[["SCT"]]@counts
        
        # 3. nUMI: Optional, a named (by pixel barcode) list of total counts or UMI’s 
        # appearing at each pixel. If not provided, nUMI will be assumed to be the total
        # counts appearing on each pixel.
        nUMI <- colSums(counts)
        
        STsample <- SpatialRNA(coords, counts, nUMI)
        
        objname <- names(listaObjST[i])
        listaObjST_RCTD[[objname]] <- STsample
        
        # RCTD
        # Create RCTD object
        myRCTD <- create.RCTD(listaObjST_RCTD[[i]], reference, max_cores = 8, UMI_min = 10, MAX_MULTI_TYPES = 5)
        # Run RCTD
        myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')
        
        # results
        assay_myRCTD <- as_AssayObject(myRCTD)
        
        listaObjST[[i]]@assays[["deconvolution.RCTD"]] <- assay_myRCTD

        assay_myRCTD_complete <- as_AssayObject_complete(myRCTD)

        listaObjST[[i]]@assays[["deconvolution.RCTD.complete"]] <- assay_myRCTD_complete
      }
      
      # Save the object
      saveRDS(listaObjST, file.path(saveDir,paste0('/ST_RCTD_YadavAnnot_patient', 
                                                        patient, "_merge_", mode, 
                                                        "_list_im", im, "_5CellsPerSpot.rds")))
      
      #create_RCTD_plots(myRCTD, paste0(saveDir,"/plots"))
    }
  }
}

