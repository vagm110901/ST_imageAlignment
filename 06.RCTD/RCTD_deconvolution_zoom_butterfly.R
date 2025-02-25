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
      swm <- t(spacexr::normalize_weights(swm))
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
      swm <- t(spacexr::normalize_weights(swm))
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

#patients <- c("19","MIX")
patients <- c("MIX")
#modes <- c("alignTR","procrustes","imageJ","alignTR_S","procrustes_S","imageJ_S")
#modes <- c("alignTR","procrustes","imageJ")
modes <- c("imageJ")
ims <- c("2","3","4")
#ims <- c("4")

for (patient in patients) {
  for (mode in modes) {
    if (mode == "imageJ" | mode == "imageJ_S") {
      dataDir <- "/home/vgaya/ST_imageAlignment/05.ImageJ" } else { 
        dataDir <- "/home/vgaya/ST_imageAlignment/04.ImageAlignment"}
    for (im in ims) {
      # Spatial Transcriptomics data 
      listaObjST <- readRDS(paste0(dataDir, "/patient", patient, "_merge_", 
                                   mode, "_list_im", im, ".rds"))
      
      # zoom
      objectST <- listaObjST[[1]]
      
      imagerow <- objectST@images[[1]]@coordinates$imagerow 
      imagecol <- objectST@images[[1]]@coordinates$imagecol

      ################################
      imagerowmin <- min(imagerow[objectST@images[[1]]@coordinates$row > 40 & 
                                  objectST@images[[1]]@coordinates$row < 63])
      imagerowmax <- max(imagerow[objectST@images[[1]]@coordinates$row > 40 & 
                                  objectST@images[[1]]@coordinates$row < 63])
      imagecolmin <- min(imagecol[objectST@images[[1]]@coordinates$col > 15 & 
                                  objectST@images[[1]]@coordinates$col < 56])
      imagecolmax <- max(imagecol[objectST@images[[1]]@coordinates$col > 15 & 
                                  objectST@images[[1]]@coordinates$col < 56])      

      rel_imagerowmin <- (imagerowmin - min(imagerow)) / (max(imagerow) - min(imagerow))
      rel_imagerowmax <- (imagerowmax - min(imagerow)) / (max(imagerow) - min(imagerow)) 
      rel_imagecolmin <- (imagecolmin - min(imagecol)) / (max(imagecol) - min(imagecol)) 
      rel_imagecolmax <- (imagecolmax - min(imagecol)) / (max(imagecol) - min(imagecol))    


      # Loop through ST objects and create subsetted, zoomed objects
      for (i in seq_along(listaObjST)) {
         objectST <- listaObjST[[i]]
  
         # Define cells within subset area
         imagerow <- objectST@images[[1]]@coordinates$imagerow
         imagecol <- objectST@images[[1]]@coordinates$imagecol
  
         if (i == 3) {
            objectST_orig <- listaObjST[[2]]
            imagerow_orig <- objectST_orig@images[[1]]@coordinates$imagerow
            imagecol_orig <- objectST_orig@images[[1]]@coordinates$imagecol
            NEWimagerowmin <- rel_imagerowmin * (max(imagerow_orig) - min(imagerow_orig)) + min(imagerow_orig)
            NEWimagerowmax <- rel_imagerowmax * (max(imagerow_orig) - min(imagerow_orig)) + min(imagerow_orig)
            NEWimagecolmin <- rel_imagecolmin * (max(imagecol_orig) - min(imagecol_orig)) + min(imagecol_orig)
            NEWimagecolmax <- rel_imagecolmax * (max(imagecol_orig) - min(imagecol_orig)) + min(imagecol_orig)
         } else {
            NEWimagerowmin <- rel_imagerowmin * (max(imagerow) - min(imagerow)) + min(imagerow)
            NEWimagerowmax <- rel_imagerowmax * (max(imagerow) - min(imagerow)) + min(imagerow)
            NEWimagecolmin <- rel_imagecolmin * (max(imagecol) - min(imagecol)) + min(imagecol)
            NEWimagecolmax <- rel_imagecolmax * (max(imagecol) - min(imagecol)) + min(imagecol)
         }

         grupo <- rownames(objectST@images[[1]]@coordinates[which(imagecol >= NEWimagecolmin & imagecol <= NEWimagecolmax &
                                                             imagerow >= NEWimagerowmin & imagerow <= NEWimagerowmax),])

          # Subset and re-normalize each ST object
          if (length(grupo) != 0) {          
          objectST.zoom <- subset(objectST, cells = grupo, invert = FALSE) %>%
                                  SCTransform(assay = "Spatial", verbose = FALSE) 
                                  # %>%
                                  # RunPCA(verbose = FALSE)
          objname <- names(listaObjST[i])
          listaObjST[[objname]] <- objectST.zoom
        } else { 
          objname <- names(listaObjST[i]) 
          listaObjST[[objname]] <- NA
        }
      }
      
      #m1 <- SpatialDimPlot(listaObjST[[1]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      #m2 <- SpatialDimPlot(listaObjST[[2]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      #m3 <- SpatialDimPlot(listaObjST[[3]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
      #print(m1 + m2 + m3)
      
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
        myRCTD <- create.RCTD(listaObjST_RCTD[[i]], reference, max_cores = 8, UMI_min = 10)
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
                                                        "_list_im", im, "_butterfly.rds")))
      
      #create_RCTD_plots(myRCTD, paste0(saveDir,"/plots"))
    }
  }
}

