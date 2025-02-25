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
library(RColorBrewer)


saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/06.RCTD"

cell.types <- c("Astrcytes", "Microglia", "OPC", "Endothelial", 
                "Neurons", "Pericytes", "Schwann", "Lymphocytes", 
                "Oligodendrocytes", "Ependymal Cells", "Meninges")

cell_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                 "#984EA3", "#FFFF33", "#A65628", "#F781BF",
                 "#999999", "#66C2A5", "#FC8D62")



#listaObjAnnot <- readRDS(paste0(saveDir, '/ST_RCTD_patient19_merge_alignTR_list_im2.rds'))
listaObjAnnot <- readRDS(paste0(saveDir, '/ST_RCTD_YadavAnnot_patientMIX_merge_alignTR_list_im4.rds'))
#listaObjAnnot <- readRDS(paste0(saveDir, '/ST_RCTD_YadavAnnot_NoNorm_patientMIX_merge_alignTR_list_im3.rds'))

#listaObjAnnot <- readRDS(paste0(saveDir, '/ST_RCTD_YadavAnnot_patient19_merge_alignTR_list_im3_5CellsPerSpot.rds'))
#listaObjAnnot <- readRDS(paste0(saveDir, '/ST_RCTD_YadavAnnot_patient19_merge_alignTR_list_im3_10CellMin.rds'))


# multi
for (i in 1:length(listaObjAnnot)) {DefaultAssay(listaObjAnnot[[i]]) <- "deconvolution.RCTD"; RCTDmode <- "deconvolution.RCTD"}
# full
for (i in 1:length(listaObjAnnot)) {DefaultAssay(listaObjAnnot[[i]]) <- "deconvolution.RCTD.complete"; RCTDmode <- "deconvolution.RCTD.complete"}

#objectST <- listaObjAnnot[[1]]

############################################################################### 
#         all_weights only with 4 more abundant cell.types                    #
############################################################################### 
library(foreach)
library(doParallel)

N <- length(cell.types)
nombres_originales <- names(listaObjAnnot)

# Define the number of cores
numCores <- detectCores() - 2
cl <- makeCluster(numCores)

clusterExport(cl, varlist = c("listaObjAnnot", "N"))
clusterEvalQ(cl, library(Matrix))  # Cargar paquetes en cada núcleo

# Paralell
listaObjAnnot <- parLapply(cl, names(listaObjAnnot), function(i) {
  data_mat <- listaObjAnnot[[i]]@assays[["deconvolution.RCTD.complete"]]@data
  
  for (j in 1:ncol(data_mat)) {
    idx <- names(sort(data_mat[, j], decreasing = FALSE)[1:(N - 4)])
    data_mat[idx, j] <- 0
    
    data_mat[, j] <- data_mat[, j] / sum(data_mat[, j])
  }
  listaObjAnnot[[i]]@assays[["deconvolution.RCTD.complete"]]@data <- data_mat
  return(listaObjAnnot[[i]])
})
stopCluster(cl)

names(listaObjAnnot) <- nombres_originales

###############################################################################
#                    Example plots                                            #
############################################################################### 
q1 <- SpatialFeaturePlot(listaObjAnnot[[1]], features = "Neurons", pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0,1)) + 
   scale_fill_gradientn(colors = c("transparent", "#984EA3"),
                        limits = c(0, 1),
                        oob = scales::squish)

q2 <- SpatialFeaturePlot(listaObjAnnot[[1]], features = "Oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE, alpha = c(0,1)) + 
  scale_fill_gradientn(colors = c("transparent", "#999999"), 
                       limits = c(0, 1), 
                       oob = scales::squish)

q3 <- SpatialFeaturePlot(listaObjAnnot[[2]], features = "Neurons", pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0,1)) + 
  scale_fill_gradientn(colors = c("transparent", "#984EA3"), 
                       limits = c(0, 1), 
                       oob = scales::squish)
q4 <- SpatialFeaturePlot(listaObjAnnot[[2]], features = "Oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE, alpha = c(0,1)) + 
  scale_fill_gradientn(colors = c("transparent", "#999999"), 
                       limits = c(0, 1), 
                       oob = scales::squish)

q5 <- SpatialFeaturePlot(listaObjAnnot[[3]], features = "Neurons", pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0,1)) + 
  scale_fill_gradientn(colors = c("transparent", "#984EA3"), 
                       limits = c(0, 1), 
                       oob = scales::squish)
q6 <- SpatialFeaturePlot(listaObjAnnot[[3]], features = "Oligodendrocytes", pt.size.factor = 1, ncol = 1, crop = FALSE, alpha = c(0,1)) + 
  scale_fill_gradientn(colors = c("transparent", "#999999"), 
                       limits = c(0, 1), 
                       oob = scales::squish)


png(paste0(saveDir, "/results/patient19_merge_imageJ_S_deconvORIG.png"))
q1 + q2
dev.off()

png(paste0(saveDir, "/results/patient19_merge_imageJ_S_deconv4NO.png"))
q3 + q4
dev.off()

png(paste0(saveDir, "/results/patient19_merge_imageJ_S_deconv4SI.png"))
q5 + q6
dev.off()

###############################################################################
#               One by one cell types visualization                           #
############################################################################### 

listaGraphs <- list()
for (i in 1:length(cell.types)) {
  for (j in 1:length(listaObjAnnot)) {
  listaGraphs[[names(listaObjAnnot[j])]][[cell.types[[i]]]] <- SpatialFeaturePlot(listaObjAnnot[[j]], 
                                                       features = cell.types[[i]], 
                                                       pt.size.factor = 1.2, 
                                                       crop = FALSE, 
                                                       alpha = c(0,1)) + 
    scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                         limits = c(0, 1), 
                         oob = scales::squish)
  }
}

referenceG <- listaGraphs[[1]]
notAlignG <- listaGraphs[[2]]
alignG <- listaGraphs[[3]]

( referenceG$Neurons + referenceG$Oligodendrocytes + 
    referenceG$Astrcytes + referenceG$Microglia +
    referenceG$`Ependymal Cells` + referenceG$Meninges + 
    referenceG$Lymphocytes + referenceG$Pericytes +
    referenceG$OPC + referenceG$Endothelial + referenceG$Schwann )

( notAlignG$Neurons + notAlignG$Oligodendrocytes + 
  notAlignG$Astrcytes + notAlignG$Microglia +
  notAlignG$`Ependymal Cells` + notAlignG$Meninges + 
  notAlignG$Lymphocytes + notAlignG$Pericytes +
  notAlignG$OPC + notAlignG$Endothelial + notAlignG$Schwann )


( referenceG$Neurons + referenceG$Oligodendrocytes + 
  referenceG$`Ependymal Cells` + referenceG$Astrcytes +
  notAlignG$Neurons + notAlignG$Oligodendrocytes + 
  notAlignG$`Ependymal Cells` + notAlignG$Astrcytes + 
  alignG$Neurons + alignG$Oligodendrocytes + 
  alignG$`Ependymal Cells` + alignG$Astrcytes ) 



###############################################################################
# Visualization with STdeconvolve: pie graphs                                 #
############################################################################### 
# To simplify the pie graphs
for (j in 1:length(listaObjAnnot)) {
  for (i in 1:length(listaObjAnnot[[j]]@assays[["deconvolution.RCTD.complete"]]@data)) {
    listaObjAnnot[[j]]@assays[["deconvolution.RCTD.complete"]]@data[i] <- 
      round(listaObjAnnot[[j]]@assays[["deconvolution.RCTD.complete"]]@data[i], 1)
  }
}



require(remotes)
remotes::install_github('JEFworks-Lab/STdeconvolve')
library(STdeconvolve)
# https://jef.works/STdeconvolve/

theta <- t(as.matrix(listaObjAnnot[[1]]@assays[["deconvolution.RCTD.complete"]]@data))
pos <- listaObjAnnot[[1]]@images[[1]]@coordinates[,c("imagecol", "imagerow")]
colnames(pos) <- c("x", "y") 
pos$y <- -pos$y
image <- listaObjAnnot[[1]]@images[[1]]@image


plt <- vizAllTopics(theta = theta,
                    pos = pos,
                    topicOrder = cell.types, 
                    topicCols = cell_colors,
                    r = 100, #130
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = NA,
                    ) +
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")

plt

###############################################################################
###############################################################################
##                                                                           ##
##                    Zoom visualization                                     ##
##                                                                           ##
############################################################################### 
############################################################################### 
# multi
for (i in 1:length(listaObjAnnot)) {DefaultAssay(listaObjAnnot[[i]]) <- "deconvolution.RCTD"; RCTDmode <- "deconvolution.RCTD"}
# full
for (i in 1:length(listaObjAnnot)) {DefaultAssay(listaObjAnnot[[i]]) <- "deconvolution.RCTD.complete"; RCTDmode <- "deconvolution.RCTD.complete"}

objectST1 <- listaObjAnnot[[1]]

imagerow1 <- objectST1@images[[1]]@coordinates$imagerow
imagecol1 <- objectST1@images[[1]]@coordinates$imagecol

min_imagerow1 <- min(imagerow1)
max_imagerow1 <- max(imagerow1)
min_imagecol1 <- min(imagecol1)
max_imagecol1 <- max(imagecol1)

################################ subset canal central #########################
imagerowmin <- min(imagerow1[objectST1@images[[1]]@coordinates$row > 18 & 
                             objectST1@images[[1]]@coordinates$row < 41])
imagerowmax <- max(imagerow1[objectST1@images[[1]]@coordinates$row > 18 & 
                             objectST1@images[[1]]@coordinates$row < 41])
imagecolmin <- min(imagecol1[objectST1@images[[1]]@coordinates$col > 53 & 
                             objectST1@images[[1]]@coordinates$col < 97])
imagecolmax <- max(imagecol1[objectST1@images[[1]]@coordinates$col > 56 & 
                             objectST1@images[[1]]@coordinates$col < 97])
###############################################################################

################################ subset extremo inferior mariposa #############
imagerowmin <- min(imagerow1[objectST1@images[[1]]@coordinates$row > 40 & 
                             objectST1@images[[1]]@coordinates$row < 63])
imagerowmax <- max(imagerow1[objectST1@images[[1]]@coordinates$row > 40 & 
                             objectST1@images[[1]]@coordinates$row < 63])
imagecolmin <- min(imagecol1[objectST1@images[[1]]@coordinates$col > 15 & 
                             objectST1@images[[1]]@coordinates$col < 56])
imagecolmax <- max(imagecol1[objectST1@images[[1]]@coordinates$col > 15 & 
                             objectST1@images[[1]]@coordinates$col < 56])
###############################################################################

rel_imagerowmin <- (imagerowmin - min_imagerow1) / (max_imagerow1 - min_imagerow1)
rel_imagerowmax <- (imagerowmax - min_imagerow1) / (max_imagerow1 - min_imagerow1) 
rel_imagecolmin <- (imagecolmin - min_imagecol1) / (max_imagecol1 - min_imagecol1) 
rel_imagecolmax <- (imagecolmax - min_imagecol1) / (max_imagecol1 - min_imagecol1) 


# Initialize lists for zoomed ST objects and images
listaObjZoomAnnot <- list()

# Loop through ST objects and create subsetted, zoomed objects
for (i in seq_along(listaObjAnnot)) {
  objectST <- listaObjAnnot[[i]]
  imagerow <- objectST@images[[1]]@coordinates$imagerow
  imagecol <- objectST@images[[1]]@coordinates$imagecol
  
  # If the image has been aligned, the coordinates to make a subset must be 
  # selected from the image before the alignment
  if (i == 3) {
    objectST_orig <- listaObjAnnot[[2]]
    imagerow_orig <- objectST_orig@images[[1]]@coordinates$imagerow
    imagecol_orig <- objectST_orig@images[[1]]@coordinates$imagecol
    
    min_imagerow_orig <- min(imagerow_orig)
    max_imagerow_orig <- max(imagerow_orig)
    min_imagecol_orig <- min(imagecol_orig)
    max_imagecol_orig <- max(imagecol_orig)
    
    NEWimagerowmin <- rel_imagerowmin * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
    NEWimagerowmax <- rel_imagerowmax * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
    NEWimagecolmin <- rel_imagecolmin * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
    NEWimagecolmax <- rel_imagecolmax * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
  
    } else {
    min_imagerow <- min(imagerow)
    max_imagerow <- max(imagerow)
    min_imagecol <- min(imagecol)
    max_imagecol <- max(imagecol)
    
    NEWimagerowmin <- rel_imagerowmin * (max_imagerow - min_imagerow) + min_imagerow
    NEWimagerowmax <- rel_imagerowmax * (max_imagerow - min_imagerow) + min_imagerow
    NEWimagecolmin <- rel_imagecolmin * (max_imagecol - min_imagecol) + min_imagecol
    NEWimagecolmax <- rel_imagecolmax * (max_imagecol - min_imagecol) + min_imagecol
  }
  
  grupo <- rownames(objectST@images[[1]]@coordinates[which(imagerow >= NEWimagerowmin & imagerow <= NEWimagerowmax &
                                                             imagecol >= NEWimagecolmin & imagecol <= NEWimagecolmax),])
  
  # Subset and re-normalize each ST object
  objectST.zoom <- subset(objectST, cells = grupo, invert = FALSE) #%>%
    #SCTransform(assay = "Spatial", verbose = FALSE) %>%
    #RunPCA(verbose = FALSE)
  objname <- names(listaObjAnnot[i])
  listaObjZoomAnnot[[objname]] <- objectST.zoom
}

# Plot the zoomed SpatialDimPlots for each ST object
m1 <- SpatialDimPlot(listaObjZoomAnnot[[1]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m2 <- SpatialDimPlot(listaObjZoomAnnot[[2]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m3 <- SpatialDimPlot(listaObjZoomAnnot[[3]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
print(m1 + m2 + m3)

############################################################################### 
#                      Only plot with spots                                   #
# Not run this code, it is better the "All code automatized for graphs"       #
############################################################################### 
listaGraphsZoom <- list()
k <- length(cell.types)
for (i in 1:length(cell.types)) {
  for (j in 1:length(listaObjZoomAnnot)) {
    listaGraphsZoom[[names(listaObjZoomAnnot[j])]][[cell.types[[i]]]] <- SpatialFeaturePlot(listaObjZoomAnnot[[j]], 
                                                                                    features = cell.types[[i]], 
                                                                                    pt.size.factor = 1.2, 
                                                                                    crop = FALSE, 
                                                                                    alpha = c(0,1)) + 
      scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                           limits = c(0, 1), 
                           oob = scales::squish)
  }
}

referenceGzoom <- listaGraphsZoom[[1]]
notAlignGzoom <- listaGraphsZoom[[2]]
alignGzoom <- listaGraphsZoom[[3]]

( referenceGzoom$Neurons + referenceGzoom$Oligodendrocytes + 
    referenceGzoom$`Ependymal Cells` + referenceGzoom$Astrcytes +
    notAlignGzoom$Neurons + notAlignGzoom$Oligodendrocytes + 
    notAlignGzoom$`Ependymal Cells` + notAlignGzoom$Astrcytes + 
    alignGzoom$Neurons + alignGzoom$Oligodendrocytes + 
    alignGzoom$`Ependymal Cells` + alignGzoom$Astrcytes ) 


###############################################################################
#                   All code automatized for graphs                           #
############################################################################### 
listaObjZoomAnnot       

      listaValuesCellType <- list()
      listaGraphsZoom <- list()
      listaStatistics <- list()
      listaStatisticsGraphs <- list()
      for (i in 1:length(cell.types)) {
        for (j in 1:length(listaObjZoomAnnot)) {
          
          # Obtain the values for each cell type
          listaValuesCellType[[cell.types[[i]]]][[names(listaObjZoomAnnot[j])]] <- 
            listaObjZoomAnnot[[j]]@assays[[RCTDmode]]@data[cell.types[[i]],]
          
          # Obtain the graphs for each cell type
          listaGraphsZoom[[names(listaObjZoomAnnot[j])]][[cell.types[[i]]]] <- 
            SpatialFeaturePlot(listaObjZoomAnnot[[j]], 
                               features = cell.types[[i]],
                               pt.size.factor = 1.2, 
                               crop = FALSE,
                               alpha = c(0,1)) + 
            scale_fill_gradientn(colors = c("transparent", cell_colors[[i]]), 
                                 limits = c(0, 1), 
                                 oob = scales::squish)
        }
        
        # Obtain the p-values for each cell type for each image to the reference one
        eachCellType <- listaValuesCellType[[cell.types[[i]]]]
        ttestvalues_before <- t.test(eachCellType[[1]],eachCellType[[2]])
        p_value_before <- ttestvalues_before$p.value
        
        ttestvalues_after <- t.test(eachCellType[[1]],eachCellType[[3]])
        p_value_after <- ttestvalues_after$p.value
        
        ttestvalues_problem <- t.test(eachCellType[[2]],eachCellType[[3]])
        p_value_problem <- ttestvalues_problem$p.value
        
        # Obtain the mean and the SD for each cell type
        stats_df <- data.frame(
          Group = c("Reference", "Original", "Transformed"),
          Mean = c(mean(eachCellType[[1]]), 
                   mean(eachCellType[[2]]), 
                   mean(eachCellType[[3]])),
          SD = c(sd(eachCellType[[1]]), 
                 sd(eachCellType[[2]]), 
                 sd(eachCellType[[3]])),
          error = c(sd(eachCellType[[1]])/sqrt(length(eachCellType[[1]])), 
                    sd(eachCellType[[2]])/sqrt(length(eachCellType[[2]])), 
                    sd(eachCellType[[3]])/sqrt(length(eachCellType[[3]])))
        )
        stats_df$Group <- factor(stats_df$Group, levels = c("Reference", "Original", "Transformed"))
        
        listaStatistics[[cell.types[[i]]]][["stats"]] <- stats_df
        listaStatistics[[cell.types[[i]]]][["p_values"]] <- c(p_value_before, p_value_after, p_value_problem)
        
        if (p_value_before <= 0.05 && !is.nan(p_value_before)) {color_pvalue_before <- "red"} else {color_pvalue_before <- "black"}
        if (p_value_after <= 0.05 && !is.nan(p_value_after)) {color_pvalue_after <- "red"} else {color_pvalue_after <- "black"}
        if (p_value_problem <= 0.05 && !is.nan(p_value_problem)) {color_pvalue_problem <- "red"} else {color_pvalue_problem <- "black"}
        
        listaStatisticsGraphs[[cell.types[[i]]]] <- ggplot(stats_df, aes(x = Group, y = Mean)) +
          geom_bar(stat = "identity",  width = 0.4, fill = c("skyblue", "green", "purple")) +  # Barras de la media
          geom_errorbar(aes(ymin = Mean - error, ymax = Mean + error), width = 0.2) +  # Barras de error
          labs(y = "Mean values", x = "", title = cell.types[[i]]) +
          theme_minimal() +
          theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),,
                axis.text.x = element_text(size = 11),
                axis.text.y = element_text(size = 11),
                axis.title.y = element_text(size = 11)
          ) +
          annotate("text", x = 1.5, y = max(stats_df[,]$Mean) + max(stats_df[,]$error) + 1/8*max(stats_df[,]$Mean),
                   #label = paste0("|-",sprintf("%.2e", p_value_before),"-|"), size = 4, color = color_pvalue_before)  +
                   label = paste0("|--",sprintf("%.4f", p_value_before),"--|"), size = 4, color = color_pvalue_before) +
                   #label = paste0("|----",p_value_before,"----|"), size = 4, color = color_pvalue_before) +
          annotate("text", x = 2, y = max(stats_df[,]$Mean) + max(stats_df[,]$error) + 1/3*max(stats_df[,]$Mean),
                   #label = paste0("|-------",sprintf("%.2e", p_value_after),"-------|"), size = 4, color = color_pvalue_after) +
                   label = paste0("|--------",sprintf("%.4f", p_value_after),"--------|"), size = 4, color = color_pvalue_after) +
                   #label = paste0("|----------",p_value_after,"----------|"), size = 4, color = color_pvalue_after) +
          annotate("text", x = 2.5, y = max(stats_df[,]$Mean) + max(stats_df[,]$error) + 1/2*max(stats_df[,]$Mean) ,
                   #label = paste0("|-",sprintf("%.2e", p_value_problem),"-|"), size = 4, color = color_pvalue_problem) 
                   label = paste0("|--",sprintf("%.4f", p_value_problem),"--|"), size = 4, color = color_pvalue_problem) 
                   #label = paste0("|----",p_value_problem,"----|"), size = 4, color = color_pvalue_problem) 
      }
      
  
      referenceGzoom <- listaGraphsZoom[[1]]
      notAlignGzoom <- listaGraphsZoom[[2]]
      alignGzoom <- listaGraphsZoom[[3]]
      
      # central canal
      ( referenceGzoom$Neurons + referenceGzoom$Oligodendrocytes + 
          referenceGzoom$`Ependymal Cells` + referenceGzoom$Astrcytes +
          notAlignGzoom$Neurons + notAlignGzoom$Oligodendrocytes + 
          notAlignGzoom$`Ependymal Cells` + notAlignGzoom$Astrcytes + 
          alignGzoom$Neurons + alignGzoom$Oligodendrocytes + 
          alignGzoom$`Ependymal Cells` + alignGzoom$Astrcytes )  
      
      # grey matter butterfly
      ( referenceGzoom$Neurons + referenceGzoom$Oligodendrocytes + 
          referenceGzoom$Endothelial + referenceGzoom$Astrcytes +
          notAlignGzoom$Neurons + notAlignGzoom$Oligodendrocytes + 
          notAlignGzoom$Endothelial + notAlignGzoom$Astrcytes + 
          alignGzoom$Neurons + alignGzoom$Oligodendrocytes + 
          alignGzoom$Endothelial + alignGzoom$Astrcytes ) 
      
      ( listaStatisticsGraphs$Neurons + listaStatisticsGraphs$Oligodendrocytes +
        listaStatisticsGraphs$`Ependymal Cells` + listaStatisticsGraphs$Astrcytes )
      
      ( listaStatisticsGraphs$Neurons + listaStatisticsGraphs$Oligodendrocytes + 
        listaStatisticsGraphs$Astrcytes + listaStatisticsGraphs$Microglia +
        listaStatisticsGraphs$`Ependymal Cells` + listaStatisticsGraphs$Meninges + 
        listaStatisticsGraphs$Lymphocytes + listaStatisticsGraphs$Pericytes +
        listaStatisticsGraphs$OPC + listaStatisticsGraphs$Endothelial + 
        listaStatisticsGraphs$Schwann )
  
difOrig <- 0
difAlign <- 0
for (cell_type in names(listaStatistics)) {
  difOrig <- difOrig + abs(listaStatistics[[cell_type]]$stats$Mean[[1]] - listaStatistics[[cell_type]]$stats$Mean[[2]])
  difAlign <- difAlign + abs(listaStatistics[[cell_type]]$stats$Mean[[1]] - listaStatistics[[cell_type]]$stats$Mean[[3]])
}
print(paste0("The difference with the reference before the alignment is ", round(difOrig,4),
      " and the difference with the reference after the alignment is ", round(difAlign,4)))
      

############################################################################### 
###############################################################################
##                                                                           ##
##      Code for the Matrix of Groups comparison (full mode)                 ##
##                                                                           ##
###############################################################################
############################################################################### 

# Matrix of the image
library(data.table)

# Build emty matrixes
imageMatrixRef <- data.table(matrix(numeric(0), ncol = 11))
imageMatrixNOT <- data.table(matrix(numeric(0), ncol = 11))
imageMatrixYES <- data.table(matrix(numeric(0), ncol = 11))
setnames(imageMatrixRef, cell.types)  
setnames(imageMatrixNOT, cell.types)  
setnames(imageMatrixYES, cell.types)  

# Coordinates from the reference image
objectST1 <- listaObjAnnot[[1]]
imagerow1 <- objectST1@images[[1]]@coordinates$imagerow
imagecol1 <- objectST1@images[[1]]@coordinates$imagecol

max_imagerow1 <- max(imagerow1)
min_imagerow1 <- min(imagerow1)
max_imagecol1 <- max(imagecol1)
min_imagecol1 <- min(imagecol1)

# Useful variables
group <- 0
listRef <- list()
listNOT <- list()
listYES <- list()

colsRef <- listaObjAnnot[[1]]@images[[1]]@coordinates$col
rowsRef <- listaObjAnnot[[1]]@images[[1]]@coordinates$row
min_colsRef <- min(colsRef)
max_colsRef <- max(colsRef)
min_rowsRef <- min(rowsRef)
max_rowsRef <- max(rowsRef)

colmin <- min_colsRef
colmax <- colmin+5
while (colmin <= max_colsRef) {
  print(paste0("columna ", colmin))
  rowmin <- min_rowsRef
  rowmax <- rowmin+3
  
  while (rowmin <= max_rowsRef) {
    group <- group + 1
    print(paste0("grupo ", group, " y fila ", rowmin))
    
    ################################ Selection of row and col values ##########
    imagerowmin <- min(imagerow1[objectST1@images[[1]]@coordinates$row >= rowmin & 
                                 objectST1@images[[1]]@coordinates$row < rowmax])
    imagerowmax <- max(imagerow1[objectST1@images[[1]]@coordinates$row >= rowmin & 
                                 objectST1@images[[1]]@coordinates$row < rowmax])
    imagecolmin <- min(imagecol1[objectST1@images[[1]]@coordinates$col >= colmin & 
                                 objectST1@images[[1]]@coordinates$col < colmax])
    imagecolmax <- max(imagecol1[objectST1@images[[1]]@coordinates$col >= colmin & 
                                 objectST1@images[[1]]@coordinates$col < colmax])
    ########################################################################### 
    
    rel_imagerowmin <- (imagerowmin - min_imagerow1) / (max_imagerow1 - min_imagerow1)
    rel_imagerowmax <- (imagerowmax - min_imagerow1) / (max_imagerow1 - min_imagerow1) 
    rel_imagecolmin <- (imagecolmin - min_imagecol1) / (max_imagecol1 - min_imagecol1) 
    rel_imagecolmax <- (imagecolmax - min_imagecol1) / (max_imagecol1 - min_imagecol1) 
  
    for (i in seq_along(listaObjAnnot)) {
      objectST <- listaObjAnnot[[i]]
      imagerow <- objectST@images[[1]]@coordinates$imagerow
      imagecol <- objectST@images[[1]]@coordinates$imagecol
      
      # If the image has been aligned, the coordinates to make a subset must be 
      # selected from the image before the alignment
      if (i == 3) {
        objectST_orig <- listaObjAnnot[[2]]
        imagerow_orig <- objectST_orig@images[[1]]@coordinates$imagerow
        imagecol_orig <- objectST_orig@images[[1]]@coordinates$imagecol
        
        min_imagerow_orig <- min(imagerow_orig)
        max_imagerow_orig <- max(imagerow_orig)
        min_imagecol_orig <- min(imagecol_orig)
        max_imagecol_orig <- max(imagecol_orig)
        
        NEWimagerowmin <- rel_imagerowmin * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
        NEWimagerowmax <- rel_imagerowmax * (max_imagerow_orig - min_imagerow_orig) + min_imagerow_orig
        NEWimagecolmin <- rel_imagecolmin * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
        NEWimagecolmax <- rel_imagecolmax * (max_imagecol_orig - min_imagecol_orig) + min_imagecol_orig
        
      } else {
        min_imagerow <- min(imagerow)
        max_imagerow <- max(imagerow)
        min_imagecol <- min(imagecol)
        max_imagecol <- max(imagecol)
        
        NEWimagerowmin <- rel_imagerowmin * (max_imagerow - min_imagerow) + min_imagerow
        NEWimagerowmax <- rel_imagerowmax * (max_imagerow - min_imagerow) + min_imagerow
        NEWimagecolmin <- rel_imagecolmin * (max_imagecol - min_imagecol) + min_imagecol
        NEWimagecolmax <- rel_imagecolmax * (max_imagecol - min_imagecol) + min_imagecol
      } # end if-else
      
      coordValues <- which(imagerow >= NEWimagerowmin & imagerow <= NEWimagerowmax &
                             imagecol >= NEWimagecolmin & imagecol <= NEWimagecolmax)
      
      if (length(coordValues) != 0) {
        grupo <- rownames(objectST@images[[1]]@coordinates[coordValues,])
        objectST.group <- subset(objectST, cells = grupo, invert = FALSE) 
        
        # We take the cell types weigths from the complete deconvolution data (full mode)
        list_rows <- as.list(sapply(cell.types, function(cell_type) {
          mean(objectST.group@assays[["deconvolution.RCTD.complete"]]@data[cell_type, ])
        }))
      } else {
        grupo <- NA
        list_rows <- as.list(sapply(cell.types, function(cell_type) { NA }))
      } # end if-else
      
      if (i == 1) {
        listRef[[group]] <- list_rows
      } else if (i == 2) {
        listNOT[[group]] <- list_rows
      } else if (i == 3) {
        listYES[[group]] <- list_rows
      }
      
    } # end for listaObjAnnot
    rowmin <- rowmax
    rowmax <- rowmin+3
  } # end for while row
  colmin <- colmax
  colmax <- colmin+5
} # end for while col


imageMatrixRef <- rbindlist(lapply(listRef, as.list), fill = TRUE)
imageMatrixNOT <- rbindlist(lapply(listNOT, as.list), fill = TRUE)
imageMatrixYES <- rbindlist(lapply(listYES, as.list), fill = TRUE)


m1 <- SpatialDimPlot(listaObjAnnot[[1]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m2 <- SpatialDimPlot(listaObjAnnot[[2]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
m3 <- SpatialDimPlot(listaObjAnnot[[3]], alpha = 0.5, crop = FALSE) & theme(legend.position = "none")
print(m1 + m2 + m3)

table(complete.cases(imageMatrixRef))
table(complete.cases(imageMatrixNOT))
table(complete.cases(imageMatrixYES))

sum(imageMatrixRef[complete.cases(imageMatrixRef) & complete.cases(imageMatrixYES)])
sum(imageMatrixNOT[complete.cases(imageMatrixNOT) & complete.cases(imageMatrixYES)])
shared_groups <- complete.cases(imageMatrixRef) & complete.cases(imageMatrixNOT) & complete.cases(imageMatrixYES)
table(shared_groups)

coeffEV_Ref_NOT <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                       imageMatrixNOT[shared_groups])

coeffEV_Ref_YES <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                       imageMatrixYES[shared_groups])

coeffEV_NOT_YES <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                       imageMatrixYES[shared_groups])


coeffEV_Ref <- FactoMineR::coeffRV(imageMatrixRef[shared_groups],
                                   imageMatrixRef[shared_groups])

coeffEV_NOT <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                   imageMatrixNOT[shared_groups])

coeffEV_YES <- FactoMineR::coeffRV(imageMatrixNOT[shared_groups],
                                   imageMatrixNOT[shared_groups])

# Visualization
rv_matrix <- matrix(c(coeffEV_Ref$rv,     coeffEV_Ref_NOT$rv, coeffEV_Ref_YES$rv,
                      coeffEV_Ref_NOT$rv, coeffEV_NOT$rv,     coeffEV_NOT_YES$rv,
                      coeffEV_Ref_YES$rv, coeffEV_NOT_YES$rv, coeffEV_YES$rv),    
                    nrow = 3, byrow = TRUE)
rownames(rv_matrix) <- c("Reference", "Original", "Transformed")
colnames(rv_matrix) <- c("Reference", "Original", "Transformed")

pval_matrix <- matrix(c(coeffEV_Ref$p.value,   coeffEV_Ref_NOT$p.value, coeffEV_Ref_YES$p.value,
                      coeffEV_Ref_NOT$p.value, coeffEV_NOT$p.value,     coeffEV_NOT_YES$p.value,
                      coeffEV_Ref_YES$p.value, coeffEV_NOT_YES$p.value, coeffEV_YES$p.value),    
                    nrow = 3, byrow = TRUE)
rownames(pval_matrix) <- c("Reference", "Original", "Transformed")
colnames(pval_matrix) <- c("Reference", "Original", "Transformed")

rvstd_matrix <- matrix(c(coeffEV_Ref$rvstd,     coeffEV_Ref_NOT$rvstd, coeffEV_Ref_YES$rvstd,
                      coeffEV_Ref_NOT$rvstd, coeffEV_NOT$rvstd,     coeffEV_NOT_YES$rvstd,
                      coeffEV_Ref_YES$rvstd, coeffEV_NOT_YES$rvstd, coeffEV_YES$rvstd),    
                    nrow = 3, byrow = TRUE)
rownames(rvstd_matrix) <- c("Reference", "Original", "Transformed")
colnames(rvstd_matrix) <- c("Reference", "Original", "Transformed")

library(reshape2)
rv_long <- melt(rv_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "RV")
pval_long <- melt(pval_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "p.value")
rvstd_long <- melt(rvstd_matrix, varnames = c("Matrix_A", "Matrix_B"), value.name = "RVstd")

rv_long <- merge(rv_long, pval_long, by = c("Matrix_A", "Matrix_B"))

library(ggplot2)
ggplot(rv_long, aes(x = Matrix_A, y = Matrix_B, fill = RV)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(sprintf("%.4f", RV), "\n(", sprintf("%.2e", p.value), ")")),
            color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5) +
  scale_y_discrete(limits = rev(levels(rv_long$Matrix_B))) +
  theme_minimal() +
  labs(x = "", y = "", title = "Heatmap of RV Coefficients and p-values") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

ggplot(rvstd_long, aes(x = Matrix_A, y = Matrix_B, fill = RVstd)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.4f", RVstd)),
            color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = coeffEV_Ref$rvstd/2) +
  scale_y_discrete(limits = rev(levels(rv_long$Matrix_B))) +
  theme_minimal() +
  labs(x = "", y = "", title = "Heatmap of standarized RV Coefficients") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
