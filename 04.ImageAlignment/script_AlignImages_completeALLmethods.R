###############################################################################
#             Main alignment code                                             #
###############################################################################
# Load required libraries
# The alignment functions in `semla` only handles rigid transformations, meaning
# that non-linear distortions are not possible to mitigate.
library(reticulate)
library(semla)
library(tibble)

# Load a custom function for alignment evaluation
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_evalAlign.R")
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_selectCoord.R")

options <- c("19", "MIX")
choice <- menu(options, title = "Select which patient do you want to align: \n(patient 19 refers to three different slices from patient 19 and \npatient MIX refers to three different samples from different patients, \nin both choices the alignment would be performed to the same slice from patient 19 \nwhich is considered as the reference)")
patient <- options[choice]

options <- c("GTEM (Geometric Transformation Estimation Model)", 
             "Procrustes (Procrustes Transformation)",
             "imageJ (ImageJ - Register Virtual Stack Slices)")
choice <- menu(options, title = "Select which mode do you want to perform the alignment.")
modesAbrs <- c("alignTR","procrustes","imageJ")
modeAbr <- modesAbrs[choice]
mode <- options[choice]

options <- c("Yes", "No")
choice <- menu(options, title = "Do you want to introduce a scale factor?")
scale <- options[choice]

# Define the files' folders 
carpetaData <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/02.Seurat/Yadav2023_adult_human_spinal_cord/ImagenesYadavMergeSlices/"
functionsDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment"
saveDir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingmentComplete"

# Load the merged Seurat object 
paciente.merge <- readRDS(file = paste0(carpetaData, "/Paciente", patient, "_merge.rds"))

if (mode == "GTEM (Geometric Transformation Estimation Model)") {
  if (scale == "No") {
    source(paste0(functionsDir, "/function_selectCoord_TR.R"))
  } else if (scale == "Yes") {
    source(paste0(functionsDir, "/function_selectCoord_TR_S.R"))
  }
} else if (mode == "Procrustes (Procrustes Transformation)") {
  source(paste0(functionsDir, "/function_selectCoord_procrustes.R"))
}

# Update Seurat object for compatibility with the semla package
pacientes.semla <- UpdateSeuratForSemla(paciente.merge)

# Load images into the Seurat object with specified image height
pacientes.semla <- LoadImages(pacientes.semla, image_height = 600)

# Plot loaded images for visual inspection
ImagePlot(pacientes.semla)

if (mode == "GTEM (Geometric Transformation Estimation Model)") {
  modeAbr <- "alignTR"
  # Select and round the coordinates from the reference image (first sample)
  coordenadas1 <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  for (j in seq_along(coordenadas1)) {
    coordenadas1[[j]] <- round(coordenadas1[[j]])
  }
  coordenadas1  # Display rounded coordinates for the reference image
  
  # Initialize lists for storing image data, coordinates, transformations, and other parameters
  aligned <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  original <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  listaCoordenadas <- list(coordenadas1)
  listaSoluciones <- list(NA)
  listaOpciones <- list(NA)
  listaOpcionesCalc <- list(NA)
  listaValores <- list(NA)
  listaTransforms <- list(NA)
  listaCoordenadasNEW <- list(coordenadas1)
  
  # Extract column and row coordinates for the first sample
  CoordinatesSemlaCol <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 2])
  CoordinatesSemlaRow <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 3])
  
  # Iterate through each image starting from the second sample
  for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
    
    # Select and round coordinates for the current image
    coordenadas2 <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    for (j in seq_along(coordenadas2)) {
      coordenadas2[[j]] <- round(coordenadas2[[j]])
    }
    print(coordenadas2)  # Display rounded coordinates for the current image
    listaCoordenadas[[i]] <- coordenadas2
    
    # Define dimensions for the current image
    xmax2 <- nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    ymax2 <- ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    
    # Solve transformations without mirroring (original orientation)
    solucionOrig <- solveCoord(coordenadas1, coordenadas2, xmax2, ymax2)
    solucionOrig[["mirrorx"]] <- 0
    solucionOrig[["mirrory"]] <- 0
    
    # Solve transformations with mirroring on the x-axis
    coordenadas2X <- coordenadas2
    for (j in seq_along(coordenadas2X$x)) {
      coordenadas2X$x[[j]] <- xmax2 - coordenadas2X$x[[j]]
    }
    solucionMirrorX <- solveCoord(coordenadas1, coordenadas2X, xmax2, ymax2)
    solucionMirrorX[["mirrorx"]] <- 10
    solucionMirrorX[["mirrory"]] <- 0
    
    # Solve transformations with mirroring on the y-axis
    coordenadas2Y <- coordenadas2
    for (j in seq_along(coordenadas2Y$y)) {
      coordenadas2Y$y[[j]] <- ymax2 - coordenadas2Y$y[[j]]
    }
    solucionMirrorY <- solveCoord(coordenadas1, coordenadas2Y, xmax2, ymax2)
    solucionMirrorY[["mirrorx"]] <- 0
    solucionMirrorY[["mirrory"]] <- 10
    
    # Solve transformations with mirroring on both x and y axes
    coordenadas2XY <- coordenadas2
    for (j in seq_along(coordenadas2XY$x)) {
      coordenadas2XY$x[[j]] <- xmax2 - coordenadas2XY$x[[j]]
    }
    for (j in seq_along(coordenadas2XY$y)) {
      coordenadas2XY$y[[j]] <- ymax2 - coordenadas2XY$y[[j]]
    }
    solucionMirrorXY <- solveCoord(coordenadas1, coordenadas2XY, xmax2, ymax2)
    solucionMirrorXY[["mirrorx"]] <- 10
    solucionMirrorXY[["mirrory"]] <- 1
    
    # Store each transformation option in a list for the current image
    option <- list()
    option[["solucionOrig"]] <- solucionOrig
    option[["solucionMirrorX"]] <- solucionMirrorX
    option[["solucionMirrorY"]] <- solucionMirrorY
    option[["solucionMirrorXY"]] <- solucionMirrorXY
    listaOpciones[[i]] <- option
  }
  
  # Iterate over each image in the dataset (starting from the second image)
  for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
    
    # Calculate the minimum dimensions between the reference image and the current image
    xmax <- min(nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    ymax <- min(ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    
    # Initialize an empty list to store calculated parameter values for each transformation option
    todosvalores <- list(NA)
    for (j in seq_along(listaOpciones[[i]])) {
      # Select the current transformation option for the image
      solucion <- listaOpciones[[i]][[j]]
      
      # Calculate parameters (e.g., angles, translations) based on the solution and image dimensions
      valores <- calcParameters(solucion, xmax, ymax)
      todosvalores[[j]] <- valores
    }
    
    # Store calculated transformation values for each option in a separate list
    listaOpcionesCalc[[i]] <- todosvalores
    
    # Extract and sum squares of selected parameter values (angle, dx, dy) across all options for comparison
    val_sum_cuad <- list(NA)
    for (j in seq_along(todosvalores)) {
      val_sum_cuad[[j]] <- todosvalores[[j]][1:3]
    }
    
    # Calculate the sum of squares for each transformation option to find the optimal alignment
    suma_de_cuadrados <- apply(do.call(rbind, listaOpcionesCalc[[i]]), 1, function(x) sum(x^2))
    indice_fila_minima <- which.min(suma_de_cuadrados)  # Index of minimum sum of squares
    
    # Determine mirroring settings based on the index of the optimal transformation
    if (indice_fila_minima == 2 | indice_fila_minima == 4) { mirrorx <- TRUE } else { mirrorx <- FALSE }
    if (indice_fila_minima == 3 | indice_fila_minima == 4) { mirrory <- TRUE } else { mirrory <- FALSE }
    
    # Update the solution with the optimal transformation option and mirror settings
    solucion <- solucionOrig
    solucion <- listaOpciones[[i]][[indice_fila_minima]]
    solucion[["mirrorx"]] <- mirrorx
    solucion[["mirrory"]] <- mirrory
    listaSoluciones[[i]] <- solucion
    
    # Store the selected transformation parameters for later use
    listaValores[[i]] <- listaOpcionesCalc[[i]][[indice_fila_minima]]
    
    # Update coordinates for the transformed image
    coord <- listaCoordenadas[[i]]
    xmax <- nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    ymax <- ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    coordNEW <- calcNewCoord(coord, solucion, xmax, ymax)
    listaCoordenadasNEW[[i]] <- coordNEW
    
    # Display the selected transformation parameters
    valores <- listaValores[[i]]
    print(valores)
    
    # Initialize a list to store all transformations applied to the image
    alltransforms <- list()
    original[[i]] <- pacientes.semla@tools$Staffli@rasterlists$raw[[i]]
    
    # Apply mirror transformation if the mirror values are non-zero
    if ( valores[['mirrorx']] != 0 || valores[['mirrory']] != 0 ) {
      transforms_mirror <- generate_rigid_transform(sampleID = i, 
                                                    mirror_x = as.logical(valores[['mirrorx']]),
                                                    mirror_y = as.logical(valores[['mirrory']]))
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_mirror)
      alltransforms[[1]] <- transforms_mirror
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply rotation if the angle is non-zero
    angulo <- valores[['angulo']]
    if ( angulo != 0 && angulo != 360 ){
      transforms_angle <- generate_rigid_transform(sampleID = i, 
                                                   angle = angulo)
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_angle)
      alltransforms[[2]] <- transforms_angle
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply translation if dx or dy is non-zero
    if ( valores[['trx']] != 0 || valores[['try']] != 0 ) {
      transforms_trans <- generate_rigid_transform(sampleID = i, 
                                                   tr_x = valores[['trx']], #round(valores[[2]], digits = 2), 
                                                   tr_y = valores[['try']]) #round(valores[[3]], digits = 2),
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_trans)
      alltransforms[[3]] <- transforms_trans
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Store aligned images and transformation lists
    aligned[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
    pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- original[[i]]
    listaTransforms[[i]] <- alltransforms
    
    # Save transformed coordinates for each sample
    CoordinatesSemlaCol[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 4]
    CoordinatesSemlaRow[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 5]
  }
  
} else if (mode == "Procrustes (Procrustes Transformation)") {
  modeAbr <- "procrustes"
  
  # Select and round the coordinates from the reference image (first sample)
  coordenadas1 <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  for (j in seq_along(coordenadas1)) {
    coordenadas1[[j]] <- round(coordenadas1[[j]])
  }
  coordenadas1  # Display rounded coordinates for the reference image
  
  # Initialize lists for storing image data, coordinates, transformations, and other parameters
  aligned <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  original <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  listaCoordenadas <- list(coordenadas1)
  listaSoluciones <- list(NA)
  listaOpciones <- list(NA)
  listaOpcionesCalc <- list(NA)
  listaValores <- list(NA)
  listaTransforms <- list(NA)
  listaCoordenadasNEW <- list(coordenadas1)
  
  # Extract column and row coordinates for the first sample
  CoordinatesSemlaCol <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 2])
  CoordinatesSemlaRow <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 3])
  
  # Loop through each sample image starting from the second one
  for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
    
    # Select and round coordinates for the current image
    coordenadas2 <- selectCoord(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    for (j in seq_along(coordenadas2)) { 
      coordenadas2[[j]] <- round(coordenadas2[[j]])
    }
    print(coordenadas2)  # Display rounded coordinates for the current image
    listaCoordenadas[[i]] <- coordenadas2
    
    ## Define dimensions for the current image
    xmax2 <- nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    ymax2 <- ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]])
    
    # Initialize lists for storing sum of squares and calculated coordinates
    val_sum_cuad <- list()
    coordCalc <- list()
    
    # Solve for the original orientation (without mirroring)
    matProb <- matrix(data = unlist(coordenadas2), ncol = 2)
    proc <- Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = TRUE, sumsq = TRUE)
    coseno <- proc$R[1,1]
    seno <- proc$R[2,1]
    dx <- proc$t[1]
    dy <- proc$t[2]
    mirrorx <- FALSE
    mirrory <- FALSE
    e <- proc$d
    solucionOrig <- c(coseno, seno, dx, dy, mirrorx, mirrory, e)
    names(solucionOrig) <- c('coseno', 'seno', 'dx', 'dy', 'mirrorx', 'mirrory', 'e')
    val_sum_cuad[["solucionOrig"]] <- proc$ss
    coordCalc[["solucionOrig"]] <- proc$X.new
    
    # Solve with mirror on x-axis
    coordenadas2X <- coordenadas2
    for (j in seq_along(coordenadas2X$x)) {
      coordenadas2X$x[[j]] <- xmax2 - coordenadas2X$x[[j]]
    }
    matProb <- matrix(data = unlist(coordenadas2X), ncol = 2)
    proc <- Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = TRUE, sumsq = TRUE)
    solucionMirrorX <- c(proc$R[1,1], proc$R[2,1], proc$t[1], proc$t[2], TRUE, FALSE, proc$d)
    names(solucionMirrorX) <- c('coseno', 'seno', 'dx', 'dy', 'mirrorx', 'mirrory', 'e')
    val_sum_cuad[["solucionMirrorX"]] <- proc$ss
    coordCalc[["solucionMirrorX"]] <- proc$X.new
    
    # Solve with mirror on y-axis
    coordenadas2Y <- coordenadas2
    for (j in seq_along(coordenadas2Y$y)) {
      coordenadas2Y$y[[j]] <- ymax2 - coordenadas2Y$y[[j]]
    }
    matProb <- matrix(data = unlist(coordenadas2Y), ncol = 2)
    proc <- Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = TRUE, sumsq = TRUE)
    solucionMirrorY <- c(proc$R[1,1], proc$R[2,1], proc$t[1], proc$t[2], FALSE, TRUE, proc$d)
    names(solucionMirrorY) <- c('coseno', 'seno', 'dx', 'dy', 'mirrorx', 'mirrory', 'e')
    val_sum_cuad[["solucionMirrorY"]] <- proc$ss
    coordCalc[["solucionMirrorY"]] <- proc$X.new
    
    # Solve with mirror on both x and y axes
    coordenadas2XY <- coordenadas2
    for (j in seq_along(coordenadas2XY$x)) {
      coordenadas2XY$x[[j]] <- xmax2 - coordenadas2XY$x[[j]]
    }
    for (j in seq_along(coordenadas2XY$y)) {
      coordenadas2XY$y[[j]] <- ymax2 - coordenadas2XY$y[[j]]
    }
    matProb <- matrix(data = unlist(coordenadas2XY), ncol = 2)
    proc <- Procrustes(X = matProb, Xstar = matrixCoord1, translate = TRUE, dilate = TRUE, sumsq = TRUE)
    solucionMirrorXY <- c(proc$R[1,1], proc$R[2,1], proc$t[1], proc$t[2], TRUE, TRUE, proc$d)
    names(solucionMirrorXY) <- c('coseno', 'seno', 'dx', 'dy', 'mirrorx', 'mirrory', 'e')
    val_sum_cuad[["solucionMirrorXY"]] <- proc$ss
    coordCalc[["solucionMirrorXY"]] <- proc$X.new
    
    # Store each transformation option in a list for the current image
    option <- list()
    option[["solucionOrig"]] <- solucionOrig
    option[["solucionMirrorX"]] <- solucionMirrorX
    option[["solucionMirrorY"]] <- solucionMirrorY
    option[["solucionMirrorXY"]] <- solucionMirrorXY
    listaOpciones[[i]] <- option
    
    # Select the transformation option with the lowest sum of squares
    indice_fila_minima <- which.min(val_sum_cuad)
    
    # Calculate final angle and translation parameters based on selected transformation
    xmax <- min(nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    ymax <- min(ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    
    valores <- calcParametersScale(option[[indice_fila_minima]], xmax, ymax)
    
    # Store chosen transformation and calculated parameters for the current image
    listaSoluciones[[i]] <- option[[indice_fila_minima]]
    listaValores[[i]] <- valores
    
    # Update coordinates based on selected transformation
    coordNEWmatrix <- coordCalc[[indice_fila_minima]]
    coordNEW <- list()
    coordNEW$x <- round(coordNEWmatrix[,1])
    coordNEW$y <- round(coordNEWmatrix[,2])
    listaCoordenadasNEW[[i]] <- coordNEW 
  }
  
  # Iterate over each image in the dataset (starting from the second image)
  for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
    
    # Display the selected transformation parameters
    valores <- listaValores[[i]]
    print(valores)
    
    # Initialize a list to store all transformations applied to the image
    alltransforms <- list()
    original[[i]] <- pacientes.semla@tools$Staffli@rasterlists$raw[[i]]
    
    # Apply mirror transformation if the mirror values are non-zero
    if ( valores[['mirrorx']] != 0 || valores[['mirrory']] != 0 ) {
      transforms_mirror <- generate_rigid_transform(sampleID = i, 
                                                    mirror_x = as.logical(valores[['mirrorx']]),
                                                    mirror_y = as.logical(valores[['mirrory']]))
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_mirror)
      alltransforms[[1]] <- transforms_mirror
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply rotation if the angle is non-zero
    angulo <- valores[['angulo']]
    if ( angulo != 0 && angulo != 360 ){
      transforms_angle <- generate_rigid_transform(sampleID = i, 
                                                   angle = angulo)
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_angle)
      alltransforms[[2]] <- transforms_angle
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply translation if dx or dy is non-zero
    if ( valores[['trx']] != 0 || valores[['try']] != 0 ) {
      transforms_trans <- generate_rigid_transform(sampleID = i, 
                                                   tr_x = valores[['trx']], #round(valores[[2]], digits = 2), 
                                                   tr_y = valores[['try']]) #round(valores[[3]], digits = 2),
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_trans)
      alltransforms[[3]] <- transforms_trans
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply scaling if the scale factor is non-zero
    if ( valores[['e']] != 1 ) {
      transforms_scale <- generate_rigid_transform(sampleID = i, 
                                                   scalefactor = valores[['e']])
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_scale)
      alltransforms[[4]] <- transforms_scale
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Store aligned images and transformation lists
    aligned[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
    pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- original[[i]]
    listaTransforms[[i]] <- alltransforms
    
    # Save transformed coordinates for each sample
    CoordinatesSemlaCol[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 4]
    CoordinatesSemlaRow[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 5]
  }
  
} else if (mode == "imageJ (ImageJ - Register Virtual Stack Slices)") {
  modeAbr <- "imageJ"
  
  # Set up lists for alignment process and coordinate retrieval
  aligned <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  original <- list(pacientes.semla@tools$Staffli@rasterlists$raw[[1]])
  listaSoluciones <- list(NA)
  listaTransforms <- list(NA)
  
  # Extract column and row coordinates for sample ID 1 from Seurat object
  CoordinatesSemlaCol <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 2])
  CoordinatesSemlaRow <- list(GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)["sampleID"] == 1), 3])
  
  
  "
  This is the point where you would save the images from the RDS file into a source_dir folder, 
  create a target_dir folder, and then use ImageJ with the Register Virtual Stack Slices plugin for alignment.
  "
  
  # Define source and target directories for ImageJ alignment
  source_dir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ/ImagenesAlinearFiji/pacienteMIX"
  target_dir <- "/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/05.ImageJ/ImagenesAlinearFiji/pacienteMIXaligned"
  #ref_name <- "original"
  ref_name <- "tissue_lowres_image_1"
  
  # Instructions for folder setup before alignment
  cat(paste("You should have a folder called: ", source_dir,
            "\nAnd a folder called: ", target_dir,
            "\nIn the first one, you should have the images from de rds document."))
  
  # Provide instructions for the manual alignment process in ImageJ
  cat(paste("It is time to open ImageJ and do the alignment of the images.",
            "\nFirst, go to File>Open and select the different images to align.",
            "The images should have only one color channel.",
            "To do this you should go to Image>Color>Split Channels and then Image>Color>Merge Channels unselecting create composite.",
            "\nTo do the align go to Plugins>Registration>Register Virtual Stack Slices.", 
            "Remember you should save the transforms.",
            "\nThe transformation parameters or Transform files should be stored in the same folder as the result images.",
            "\nAnd finally you have to select the reference image."))
  
  # Wait for user confirmation before proceeding
  answer <- "no"
  while (answer != 'yes') {
    cat(paste0("When you have all this completed, write: 'yes': "))
    answer <- readline()
  }
  if (answer == 'yes') {
    cat(paste0("The alignment is done and the transformation parameters are saved."))
  }
  
  # Extract translation and rotation parameters from XML files in the target directory
  filesTarget <- list.files(target_dir)
  files_xml <- filesTarget[grep("\\.xml$", filesTarget)]
  
  Nimage <- 1
  for (xml in files_xml) {
    # Skip the reference image file and process each alignment XML file
    if (!startsWith(xml, ref_name)) {
      Nimage <- Nimage + 1
      lines <- readLines(paste0(target_dir, "/", xml))
      lines <- lines[startsWith(lines, "\t<iict_transform")]
      lines <- sapply(lines, function(line) substr(line, start = nchar("\\t<iict_transform"), stop = nchar(line) - 3))
      lines <- sapply(lines, function(line) strsplit(line, "\""))
      
      # Extract rotation and translation parameters based on transform class type
      for (line in lines) {
        classN <- which(startsWith(line, " class="))
        if (endsWith(line[[classN + 1]], "transform.RigidModel2D")) { 
          dataN <- which(startsWith(line, " data="))
          transformsR <- line[[dataN + 1]]
          transformsR <- strsplit(transformsR, " ")
        }
        if (endsWith(line[[classN + 1]], "transform.TranslationModel2D")) { 
          dataN <- which(startsWith(line, " data="))
          transformsT <- line[[dataN + 1]]
          transformsT <- strsplit(transformsT, " ")
        }
      }
      # Save transformation parameters: rotation angle, dx, and dy
      transformsParams <- list()
      transformsParams["angle"] <- as.numeric(transformsR[[1]][1])  # radians
      transformsParams["dx"] <- as.numeric(transformsR[[1]][2]) + as.numeric(transformsT[[1]][1])  # without normalization (-1,1)
      transformsParams["dy"] <- as.numeric(transformsR[[1]][3]) + as.numeric(transformsT[[1]][2])  # without normalization (-1,1)
      listaSoluciones[[Nimage]] <- transformsParams
    }
  }
  
  # Apply transformations for each image in the Seurat object
  for (i in 2:length(pacientes.semla@tools$Staffli@rasterlists$raw)) {
    # Set image dimensions based on minimum height and width across images
    xmax <- min(nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                nrow(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    ymax <- min(ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[1]]),
                ncol(pacientes.semla@tools$Staffli@rasterlists$raw[[i]]))
    
    # Retrieve and convert transformation parameters
    solucion <- listaSoluciones[[i]]
    valores <- solucion
    valores$angle <- solucion$angle * (180/pi)  # Convert angle to degrees
    
    # Adjust dx and dy based on image center for transformation normalization
    valores$dx <-  (solucion$dx - (xmax/2 - xmax/2 * cos(solucion$angle) + ymax/2 * sin(solucion$angle))) / xmax
    valores$dy <-  (solucion$dy - (ymax/2 - ymax/2 * cos(solucion$angle) - xmax/2 * sin(solucion$angle))) / ymax
    
    alltransforms <- list()
    original[[i]] <- pacientes.semla@tools$Staffli@rasterlists$raw[[i]]
    
    
    # Apply rotation if the angle is non-zero
    if (valores$angle != 0 && valores$angle != 360) {
      transforms_angle <- generate_rigid_transform(sampleID = i, 
                                                   angle = valores$angle)
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_angle)
      alltransforms[[2]] <- transforms_angle
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Apply translation if dx or dy is non-zero
    if (valores$dx != 0 || valores$dy != 0) {
      transforms_trans <- generate_rigid_transform(sampleID = i, 
                                                   tr_x = valores$dx, 
                                                   tr_y = valores$dy)
      pacientes.semla <- RigidTransformImages(pacientes.semla, transforms = transforms_trans)
      alltransforms[[3]] <- transforms_trans
      pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
      
      # Update meta_data with new transformations
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 2] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 4]
      pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 3] <- pacientes.semla@tools[["Staffli"]]@meta_data[which(pacientes.semla@tools[["Staffli"]]@meta_data[6] == i), 5]
    }
    
    # Store aligned images and transformation lists
    aligned[[i]] <- pacientes.semla@tools$Staffli@rasterlists$transformed[[i]]
    pacientes.semla@tools$Staffli@rasterlists$raw[[i]] <- original[[i]]
    listaTransforms[[i]] <- alltransforms
    
    # Save transformed coordinates for each sample
    CoordinatesSemlaCol[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 4]
    CoordinatesSemlaRow[[i]] <- GetCoordinates(pacientes.semla)[which(GetCoordinates(pacientes.semla)[6] == i), 5]
  }
}


# Update transformed image list in Seurat object
pacientes.semla@tools$Staffli@rasterlists[["transformed"]] <- aligned

"Since I do the transformations sequentially, every time you do one 
new, only the modification of the last one you make is saved. That's why I'm 
saving at the end of the transformations of each image the new values of 
the positions to, once done on all the images, save the 
definitive versions of each of the images in the Seurat object in Staffli by semla.
"

# Aggregate transformed column and row coordinates
CoordCol <- integer()
CoordRow <- integer()
for (i in seq_along(CoordinatesSemlaCol)) {
  CoordCol <- c(CoordCol, CoordinatesSemlaCol[[i]][[1]])
  CoordRow <- c(CoordRow, CoordinatesSemlaRow[[i]][[1]])
}
pacientes.semla@tools$Staffli@meta_data$pxl_col_in_fullres_transformed <- CoordCol
pacientes.semla@tools$Staffli@meta_data$pxl_row_in_fullres_transformed <- CoordRow

# Plot the original and transformed images
png(paste0(saveDir, "/results/patient", patient, "_merge.png"))
ImagePlot(pacientes.semla)
title("Original H&E image", line = -12.8)
dev.off()
png(paste0(saveDir, "/results/patient", patient, "_merge_", modeAbr,".png"))
ImagePlot(pacientes.semla, image_use = "transformed")
title("Transformed H&E image", line = -12.8)
dev.off()

###############################################################################
#           Evaluation of the alignment                                       #
###############################################################################
listaRawImages <- list()
listaTransImages <- list()
for (i in 1:length(listaCoordenadasNEW)) {
  listaRawImages[[i]] <- pacientes.semla@tools[["Staffli"]]@rasterlists[["raw"]][[i]]
  listaTransImages[[i]] <- pacientes.semla@tools[["Staffli"]]@rasterlists[["transformed"]][[i]]
}

Evaluation <- evaluationComplete(listaCoordenadasNEW, listaCoordenadas,
                                 listaRawImages, listaTransImages)
  
# Update the transformed images in the original Seurat object
pacientes.seurat <- UpdateSeuratFromSemla(pacientes.semla, image_use = "transformed")

############################################################################### 
#                     Visualization 1                                         #
###############################################################################

# Visualize transformed images with spatial dimensions
palign <- SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
                         images = c("slice1", "slice2", "slice3", "slice4")) & 
  theme(legend.position = "none",)

# Visualize original images
porig <- SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
                        images = c("slice1.1", "slice1.2", "slice1.3", "slice1.4")) & 
  theme(legend.position = "none")

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_region.png"), plot = palign)
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_region.png"), plot = porig)

SpatialDimPlot(pacientes.seurat, alpha = 0.5, crop = FALSE, ncol = 2,
               images = c("slice1.3", "slice3")) & 
  theme(legend.position = "none",)

############################################################################### 
#                     Visualization 1                                         #
###############################################################################

pmt <- SpatialPlot(object = pacientes.seurat, images = c("slice1", "slice2", "slice3", "slice4"), 
                   features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 2) & 
  theme(legend.position = "top")
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_genesButterfly_trans.png"), plot = pmt)

pmt_orig <- SpatialPlot(object = pacientes.seurat, images = c("slice1.1", "slice1.2", "slice1.3", "slice1.4"), 
                        features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 2) & 
  theme(legend.position = "top")
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_genesButterfly.png"), plot = pmt_orig)

############################################################################### 
#                     Save Objects                                            #
###############################################################################

saveRDS(pacientes.seurat, paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,".rds"))
pacientes.seurat <- readRDS(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,".rds"))
pacientes.semla <- pacientes.seurat
saveRDS(Evaluation, paste0(saveDir, "/results/patient",patient,"_EVALUATION_merge_",modeAbr,".rds"))

############################################################################### 
#                Create objects for deconvolution                             #
############################################################################### 

# Obtain the individual object of the slide
# Splits the Seurat object 'paciente19.merge' into a list of Seurat objects based on the 'name' metadata
pacientes.seurat.split.orig <- SplitObject(paciente.merge, split.by = "name")

# Generate a spatial dimensional plot for the second object in the split list
# The alpha parameter controls the transparency of the plot, and crop determines whether to crop the image
SpatialDimPlot(pacientes.seurat.split.orig[[2]], alpha = 0.5, crop = FALSE) & 
  theme(legend.position = "none",)

# Create a copy of the split objects for transformed images
pacientes.seurat.split.trans <- pacientes.seurat.split.orig

# Assign transformed images to the corresponding objects in the new split list
pacientes.seurat.split.trans[[2]]@images$slice1.2 <- pacientes.seurat@images$slice2
pacientes.seurat.split.trans[[3]]@images$slice1.3 <- pacientes.seurat@images$slice3
pacientes.seurat.split.trans[[4]]@images$slice1.4 <- pacientes.seurat@images$slice4

# Generate spatial plots for the "percent.butterfly" feature for the transformed and original objects
p1 <- SpatialPlot(object = pacientes.seurat.split.trans[[3]], 
                  features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 1) & 
  theme(legend.position = "none")

p2 <- SpatialPlot(object = pacientes.seurat.split.orig[[3]], 
                  features = "percent.butterfly", alpha = 1, crop = FALSE, ncol = 1) & 
  theme(legend.position = "none")

# Print both plots side by side
print(p1 + p2)

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
  
  savename <- paste0("/results/patient",patient,"_merge_",modeAbr,"_list_im", i, ".rds")
  saveRDS(listaObjDeconv, paste0(saveDir, savename))
}

############################################################################### 
#                    Alignment Evaluation metrics                             #
#                                                                             #    
# This graphs are better performed in the script: graphs.R, where they are    #
# more complete in terms of error bar and p-values, and are visualized as a   #
# single graph for the three methods.                                         #
###############################################################################

###############################################################################
#      file graphs.R has automatized the graphs combining the different       # 
#                      alignment methodologies used                           #
###############################################################################

Evaluation <- readRDS(paste0(saveDir, "/results/patient",patient,"_EVALUATION_merge_",modeAbr,".rds"))

# Initialize vectors to store MSE, SSIM, and Euclidean values for original and transformed images
MSE_controlPos <- c()
MSE_controlNeg <- c()
MSE_controlMov <- c()
MSE_original <- c()
MSE_transformado <- c()
MSE_gray_controlPos <- c()
MSE_gray_controlNeg <- c()
MSE_gray_controlMov <- c()
MSE_gray_original <- c()
MSE_gray_transformado <- c()
SSIM_controlPos <- c()
SSIM_controlNeg <- c()
SSIM_controlMov <- c()
SSIM_original <- c()
SSIM_transformado <- c()
Eucl_original <- c()
Eucl_transformado <- c()

# Loop through the first three images to extract MSE, SSIM, and Euclidean values
for (i in 1:3) {
  MSE_original <-          c(MSE_original, Evaluation$original$sameRegion_samePoint[[i]]$mse_value)
  MSE_transformado <-      c(MSE_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$mse_value)
  MSE_gray_original <-     c(MSE_gray_original, Evaluation$original$sameRegion_samePoint[[i]]$mse_gray_value)
  MSE_gray_transformado <- c(MSE_gray_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$mse_gray_value)
  SSIM_original <-         c(SSIM_original, Evaluation$original$sameRegion_samePoint[[i]]$ssim_value)
  SSIM_transformado <-     c(SSIM_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$ssim_value)
  Eucl_original <-         c(Eucl_original, Evaluation$original$sameRegion_samePoint[[i]]$Eucl_value)
  Eucl_transformado <-     c(Eucl_transformado, Evaluation$transformed$sameRegion_samePoint[[i]]$Eucl_value)
}

# Loop through the control data to extract MSE and SSIM values
for (i in 2:4) {
  MSE_controlPos <-      c(MSE_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_value)
  MSE_controlNeg <-      c(MSE_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_value)
  MSE_controlMov <-      c(MSE_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_value)
  MSE_gray_controlPos <- c(MSE_gray_controlPos, Evaluation$control[[i]]$`Positive control solution`$mse_gray_value)
  MSE_gray_controlNeg <- c(MSE_gray_controlNeg, Evaluation$control[[i]]$`Negative control solution`$mse_gray_value)
  MSE_gray_controlMov <- c(MSE_gray_controlMov, Evaluation$control[[i]]$`Movement control solution`$mse_gray_value)
  SSIM_controlPos <-     c(SSIM_controlPos, Evaluation$control[[i]]$`Positive control solution`$ssim_value)
  SSIM_controlNeg <-     c(SSIM_controlNeg, Evaluation$control[[i]]$`Negative control solution`$ssim_value)
  SSIM_controlMov <-     c(SSIM_controlMov, Evaluation$control[[i]]$`Movement control solution`$ssim_value)
}

# Create a data frame to hold the values from the same region and same point
data_sameRegion_samePoint <- data.frame(
  Image = c("Image2", "Image3", "Image4"),
#  MSE_controlPos = MSE_controlPos,
#  MSE_controlNeg = MSE_controlNeg,
  MSE_controlMov = MSE_controlMov,
  MSE_original = MSE_original,
  MSE_transformado = MSE_transformado,
#  MSE_gray_controlPos = MSE_gray_controlPos,
#  MSE_gray_controlNeg = MSE_gray_controlNeg,
  MSE_gray_controlMov = MSE_gray_controlMov,
  MSE_gray_original = MSE_gray_original,
  MSE_gray_transformado = MSE_gray_transformado,
#  SSIM_controlPos = SSIM_controlPos,
#  SSIM_controlNeg = SSIM_controlNeg,
  SSIM_controlMov = SSIM_controlMov,
  SSIM_original = SSIM_original,
  SSIM_transformado = SSIM_transformado,
  Eucl_original = Eucl_original,
  Eucl_transformado = Eucl_transformado
)

# Set row names for the data frame
rownames(data_sameRegion_samePoint) <- c("Image2", "Image3", "Image4")
data_sameRegion_samePoint

library(ggplot2)
# Create a scatter plot for the differnet parameters
pMSE <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = MSE_controlPos, color = "MSE Positive Control"), size = 5) +
  #  geom_point(aes(y = MSE_controlNeg, color = "MSE Negative Control"), size = 5) +
  geom_point(aes(y = MSE_controlMov, color = "MSE Movement Control"), size = 5) +
  geom_point(aes(y = MSE_original, color = "MSE original"), size = 5) +
  geom_point(aes(y = MSE_transformado, color = "MSE transformed"), size = 5) +
  labs(x = "", y = "MSE value") +
  scale_color_manual(name = "", values = c("MSE Positive Control" = "orange","MSE Negative Control" = "red", "MSE Movement Control" = "blue", "MSE original" = "green", "MSE transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pMSEgray <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = MSE_gray_controlPos, color = "MSE gray Positive Control"), size = 5) +
  #  geom_point(aes(y = MSE_gray_controlNeg, color = "MSE gray Negative Control"), size = 5) +
  geom_point(aes(y = MSE_gray_controlMov, color = "MSE gray Movement Control"), size = 5) +
  geom_point(aes(y = MSE_gray_original, color = "MSE gray Original"), size = 5) +
  geom_point(aes(y = MSE_gray_transformado, color = "MSE gray Transformed"), size = 5) +
  labs(x = "", y = "MSE gray value") +
  scale_color_manual(name = "", values = c("MSE gray Positive Control" = "orange", "MSE gray Negative Control" = "red", "MSE gray Movement Control" = "blue", "MSE gray Original" = "green", "MSE gray Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pSSIM <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = SSIM_controlPos, color = "SSIM Positive Control"), size = 5) +
  geom_point(aes(y = SSIM_controlNeg, color = "SSIM Negative Control"), size = 5) +
  geom_point(aes(y = SSIM_controlMov, color = "SSIM Movement Control"), size = 5) +
  geom_point(aes(y = SSIM_original, color = "SSIM Original"), size = 5) +
  geom_point(aes(y = SSIM_transformado, color = "SSIM Transformed"), size = 5) +
  labs(x = "", y = "SSIM value") +
  scale_color_manual(name = "", values = c("SSIM Positive Control" = "orange", "SSIM Negative Control" = "red", "SSIM Movement Control" = "blue", "SSIM Original" = "green", "SSIM Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

pEucl <- ggplot(data_sameRegion_samePoint, aes(x = Image)) +
  geom_point(aes(y = Eucl_original, color = "Euclidean Original"), size = 5) +
  geom_point(aes(y = Eucl_transformado, color = "Euclidean Transformed"), size = 5) +
  labs(x = "", y = "Euclidean value") +
  scale_color_manual(name = "", values = c("Euclidean Original" = "green", "Euclidean Transformed" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom"
  )

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_evalMSE.png"), plot = pMSE)
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_evalMSEgray.png"), plot = pMSEgray)
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_evalSSIM.png"), plot = pSSIM)
ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_evalEucl.png"), plot = pEucl)

# Calculate the mean and standard deviation
mean_values <- colMeans(data_sameRegion_samePoint[, -1])  # Ignorar la columna "Image"
sd_values <- apply(data_sameRegion_samePoint[, -1], 2, sd)  # Ignorar la columna "Image"

# Create a new dataframe for statistics
stats_df <- data.frame(
  Parameter = c("MSE Movement Control","MSE Original","MSE Tranformed",
                "MSE gray Movement Control","MSE gray Original","MSE gray Transformed",
                "SSIM Movement Control","SSIM Original","SSIM Transformed",
                "Euclidean Original","Euclidean Transformed"),
  Mean = mean_values,
  SD = sd_values
)

# Perform t-tests
mse_test <- t.test(data_sameRegion_samePoint$MSE_original,
                   data_sameRegion_samePoint$MSE_transformado)
mse_p_value <- mse_test$p.value

ssim_test <- t.test(data_sameRegion_samePoint$SSIM_original,
                    data_sameRegion_samePoint$SSIM_transformado)
ssim_p_value <- ssim_test$p.value

eucl_test <- t.test(data_sameRegion_samePoint$Eucl_original,
                    data_sameRegion_samePoint$Eucl_transformado)
eucl_p_value <- eucl_test$p.value

mse_gray_test <- t.test(data_sameRegion_samePoint$MSE_gray_original,
                        data_sameRegion_samePoint$MSE_gray_transformado)
mse_gray_p_value <- mse_gray_test$p.value

library(tidyr)

# Transform the data to long format
stats_long <- stats_df %>%
  pivot_longer(cols = c(Mean, SD), names_to = "Statistic", values_to = "Value")

# Create bar plots for evaluations
if (mse_p_value <= 0.05 && !is.nan(mse_p_value)) {color_mse_p_value <- "red"} else {color_mse_p_value <- "black"}
peval1 <- ggplot(stats_df[2:3,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "MSE") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 1.5, y = stats_df[3,]$Mean + 150,
           label = paste0("|--",sprintf("%.4f", mse_p_value),"--|"), size = 8, color = color_mse_p_value)

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_eval1.png"), plot = peval1)

if (mse_gray_p_value <= 0.05 && !is.nan(mse_gray_p_value)) {color_mse_gray_p_value <- "red"} else {color_mse_gray_p_value <- "black"}
peval2 <- ggplot(stats_df[5:6,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "MSE gray") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 1.5, y = stats_df[6,]$Mean + 150,
           label = paste0("|--",sprintf("%.4f", mse_gray_p_value),"--|"), size = 8, color = color_mse_gray_p_value)

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_eval2.png"), plot = peval2)

if (ssim_p_value <= 0.05 && !is.nan(ssim_p_value)) {color_ssim_p_value <- "red"} else {color_ssim_p_value <- "black"}
peval3 <- ggplot(stats_df[8:9,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "SSIM") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 1.5, y = stats_df[9,]$Mean + 0.5,
           label = paste0("|--",sprintf("%.4f", ssim_p_value),"--|"), size = 8, color = color_ssim_p_value)

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_eval3.png"), plot = peval3)

if (eucl_p_value <= 0.05 && !is.nan(eucl_p_value)) {color_eucl_p_value <- "red"} else {color_eucl_p_value <- "black"}
peval4 <- ggplot(stats_df[10:11,], aes(x = Parameter, y = Mean)) +
  geom_bar(stat = "identity", fill = c("green","purple"), width = 0.4) +  # Barras de la media
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Barras de error
  labs(y = "Mean value", x = "", title = "Euclidean distance") +
  scale_x_discrete(labels = c("Original", "Transformed")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, face = "bold", hjust = 0.5),,
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 28)  # Tamaño del título del eje y
  ) +
  annotate("text", x = 1.5, y = stats_df[11,]$Mean + 150,
           label = paste0("|--",sprintf("%.4f", eucl_p_value),"--|"), size = 8, color = color_eucl_p_value)

ggsave(paste0(saveDir, "/results/patient",patient,"_merge_",modeAbr,"_eval4.png"), plot = peval4)
