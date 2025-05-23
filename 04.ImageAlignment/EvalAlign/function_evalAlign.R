library(SpatialPack)
source("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/function_selectCoord.R")

# Convert hexadecimal color data into an RGB array with dimensions (x, y, 4).
# The function organizes color values (R, G, B, Alpha) into a 3D array.
hex2rgb_table <- function(objetHex) {
  rgb_values <- col2rgb(objetHex, alpha = TRUE)  # Convert hex to RGB and alpha
  
  # Create an empty array for storing RGBA values and populate with color channels.
  x <- ncol(objetHex)
  y <- nrow(objetHex)
  rgba_array <- array(0, dim = c(y, x, 4))
  
  # Assign RGB and alpha values to the array for each color channel.
  rgba_array[,,1] <- matrix(rgb_values[1, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,2] <- matrix(rgb_values[2, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,3] <- matrix(rgb_values[3, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,4] <- matrix(rgb_values[4, ], nrow = y, ncol = x, byrow = TRUE)
  
  return(rgba_array)
}

# Similar to the previous function, but normalizes RGB values to a 0-1 range.
hex2rgb_table_norm <- function(objetHex) {
  rgb_values <- col2rgb(objetHex, alpha = TRUE) / 255
  
  # Prepare array to store normalized RGBA values.
  x <- ncol(objetHex)
  y <- nrow(objetHex)
  rgba_array <- array(0, dim = c(y, x, 4))
  
  # Populate array with normalized color channel data.
  rgba_array[,,1] <- matrix(rgb_values[1, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,2] <- matrix(rgb_values[2, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,3] <- matrix(rgb_values[3, ], nrow = y, ncol = x, byrow = TRUE)
  rgba_array[,,4] <- matrix(rgb_values[4, ], nrow = y, ncol = x, byrow = TRUE)
  
  return(rgba_array)
}

# Calculate Mean Squared Error (MSE) between two images by comparing pixel values.
mse <- function(image1, image2) {
  sum_error <- 0  # Initialize total error
  
  # Loop through each color channel to calculate MSE.
  for (channel in seq_along(image1[1,1,]) - 1) {
    error <- sum((image1[,,channel] - image2[,,channel])^2)
    error <- error / (nrow(image1) * ncol(image2))
    sum_error <- sum_error + error
  }
  
  # Average the total error across color channels.
  total_error <- sum_error / 3
  return(total_error)
}

# Compute the grayscale MSE between two images by first converting them to grayscale.
mseGS <- function(image1, image2) {
  im1 <- RGB2gray(image1)  # Convert image 1 to grayscale
  im2 <- RGB2gray(image2)  # Convert image 2 to grayscale
  
  # Calculate grayscale MSE by summing squared differences.
  error <- sum((im1[,] - im2[,])^2) / (nrow(im1) * ncol(im2))
  return(error)
}

# Calculate Euclidean distance between two sets of coordinates.
EuclDist <- function(listaCoordenadas, nIm) {
  # Retrieve coordinates for images i and j
  i <- nIm[[1]]
  j <- nIm[[2]]
  xi <- listaCoordenadas[[i]]$x
  yi <- listaCoordenadas[[i]]$y
  xj <- listaCoordenadas[[j]]$x
  yj <- listaCoordenadas[[j]]$y
  
  # Calculate the average Euclidean distance across all points.
  distance <- mean(sqrt((xj - xi)^2 + (yj - yi)^2))
  return(distance)
}

# Evaluate alignment between two images using various metrics (MSE, SSIM, etc.).
evalAlign <- function(image1, image2, listaCoordenadas, nIm) {
  # Ensure images have the same dimensions, converting to RGB if necessary.
  if (!all(dim(image1) == dim(image2))) stop("Images must have the same dimensions.")
  if (is.na(dim(image1)[3])) image1 <- hex2rgb_table(image1)
  if (is.na(dim(image2)[3])) image2 <- hex2rgb_table(image2)
  
  # Calculate alignment metrics: MSE, grayscale MSE, RMSE, SSIM, and Euclidean distance.
  # Calculate the MSE and RMSE
  mse_value <- mse(image1, image2)
  mse_gray_value <- mseGS(image1, image2)
  rmse_value <- sqrt(mse_value)
  # Calculate the SSIM
  im1 <- RGB2gray(image1)
  im2 <- RGB2gray(image2)
  ssim_value <- SpatialPack::SSIM(im1, im2)[["SSIM"]]
  # Calculate the Euclidean distance
  Eucl_value <- EuclDist(listaCoordenadas, nIm)
  
  # Compile metrics into a solution list.
  solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                   rmse_value = rmse_value, ssim_value = ssim_value,
                   Eucl_value = Eucl_value)
  return(solution)
}

# Evaluate alignment metrics for a reference image with various controls.
controlAlign <- function(image1) { 
  # If needed, convert the image to RGB format.
  if (is.na(dim(image1)[3])) image1 <- hex2rgb_table(image1)
  
  # Define the central region for analysis.
  x1 <- round(dim(image1)[1] / 2)
  y1 <- round(dim(image1)[2] / 2)
  image1rec <- image1[(x1-200):(x1+200), (y1-200):(y1+200), ]
  
  # Calculate alignment metrics for positive control (image vs itself).
  mse_value <- mse(image1rec, image1rec)
  mse_gray_value <- mseGS(image1rec, image1rec)
  rmse_value <- sqrt(mse_value)
  im1 <- RGB2gray(image1rec)
  ssim_value <- SSIM(im1, im1)[["SSIM"]]
  positive_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Calculate alignment metrics for a negative control (image vs random image).
  library(jpeg)
  image2 <- readJPEG("/Users/vagm_110901/Documents/Universitat/ConesaLab/Lab/04.ImageAlingment/EvalAlign/image_negative.JPG")
  if (is.na(dim(image2)[3])) image2 <- hex2rgb_table(image2)
  image2 <- image2[1:dim(image1rec)[[1]], 1:dim(image1rec)[[2]], ]
  mse_value <- mse(image1rec, image2)
  mse_gray_value <- mseGS(image1rec, image2)
  rmse_value <- sqrt(mse_value)
  im2 <- RGB2gray(image2)
  ssim_value <- SSIM(im1, im2)[["SSIM"]]
  negative_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Calculate alignment metrics for movement control (image vs shifted version 10 positions to the rigth).
  image3 <- image1[(x1-190):(x1+210), (y1-200):(y1+200), ]
  mse_value <- mse(image1rec, image3)
  mse_gray_value <- mseGS(image1rec, image3)
  rmse_value <- sqrt(mse_value)
  im3 <- RGB2gray(image3)
  ssim_value <- SSIM(im1, im3)[["SSIM"]]
  movement_solution <- list(mse_value = mse_value, mse_gray_value = mse_gray_value,
                            rmse_value = rmse_value, ssim_value = ssim_value)
  
  # Compile all control solutions into one list.
  solution <- list("Positive control solution" = positive_solution,
                   "Negative control solution" = negative_solution,
                   "Movement control solution" = movement_solution)
  return(solution)
}


evaluationComplete <- function(listaCoordenadasNEW, listaCoordenadas,
                       listaRawImages, listaTransImages, patient = c('unique','multiple')) {
  
  patient <- match.arg(patient)
  
  # Initialize an evaluation list
  Evaluation <- list()
  # Define coordinates for comparison
  x1 <- listaCoordenadasNEW[[1]]$x[[5]]
  y1 <- listaCoordenadasNEW[[1]]$y[[5]]
  
  # Create a list to store control alignment parameters
  control <- list()
  
  # Loop through the new coordinates and evaluate the alignment
  for (i in 1:length(listaCoordenadasNEW)) {
    control[[i]] <- controlAlign(listaTransImages[[i]])
    print(control)
  }
  
  # Store the control evaluations in the Evaluation list
  Evaluation$control <- control
  
  # Evaluate raw images (same region, same point)
  for ( i in 1:length(listaCoordenadas)) {
    if (i != length(listaCoordenadas)) {
      for ( j in i:length(listaCoordenadas))  {
        if (i != j) {
          x1 <- listaCoordenadas[[i]]$x[[5]]
          y1 <- listaCoordenadas[[i]]$y[[5]]
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patient == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            
          } else if (patient == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
          }
          print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " comparing the same region without selecting a different point."))
          
          # Evaluate the alignment and store parameters
          parameters <- 
            evalAlign(
              listaRawImages[[i]][coordRange[[1]],coordRange[[2]],], 
              listaRawImages[[j]][coordRange[[1]],coordRange[[2]],], 
              listaCoordenadas, c(i,j))
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$original$sameRegion_samePoint[[imagescompare]] <- parameters
          
        }}}}
  
  # Evaluate raw images (common region, different point)
  for ( i in 1:length(listaCoordenadas)) {
    if (i != length(listaCoordenadas)) {
      for ( j in i:length(listaCoordenadas))  {
        if (i != j) {
          x1 <- listaCoordenadas[[i]]$x[[5]]
          y1 <- listaCoordenadas[[i]]$y[[5]]
          x2 <- listaCoordenadas[[j]]$x[[5]]
          y2 <- listaCoordenadas[[j]]$y[[5]]
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patient == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 200)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2 + 200)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
            
          } else if (patient == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 310)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
          }
          print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " selecting an area from a reference point in each image."))
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaRawImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaRawImages[[j]][coordRange2[[1]],coordRange2[[2]],],
              listaCoordenadas, c(i,j))
          }, error = function(e) {
            message(sprintf("Error en comparación %d-%d: %s", i, j, e$message))
            return(NULL)
            })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$original$commonRegion_differentPoint[[imagescompare]] <- parameters
        }}}}
  
  # Evaluate transformed images (same region, same point)
  for ( i in 1:length(listaCoordenadasNEW)) {
    if (i != length(listaCoordenadasNEW)) {
      for ( j in i:length(listaCoordenadasNEW))  {
        if (i != j) {
          x1 <- listaCoordenadasNEW[[i]]$x[[5]]
          y1 <- listaCoordenadasNEW[[i]]$y[[5]]
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patient == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            
          } else if (patient == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
          }
          print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " comparing the same region without selecting a different point."))
          
          # Evaluate the alignment and store parameters
          parameters <- 
            evalAlign(
              listaTransImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaTransImages[[j]][coordRange[[1]],coordRange[[2]],],
              listaCoordenadasNEW, c(i,j))
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$transformed$sameRegion_samePoint[[imagescompare]] <- parameters
          
        }}}}
  
  # Evaluate transformed images (common region, different point)
  for ( i in 1:length(listaCoordenadasNEW)) {
    if (i != length(listaCoordenadasNEW)) {
      for ( j in i:length(listaCoordenadasNEW))  {
        if (i != j) {
          x1 <- listaCoordenadasNEW[[i]]$x[[5]]
          y1 <- listaCoordenadasNEW[[i]]$y[[5]]
          x2 <- listaCoordenadasNEW[[j]]$x[[5]]
          y2 <- listaCoordenadasNEW[[j]]$y[[5]]
          dim1 <- dim(listaTransImages[[i]])
          dim2 <- dim(listaTransImages[[j]])
          
          if (patient == 'unique') { 
            x1_min <- max(0, x1 - 200)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1 + 200)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 200)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2 + 200)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
            
          } else if (patient == 'multiple') {
            x1_min <- max(0, x1 - 310)
            y1_min <- max(0, y1 - 200)
            x1_max <- min(dim1[1], dim2[1], x1)
            y1_max <- min(dim1[2], dim2[2], y1 + 200)
            x2_min <- max(0, x2 - 310)
            y2_min <- max(0, y2 - 200)
            x2_max <- min(dim1[1], dim2[1], x2)
            y2_max <- min(dim1[2], dim2[2], y2 + 200)
            
            coordRange  <- list(x1_min:x1_max, y1_min:y1_max)
            coordRange2 <- list(x2_min:x2_max, y2_min:y2_max)
          }
          print(paste0("Evaluation of the alignment between images ", as.character(i), " and ", as.character(j), " selecting an area from a reference point in each image."))
          
          # Evaluate the alignment and store parameters
          parameters <- tryCatch({
            evalAlign(
              listaTransImages[[i]][coordRange[[1]],coordRange[[2]],],
              listaTransImages[[j]][coordRange2[[1]],coordRange2[[2]],],
              listaCoordenadasNEW, c(i,j)) 
          }, error = function(e) {
            message(sprintf("Error en comparación %d-%d: %s", i, j, e$message))
            return(NULL)
          })
          print(parameters)
          imagescompare <- paste0("comparing_", as.character(i), '_', as.character(j))
          Evaluation$transformed$commonRegion_differentPoint[[imagescompare]] <- parameters
        }}}}
  
  return(Evaluation)
}
