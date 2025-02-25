### function selectCoord(image)
# given an image, the user can select as point as wanted directly from de image by clicking on it
# it returns the coordinates points selected on a list with two elements: x and y
selectCoord <- function(image) {
  library(imager)
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}

### function solveCoord(coordinates)
# given the coordinates result of locator(type = "p"), 
# the values cos(ang), sen(ang), dx, dy are calculated by resolving an equation system
solveCoord <- function(coord1, coord2, xmax, ymax) {
  cx <- xmax/2
  cy <- ymax/2
  
  coord2MOD <- coord2
  for (j in seq_along(coord2$x)) {
    coord2MOD$x[[j]] <- coord2$x[[j]] - cx
    coord2MOD$y[[j]] <- coord2$y[[j]] - cy
  }
  
  coord1MOD <- coord1
  for (j in seq_along(coord2$x)) {
    coord1MOD$x[[j]] <- coord1$x[[j]] - cx
    coord1MOD$y[[j]] <- coord1$y[[j]] - cy
  }
  
  
  # optimization overdeterminated equation system
  
  fn <- function(x){
    y <- numeric(0)
    for (i in seq_along(coord1$x)) {
      y[length(y) + 1] <- coord2MOD$x[[i]]*x[1] + coord2MOD$y[[i]]*x[2] + 1*x[3] + 0*x[4] + x[1]*x[3] + x[2]*x[4] - coord1MOD$x[[i]]
      y[length(y) + 1] <- coord2MOD$y[[i]]*x[1] - coord2MOD$x[[i]]*x[2] + 0*x[3] + 1*x[4] - x[2]*x[3] + x[1]*x[4] - coord1MOD$y[[i]]
    }
    y[length(y) + 1] <- (x[1])^2 + (x[2])^2 - 1  
    return(sum(y^2))
  }
  
  xstart <- c(cos(pi/4), sin(pi/4), 2, -2)
  solution <- optim(xstart, fn)
  solution <- solution$par
  names(solution) <- c('coseno', 'seno', 'dx', 'dy')
  solution[['dx']] <- solution[['dx']]
  solution[['dy']] <- solution[['dy']]
  
 
  # non linear equations solver
  " 
  library(nleqslv)
  fn <- function(x){
    y <- numeric(0)
    for (i in seq_along(coord1$x)) {
      y[length(y) + 1] <- coord2MOD$x[[i]]*x[1] + coord2MOD$y[[i]]*x[2] + 1*x[3] + 0*x[4] - coord1MOD$x[[i]]
      y[length(y) + 1] <- coord2MOD$y[[i]]*x[1] - coord2MOD$x[[i]]*x[2] + 0*x[3] + 1*x[4] - coord1MOD$y[[i]]
    }
    y[length(y)] <- (x[1])^2 + (x[2])^2 - 1  
    y
  }

  xstart <- c(cos(pi/4), sin(pi/4), 2, -2)
  solution <- nleqslv(xstart, fn) 
  solution <- solution$x
  names(solution) <- c('coseno', 'seno', 'dx', 'dy')
  solution[['dx']] <- solution[['dx']]
  solution[['dy']] <- solution[['dy']]
  " 
  
  # linear equations solver
  "
  coseno <- c(coord2MOD$x, coord2MOD$y)
  seno <- c(coord2MOD$y, -coord2MOD$x)
  dx <- c(rep(1, length(coord2$x)), rep(0, length(coord2$y))) 
  dy <- c(rep(0, length(coord2$x)), rep(1, length(coord2$y)))
  
  A <- data.frame(coseno, seno, dx, dy)
  b <- c(coord1MOD$x, coord1MOD$y)

  solution <- solve(A,b) 
  
  solution[['dx']] <- solution[['dx']] 
  solution[['dy']] <- solution[['dy']]
  "
  
  # modelo lineal
  "
  coseno <- c(coord2MOD$x, coord2MOD$y)
  seno <- c(coord2MOD$y, -coord2MOD$x)
  dx <- c(rep(1, length(coord2$x)), rep(0, length(coord2$y))) 
  dy <- c(rep(0, length(coord2$x)), rep(1, length(coord2$y)))
  
  A <- data.frame(coseno, seno, dx, dy)
  b <- c(coord1MOD$x, coord1MOD$y)
  
  solution <- lm(b ~ coseno + seno + dx + dy - 1, data = A)
  solution <- solution$coefficients
  solution[['dx']] <- solution[['dx']] 
  solution[['dy']] <- solution[['dy']]
  "
  return(solution)
}

### function calcParameters(solution, xmax, ymax)
# given the solution, the xmax and ymax parameters 
# all the necessary parameters for semla functions are calculated
calcParameters <- function(solution, xmax, ymax) {
  # angulos
  coseno <- solution[["coseno"]]
  seno <- solution[["seno"]]
  if (coseno > 1) {coseno <- 1} else {if (coseno < -1) {coseno <- -1}}
  if (seno > 1) {seno <- 1} else {if (seno < -1) {seno <- -1}}

  tangente <- seno/coseno

  angulo <- atan2(seno, coseno) * (180/pi)

  # desplazamientos en x e y 
  trx <-  solution[["dx"]] / xmax
  try <- - solution[["dy"]] / ymax
  
  mirrorx <- solution[["mirrorx"]]
  mirrory <- solution[["mirrory"]]
  
  valores <- c(angulo, trx, try, mirrorx, mirrory)
  names(valores) <- c('angulo', 'trx', 'try', 'mirrorx', 'mirrory')

  for ( i in seq_along(valores) ) {
    if (is.na(valores[[i]])) {valores[[i]] <- 0}
  }

  return(valores)
}


