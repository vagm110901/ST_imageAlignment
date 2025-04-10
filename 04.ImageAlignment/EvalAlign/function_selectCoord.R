###############################################################################
#             function selectCoord(image)                                     #
#                                                                             #
# Allows the user to select multiple points on the provided image by clicking #
# on it.                                                                      # 
# Returns a list containing the coordinates of the selected points with two   #      
# elements: x and y.                                                          #
###############################################################################
library(imager)
selectCoord <- function(image) {
  plot(image)
  coordinates <- locator(type = "p")
  return(coordinates)
}

###############################################################################
#                     function selectCoordSTalign(image)                      #        
#                                                                             #
# Allows the user to select multiple points on the provided image by clicking # 
# on it.                                                                      #
# Returns a list containing the coordinates of the selected points with two   #
# elements: x and y.                                                          #
###############################################################################
selectCoordSTalign <- function(image) {
  par(mfrow=c(1,1))
  plot(c(0,0), xlim=c(0,ncol(image)), ylim=c(0,nrow(image)),
       xlab=NA, ylab=NA)
  rasterImage(image, xleft = 0, xright = ncol(image),
              ytop = 0, ybottom = nrow(image), interpolate = FALSE)
  coordinates <- locator(type = "p")
  return(coordinates)
}