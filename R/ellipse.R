#' 
#' Minimum Volume Enclosing Ellipsoid
#' 
#' @description Uses the Khachiyan Algorithm to find the the minimum volume enclosing
#' ellipsoid (MVEE) of a set of points. In two dimensions, this is just
#' the bounding ellipse (this function is limited to two dimensions).
#' Adapted by Andy Lyons from Matlab code by Nima Moshtagh.
#' Copyright (c) 2009, Nima Moshtagh
#' 
#' @param xy a two-column data frame containing x and y coordinates if NULL then a random sample
#' set of 10 points will be generated
#' @param tolerance a tolerance value (default = 0.005)
#' @param plotme FALSE/TRUE. Plots the points and ellipse. Default TRUE.
#' @param max.iter The maximum number of iterations. If the script tries this number of iterations
#' but still can't get to the tolerance value, it displays an error message and returns NULL
#' @param shiftxy TRUE/FALSE. If True, will apply a shift to the coordinates to make them
#' smaller and speed up the matrix calculations, then reverse the shift to the center point of the
#' resulting ellipoid. Default TRUE
#' 
#' @return Output A list containing the "center form" matrix equation of the ellipse. i.e. a 2x2
#' matrix "A" and a 2x1 vector "C" representing the center of the ellipse such that:
#' (x - C)' A (x - C) <= 1. Also in the list is a 2x1 vector elps.axes.lngth whose elements
#' are  one-half  the lengths of the major and minor axes (variables a and b). Also in list is
#' alpha, the angle of rotation.
#' 
#
# @references https://stat.ethz.ch/pipermail/r-help/2011-March/272997.html
#             http://www.mathworks.com/matlabcentral/fileexchange/9542
#             http://www.mathworks.com/matlabcentral/fileexchange/13844
#             http://stackoverflow.com/questions/1768197/bounding-ellipse
#             http://stackoverflow.com/questions/1768197/bounding-ellipse
#             https://stat.ethz.ch/pipermail/r-help/2006-October/114652.html
#
# TODO: remove copied code and adapt
#
mvee <- function(xy = NULL, tolerance = 0.005, plotme = TRUE, max.iter = 500, shiftxy = TRUE) {
  if (is.null(xy)) {
    xy <- data.frame(x = runif(10, 100, 200), y = runif(10, 100, 200))
  } else if (ncol(xy) != 2) {
    warning("xy must be a two-column data frame")
    return(NULL)
  }
  ## Number of points
  n <- nrow(xy)
  ## Dimension of the points (2)
  d <- ncol(xy)

  if (n <= d) return(NULL)

  ##  Apply  a  uniform shift to the x&y coordinates to make matrix calculations computationally
  ## simpler (if x and y are very large, for example UTM coordinates, this may be necessary
  ## to prevent a 'compuationally singular matrix' error
  if (shiftxy) {
    xy.min <- sapply(xy, FUN = "min")
  } else {
    xy.min <- c(0, 0)
  }
  xy.use <- xy - rep(xy.min, each = n)
  ## Add a column of 1s to the (n x 2) matrix xy - so it is now (n x 3)
  Q <- t(cbind(xy.use, rep(1, n)))
  ## Initialize
  count <- 1
  err <- 1
  u <- rep(1 / n, n)
  ## Khachiyan Algorithm
  while (err > tolerance) {
    # see [4]http://stackoverflow.com/questions/1768197/bounding-ellipse
    # for commented code
    X <- Q %*% diag(u) %*% t(Q)
    M <- diag(t(Q) %*% solve(X) %*% Q)
    maximum <- max(M)
    j <- which(M == maximum)
    step_size <- (maximum - d - 1) / ((d + 1) * (maximum - 1))
    new_u <- (1 - step_size) * u
    new_u[ j ] <- new_u[ j ] + step_size
    err <- sqrt(sum((new_u - u)^2))
    count <- count + 1
    if (count > max.iter) {
      warning(paste("Iterated", max.iter, "times and still can't find the bounding ellipse. \n", sep = ""))
      warning("Either increase the tolerance or the maximum number of iterations")
      return(NULL)
    }
    u <- new_u
  }
  ## Put the elements of the vector u into the diagonal of a matrix
  U <- diag(u)
  ## Take the transpose of xy
  P <- t(xy.use)
  ## Compute the center, adding back the shifted values
  c <- as.vector((P %*% u) + xy.min)
  ## Compute the A-matrix
  A <- (1 / d) * solve(P %*% U %*% t(P) - (P %*% u) %*% t(P %*% u))
  ## Find the Eigenvalues of matrix A which will be used to get the major and minor axes
  A.eigen <- eigen(A)
  ## Calculate the length of the semi-major and semi-minor axes
  ## from the Eigenvalues of A.
  semi.axes <- sort(1 / sqrt(A.eigen$values), decreasing = TRUE)
  ## Calculate the rotation angle from the first Eigenvector
  alpha <- atan2(A.eigen$vectors[ 2, 1 ], A.eigen$vectors[ 1, 1 ]) - pi / 2
  if (plotme) {
    # Plotting commands adapted from code by Alberto Monteiro
    # [5]https://stat.ethz.ch/pipermail/r-help/2006-October/114652.html
    # Create the points for the ellipse
    theta <- seq(0, 2 * pi, length = 72)
    a <- semi.axes[ 1 ]
    b <- semi.axes[ 2 ]
    elp.plot.xs <- c[ 1 ] + a * cos(theta) * cos(alpha) - b * sin(theta) *
      sin(alpha)
    elp.plot.ys <- c[ 2 ] +
      a * cos(theta) * sin(alpha) +
      b * sin(theta) *
        cos(alpha)
    # Plot the ellipse with the same scale on each axis
    plot(elp.plot.xs, elp.plot.ys, type = "l", lty = "dotted", col = "blue",
         asp = 1,
         main = "minimum volume enclosing ellipsoid", xlab = names(xy)[ 1 ],
         ylab = names(xy)[ 2 ])
    # Plot the original points
    points(xy[ , 1 ], xy[ , 2 ], type = "p", pch = 16)
    # add the center of the ellipse using a triangle symbol
    points(c[ 1 ], c[ 2 ], pch = 2, col = "blue")
  }
  ellipse.params <- list("A" = A, "c" = c, "a" = semi.axes[ 1 ], "b" = semi.axes[ 2 ], alpha = alpha)

  ellipse.params
}
