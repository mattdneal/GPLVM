
#' Generate structured synthetic data
#'
#' @param num.samples
#' @param nrow
#' @param ncol
#' @param num.latent.variables
#' @param num.structures
#' @param function.noise
#' @param nonlinear
#' @param Z
#'
#' @return
#' @export
generate.structured.dataset <- function(num.samples=NULL,
                                        nrow,
                                        ncol,
                                        num.latent.variables=NULL,
                                        num.structures,
                                        function.noise=0,
                                        nonlinear=FALSE,
                                        Z=NULL) {
  if (is.null(Z)) {
    Z <- matrix(rnorm(num.samples * num.latent.variables), ncol=num.latent.variables)
  } else {
    num.samples <- nrow(Z)
    num.latent.variables <- ncol(Z)
  }

  if (nonlinear) {
    l <- quantile(as.matrix(dist(Z)), .25)
    K <- gplvm.SE(Z=Z, l, 1)
  } else {
    K <- Z %*% t(Z)
  }

  diag(K) <- diag(K) + function.noise

  A <- matrix(rnorm(num.structures * num.samples), ncol=num.structures)

  structure.heights <- K %*% A

  structures <- list()
  structures$x <- sample(seq(ncol), num.structures)
  structures$y <- sample(seq(nrow), num.structures)

  structures$cov.matrices <- list()
  structures$templates <- list()

  coords <- matrix(0, ncol=2, nrow=nrow * ncol)

  counter <- 0
  for (j in 1:ncol) {
    for (i in 1:nrow) {
      counter <- counter + 1
      coords[counter, ] <- c(i, j)
    }
  }

  for (i in seq(num.structures)) {
    R <- matrix(0, 2, 2)
    diag(R) <- rgamma(2, max(nrow, ncol) / 3)
    R[1, 2] <- rnorm(1, sd=max(nrow, ncol) / 3)
    structures$cov.matrix.R[[i]] <- R

    mu <- c(structures$y[i], structures$x[i])

    dists <- t(t(coords) - mu)

    structures$templates[[i]] <- (  1 / (2 * pi)
                                    * prod(diag(R))^-1
                                    * exp(  -1/2
                                            * diag(dists %*% chol2inv(R) %*% t(dists))
                                    )
    )
  }

  data <- matrix(0, nrow=num.samples, ncol=nrow * ncol)
  for (sample in seq(num.samples)) {
    for (structure in seq(num.structures)) {
      data[sample,] <-
        data[sample, ] +
        structures$templates[[structure]] * structure.heights[sample, structure]
    }
  }

  out <- list()
  out$data <- data
  out$structures <- structures
  out$structure.heights <- structure.heights
  out$A <- A
  out$K <- K
  out$Z <- Z
  return(out)
}
