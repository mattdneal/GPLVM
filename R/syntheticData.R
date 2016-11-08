
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
                                        nonlinear.lengthscale.quantile=0.25,
                                        Z=NULL) {
  if (is.null(Z)) {
    Z <- matrix(rnorm(num.samples * num.latent.variables), ncol=num.latent.variables)
  } else {
    num.samples <- nrow(Z)
    num.latent.variables <- ncol(Z)
  }

  if (nonlinear) {
    l <- quantile(as.matrix(dist(Z)), nonlinear.lengthscale.quantile)
    K <- gplvm.SE(Z=Z, l, 1)
  } else {
    K <- Z %*% t(Z)
  }

  diag(K) <- diag(K) + function.noise

  A <- matrix(rnorm(num.structures * num.samples), ncol=num.structures)

  structure.heights <- K %*% A

  structures <- list()
  structures$x <- sample(seq(ncol), num.structures, replace = TRUE)
  structures$y <- sample(seq(nrow), num.structures, replace = TRUE)

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
    structures$templates[[i]] <- structures$templates[[i]] / max(as.numeric(structures$templates[[i]]))
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


#' Sample from a structured GPLVM model
#'
#' @param Z
#' @param nrow.X
#' @param K_S
#' @param K_Z
#'
#' @return
#' @export
#'
#' @examples
#' Z <- rnorm(50, 2)
#' Z <- c(Z, rnorm(50, -2))
#' Z <- cbind(Z, rnorm(100))
#' K_Z <- GPLVM:::gplvm.SE(Z, 1, 1, 0)
#' K_S <- GPLVM:::structured.kernel.Matrix(50, 50, c(10,10), 1)
#' K_S <- K_S - Diagonal(nrow(K_S))
#' X <- sample.from.model(Z, 50, K_S, K_Z)
sample.from.model <- function(Z, nrow.X, K_S, K_Z) {
  K_S.chol <- chol(K_S, pivot=FALSE)
  K_Z.chol <- chol(K_Z, pivot=FALSE)
  X <- matrix(rnorm(nrow(K_S.chol) * nrow(K_Z.chol)), ncol=ncol(K_S.chol))
  X <- t(K_Z.chol) %*% X %*% K_S.chol
  layout(matrix(1:8, 2))
  for (i in sample(nrow(Z), 8)) image(matrix(X[i,], nrow.X), zlim=range(X))
  return(X)
}

call.sample.from.model <- function() {

  Z <- rnorm(50, 2)
  Z <- c(Z, rnorm(50, -2))
  Z <- cbind(Z, rnorm(100))
  plot(Z)
}
