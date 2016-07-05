gplvm.SE.dist <- function(dist.matrix, l, alpha, sigma=0) {
  K <- alpha^2 * exp(-dist.matrix^2/(2*l^2))
  diag(K) <- diag(K) + sigma^2
  return(K)
}


gplvm.SE <- function(Z, l, alpha, sigma=0) {
  gplvm.SE.dist(dist.matrix=as.matrix(dist(Z)),
                l=l,
                alpha=alpha,
                sigma=sigma)
}

dK.dZij <- function(Z, K, i, j, l) {
  out <- K
  out[] <- 0
  out[i, ] <- out[, i] <- (Z[, j] - Z[i, j]) / l^2 * K[i, ]
  return(out)
}

dK.dl <- function(Z, K, l) {
  d <- as.matrix(dist(Z)^2)
  out <- d / l^3 * K
  return(out)
}

dK.dalpha <- function(Z, K, alpha, sigma) {
  if (alpha == 0) {
    out <- 0
  } else {
    out <- 2 * (K - diag(sigma^2, nrow(K))) / alpha
  }
  return(out)
}

###############################################################################
## Structured GPLVM Kernels
###############################################################################

#q=2 from the list of compact kernels in rasmussen
compact.kernel <- function(x1, x2, C) {
  r <- sqrt(sum(((x1-x2) / C)^2))
  j <- floor(length(x1) / 2) + 3
  (max(0, 1 - r)^(j+2) * ((j^2 + 4*j + 3) * r^2 + (3 * j + 6) * r + 3) / 3)
}

#q=2 from the list of compact kernels in rasmussen
compact.kernel.grad <- function(x1, x2, C) {
  r <- sqrt(sum(((x1-x2) / C)^2))
  j <- floor(length(x1) / 2) + 3

  K <- (max(0, 1 - r)^(j+2) * ((j^2 + 4*j + 3) * r^2 + (3 * j + 6) * r + 3) / 3)

  if (r == 0) {
    dr.dc <- rep(0, length(C))
  } else {
    dr.dc <- -1 / r * (x1 - x2)^2 / C^3
  }

  dK.dr <- (-(j + 2) * max(0, 1 - r)^(j+1) * (  (j^2 + 4 * j + 3) * r^2
                                                + (3 * j + 6) * r
                                                + 3) / 3
            + max(0, 1 - r)^(j+2) * (  2 * (j^2 + 4 * j + 3) * r
                                       + 3 * j + 6) / 3
  )
  dK.dC <- dK.dr * dr.dc

  names(dK.dC) <- paste("C", seq_along(C), sep="_")

  return(dK.dC)
}

#' Structured data kernel
#'
#' @param nrows
#' @param ncols
#' @param C
#' @param alpha
#' @param grad
#'
#' @return
#'
#' @import Matrix
structured.kernel.Matrix <- function(nrows, ncols, C, alpha, grad=0) {

  cutoff <- ceiling(abs(C))
  proto.cov.mat <- matrix(0, nrow=cutoff[1], ncol=cutoff[2])

  for (i in 1:cutoff[1]) {
    for (j in 1:cutoff[2]) {
      if (!grad) {
        proto.cov.mat[i, j] <- compact.kernel(c(1, 1), c(i, j), C)
      } else {
        proto.cov.mat[i, j] <- compact.kernel.grad(c(1, 1), c(i, j), C)[grad]
      }
    }
  }

  proto.cov.mat <- alpha^2 * proto.cov.mat

  if (!grad) {
    proto.cov.mat[1,1] <- proto.cov.mat[1,1] + 1
  }

  Kdim <- nrows * ncols
  ra.diag <- rep(proto.cov.mat[1,1], Kdim)
  ja.diag <- 1:Kdim
  ia.diag <- 1:Kdim

  ra.offdiag <- c()
  ja.offdiag <- c()
  ia.offdiag <- c()

  for (i in 1:cutoff[1]) {
    for (j in 1:cutoff[2]) {
      if (i != 1 | j != 1) {
        starti <- i + nrows * (j - 1)
        endi <- nrows * j
        ia.temp <- starti:endi
        ja.temp <- seq_along(ia.temp)
        blocks <- rep(0:(ncols-j) * nrows, each=length(ia.temp))
        ia.temp <- ia.temp + blocks
        ja.temp <- ja.temp + blocks

        ra.temp <- rep(proto.cov.mat[i, j], length(ia.temp))
        ra.offdiag <- c(ra.offdiag, ra.temp)
        ia.offdiag <- c(ia.offdiag, ia.temp)
        ja.offdiag <- c(ja.offdiag, ja.temp)

        if (j != 1 & i != 1) {
          startj <- i
          endj <- nrows
          ja.temp <- startj:endj
          ia.temp <- seq_along(ja.temp) + nrows * (j - 1)
          blocks <- rep(0:(ncols-j) * nrows, each=length(ia.temp))
          ia.temp <- ia.temp + blocks
          ja.temp <- ja.temp + blocks

          ra.temp <- rep(proto.cov.mat[i, j], length(ia.temp))
          ra.offdiag <- c(ra.offdiag, ra.temp)
          ia.offdiag <- c(ia.offdiag, ia.temp)
          ja.offdiag <- c(ja.offdiag, ja.temp)

        }
      }
    }
  }

  ra <- c(ra.diag, ra.offdiag)
  ja <- c(ja.diag, ja.offdiag)
  ia <- c(ia.diag, ia.offdiag)

  K <- sparseMatrix(i=as.integer(ia),
                    j=as.integer(ja),
                    x=ra,
                    dims=c(Kdim, Kdim),
                    symmetric = T)

  return(K)
}
