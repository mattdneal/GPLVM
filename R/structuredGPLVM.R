
getNeighbourIndices <- function(index, width, height, matrixDim) {
  rowNeighbours <- (index - width):(index + width)
  keepRowNeighbours <- ceiling(rowNeighbours / 512) == ceiling(index / 512)
  rowNeighbours <- rowNeighbours[keepRowNeighbours]

  neighbours <- c()
  for (i in 0:height) {
    neighbours <- c(neighbours, rowNeighbours + matrixDim[1] * i)
    if (i != 0) {
      neighbours <- c(neighbours, rowNeighbours - matrixDim[1] * i)
    }
  }

  numPixels <- prod(matrixDim)
  matrixNum <- ceiling(index / numPixels)

  keep <- which(neighbours > numPixels * (matrixNum - 1)
                & neighbours <= numPixels * matrixNum
                & neighbours != index)

  neighbours <- neighbours[keep]
  return(neighbours)
}

structured.likelihood.calc <- function(X, K_Z, K_S, K_Z.chol, K_S.chol) {
  n <- nrow(K_Z)
  p <- nrow(K_S)
  A_S <- solve(K_S.chol, t(X))
  A_Z <- backsolve(K_Z.chol, forwardsolve(t(K_Z.chol), X))
  log.det.K_S <- as.numeric(determinant(K_S.chol, logarithm=T)$modulus) * 2
  log.det.K_Z <- 2 * sum(log(diag(K_Z.chol)))
  L <- -1 / 2 * (   n * p * log(2 * pi)
                    + n * log.det.K_S
                    + p * log.det.K_Z
                    + sum(t(A_S) * A_Z)
  )
  print("likelihood:")
  print(L)
  return(L)
}

structured.likelihood <- function(X, nrows,
                               se.l, se.alpha, se.sigma,
                               structured.C, structured.alpha,
                               Z) {
  print("pars:")
  print(c(se.l, se.alpha, se.sigma,
          structured.C, structured.alpha))
  ncols <- ncol(X) / nrows
  if (ncols != floor(ncols)) {
    stop("Incorrect nrows specified")
  }

  K_Z <- gplvm.SE(Z=Z,
            l=se.l,
            alpha=se.alpha,
            sigma=se.sigma)

  K_Z.chol <- chol(K_Z)

  K_S <- structured.kernel.Matrix(nrows=nrows,
                               ncols=ncols,
                               C=structured.C,
                               alpha=structured.alpha)

  K_S.chol <- Cholesky(K_S, perm=T)

  return(structured.likelihood.calc(X=X, K_Z=K_Z, K_S=K_S, K_Z.chol=K_Z.chol, K_S.chol=K_S.chol))
}

structured.likelihood.optimx <- function(par, X, nrows.X,
                                      probabilistic.trace.estimate=TRUE) {
  structured.likelihood(X, nrows.X,
                     se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                     structured.C=par[4:5], structured.alpha=par[6],
                     Z=matrix(par[-(1:6)], nrow=nrow(X)))
}


structured.likelihood.fixedZ.optimx <- function(par, X, nrows.X, Z,
                                             probabilistic.trace.estimate=TRUE) {
  structured.likelihood(X, nrows.X,
                     se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                     structured.C=par[4:5], structured.alpha=par[6],
                     Z=Z)
}

structured.likelihood.fixedZC.optimx <- function(par, X, nrows.X, Z, structured.C,
                                              probabilistic.trace.estimate=TRUE) {
  structured.likelihood(X, nrows.X,
                     se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                     structured.C=structured.C, structured.alpha=par[4],
                     Z=Z)
}

structured.likelihood.fixedC.optimx <- function(par, X, nrows.X, structured.C,
                                             probabilistic.trace.estimate=TRUE) {
  structured.likelihood(X, nrows.X,
                     se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                     structured.C=structured.C, structured.alpha=par[4],
                     Z=matrix(par[-(1:4)], nrow=nrow(X)))
}

structured.likelihood.fixed.K_S.optimx <- function(par, X, nrows, K_S, K_S.chol) {
  print(par[1:3])
  Z <- matrix(par[-(1:3)], nrow=nrow(X))
  se.l <- par[1]
  se.alpha <- par[2]
  se.sigma <- par[3]
  K_Z <- gplvm.SE(Z=Z,
            l=se.l,
            alpha=se.alpha,
            sigma=se.sigma)

  K_Z.chol <- chol(K_Z)

  structured.likelihood.calc(X, K_Z, K_S, K_Z.chol, K_S.chol)
}



estimate.trace <- function(FUN, n, d, ...) {
  # From Hutchinson 1990
  V <- Matrix(sample(c(-1,1), n * d, replace = T), nrow=d)
  trace.estimates <- FUN(V, ...)
  trace <- mean(trace.estimates)
  se <- sd(trace.estimates) / sqrt(n)
  return(list(trace=trace, se=se))
}


A.inv..B.trace.estimator <- function(V, A.chol, B) {
  colSums(V * solve(A.chol, B %*% V))
}



structured.likelihood.grad.calc <- function(K_S.chol, K_Z.chol,
                                         dK_S.dtheta, dK_Z.dtheta,
                                         A_S, A_Z,
                                         trace.estimate.error.perc.thresh=0.1,
                                         probabilistic.trace.estimate=TRUE) {

  n <- nrow(K_Z.chol)
  p <- nrow(K_S.chol)

  trace.1 <- 0
  trace.2 <- 0

  final.trace.inner.value <- 0
  if (!identical(0, dK_S.dtheta)) {
    if (probabilistic.trace.estimate) {
      trace.estimate.error.perc <- Inf
      num.random.vectors <- 128
      counter <- 0
      trace.estimate.mean <- 0
      trace.estimate.var <- 0
      while (trace.estimate.error.perc > trace.estimate.error.perc.thresh) {
        trace.estimate <- estimate.trace(A.inv..B.trace.estimator,
                                         num.random.vectors,
                                         nrow(dK_S.dtheta),
                                         A.chol=K_S.chol,
                                         B=dK_S.dtheta)
        trace.estimate.mean <- (  (  trace.estimate.mean * counter
                                     + trace.estimate$trace * num.random.vectors)
                                  / (counter + num.random.vectors))

        trace.estimate.var <- (   (  trace.estimate.var * (counter - 1)
                                     + trace.estimate$se^2 * num.random.vectors * (num.random.vectors - 1))
                                  / (counter + num.random.vectors - 1))

        counter <- counter + num.random.vectors

        trace.estimate.se <- sqrt(trace.estimate.var / counter)
        if (trace.estimate.mean == 0 & trace.estimate.se == 0) {
          trace.estimate.error.perc <- 0
        } else {
          trace.estimate.error.perc <- trace.estimate.se / abs(trace.estimate.mean)
        }
      }

      cat("Trace estimation SE: ", signif(trace.estimate.se, 3), "; Num RVs: ", counter, fill=T)

      trace.1 <- trace.1 + n * trace.estimate.mean

    } else {
      # We need to do it like this because it starts paging and slows right down otherwise (not enough memory)
      min.col <- 1
      increment <- 5000
      while (min.col < ncol(dK_S.dtheta)) {
        max.col <- min(ncol(dK_S.dtheta), min.col + increment - 1)
        K_S.inv..dK_S.dtheta.temp <- solve(K_S.chol,
                                           dK_S.dtheta[,min.col:max.col])
        trace.1 <- trace.1 + n * sum(diag(K_S.inv..dK_S.dtheta.temp[min.col:max.col, ]))
        min.col <- min.col + increment
      }
    }

    final.trace.inner.value <- final.trace.inner.value + dK_S.dtheta %*% A_S

  }

  if (!identical(0, dK_Z.dtheta)) {
    K_Z.inv..dK_Z.dtheta <- backsolve(K_Z.chol, forwardsolve(t(K_Z.chol), dK_Z.dtheta))

    trace.2 <- p * sum(diag(K_Z.inv..dK_Z.dtheta))

    final.trace.inner.value <- final.trace.inner.value + t(A_Z) %*% dK_Z.dtheta
  }

  last.term.left <- solve(K_S.chol, final.trace.inner.value)



  trace.3 <- sum(A_Z * t(last.term.left))

  return(1/2 * (trace.3 - trace.2 - trace.1))
}




dL.dK_Z <- function(K_S.chol, K_Z.chol,
                    A_S, A_Z) {

  n <- nrow(K_Z.chol)
  p <- nrow(K_S.chol)

  K_Z.inv <- chol2inv(K_Z.chol)

  temp.summand <- t(backsolve(K_Z.chol, forwardsolve(t(K_Z.chol), t(A_Z %*% A_S))))

  return(-1/2 * (p * K_Z.inv - temp.summand))
}


dK_Z.dZij <- function(Z, K, i, j, l) {
  out <- K
  out[] <- 0
  for (n in 1:nrow(K)) {
    out[i, n] <- (Z[n, j] - Z[i, j]) / l^2 * K[i, n]
    out[n, i] <- (Z[n, j] - Z[i, j]) / l^2 * K[i, n]
  }

  return(out)
}


# Grad w.r.t. all parameters bar X and nrows returned in same order as params
structured.likelihood.grad <- function(X, nrows,
                                    se.l, se.alpha, se.sigma,
                                    structured.C, structured.alpha,
                                    Z,
                                    probabilistic.trace.estimate=TRUE) {
  ncols <- ncol(X) / nrows
  if (ncols != floor(ncols)) {
    stop("Incorrect nrows specified")
  }

  K_S <- structured.kernel.Matrix(nrows=nrows,
                               ncols=ncols,
                               C=structured.C,
                               alpha=structured.alpha)

  K_S.noisefree <- K_S - Diagonal(n=nrow(K_S), x=1)

  K_S.chol <- Cholesky(K_S, perm=T)


  Z.dist <- Matrix(as.matrix(dist(Z)), sparse = F, doDiag = F)

  K_Z.noisefree <- gplvm.SE.dist(dist.matrix=Z.dist,
                                 l=se.l,
                                 alpha=se.alpha,
                                 sigma=0)

  K_Z <- K_Z.noisefree + Diagonal(nrow(K_Z.noisefree), x=se.sigma^2)

  K_Z.chol <- chol(K_Z)

  A_S <- solve(K_S.chol, t(X))
  A_Z <- backsolve(K_Z.chol, forwardsolve(t(K_Z.chol), X))

  dL.dK_Z <- dL.dK_Z(K_S.chol, K_Z.chol,
                     A_S, A_Z)

  out <- c()

  #theta == se.l
  dK_Z.dtheta <- Z.dist^2 / se.l^3 * K_Z.noisefree
  out[1] <- sum(dL.dK_Z * dK_Z.dtheta)

  #theta == se.alpha
  if (se.alpha==0) {
    dK_Z.dtheta <- 0
  } else {
    dK_Z.dtheta <- 2 * K_Z.noisefree / se.alpha
  }
  out[2] <- sum(dL.dK_Z * dK_Z.dtheta)

  #theta == se.sigma
  dK_Z.dtheta <- Diagonal(n=nrow(K_Z), x=(2 * se.sigma))
  out[3] <- sum(dL.dK_Z * dK_Z.dtheta)

  #theta == structured.C[1]
  print("C1")
  print(structured.C[1])
  dk_S.dtheta <- structured.kernel.Matrix(nrows=nrows,
                                       ncols=ncols,
                                       C=structured.C,
                                       alpha=structured.alpha,
                                       grad=1)
  out[4] <- structured.likelihood.grad.calc(K_S.chol, K_Z.chol,
                                         dk_S.dtheta, 0,
                                         A_S, A_Z,
                                         probabilistic.trace.estimate=probabilistic.trace.estimate)

  #theta == structured.C[2]
  print("C2")
  print(structured.C[2])
  dk_S.dtheta <- structured.kernel.Matrix(nrows=nrows,
                                       ncols=ncols,
                                       C=structured.C,
                                       alpha=structured.alpha,
                                       grad=2)
  out[5] <- structured.likelihood.grad.calc(K_S.chol, K_Z.chol,
                                         dk_S.dtheta, 0,
                                         A_S, A_Z,
                                         probabilistic.trace.estimate=probabilistic.trace.estimate)

  #theta == structured.alpha (we can do this using structured.kernel.Matrix but faster to just calculate
  # directly from K_S.noisefree)
  print("alpha")
  print(structured.alpha)
  if (structured.alpha == 0) {
    dk_S.dtheta <- Matrix(0, nrow=nrow(K_S), ncol=ncol(K_S))
  } else {
    dk_S.dtheta <- 2 * K_S.noisefree / structured.alpha
  }
  out[6] <- structured.likelihood.grad.calc(K_S.chol, K_Z.chol,
                                         dk_S.dtheta, 0,
                                         A_S, A_Z,
                                         probabilistic.trace.estimate=probabilistic.trace.estimate)

  #theta == Z
  dL.dZ <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  K_Z.noisefree <- as.matrix(K_Z.noisefree)
  dL.dK_Z <- as.matrix(dL.dK_Z)
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      dK_Z.dtheta <- dK_Z.dZij(Z=Z, K=K_Z.noisefree, i=i, j=j, l=se.l)
      dL.dZ[i, j] <- sum(dL.dK_Z * dK_Z.dtheta)
    }
  }
  out <- c(out, as.numeric(dL.dZ))
  return(out)
}


structured.likelihood.grad.optimx <- function(par, X, nrows.X) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                                 structured.C=par[4:5], structured.alpha=par[6],
                                 Z=matrix(par[-(1:6)], nrow=nrow(X)))
  return(out)
}


structured.likelihood.grad.fixedZ.optimx <- function(par, X, nrows.X, Z, probabilistic.trace.estimate=TRUE) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                                 structured.C=par[4:5], structured.alpha=par[6],
                                 Z=Z,
                                 probabilistic.trace.estimate=probabilistic.trace.estimate)[1:6]
  return(out)
}


structured.likelihood.grad.fixedZC.optimx <- function(par, X, nrows.X, Z, structured.C, probabilistic.trace.estimate=TRUE) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                                 structured.C=structured.C, structured.alpha=par[4],
                                 Z=Z,
                                 probabilistic.trace.estimate=probabilistic.trace.estimate)[c(1:3,6)]
  print(out)
  return(out)
}


structured.likelihood.grad.fixedC.optimx <- function(par, X, nrows.X, Z, structured.C, probabilistic.trace.estimate=TRUE) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 se.l=par[1], se.alpha=par[2], se.sigma=par[3],
                                 structured.C=structured.C, structured.alpha=par[4],
                                 Z=matrix(par[-(1:4)], nrow=nrow(X)),
                                 probabilistic.trace.estimate=probabilistic.trace.estimate)[-c(4,5)]
  return(out)
}


K_S.param.L.optimx <- function(par, X, nrows.X) {
  structured.likelihood(X, nrows.X,
                     se.l=1, se.alpha=0, se.sigma=1,
                     structured.C=par[1:2], structured.alpha=par[3],
                     Z=matrix(0, nrow=nrow(X), ncol=1))
}


K_S.param.L.grad.optimx <- function(par, X, nrows.X) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 1, 0, 1,
                                 par[1:2], par[3],
                                 matrix(0, nrow=nrow(X), ncol=1))[3:6]
  return(out)
}


K_S.param.fixedC.L.optimx <- function(par, X, nrows.X, structured.C) {
  out <- structured.likelihood(X, nrows.X,
                            se.l=1, se.alpha=0, se.sigma=1,
                            structured.C=structured.C, structured.alpha=par[1],
                            Z=matrix(0, nrow=nrow(X), ncol=1))
  return(out)
}


K_S.param.fixedC.L.grad.optimx <- function(par, X, nrows.X, structured.C) {
  out <- structured.likelihood.grad(X, nrows.X,
                                 1, 0, 1,
                                 structured.C, par[1],
                                 matrix(0, nrow=nrow(X), ncol=1))[6]
  return(out)
}


optimize.structured.params <- function(X, nrows.X, starting.par=c(c1=3,c2=3,alpha=1)) {
  optimx(par=starting.par, fn=K_S.param.L.optimx, gr=K_S.param.L.grad.optimx,
         method="L-BFGS-B",
         lower=c(1,1,0.01), upper=c(20, 5, Inf),
         control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE),
         X=X, nrows.X=nrows.X)
}

optimize.structured.params.fixedC <- function(X, nrows.X, structured.C=c(10, 1), starting.par=c(alpha=1)) {
  optimx(par=starting.par, fn=K_S.param.fixedC.L.optimx, gr=K_S.param.fixedC.L.grad.optimx,
         method="L-BFGS-B",
         lower=c(0.01),
         control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE),
         X=X, nrows.X=nrows.X, structured.C=structured.C)
}


structured.gplvm <- function(X, nrows.X,
                          se.par=c(3, 1, 1),
                          structured.par=c(20, 1, 1),
                          Z=NULL,
                          q=2) {
  if (is.null(Z)) {
    Z <- matrix(rnorm(q * nrow(X)), nrow=nrow(X))
  }

  par <- c(se.par, structured.par)
  lower <- rep(-Inf, length(par))
  upper <- rep(Inf, length(par))

  lower[4:6] <- c(1, 1, 0.01)
  upper[4:6] <- c(20, 5, Inf)


  optimx.obj <- optimx(par=par,
                       fn=structured.likelihood.fixedZ.optimx,
                       gr=structured.likelihood.grad.fixedZ.optimx,
                       method="L-BFGS-B",
                       lower=lower, upper=upper,
                       control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE),
                       X=X, nrows.X=nrows.X, Z=Z)

  par <- c(optimx.obj[1:length(par)], as.numeric(Z))
  lower <- rep(-Inf, length(par))
  upper <- rep(Inf, length(par))

  lower[4:6] <- c(1, 1, 0.01)
  upper[4:6] <- c(20, 5, Inf)

  optimx(par=par,
         fn=structured.likelihood.optimx,
         gr=structured.likelihood.grad.optimx,
         method="L-BFGS-B",
         lower=lower, upper=upper,
         control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE),
         X=X, nrows.X=nrows.X)
}



#' Fit a structured GPLVM model with fixed structure covariance
#'
#' @param X
#' @param nrows.X
#' @param se.par
#' @param structured.par
#' @param Z
#' @param q
#' @param probabilistic.trace.estimate
#' @param maxit
#'
#' @return
#' @export
#'
#' @examples
#' Z.true <- cbind(rnorm(50, mean=-1), rnorm(50, mean=-1), rnorm(50, mean=-1))
#' Z.true <- rbind(Z.true, cbind(rnorm(50, mean=1), rnorm(50, mean=1), rnorm(50, mean=1)))
#' classes <- rep(1:2, each=50)
#' pairs(Z.true, col=classes, pty=".")
#' X <- generate.structured.dataset(100,20,20,num.latent.variables=NULL,5,0.1,nonlinear=TRUE, Z=Z.true)
#' X$data <- scale(X$data, center=T, scale=F)
#' X.sd <- sd(as.numeric(X$data))
#' X$data <- X$data / X.sd
#' X.noisy <- X$data + rnorm(length(X$data))
#' Z.init <- prcomp(X.noisy)$x[,1:5]
#' sgplvm.model <- structured.gplvm.fixedC(X=X$data, nrows.X=20, probabilistic.trace.estimate=F, structured.par=c(2, 2, 1), Z=Z.init, maxit=1000, q=5)
structured.gplvm.fixedC <- function(X, nrows.X,
                                 se.par=c(3, 1, 1),
                                 structured.par=c(20, 1, 1),
                                 Z=NULL,
                                 q=2,
                                 probabilistic.trace.estimate=TRUE,
                                 maxit=100) {
  if (is.null(Z)) {
    Z <- matrix(rnorm(q * nrow(X)), nrow=nrow(X))
  }

  par <- c(se.par, structured.par[3])
  structured.C <- structured.par[1:2]

  optimx.obj <- optimx(par=par,
                       fn=structured.likelihood.fixedZC.optimx,
                       gr=structured.likelihood.grad.fixedZC.optimx,
                       method="BFGS",
                       control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE, type=2, maxit=maxit),
                       X=X, nrows.X=nrows.X, Z=Z, structured.C=structured.C,
                       probabilistic.trace.estimate=probabilistic.trace.estimate)

  par <- c(as.numeric(optimx.obj[1:length(par)]), as.numeric(Z))

  optimx.obj <- optimx(par=par,
                       fn=structured.likelihood.fixedC.optimx,
                       gr=structured.likelihood.grad.fixedC.optimx,
                       method="L-BFGS-B",
                       control=list(maximize=TRUE, trace=1, kkt=FALSE, starttests=FALSE, type=2, maxit=maxit),
                       X=X, nrows.X=nrows.X, structured.C=structured.C,
                       probabilistic.trace.estimate=probabilistic.trace.estimate)

  par <- head(as.numeric(optimx.obj), length(par))
  Z <- matrix(par[seq_along(Z) + 4], ncol=q)
  se.par <- par[1:3]
  names(se.par) <- c("l", "alpha", "sigma")
  structured.par <- c(structured.C, par[4])
  names(structured.par) <- c("C1", "C2", "alpha")

  out <- list(Z=Z, se.par=se.par, structured.par=structured.par)
  return(out)
}


L.Ztest <- function(Z, X, nrows,
                    se.l, se.alpha, se.sigma,
                    structured.C, structured.alpha) {
  structured.likelihood(X, nrows,
                     se.l, se.alpha, se.sigma,
                     structured.C, structured.alpha,
                     Z)
}
