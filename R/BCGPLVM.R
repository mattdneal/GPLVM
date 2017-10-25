
bc.Z <- function(K.bc, A) {
  Z <- K.bc %*% A
  return(Z)
}

bc.L <- function(A, X, l, alpha, sigma, K.bc, Z.normal.prior) {
  Z <- bc.Z(K.bc, A)
  gplvm.L(Z, X, l, alpha, sigma, Z.normal.prior=Z.normal.prior)
}


dZ.dAij <- function(Z, K.bc, i, j) {
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  out[, j] <- K.bc[, i]
  return(out)
}

dL.dA <- function(X, Z, l, alpha, sigma, K.bc, K, dL.dK, Z.normal.prior) {
  dL.dZ <- dL.dZ(X, Z, l, alpha, sigma, K, dL.dK, Z.normal.prior=Z.normal.prior)
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dZ * dZ.dAij(Z, K.bc, i, j))
    }
  }
  return(out)
}

dL.dA.check <- function(X, A, l, alpha, sigma, K.bc, Z.normal.prior) {
  Z <- bc.Z(K.bc, A)
  K <- gplvm.SE(Z, l, alpha, sigma)
  dL.dK <- dL.dK(X, K)
  testfun <- function(A, X, l, alpha, sigma, K.bc, Z.normal.prior) {
    A <- matrix(A, nrow=nrow(X))
    bc.L(A, X, l, alpha, sigma, K.bc, Z.normal.prior=Z.normal.prior)
  }
  library(numDeriv)
  out <- list()
  out$anGrad <- dL.dA(X, Z, l, alpha, sigma, K.bc, K=K, dL.dK=dL.dK,  Z.normal.prior=Z.normal.prior)
  out$numGrad <-  matrix(grad(testfun,
                              as.numeric(A),
                              X=X,
                              l=l,
                              alpha=alpha,
                              sigma=sigma,
                              K.bc=K.bc,
                              Z.normal.prior=Z.normal.prior),
                         nrow=nrow(Z))
  out$maxDiff <- max(abs(out$numGrad-out$anGrad))
  return(out)
}


bcgplvm.f <- function(par, X, K.bc, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  A <- matrix(par[-(1:3)], nrow=nrow(X))
  bc.L(A, X, l, alpha, sigma, K.bc, Z.normal.prior=Z.normal.prior)
}

bcgplvm.gr <- function(par, X, K.bc, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  A <- matrix(par[-(1:3)], nrow=nrow(X))
  Z <- bc.Z(K.bc=K.bc, A=A)
  K <- gplvm.SE(Z, l, alpha, sigma)
  dL.dK <- dL.dK(X, K)
  c(dL.dl(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dalpha(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dsigma(X, Z, l, alpha, sigma, K, dL.dK),
    as.numeric(dL.dA(X, Z, l, alpha, sigma, K.bc, K, dL.dK, Z.normal.prior=Z.normal.prior))
  )
}

bc.S <- function(target, Z) {
  sum((Z - target)^2)
}

dS.dZ <- function(target, Z) {
  2 * (Z - target)
}

dS.dA <- function(target, Z, K.bc) {
  dS.dZ <- dS.dZ(target, Z)
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dS.dZ * dZ.dAij(Z, K.bc, i, j))
    }
  }
  return(out)
}

optimx.S <- function(A, target, K.bc) {
  A <- matrix(A, nrow=nrow(target), ncol=ncol(target))
  Z <- bc.Z(K.bc=K.bc, A=A)
  return(bc.S(target, Z))
}

optimx.S.grad <- function(A, target, K.bc) {
  A <- matrix(A, nrow=nrow(target), ncol=ncol(target))
  Z <- bc.Z(K.bc=K.bc, A=A)
  return(as.numeric(dS.dA(target, Z, K.bc)))
}

#' Fit function coefficients for target latent variable values.
#'
#' @param target
#' @param K.bc
#'
#' @return
#'
#' @importFrom optimx optimx
fit.A <- function(target, K.bc) {
  A <- tryCatch({
      K.bc.chol <- chol(K.bc)
      return(backsolve(K.bc.chol, forwardsolve(t(K.bc.chol), target)))
    },
    error=function(e) {
      par <- rnorm(length(target))
      optimx.out <- optimx(par,
                           optimx.S,
                           optimx.S.grad,
                           target=target,
                           K.bc=K.bc,
                           method="L-BFGS-B",
                           control=list(trace=T,
                                        maximize=F,
                                        kkt=FALSE,
                                        maxit=1000,
                                        starttests=TRUE)
      )
      return(matrix(as.numeric(optimx.out)[1:length(target)], nrow=nrow(target), ncol=ncol(target)))
    }
  )
  return(A)
}


#' Fit a back-constrained GPLVM model
#'
#' @param X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param Z.init
#' @param A.init
#' @param num.init.params
#' @param K.bc.l
#' @param Z.normal.prior
#'
#' @return
#' @export
#'
#' @importFrom optimx optimx
fit.bcgplvm <- function(X,
                        q,
                        iterations=1000,
                        plot.freq=100,
                        classes=1,
                        Z.init=NULL,
                        A.init=NULL,
                        num.init.params=100,
                        K.bc.l=0.1,
                        Z.normal.prior=TRUE) {

  K.bc <- gplvm.SE(Z=X, l=K.bc.l, alpha=1, sigma=0)

  if (!(is.null(Z.init)) & !(is.null(A.init))) {
    stop("Can't pass both Z.init and A.init")
  }

  if (is.null(Z.init)) {
    if (is.null(A.init)) {
      A.init <- matrix(rnorm(nrow(X)*q), ncol=q)
    } else {
      if (ncol(as.matrix(A.init)) != q) stop("Mismatch between A.init and q")
      A.init <- as.matrix(A.init)
    }
  } else if (identical(Z.init, "PCA")) {
    X.pca <- prcomp(X)
    target <- X.pca$x[,1:q]
    target.sd <- sd(as.numeric(target))
    target <- target / target.sd
    A.init <- fit.A(target = target, K.bc = K.bc)
  } else {
    if (ncol(as.matrix(Z.init)) != q) stop("Mismatch between Z.init and q")
    A.init <- fit.A(target = as.matrix(Z.init)[,1:q], K.bc = K.bc)
  }

  Z <- bc.Z(K.bc=K.bc, A=A.init)

  if (plot.freq == 0) {
    plot.freq <- iterations
  }

  if (q > 2) {
    pairs(Z, col=classes)
  } else if (q == 2) {
    plot(Z, col=classes)
  } else {
    plot(as.numeric(Z), rep(0, length(Z)), col=classes)
  }

  Z.dist <- as.matrix(dist(Z))

  dist.mean <- mean(as.numeric(Z.dist))
  dist.var <- var(as.numeric(Z.dist))

  l.rate <- dist.mean / dist.var
  l.shape <- dist.mean * l.rate

  alpha.shape <- var(as.numeric(X))

  par.init <- c(median(Z.dist), alpha.shape, 1)
  #par.init.L <- gplvm.hp.f(par.init, Z=Z, X=X)
  ##first optimize the hyperparams for the initial Z
  #for (i in seq(num.init.params)) {
  #  test.par.init <- c(rgamma(1, shape=l.shape, rate=l.rate), rgamma(1, alpha.shape), rgamma(1, 1))
  #  test.par.init.L <- gplvm.hp.f(test.par.init, Z=Z, X=X)
  #  if (test.par.init.L > par.init.L) {
  #    par.init <- test.par.init
  #    par.init.L <- test.par.init.L
  #  }
  #}


  optout <- optimx(par.init,
                   gplvm.hp.f,
                   gplvm.hp.gr,
                   Z=Z,
                   X=X,
                   Z.normal.prior=Z.normal.prior,
                   method="L-BFGS-B",
                   control=list(trace=T,
                                maximize=T,
                                kkt=FALSE,
                                maxit=1000,
                                starttests=FALSE)
  )

  init.l <- as.numeric(optout[1])
  init.alpha <- as.numeric(optout[2])
  init.sigma <- as.numeric(optout[3])

  cat("Starting params: l=", init.l, "; alpha=", init.alpha, "; sigma=", init.sigma, fill=TRUE)

  par <- c(init.l, init.alpha, init.sigma, as.numeric(A.init))

  its <- 0
  convcode=1
  while (its < iterations & convcode != 0) {
    out <- optimx(par,
                  bcgplvm.f,
                  bcgplvm.gr,
                  X=X,
                  K.bc=K.bc,
                  Z.normal.prior=Z.normal.prior,
                  method="L-BFGS-B",
                  control=list(trace=T,
                               maximize=T,
                               kkt=FALSE,
                               maxit=plot.freq,
                               starttests=FALSE)
    )
    convcode <- out$convcode
    par <- head(as.numeric(out), length(A.init) + 3)
    l <- par[1]
    alpha <- par[2]
    sigma <- par[3]
    A <- matrix(par[-(1:3)], ncol=q)
    Z <- bc.Z(K.bc=K.bc, A=A)
    if (q > 2) {
      pairs(Z,
            col=classes,
            main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else if (q == 2) {
      plot(Z,
           col=classes,
           main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    } else {
      plot(as.numeric(Z),
           rep(0, length(Z)),
           col=classes,
           main=paste("l=", signif(l, 3), "; alpha=", signif(alpha, 3), "; sigma=", signif(sigma, 3), sep=""))
    }

    its <- its + plot.freq
  }

  return(list(Z=Z, A=A, l=l, alpha=alpha, sigma=sigma, convcode=convcode))
}

#' Fit a backconstrained GPLVM model sequentially
#'
#' @param X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param num.init.params
#' @param K.bc.l
#' @param Z.normal.prior
#'
#' @return
#' @export
#'
#' @importFrom optimx optimx
fit.bcgplvm.sequential <- function(X,
                                   q,
                                   iterations=1000,
                                   plot.freq=100,
                                   classes=NULL,
                                   num.init.params=100,
                                   K.bc.l=0.1,
                                   Z.normal.prior=TRUE) {
  current.q <- 1
  while (current.q < q) {
    current.q <- current.q + 1

    if (current.q == 2) {
      A.init <- matrix(rnorm(nrow(X) * 2), ncol=2)
    } else {
      A.init <- cbind(temp.fit$A, as.matrix(rnorm(nrow(X))))
    }

    temp.fit <- fit.bcgplvm(X, q=current.q,
                            iterations=iterations,
                            plot.freq=plot.freq,
                            classes=classes,
                            Z.init=NULL,
                            A.init=A.init,
                            num.init.params=num.init.params,
                            K.bc.l=K.bc.l,
                            Z.normal.prior=Z.normal.prior)
  }

  return(temp.fit)
}

#' Select a lengthscale for constrained optimization by min and max
#'
#' @param X
#' @param target.max.min
#' @param target.min.max
#'
#' @return
#' @export
select.bc.l.minmax <- function(X, target.max.min, target.min.max) {
  cost.fun <- function(l, X.dist, target.max.min, target.min.max) {
    K.bc <- gplvm.SE.dist(X.dist, l, 1)
    max.min <- max(apply(K.bc, 1, min))
    min.max <- min(apply(K.bc - diag(1, nrow(K.bc)), 1, max))
    max(max.min - target.max.min, 0)^2 + max(target.min.max - min.max, 0)^2
  }
  X.dist <- as.matrix(dist(X))
  optimize(cost.fun, c(0, max(X.dist)*10),
           X.dist=X.dist, target.max.min=target.max.min, target.min.max=target.min.max)$minimum
}

#' Select a lengthscale for constrained optimization by specifying the n'th
#' centile of the off-diagonal entries of the correlation matrix
#'
#' @param params a vector with named entries "centile" and "target"
#' @param X
#'
#' @return
#' @export
select.bc.l.centile <- function(X, params=NULL) {
  if (is.null(params)) {
    params <- c(centile=0.9, target=0.5)
  }
  if (any(!(c("centile", "target") %in% names(params)))) {
    stop("params does not contain named entries for 'centile' and 'target'")
  }
  X.dist <- as.matrix(dist(X))
  X.dist.centile <- quantile(X.dist[X.dist!=0], 1-params["centile"])
  l <- as.numeric(X.dist.centile / sqrt(-2*log(params["target"])))
  return(l)
}

#' Select a lengthscale for constrained optimization by median
#'
#' @param X
#' @param target.median.cor
#'
#' @return
#' @export
select.bc.l.median <- function(X, target.median.cor) {
  return(select.bc.l.centile(X, c(centile=0.5, target=target.median.cor)))
}

bc.l.selection.plots <- function(X, steps=1024, llim=NULL, chosen.lengthscale=NULL) {
  X.dist <- as.matrix(dist(X))
  X.dist.vec <- X.dist
  diag(X.dist.vec) <- NA
  X.dist.vec <- X.dist.vec[!is.na(X.dist.vec)]
  plot(density(X.dist.vec), main="Density of dist(X)")
  dispMatrix <- matrix(0, steps, steps)
  ecdfMatrix <- matrix(0, steps, steps)
  if (is.null(llim)) {
    llim <- c(min(X.dist.vec), max(X.dist.vec))
  }
  #rescale chosen.lengthscale to w.r.t. llim sit in [0,1]
  if(!is.null(chosen.lengthscale)) {
    chosen.lengthscale.scaled <- chosen.lengthscale - llim[1] / (llim[2] - llim[1])
  }
  lValues <- seq(llim[1], llim[2], length.out = steps)
  corValues <- seq(0, 1, length.out = steps)
  rownames(dispMatrix) <- lValues
  colnames(dispMatrix) <- corValues
  for (i in 1:steps) {
    tempCor <- GPLVM:::gplvm.SE.dist(X.dist, lValues[i], 1)
    diag(tempCor) <- NA
    tempCor <- tempCor[!is.na(tempCor)]
    tempDens <- density(tempCor, n=steps, from=0, to=1, bw=10^-2)
    #plot(tempDens)
    #temphist <- hist(tempCor, plot=F, breaks=corValues)$counts
    #temphist <- (temphist - min(temphist)) / (max(temphist) - min(temphist))
    dispMatrix[i, ] <- tempDens$y#(tempDens$y - min(tempDens$y))/(max(tempDens$y) - min(tempDens$y)) #temphist
    tempECDF <- ecdf(tempCor)
    ecdfMatrix[i, ] <- tempECDF(corValues)
  }
  image(dispMatrix, col=viridis::viridis(256), useRaster=T, axes=F, zlim=c(0, quantile(dispMatrix, 0.95)), main="Correlation Matrix Density Smear", xlab="Lengthscale", ylab="Correlation")
  ticks <- round(seq(1, steps, length.out=10))
  axis(1, at=corValues[ticks], labels = signif(lValues[ticks], 2))
  axis(2)
  if(!is.null(chosen.lengthscale)) {
    abline(v=chosen.lengthscale.scaled, lty=2)
  }

  image(ecdfMatrix, col=viridis::viridis(256), useRaster=T, axes=F,main="eCDF smear", xlab="Lengthscale", ylab="Correlation")
  ticks <- round(seq(1, steps, length.out=10))
  axis(1, at=corValues[ticks], labels = signif(lValues[ticks], 2))
  axis(2)
  if(!is.null(chosen.lengthscale)) {
    abline(v=chosen.lengthscale.scaled, lty=2)
  }


  contour(ecdfMatrix, axes=F, add=T)


}
