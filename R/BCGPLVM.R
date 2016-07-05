
bc.Z <- function(K.bc, A) {
  Z <- K.bc %*% A
  return(Z)
}

bc.L <- function(A, X, l, alpha, sigma, K.bc, Z.normal.prior=FALSE) {
  Z <- bc.Z(K.bc, A)
  gplvm.L(Z, X, l, alpha, sigma, Z.normal.prior=Z.normal.prior)
}


dZ.dAij <- function(Z, K.bc, i, j) {
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  out[, j] <- K.bc[, i]
  return(out)
}

dL.dA <- function(X, Z, l, alpha, sigma, K.bc, K, dL.dK, Z.normal.prior=FALSE) {
  dL.dZ <- dL.dZ(X, Z, l, alpha, sigma, K, dL.dK, Z.normal.prior=Z.normal.prior)
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dZ * dZ.dAij(Z, K.bc, i, j))
    }
  }
  return(out)
}

dL.dA.check <- function(X, A, l, alpha, sigma, K.bc, Z.normal.prior=FALSE) {
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


bcgplvm.f <- function(par, X, K.bc, Z.normal.prior=FALSE) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  A <- matrix(par[-(1:3)], nrow=nrow(X))
  bc.L(A, X, l, alpha, sigma, K.bc, Z.normal.prior=Z.normal.prior)
}

bcgplvm.gr <- function(par, X, K.bc, Z.normal.prior=FALSE) {
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
  par <- rnorm(length(target))
  optimx.out <- optimx(par,
                       optimx.S,
                       optimx.S.grad,
                       target=target,
                       K.bc=K.bc,
                       method="BFGS",
                       control=list(trace=T,
                                    maximize=F,
                                    kkt=FALSE,
                                    maxit=1000,
                                    starttests=TRUE)
  )
  A <- matrix(as.numeric(optimx.out)[1:length(target)], nrow=nrow(target), ncol=ncol(target))
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
                        classes=NULL,
                        Z.init=NULL,
                        A.init=NULL,
                        num.init.params=100,
                        K.bc.l=0.1,
                        Z.normal.prior=FALSE) {

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
                   method="BFGS",
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
                  method="BFGS",
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
                                   Z.normal.prior=FALSE) {
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
