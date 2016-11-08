
dL.dK <- function(X, K) {
  K.chol <- chol(K)
  K.inv <- chol2inv(K.chol)
  temp <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  D <- ncol(X)
  (temp %*% t(temp) - D * K.inv) / 2
}

dL.dZ <- function(X, Z, l, alpha, sigma, K, dL.dK, Z.normal.prior) {
  out <- Z
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dK * dK.dZij(Z, K, i, j, l))
    }
  }
  if (Z.normal.prior) {
    out <- out - Z / 10^2
  }
  return(out)
}

dL.dl <- function(X, Z, l, alpha, sigma, K, dL.dK) {
  dK.dl <- dK.dl(Z, K, l)
  out <- sum(dL.dK * dK.dl)
  return(out)
}

dL.dalpha <- function(X, Z, l, alpha, sigma, K, dL.dK) {
  dK.dalpha <- dK.dalpha(Z, K, alpha, sigma)
  out <- sum(dL.dK * dK.dalpha)
  return(out)
}

dL.dsigma <- function(X, Z, l, alpha, sigma, K, dL.dK) {
  dK.dsigma <- diag(2 * sigma, nrow=nrow(Z))
  out <- sum(dL.dK * dK.dsigma)
  return(out)
}


gplvm.L.from.K <- function(Z, X, K, Z.normal.prior=TRUE) {
  X <- as.matrix(X)
  D <- ncol(X)
  N <- nrow(X)
  K.chol <- chol(K)
  K.inv.X <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  K.log.det <- 2 * sum(log(diag(K.chol)))
  Z.prior.term <- 0
  if (Z.normal.prior) {
    Z.prior.term <- - length(Z)/2 * log(2 * 10^2 * pi) - 1 / (2 * 10 ^2) * sum(as.numeric(Z)^2)
  }
  return(-D * N * log(2 * pi) / 2 - D * K.log.det / 2 - 1/2 * sum(K.inv.X * X) + Z.prior.term)
}

gplvm.L <- function(Z, X, l, alpha, sigma, Z.normal.prior=TRUE) {
  K <- gplvm.SE(Z, l, alpha, sigma)
  return(gplvm.L.from.K(Z, X, K, Z.normal.prior))
}

gplvm.f <- function(par, X, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  Z <- matrix(par[-(1:3)], nrow=nrow(X))
  gplvm.L(Z, X, l, alpha, sigma, Z.normal.prior)
}

gplvm.gr <- function(par, X, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  Z <- matrix(par[-(1:3)], nrow=nrow(X))
  K <- gplvm.SE(Z, l, alpha, sigma)
  dL.dK <- dL.dK(X, K)
  c(dL.dl(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dalpha(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dsigma(X, Z, l, alpha, sigma, K, dL.dK),
    as.numeric(dL.dZ(X, Z, l, alpha, sigma, K, dL.dK, Z.normal.prior))
  )
}

gplvm.hp.f <- function(par, Z, X, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  gplvm.L(Z, X, l, alpha, sigma, Z.normal.prior=Z.normal.prior)
}

gplvm.hp.gr <- function(par, Z, X, Z.normal.prior) {
  l <- par[1]
  alpha <- par[2]
  sigma <- par[3]
  K <- gplvm.SE(Z, l, alpha, sigma)
  dL.dK <- dL.dK(X, K)
  c(dL.dl(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dalpha(X, Z, l, alpha, sigma, K, dL.dK),
    dL.dsigma(X, Z, l, alpha, sigma, K, dL.dK)
  )
}


#' Fit an unconstrained GPLVM model
#'
#' @param X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param Z.init
#' @param num.init.params
#' @param Z.normal.prior
#'
#' @return
#' @export
#'
#' @importFrom optimx optimx
fit.gplvm <- function(X,
                      q,
                      iterations=1000,
                      plot.freq=100,
                      classes=NULL,
                      Z.init=NULL,
                      num.init.params=100,
                      Z.normal.prior) {
  if (is.null(Z.init)) {
    Z <- matrix(rnorm(nrow(X)*q, sd=0.2), ncol=q)
  } else if (identical(Z.init, "PCA")) {
    X.pca <- prcomp(X)
    Z <- X.pca$x[,1:q]
  } else {
    if (ncol(as.matrix(Z.init)) != q) warning("Mismatch between Z.init and q")
    Z <- as.matrix(Z.init)[,1:q]
  }

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


  par.init <- c(1, 1, 1)
  par.init.L <- gplvm.hp.f(par.init, Z=Z, X=X)
  #first optimize the hyperparams for the initial Z
  for (i in seq(num.init.params)) {
    test.par.init <- rgamma(3, 1)
    test.par.init.L <- gplvm.hp.f(test.par.init, Z=Z, X=X)
    if (test.par.init.L > par.init.L) {
      par.init <- test.par.init
      par.init.L <- test.par.init.L
    }
  }


  optout <- optimx(par.init,
                   gplvm.hp.f,
                   gplvm.hp.gr,
                   Z=Z,
                   X=X,
                   method="L-BFGS-B",
                   control=list(trace=T,
                                maximize=T,
                                kkt=FALSE,
                                maxit=1000,
                                starttests=TRUE)
  )

  init.l <- as.numeric(optout[1])
  init.alpha <- as.numeric(optout[2])
  init.sigma <- as.numeric(optout[3])

  cat("Starting params: l=", init.l, "; alpha=", init.alpha, "; sigma=", init.sigma, fill=TRUE)

  par <- c(init.l, init.alpha, init.sigma, as.numeric(Z))

  its <- 0
  convcode=1
  while (its < iterations & convcode != 0) {
    out <- optimx(par,
                  gplvm.f,
                  gplvm.gr,
                  X=X,
                  Z.normal.prior=Z.normal.prior,
                  method="L-BFGS-B",
                  control=list(trace=T,
                               maximize=T,
                               kkt=FALSE,
                               maxit=plot.freq,
                               starttests=FALSE)
    )
    convcode <- out$convcode
    par <- head(as.numeric(out), length(Z) + 3)
    l <- par[1]
    alpha <- par[2]
    sigma <- par[3]
    Z <- matrix(par[-(1:3)], ncol=q)
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

  return(list(Z=Z, l=l, alpha=alpha, sigma=sigma, convcode=convcode))

}
