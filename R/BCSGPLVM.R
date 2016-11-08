
bcs.L <- function(X, nrows,
                  se.l, se.alpha, se.sigma,
                  structured.C, structured.alpha,
                  A, Z.normal.prior,
                  K.bc) {
  Z <- bc.Z(K.bc, A)
  structured.likelihood(X, nrows,
                        se.l, se.alpha, se.sigma,
                        structured.C, structured.alpha,
                        Z, Z.normal.prior)
}

bcs.L.grad <- function(X, nrows,
                       se.l, se.alpha, se.sigma,
                       structured.C, structured.alpha,
                       A,
                       probabilistic.trace.estimate,
                       Z.normal.prior,
                       K.bc) {
  Z <- bc.Z(K.bc, A)

  out <- structured.likelihood.grad(X, nrows,
                                    se.l, se.alpha, se.sigma,
                                    structured.C, structured.alpha,
                                    Z,
                                    probabilistic.trace.estimate,
                                    Z.normal.prior)

  dL.dZ <- matrix(out[-(1:6)], nrow=nrow(Z), ncol=ncol(Z))

  dL.dA <- matrix(0, nrow=nrow(A), ncol=ncol(A))
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      dL.dA[i, j] <- sum(dL.dZ * dZ.dAij(Z, K.bc, i, j))
    }
  }

  out[-(1:6)] <- as.numeric(dL.dA)

  return(out)
}

bcs.L.optimx <- function(par,
                         X, nrows.X,
                         probabilistic.trace.estimate,
                         Z.normal.prior,
                         K.bc) {
  bcs.L(X, nrows.X,
        se.l=par[1], se.alpha=par[2], se.sigma=par[3],
        structured.C=par[4:5], structured.alpha=par[6],
        A=matrix(par[-(1:6)], nrow=nrow(X)),
        Z.normal.prior=Z.normal.prior,
        K.bc=K.bc)
}


bcs.L.grad.optimx <- function(par,
                              X, nrows.X,
                              probabilistic.trace.estimate,
                              Z.normal.prior,
                              K.bc) {
  bcs.L.grad(X, nrows.X,
             se.l=par[1], se.alpha=par[2], se.sigma=par[3],
             structured.C=par[4:5], structured.alpha=par[6],
             A=matrix(par[-(1:6)], nrow=nrow(X)),
             probabilistic.trace.estimate=probabilistic.trace.estimate,
             Z.normal.prior=Z.normal.prior,
             K.bc=K.bc)
}


## Fixed A Optimx functions

bcs.L.fixedA.optimx <- function(par,
                                X, nrows.X,
                                A,
                                probabilistic.trace.estimate,
                                Z.normal.prior,
                                K.bc) {
  par <- c(par, as.numeric(A))
  bcs.L.optimx(par,
               X, nrows.X,
               probabilistic.trace.estimate,
               Z.normal.prior,
               K.bc)
}

bcs.L.fixedA.grad.optimx <- function(par,
                                     X, nrows.X,
                                     A,
                                     probabilistic.trace.estimate,
                                     Z.normal.prior,
                                     K.bc) {
  par <- c(par, as.numeric(A))
  bcs.L.grad.optimx(par,
                    X, nrows.X,
                    probabilistic.trace.estimate,
                    Z.normal.prior,
                    K.bc)[1:6]
}


## Fixed C Optimx functions

bcs.L.fixedC.optimx <- function(par,
                                X, nrows.X,
                                structured.C,
                                probabilistic.trace.estimate,
                                Z.normal.prior,
                                K.bc) {
  new.par <- numeric(length(par) + 2)
  new.par[1:3] <- par[1:3]
  new.par[4:5] <- structured.C
  new.par[6] <- par[4]
  new.par[-(1:6)] <- par[-(1:4)]
  bcs.L.optimx(new.par,
               X, nrows.X,
               probabilistic.trace.estimate,
               Z.normal.prior,
               K.bc)
}

bcs.L.fixedC.grad.optimx <- function(par,
                                     X, nrows.X,
                                     structured.C,
                                     probabilistic.trace.estimate,
                                     Z.normal.prior,
                                     K.bc) {
  new.par <- numeric(length(par) + 2)
  new.par[1:3] <- par[1:3]
  new.par[4:5] <- structured.C
  new.par[6] <- par[4]
  new.par[-(1:6)] <- par[-(1:4)]
  bcs.L.grad.optimx(new.par,
                    X, nrows.X,
                    probabilistic.trace.estimate,
                    Z.normal.prior,
                    K.bc)[-(4:5)]
}

## Fixed A and C Optimx functions

bcs.L.fixedAC.optimx <- function(par,
                                 X, nrows.X,
                                 A,
                                 structured.C,
                                 probabilistic.trace.estimate,
                                 Z.normal.prior,
                                 K.bc) {
  new.par <- numeric(length(par) + 2 + length(A))
  new.par[1:3] <- par[1:3]
  new.par[4:5] <- structured.C
  new.par[6] <- par[4]
  new.par[-(1:6)] <- par[-(1:4)]
  bcs.L.optimx(new.par,
               X, nrows.X,
               probabilistic.trace.estimate,
               Z.normal.prior,
               K.bc)
}

bcs.L.fixedAC.grad.optimx <- function(par,
                                      X, nrows.X,
                                      A,
                                      structured.C,
                                      probabilistic.trace.estimate,
                                      Z.normal.prior,
                                      K.bc) {
  new.par <- numeric(length(par) + 2 + length(A))
  new.par[1:3] <- par[1:3]
  new.par[4:5] <- structured.C
  new.par[6] <- par[4]
  new.par[-(1:6)] <- as.numeric(A)
  bcs.L.grad.optimx(new.par,
                    X, nrows.X,
                    probabilistic.trace.estimate,
                    Z.normal.prior,
                    K.bc)[c(1:3, 6)]
}


#' Fit a backconstrained structured GPLVM model
#'
#' @param X
#' @param nrows.X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param Z.init
#' @param A.init
#' @param K.bc.l
#' @param Z.normal.prior
#' @param structured.C.init
#' @param structured.alpha.init
#' @param se.par.init
#' @param fixed.C
#' @param structured.C.bounds
#' @param noise.lower.bound
#' @param probabilistic.trace.estimate
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
#' Z.init <- prcomp(X.noisy)$x[,1:3]
#' bcsgplvm.model <- fit.bcsgplvm(X=X.noisy, nrows.X=20,
#'                                q=3, iterations=1000, plot.freq=100,
#'                                classes=classes, Z.init="PCA", K.bc.l=50,
#'                                probabilistic.trace.estimate=FALSE)
fit.bcsgplvm <- function(X,
                         nrows.X,
                         q=2,
                         iterations=1000,
                         plot.freq=100,
                         classes=1,
                         Z.init=NULL,
                         A.init=NULL,
                         K.bc.l=0.1,
                         Z.normal.prior=TRUE,
                         structured.C.init=c(2, 2),
                         structured.alpha.init=1,
                         se.par.init=NULL,
                         fixed.C=FALSE,
                         structured.C.bounds=c(1, 1, 5, 5),
                         noise.lower.bound=10^-4,
                         probabilistic.trace.estimate=TRUE) {
  if (!(is.null(Z.init) | is.null(A.init))) {
    stop("Can't specify both Z.init and A.init")
  }

  K.bc <- gplvm.SE(Z=X, l=K.bc.l, alpha=1, sigma=0)


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

  if (is.null(se.par.init)) {
    se.par.init <- c(median(Z.dist), alpha.shape / (structured.alpha.init + 1), 1)
  }


  if (fixed.C) {
    par.init <- c(se.par.init, structured.alpha.init)
    par.lower <- rep(-Inf, length(par.init))
    par.lower[3] <- noise.lower.bound
    par.upper <- rep(Inf, length(par.init))
    optout <- optimx(par.init,
                     bcs.L.fixedAC.optimx,
                     bcs.L.fixedAC.grad.optimx,
                     X=X, nrows.X=nrows.X,
                     A=A.init,
                     structured.C=structured.C.init,
                     probabilistic.trace.estimate=probabilistic.trace.estimate,
                     Z.normal.prior=Z.normal.prior,
                     K.bc=K.bc,
                     method="L-BFGS-B",
                     lower=par.lower, upper=par.upper,
                     control=list(trace=T,
                                  maximize=T,
                                  kkt=FALSE,
                                  maxit=iterations,
                                  starttests=FALSE)
    )
    par <- numeric(4 + length(A.init))
    par[1:4] <- as.numeric(optout)[1:4]
    par[-(1:4)] <- as.numeric(A.init)
  } else {
    par.init <- c(se.par.init, structured.C.init, structured.alpha.init)
    par.lower <- rep(-Inf, length(par.init))
    par.lower[4:5] <- structured.C.bounds[1:2]
    par.lower[3] <- noise.lower.bound
    par.upper <- rep(Inf, length(par.init))
    par.upper[4:5] <- structured.C.bounds[3:4]

    optout <- optimx(par.init,
                     bcs.L.fixedA.optimx,
                     bcs.L.fixedA.grad.optimx,
                     X=X, nrows.X=nrows.X,
                     A=A.init,
                     probabilistic.trace.estimate=probabilistic.trace.estimate,
                     Z.normal.prior=Z.normal.prior,
                     K.bc=K.bc,
                     method="L-BFGS-B",
                     lower=par.lower, upper=par.upper,
                     control=list(trace=T,
                                  maximize=T,
                                  kkt=FALSE,
                                  maxit=iterations,
                                  starttests=FALSE)
    )
    par <- numeric(6 + length(A.init))
    par[1:6] <- as.numeric(optout)[1:6]
    par[-(1:6)] <- as.numeric(A.init)
  }

  its <- 0
  convcode=1
  while (its < iterations & convcode != 0) {

    par.lower <- c(par.lower, rep(-Inf, length(A.init)))
    par.upper <- c(par.upper, rep(Inf, length(A.init)))

    if (fixed.C) {
      out <- optimx(par,
                    bcs.L.fixedC.optimx,
                    bcs.L.fixedC.grad.optimx,
                    X=X, nrows.X=nrows.X,
                    structured.C=structured.C.init,
                    probabilistic.trace.estimate=probabilistic.trace.estimate,
                    Z.normal.prior=Z.normal.prior,
                    K.bc=K.bc,
                    method="L-BFGS-B",
                    lower=par.lower, upper=par.upper,
                    control=list(trace=T,
                                 maximize=T,
                                 kkt=FALSE,
                                 maxit=plot.freq,
                                 starttests=FALSE)
      )
    } else {
      out <- optimx(par,
                    bcs.L.optimx,
                    bcs.L.grad.optimx,
                    X=X, nrows.X=nrows.X,
                    probabilistic.trace.estimate=probabilistic.trace.estimate,
                    Z.normal.prior=Z.normal.prior,
                    K.bc=K.bc,
                    method="L-BFGS-B",
                    lower=par.lower, upper=par.upper,
                    control=list(trace=T,
                                 maximize=T,
                                 kkt=FALSE,
                                 maxit=plot.freq,
                                 starttests=FALSE)
      )
    }
    convcode <- out$convcode
    par <- head(as.numeric(out), length(par))
    se.l <- par[1]
    se.alpha <- par[2]
    se.sigma <- par[3]

    if (fixed.C) {
      structured.C <- structured.C.init
      structured.alpha <- par[4]
      kern.par.inds <- 1:4
    } else {
      structured.C <- par[4:5]
      structured.alpha <- par[6]
      kern.par.inds <- 1:6
    }

    A <- matrix(par[-kern.par.inds], ncol=q)
    Z <- bc.Z(K.bc=K.bc, A=A)
    title.text <- paste("l=", signif(se.l, 3),
                        "; se.a=", signif(se.alpha, 3),
                        "; sigma=", signif(se.sigma, 3),
                        "; str.C=", paste(signif(structured.C, 3), sep=","),
                        "; str.a=", signif(structured.alpha, 3),
                        sep="")
    if (q > 2) {
      pairs(Z,
            col=classes,
            main=title.text)
    } else if (q == 2) {
      plot(Z,
           col=classes,
           main=title.text)
    } else {
      plot(as.numeric(Z),
           rep(0, length(Z)),
           col=classes,
           main=title.text)
    }

    its <- its + plot.freq
  }

  return(list(Z=Z, A=A,
              se.l=se.l, se.alpha=se.alpha, se.sigma=se.sigma,
              structured.C=structured.C, structured.alpha=structured.alpha,
              convcode=convcode))
}
