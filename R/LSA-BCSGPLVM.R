ZPriorUnif_ <- "uniform"
ZPriorNorm_ <- "normal"
ZPriorDisc_ <- "discriminative"
ZPriorOptions_ <- c(ZPriorNorm_, ZPriorUnif_, ZPriorDisc_)


arr_ind_conv_multipliers <- function(limits) {
  dim <- length(limits)
  multipliers <- numeric(dim) + 1
  if (dim > 1) {
    for (i in 2:dim) {
      multipliers[i:dim] <- multipliers[i:dim] * limits[i-1]
    }
  }
  return(multipliers)
}

arr_to_vec_index <- function(index, limits) {
  if (any(index > limits | index < 1)) {
    stop(paste("Index is outside limits of array:", index))
  }
  if (length(limits) != length(index)) {
    stop("Index dimension does not match limits dimension.")
  }
  multipliers <- arr_ind_conv_multipliers(limits)
  out <- sum((index - 1) * multipliers) + 1
  return(out)
}

vec_to_arr_index <- function(index, limits) {
  if (index > prod(limits) | index < 1) {
    stop(paste("Index is outside limits of array:", index))
  }
  multipliers <- arr_ind_conv_multipliers(limits)
  dim <- length(limits)
  out <- numeric(dim)
  if (dim > 1) {
    for (i in 1:(dim - 1)) {
      out[i] <- ((index - 1) %% multipliers[i + 1]) / multipliers[i] + 1
      index <- index - (out[i] - 1) * multipliers[i]
    }
  }
  out[dim] <- (index - 1) / multipliers[dim] + 1
  return(out)
}

LSA_BCSGPLVM.kernel <- function(W, l, alpha, sigma) {
  structure.dim <- length(l) - 1
  z.dim <- ncol(W) - structure.dim
  l <- c(rep(l[1], z.dim), l[-1])
  W <- t(t(W) / l)
  W.dist <- as.matrix(dist(W))
  K <- alpha^2 * exp(-W.dist^2 / 2)
  diag(K) <- diag(K) + sigma^2
  return(K)
}

LSA_BCSGPLVM.Z.from.W <- function(W, l) {
  Z.ncol <- ncol(W) - length(l) + 1
  Z <- W[, 1:Z.ncol]
  return(Z)
}

class_var_matrices <- function(Z, classes, grad=FALSE) {
  # There is a typo in the Discriminative GPLVM paper (Urtasun 2007) (S_w and
  # S_b are switched) so the definitions of S_w and S_b used here are from
  # Sugiyama 2006
  out <- list()
  if (!is.factor(classes)) stop("classes parameter must be a factor")
  if (!is.numeric(Z)) {
    stop("Z must be numeric")
  } else {
    Z <- as.matrix(Z)
  }
  M_0 <- colMeans(Z)
  nclasses <- nlevels(classes)
  M <- matrix(0, nrow=nclasses, ncol=ncol(Z))
  S_w <- S_b <- matrix(0, ncol(Z), ncol(Z))
  if (grad) {
    dS_w.dZ <- dS_b.dZ <- array(0, dim=c(dim(Z), dim(S_w)))
  }
  for (i in 1:nclasses) {
    index <- which(classes==levels(classes)[i])
    M[i, ] <- colMeans(Z[index, ])
    temp.sum <- M[i, ] - M_0
    S_b <- S_b + length(index) * outer(temp.sum, temp.sum)
    Z.centered.class <- scale(Z[index, ], center=M[i, ], scale=F)
    S_w <- S_w + t(Z.centered.class) %*% Z.centered.class
    if (grad) {
      for (l in index) {
        for (m in 1:ncol(Z)) {
          dS_b.dZ[l, m, m, ] <- dS_b.dZ[l, m, , m] <- temp.sum
          dS_b.dZ[l, m, m, m] <- 2 * dS_b.dZ[l, m, m, m]

          dS_w.dZ[l, m, m, ] <- dS_w.dZ[l, m, , m] <- Z[l, ] - M[i, ]
          dS_w.dZ[l, m, m, m] <- 2 * dS_w.dZ[l, m, m, m]
        }
      }
    }
  }
  out$S_b <- S_b
  out$S_w <- S_w

  if (grad) {
    out$dS_b.dZ <- dS_b.dZ
    out$dS_w.dZ <- dS_w.dZ
  }

  return(out)
}

discriminative.prior <- function(Z, Z.prior.params, grad=FALSE) {
  classes <- Z.prior.params$classes
  sigma_d <- Z.prior.params$sigma_d

  # Check for missing parameters
  if (is.null(classes) | is.null(sigma_d)) stop("Z.prior.params incorrectly specified for discriminative prior")

  class.var <- class_var_matrices(Z = Z, classes = classes, grad=grad)

  S_w.inv.S_b <- solve(class.var$S_w, class.var$S_b)

  J <- sum(diag(S_w.inv.S_b))

  out <- list()
  out$J <- J
  out$L <- -1 / (sigma_d^2 * J)
  if (grad) {
    dJ.dZ <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))

    dS_b.dZ <- class.var$dS_b.dZ
    dS_w.dZ <- class.var$dS_w.dZ

    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        S_w.inv.dS_b <- solve(class.var$S_w, dS_b.dZ[i, j, ,])
        S_w.inv.dS_w <- solve(class.var$S_w, dS_w.dZ[i, j, ,])
        dJ.dZ[i, j] <- sum(diag(S_w.inv.dS_b)) - sum(S_w.inv.dS_w * t(S_w.inv.S_b))
      }
    }
    out$dJ.dZ <- dJ.dZ
    out$dL.dZ <- 1/(sigma_d^2 * J^2) * dJ.dZ
  }
  return(out)
}

LSA_BCSGPLVM.L <- function(W, X, l, alpha, sigma, Z.prior=c("normal", "uniform", "discriminative"), K=NULL, Z.prior.params=list()) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  if (missing(K)) {
    K <- LSA_BCSGPLVM.kernel(W=W, l=l, alpha=alpha, sigma=sigma)
  }

  X <- as.matrix(X)
  N <- nrow(X)
  K.chol <- chol(K)
  K.inv.X <- backsolve(K.chol, forwardsolve(t(K.chol), X))
  K.log.det <- 2 * sum(log(diag(K.chol)))

  Z.prior.term <- 0
  Z <- LSA_BCSGPLVM.Z.from.W(W, l)
  if (Z.prior == ZPriorNorm_) {
    Z.prior.term <- -length(Z)/2 * log(2 * 10^2 * pi) - 1 / (2 * 10 ^2) * sum(as.numeric(Z)^2)
  } else if (Z.prior == ZPriorDisc_) {
    Z.prior.term <- discriminative.prior(Z, Z.prior.params)$L
  }

  return(N * log(2 * pi) / 2 - K.log.det / 2 - 1/2 * sum(K.inv.X * X) + Z.prior.term)
}

LSA_BCSGPLVM.dK.dL_Z <- function(W, l, K) {
  Z <- LSA_BCSGPLVM.Z.from.W(W, l)
  Z.dist <- as.matrix(dist(Z)^2)
  out <- Z.dist / l[1]^3 * K
  return(out)
}

LSA_BCSGPLVM.dK.dL_S_i <- function(W, l, K, i) {
  Z.ncol <- ncol(W) - length(l) + 1
  S_i <- W[, i + Z.ncol]
  S_i.dist <- as.matrix(dist(S_i)^2)
  out <- S_i.dist / l[i + 1]^3 * K
  return(out)
}

LSA_BCSGPLVM.dL.dpar <- function(W, X, l, alpha, sigma, K=NULL, dL.dK=NULL) {
  X <- as.matrix(X)
  if (missing(K)) {
    K <- LSA_BCSGPLVM.kernel(W=W, l=l, alpha=alpha, sigma=sigma)
  }
  if (missing(dL.dK)) {
    dL.dK <- dL.dK(X, K)
  }
  dL.dpar <- numeric(length(l) + 2)
  #alpha
  dK.dalpha <- dK.dalpha(Z=W, K=K, alpha=alpha, sigma=sigma)
  dL.dpar[1] <- sum(dL.dK * dK.dalpha)

  #sigma
  dL.dpar[2] <- sum(diag(dL.dK)) * 2 * sigma

  #l_Z
  dK.dl_Z <-LSA_BCSGPLVM.dK.dL_Z(W, l, K)
  dL.dpar[3] <- sum(dL.dK * dK.dl_Z)

  num.structural.dim <- length(l) - 1

  #l_S_i
  for (i in 1:num.structural.dim) {
    dK.dl_S_i <- LSA_BCSGPLVM.dK.dL_S_i(W, l, K, i)
    dL.dpar[3 + i] <- sum(dL.dK * dK.dl_S_i)
  }

  return(dL.dpar)
}

LSA_BCSGPLVM.dK.dZij <- function(Z, K, i, j, l_Z, W.pre) {
  out <- matrix(0, nrow=nrow(K), ncol=ncol(K))
  index <- which(W.pre[, 1] == i)
  if (length(index) != 0) {
    out[, index] <- (Z[W.pre[, 1], j] - Z[i, j]) / l_Z^2 * K[, index]
    out[index, ] <- t(out[, index])
  }
  return(out)
}

LSA_BCSGPLVM.dL.dZ <- function(W, W.pre, X, Z,
                               l, alpha, sigma,
                               Z.prior=c("normal", "uniform", "discriminative"),
                               K=NULL, dL.dK=NULL,
                               Z.prior.params=list()) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  X <- as.matrix(X)
  if (missing(K)) {
    K <- LSA_BCSGPLVM.kernel(W=W, l=l, alpha=alpha, sigma=sigma)
  }
  if (missing(dL.dK)) {
    dL.dK <- dL.dK(X, K)
  }
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      index <- which(W.pre[, 1] == i)
      dK.dZij <-  LSA_BCSGPLVM.dK.dZij(Z, K, i, j, l[1], W.pre)
      out[i, j] <- (2 * sum(dL.dK[index, ] * dK.dZij[index, ]) -
        sum(dL.dK[index, index] * dK.dZij[index, index]))
    }
  }
  if (Z.prior==ZPriorNorm_) {
    out <- out - Z / 10^2
  } else if (Z.prior==ZPriorDisc_) {
    out <- out + discriminative.prior(Z, Z.prior.params, grad=TRUE)$dL.dZ
  }
}

LSA_BCSGPLVM.dL.dA <- function(W, W.pre, X, Z,
                               l, alpha, sigma,
                               K.bc,
                               Z.prior=c("normal", "uniform", "discriminative"),
                               K=NULL, dL.dK=NULL,
                               Z.prior.params=list()) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  X <- as.matrix(X)
  if (missing(K)) {
    K <- LSA_BCSGPLVM.kernel(W=W, l=l, alpha=alpha, sigma=sigma)
  }
  if (missing(dL.dK)) {
    dL.dK <- dL.dK(X, K)
  }
  dL.dZ <- LSA_BCSGPLVM.dL.dZ(W=W, W.pre=W.pre, X=X, Z=Z,
                              l=l, alpha=alpha, sigma=sigma,
                              Z.prior=Z.prior, K=K, dL.dK=dL.dK,
                              Z.prior.params=Z.prior.params)
  out <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
  for (i in 1:nrow(Z)) {
    for (j in 1:ncol(Z)) {
      out[i, j] <- sum(dL.dZ[, j] * dZ.dAij(Z, K.bc, i, j)[, j])
    }
  }
  return(out)
}


array_to_flat_matrix <- function(data.array) {
  array.dim <- dim(data.array)
  out <- matrix(NA, nrow=prod(array.dim), ncol=length(array.dim) + 1)
  each <- 1
  for (i in 1:length(array.dim)) {
    out[,i] <- rep(1:array.dim[i], length.out=prod(array.dim), each=each)
    each <- each * array.dim[i]
  }
  out[, length(array.dim) + 1] <- as.numeric(data.array)
  if (!is.null(dimnames(data.array))) {
    colnames(out) <- c(dimnames(data.array), "data")
  }
  return(out)
}

LSA_BCSGPLVM.plot_iteration <- function(Z, par.hist, plot.A, iteration, classes) {
  close.screen(all.screens=TRUE)
  num.pars <- ncol(par.hist)
  if (plot.A) {
    split.screen(c(1, 2))
    par.screens <- split.screen(c(num.pars + 1, 1), 1)
    q <- ncol(Z)
    if (q > 2) {
      plots <-q * (q - 1) / 2
      ncols <- floor(sqrt(plots))
      nrows <- ceiling(plots / ncols)
      plot.screens <- split.screen(c(nrows, ncols), 2)
    }
  } else {
    par.screens <- split.screen(c(num.pars + 1, 1))
  }
  smoothed.L <- numeric(iteration)
  smoothed.L[1] <- par.hist[1, ncol(par.hist)]
  for (i in 2:iteration) {
    smoothed.L[i] <- 19/20 * smoothed.L[i-1] + 1/20 * par.hist[i, ncol(par.hist)]
  }
  par.hist <- par.hist[1:iteration, ]
  par.hist[, ncol(par.hist)] <- smoothed.L
  for (i in 1:ncol(par.hist)) {
    screen(par.screens[i])
    par(mar=c(0,4,0,4)+0.1)
    plot(1:iteration, par.hist[1:iteration, i],
         type="l", xlim=c(0, iteration),
         ylim=range(par.hist[1:iteration, i]),
         axes=F, ann=F, col=rainbow(ncol(par.hist))[i])
    box()
    range <- range(par.hist[1:iteration, i])
    at <- mean(range)
    at <- c(at - diff(range)/3, at, at + diff(range)/3)
    labels <- format(at, digits=2)
    axis(2, at=at, labels = labels, las=1)
    at <- mean(range)
    labels <- format(par.hist[iteration, i], digits=2)
    axis(4, at=at, labels=labels, las=1, tick=FALSE)
    par(mar=c(5.1,4.1,4.1,2.1))
  }
  screen(par.screens[length(par.screens)])
  par(mar=c(0,0,0,0))
  legend("center", legend=c("alpha", "sigma", "l_Z", paste("l_S", 1:(ncol(par.hist) - 4), sep=""), "L"),
         col=rainbow(ncol(par.hist)), lty=1, ncol=4, bty="n")

  if (plot.A) {
    if (q > 2) {
      screen.num <- 0
      for (i in 1:(q-1)) {
        for (j in (i+1):q) {
          screen.num <- screen.num + 1
          screen(plot.screens[screen.num])
          par(mar=c(0,0,0,0) + 0.1)
          plot(Z[, c(i,j)], col=classes, axes=FALSE, ann=FALSE)
          box()
        }
      }
      par(mar=c(5.1,4.1,4.1,2.1))
    } else if (q == 2) {
      screen(2)
      plot(Z, col=classes)
    } else {
      screen(2)
      plot(as.numeric(Z), rep(0, length(Z)), col=classes)
    }
  }
  close.screen(all.screens=TRUE)
}

LSA_BCSGPLVM.sgdopt <- function(X, A.init, par.init, K.bc, points.in.approximation,
                                iterations, step.size.range,
                                optimize.A, classes, plot.freq, learning.rate,
                                momentum.rate=momentum.rate, adam.epsilon=adam.epsilon,
                                par.step.size.range=NULL, par.fixed=NULL,
                                verbose=FALSE, Z.prior=c("normal", "uniform", "discriminative"),
                                Z.prior.params=list()) {

  ## TODO: implement Stochastic Meta-Descent


  if (!is.null(par.fixed)) {
    if (length(par.fixed) != length(par.init) | class(par.fixed) != "logical") {
      stop("par.fixed specified incorrectly")
    }
  } else {
    par.fixed <- rep(FALSE, length(par.init))
  }
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  iteration <- 0
  A <- A.init
  par <- par.init
  Z <- bc.Z(K.bc, A)

  delta.A <- matrix(0, nrow=nrow(A), ncol=ncol(A))
  delta.par <- numeric(length(par))

  # Variables for ADAM update
  m.par <- numeric(length(par))
  m.A <- matrix(0, nrow=nrow(A), ncol=ncol(A))
  v.par <- numeric(length(par))
  v.A <- matrix(0, nrow=nrow(A), ncol=ncol(A))

  par.hist <- matrix(0, ncol=length(par)+1, nrow=iterations + 1)
  if (optimize.A) {
    A.hist <- array(0, dim=c(iterations+1, dim(A)))
  }

  step.size <- step.size.range[1]
  step.size.change <- (step.size.range[1] - step.size.range[2]) / iterations

  if (is.null(par.step.size.range)) {
    par.step.size <- step.size
    par.step.size.change <- step.size.change
  } else {
    par.step.size <- par.step.size.range[1]
    par.step.size.change <- (par.step.size.range[1] - par.step.size.range[2]) / iterations
  }

  while (iteration < iterations) {
    iteration <- iteration + 1
    if (optimize.A) {
      Z <- bc.Z(K.bc, A)
    }

    # Select a subset of points to use in this step
    num.points <- prod(dim(X))
    sample.points <- sample.int(n=num.points, size=points.in.approximation, replace=FALSE)
    W.pre <- t(sapply(sample.points, vec_to_arr_index, limits=dim(X)))
    # Order it now to make computing dL.dA_ij more efficient.
    W.pre <- W.pre[order(W.pre[,1]),]

    # Replace the first column with the Z values for those samples
    W <- cbind(Z[W.pre[, 1], ], W.pre[, -1])

    X.sample <- as.matrix(X[W.pre])

    # Calculate the gradient w.r.t. the parameters
    K <- LSA_BCSGPLVM.kernel(W, par[-(1:2)], par[1], par[2])
    if (verbose) print(paste("Squared sum of K off-diag:", sum((K^2)) - sum(diag(K)^2)))
    dL.dK <- dL.dK(X.sample, K)
    L <- LSA_BCSGPLVM.L(W=W, X=X.sample, l=par[-(1:2)], alpha=par[1], sigma=par[2], Z.prior=Z.prior, K=K, Z.prior.params=Z.prior.params)

    if (optimize.A) {
      A.hist[iteration, , ] <- A
    }
    par.hist[iteration, ] <- c(par, L)

    if (iteration > 1 & plot.freq != 0 & iteration %% plot.freq == 0) {
      LSA_BCSGPLVM.plot_iteration(Z, par.hist, optimize.A, iteration, classes)
    }

    if (optimize.A) {
      dL.dA <- LSA_BCSGPLVM.dL.dA(W=W, W.pre=W.pre, X=X.sample, Z=Z,
                                  l=par[-(1:2)], alpha=par[1], sigma=par[2],
                                  K.bc=K.bc, Z.prior=Z.prior, K=K, dL.dK=dL.dK,
                                  Z.prior.params=Z.prior.params)
      m.A <- momentum.rate * m.A + (1 - momentum.rate) * dL.dA
      v.A <- learning.rate * v.A + (1 - learning.rate) * dL.dA^2

      m.A.hat <- m.A / (1 - momentum.rate^iteration)
      v.A.hat <- v.A / (1 - learning.rate^iteration)
      delta.A <- step.size / (sqrt(v.A.hat) + adam.epsilon) * m.A.hat
    }

    dL.dpar <- LSA_BCSGPLVM.dL.dpar(W, X.sample, par[-(1:2)], par[1], par[2], K=K, dL.dK=dL.dK)

    # Update the Adam variables
    m.par <- momentum.rate * m.par + (1 - momentum.rate) * dL.dpar
    v.par <- learning.rate * v.par + (1 - learning.rate) * dL.dpar^2

    m.par.hat <- m.par / (1 - momentum.rate^iteration)
    v.par.hat <- v.par / (1 - learning.rate^iteration)

    # Calculate the change (including momentum)
    delta.par <- par.step.size / (sqrt(v.par.hat) + adam.epsilon) * m.par.hat

    # Update par and step size
    par[!par.fixed] <- par[!par.fixed] + delta.par[!par.fixed]
    if (optimize.A) {
      A <- A + delta.A
    }
    if (verbose) print(par)
    par.step.size <- par.step.size - par.step.size.change
    step.size <- step.size - step.size.change
  }

  L <- LSA_BCSGPLVM.L(W=W, X=X.sample, l=par[-(1:2)], alpha=par[1], sigma=par[2], Z.prior=Z.prior, Z.prior.params=Z.prior.params)

  par.hist[iteration + 1,] <- c(par, L)

  if (optimize.A) {
    A.hist[iteration + 1, , ] <- A
  }

  out <- list(par=par, par.hist=par.hist)

  if (optimize.A) {
    out$A <- A
    out$A.hist <- A.hist
    out$Z <- bc.Z(K.bc, A)
    out$K.bc <- K.bc
    out$classes <- classes
  }

  return(out)

}


###############################################################################
#
# Obsolete. Implemented for implicit only for now:
#
# X.structure.type = "implicit" implies that X is an array and the structure
# coordinates are defined by the array index of each data point, with the first
# dimension of the array specifying the sample ID.
# X.structure.type = "explicit" implies that the structure coordinates of X are
# explicitly given in the leading columns of X, with the first column being the
# sample ID, the following columns being the structure coordinates for that
# point, and the final column being the data point value.

# Consider seperating into X.implicit (an array) and X.explicit (any additional
# non-implicit structure, with possibly multiple arrays per sample and an
# equivalence between row numbers of the two)
###############################################################################


#' Fit a large scale approximation back constrained structured GPLVM model
#'
#' @param X
#' @param q
#' @param iterations
#' @param plot.freq
#' @param classes
#' @param Z.init Either a matrix of initial latent values, or "PCA" for PCA start values, or ISOMAP for ISOMAP start values.
#' @param A.init
#' @param K.bc.l
#' @param K.bc.target.median
#' @param par.init Vector of parameters: alpha, sigma, l_Z, followed by the lengthscales for the structural dimensions
#' @param points.in.approximation
#' @param initial.step.size
#' @param final.step.size
#' @param learning.rate
#' @param Z.prior
#' @param initial.step.size.par
#' @param final.step.size.par
#' @param momentum.rate
#' @param adam.epsilon
#' @param parameter.opt.iterations
#' @param par.fixed.A.opt
#' @param verbose
#' @param subsample.flat.X
#' @param Z.prior.params
#'
#' @return
#' @export
#' @importFrom vegan isomap
#'
#' @examples
fit.lsa_bcsgplvm <- function(X,
                             q=2,
                             iterations=1000,
                             plot.freq=100,
                             classes=1,
                             Z.init=NULL,
                             A.init=NULL,
                             K.bc.l="auto",
                             K.bc.target.median=0.3,
                             Z.prior=c("normal", "uniform", "discriminative"),
                             par.init=NULL,
                             points.in.approximation=1024,
                             initial.step.size.par=10^-1,
                             final.step.size.par=10^-3,
                             initial.step.size=10^-3,
                             final.step.size=10^-3,
                             momentum.rate=0.5,
                             learning.rate=0.9,
                             adam.epsilon=10^-8,
                             parameter.opt.iterations=300,
                             par.fixed.A.opt=NULL,
                             verbose=FALSE,
                             subsample.flat.X=NULL,
                             Z.prior.params=list()
) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  out <- list()
  out$Z.prior <- Z.prior
  out$Z.prior.params <- Z.prior.params
  step.size.range <- c(initial.step.size, final.step.size)
  par.step.size.range <- c(initial.step.size.par, final.step.size.par)
  if (!(is.null(Z.init) | is.null(A.init))) {
    stop("Can't specify both Z.init and A.init")
  }

  if (is.null(subsample.flat.X)) {
    X.unstructured <- t(apply(X, 1, as.numeric))
  } else {
    if (!is.numeric(subsample.flat.X)) stop("subsample.flat.X should be NULL (for no subsampling), or an integer.")
    ncols <- prod(dim(X)[-1])
    sampled.cols <- sample(ncols, subsample.flat.X)
    ind <- t(sapply(sampled.cols, vec_to_arr_index, limits=dim(X)[-1]))
    X.unstructured <- matrix(0, nrow=dim(X)[1], ncol=subsample.flat.X)
    for (i in 1:dim(X)[1]) {
      X.unstructured[i,] <- X[cbind(i, ind)]
    }
  }

  if (is.null(K.bc.l)) {
    #Turn off constraint
    K.bc <- diag(nrow(X.unstructured))
  } else {
    if (K.bc.l == "auto") {
      K.bc.l <- select.bc.l.median(X.unstructured, target.median.cor = K.bc.target.median)
    }
    K.bc <- gplvm.SE(Z=X.unstructured, l=K.bc.l, alpha=1, sigma=0)
  }



  if (is.null(Z.init)) {
    if (is.null(A.init)) {
      A.init <- matrix(rnorm(nrow(X)*q), ncol=q)

      # Center and normalise the Z values
      Z <- bc.Z(K.bc=K.bc, A=A.init)
      m <- colMeans(Z)
      K.bc.chol <- chol(K.bc)
      M.1 <- backsolve(K.bc.chol, forwardsolve(t(K.bc.chol), rep(1, nrow(K.bc))))
      M <- matrix(M.1, nrow=nrow(A.init), ncol=ncol(A.init)) %*% diag(m)
      A.init <- A.init - M

      Z.sd <- apply(t(Z) - colMeans(Z), 1, sd)
      A.init <- A.init %*% diag(1 / Z.sd)
      Z <- bc.Z(K.bc=K.bc, A=A.init)
      Z.svd <- svd(Z)
      A.init <- A.init %*% Z.svd$v
    } else {
      if (ncol(as.matrix(A.init)) != q) stop("Mismatch between A.init and q")
      A.init <- as.matrix(A.init)
    }
  } else if (identical(Z.init, "PCA")) {
    X.pca <- prcomp(X.unstructured)
    target <- X.pca$x[,1:q]
    target.sd <- sd(as.numeric(target))
    target <- target / target.sd
    A.init <- fit.A(target = target, K.bc = K.bc)
  } else if (identical(Z.init, "ISOMAP")) {
    X.isomap <- isomap(dist(X.unstructured), k=25, ndim=q, fragmentedOK=FALSE)
    target <- X.isomap$points
    target.sd <- sd(as.numeric(target))
    target <- target / target.sd
    A.init <- fit.A(target = target, K.bc = K.bc)
  } else {
    if (ncol(as.matrix(Z.init)) != q) stop("Mismatch between Z.init and q")
    A.init <- fit.A(target = as.matrix(Z.init)[,1:q], K.bc = K.bc)
  }

  Z <- bc.Z(K.bc=K.bc, A=A.init)
  out$A.init <- A.init
  out$Z.init <- Z
  out$K.bc.l <- K.bc.l
  out$K.bc <- K.bc
  if (plot.freq == 0) {
    plot.freq <- iterations+1
  }

  if (q > 2) {
    pairs(Z, col=classes)
  } else if (q == 2) {
    plot(Z, col=classes)
  } else {
    plot(as.numeric(Z), rep(0, length(Z)), col=classes)
  }

  # Set initial parameters
  if (is.null(par.init)) {
    l_Z <- as.numeric(quantile(dist(Z), 0.1))
    alpha <- sd(as.numeric(X.unstructured)) / sqrt(2)
    sigma <- alpha
    l <- c(l_Z, rep(1, length(dim(X)) - 1))
    par <- c(alpha=alpha, sigma=sigma, l)
  } else {
    if(length(par.init) != length(dim(X)) + 2) stop("par.init is not the correct length")
    par <- par.init
  }

  #optimize with fixed A
  opt <- LSA_BCSGPLVM.sgdopt(X=X, A.init=A.init, par.init=par, K.bc=K.bc,
                             points.in.approximation=points.in.approximation,
                             iterations=parameter.opt.iterations,
                             step.size.range=step.size.range,
                             optimize.A=FALSE, classes=classes, plot.freq=plot.freq,
                             learning.rate=learning.rate, momentum.rate=momentum.rate,
                             adam.epsilon=adam.epsilon,
                             par.step.size.range=par.step.size.range,
                             verbose=verbose, Z.prior=Z.prior,
                             Z.prior.params=Z.prior.params)

  par <- opt$par

  out$par.opt <- opt

  # Optimize pars and A
  opt <- LSA_BCSGPLVM.sgdopt(X=X, A.init=A.init, par.init=par, K.bc=K.bc,
                             points.in.approximation=points.in.approximation,
                             iterations=iterations,
                             step.size.range=step.size.range,
                             optimize.A=TRUE, classes=classes, plot.freq=plot.freq,
                             learning.rate=learning.rate, momentum.rate=momentum.rate,
                             adam.epsilon=adam.epsilon,
                             par.step.size.range=par.step.size.range,
                             par.fixed=par.fixed.A.opt,
                             verbose=verbose, Z.prior=Z.prior,
                             Z.prior.params=Z.prior.params)


  out$A.opt <- opt

  return(out)
}

#' Replay the plots from fitting an LSA_BCSGPLVM model
#'
#' @param fit.lsa_bcsgplvm.output optimization to replay
#' @param time runtime
#' @param pps plots per second
#'
#' @return
#' @export
replay.plots <- function(fit.lsa_bcsgplvm.output, time=30, pps=5) {
  data <- fit.lsa_bcsgplvm.output
  par.hist <- data$par.hist
  A.hist <- data$A.hist
  K.bc <- data$K.bc
  classes <- data$classes

  plots <- min(time * pps, nrow(par.hist))
  interval <- time / plots
  plot.freq <- max(floor(nrow(par.hist) / plots), 1)
  start.time <- proc.time()[3]
  num.it <- 0
  time.mean <- 0
  it.start.time <- proc.time()[3]
  for (iteration in 2:nrow(par.hist)) {
    if (iteration%%plot.freq==0) {
      Z <- bc.Z(K.bc, A.hist[iteration, , ])
      LSA_BCSGPLVM.plot_iteration(Z, par.hist, TRUE, iteration, classes)
      it.run.time <- proc.time()[3] - it.start.time
      it.start.time <- proc.time()[3]
      time.mean <- (time.mean * num.it + it.run.time) / (num.it + 1)
      num.it <- num.it + 1
      if (time.mean > interval) {
        plots <- floor(time / time.mean)
        plot.freq <- max(floor(nrow(par.hist) / plots), 1)
      }
      if (it.run.time < interval) Sys.sleep(interval - it.run.time)
    }
  }
  layout(1)
}
