# Priors on Z
ZPriorUnif_ <- "uniform"
ZPriorNorm_ <- "normal"
ZPriorDisc_ <- "discriminative"
ZPriorOptions_ <- c(ZPriorNorm_, ZPriorUnif_, ZPriorDisc_)

# Optimization methods
ADAM_ <- "ADAM"
SMD_ <- "SMD"

#' Title
#'
#' @param g
#' @param x
#' @param v
#' @param r
#'
#' @return
#' @import numDeriv
#'
#' @examples
approx.hessian.vector.product <- function(g, x, v) {
  num <- as.numeric(jacobian(function(r) g(x + r*v), 0))
  return(num)
}

calc.W <- function(W.pre, Z) {
  return(cbind(Z[W.pre[, 1], ], W.pre[, -1]))
}

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

sample.cols.from.array <- function(dataArray, ind) {
  X.unstructured <- matrix(0, nrow=dim(dataArray)[1], ncol=nrow(ind))
  for (i in 1:dim(dataArray)[1]) {
    X.unstructured[i,] <- dataArray[cbind(i, ind)]
  }
  return(X.unstructured)
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
    M[i, ] <- colMeans(Z[index, , drop=F])
    temp.sum <- M[i, ] - M_0
    S_b <- S_b + length(index) * outer(temp.sum, temp.sum)
    Z.centered.class <- scale(Z[index, , drop=F], center=M[i, ], scale=F)
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

LSA_BCSGPLVM.dL.dZ <- function(W.pre, X, Z,
                               l, alpha, sigma,
                               Z.prior=c("normal", "uniform", "discriminative"),
                               K=NULL, dL.dK=NULL,
                               Z.prior.params=list()) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  X <- as.matrix(X)
  if (missing(K)) {
    W <- calc.W(W.pre=W.pre, Z=Z)
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
  return(out)
}

LSA_BCSGPLVM.dL.dA <- function(W.pre, X, Z,
                               l, alpha, sigma,
                               K.bc,
                               Z.prior=c("normal", "uniform", "discriminative"),
                               K=NULL, dL.dK=NULL,
                               Z.prior.params=list()) {
  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  X <- as.matrix(X)
  if (missing(K)) {
    W <- calc.W(W.pre=W.pre, Z=Z)
    K <- LSA_BCSGPLVM.kernel(W=W, l=l, alpha=alpha, sigma=sigma)
  }
  if (missing(dL.dK)) {
    dL.dK <- dL.dK(X, K)
  }
  dL.dZ <- LSA_BCSGPLVM.dL.dZ(W.pre=W.pre, X=X, Z=Z,
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

LSA_BCSGPLVM.plot_iteration <- function(A.hist, par.hist, plot.Z, iteration, classes, K.bc, max.iterations=1000) {
  close.screen(all.screens=TRUE)
  num.pars <- ncol(par.hist)
  iteration.range <- max(1, iteration - max.iterations + 1):iteration
  iteration.lims <- range(iteration.range)
  iteration.lims[1] <- iteration.lims[1] - 1
  if (plot.Z) {
    split.screen(c(1, 2))
    par.screens <- split.screen(c(num.pars + 1, 1), 1)
    Z <- bc.Z(K.bc, A.hist[iteration, ,])
    if (iteration > 1) {
      prev.Z <- bc.Z(K.bc, A.hist[iteration - 1, ,])
    }
    q <- ncol(Z)
    if (q > 2) {
      plots <- q * (q - 1) / 2
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

  # Do the legend first so we can plot the bottom axis over it.
  screen(par.screens[length(par.screens)])
  par(mar=c(0,0,0,0))
  legend("center", legend=c(paste("It.", iteration), "alpha", "tau", "l_Z", paste("l_S", 1:(ncol(par.hist) - 4), sep=""), "L"),
         col=c(0, viridis::viridis(ncol(par.hist), end=0.8)), lty=1, ncol=4, bty="n")

  for (i in 1:ncol(par.hist)) {
    screen(par.screens[i])
    par(mar=c(0,4,0,4)+0.1)
    plot(iteration.range, par.hist[iteration.range, i],
         type="l", xlim=iteration.lims,
         ylim=range(par.hist[iteration.range, i]),
         axes=F, ann=F, col=viridis::viridis(ncol(par.hist), end=0.8)[i])
    box()
    range <- range(par.hist[iteration.range, i])
    at <- mean(range)
    at <- c(at - diff(range)/3, at, at + diff(range)/3)
    labels <- format(at, digits=2)
    axis(2, at=at, labels = labels, las=1)
    at <- mean(range)
    labels <- format(par.hist[iteration, i], digits=2)
    axis(4, at=at, labels=labels, las=1, tick=FALSE)
    if (i == ncol(par.hist)) axis(1, at=iteration.lims)
    par(mar=c(5.1,4.1,4.1,2.1))
  }
  if (plot.Z) {
    class.colours <- viridis::viridis(max(classes), end=0.8)[classes]
    if (q > 2) {
      screen.num <- 0
      for (i in 1:(q-1)) {
        for (j in (i+1):q) {
          screen.num <- screen.num + 1
          screen(plot.screens[screen.num])
          par(mar=c(0,0,0,0) + 0.1)
          plot(Z[, c(i,j)], col=class.colours, axes=FALSE, ann=FALSE)
          arrows(prev.Z[, i], prev.Z[, j], Z[, i], Z[, j], length=0.05, col=class.colours)
          box()
        }
      }
      par(mar=c(5.1,4.1,4.1,2.1))
    } else if (q == 2) {
      screen(2)
      plot(Z, col=class.colours)
      arrows(prev.Z[, 1], prev.Z[, 2], Z[, 1], Z[, 2], length=0.05, col=class.colours)
    } else {
      screen(2)
      plot(as.numeric(Z), rep(0, length(Z)), col=class.colours)
      arrows(as.numeric(prev.Z), rep(0, length(Z)), as.numeric(Z), rep(0, length(Z)), length=0.05, col=class.colours)
    }
  }
  close.screen(all.screens=TRUE)
}

LSA_BCSGPLVM.sgdopt <- function(X, A.init, par.init, K.bc, points.in.approximation,
                                iterations,
                                optimize.A, classes, plot.freq, par.fixed=NULL,
                                verbose=FALSE, Z.prior=c("normal", "uniform", "discriminative"),
                                Z.prior.params=list(),
                                optimization.method=c("SMD", "ADAM"),
                                optimization.method.pars,
                                ivm=FALSE, ivm.selection.size=NULL,
                                optimizing.structure=FALSE) {

  optimization.method <- match.arg(optimization.method)

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

  if (optimization.method == ADAM_) {
    learning.rate <- optimization.method.pars$learning.rate
    momentum.rate <- optimization.method.pars$momentum.rate
    adam.epsilon <- optimization.method.pars$adam.epsilon
    step.size.range <- optimization.method.pars$step.size.range
    par.step.size.range <- optimization.method.pars$par.step.size.range

    step.size <- step.size.range[1]
    step.size.change <- (step.size.range[1] - step.size.range[2]) / iterations

    if (is.null(par.step.size.range)) {
      par.step.size <- step.size
      par.step.size.change <- step.size.change
    } else {
      par.step.size <- par.step.size.range[1]
      par.step.size.change <- (par.step.size.range[1] - par.step.size.range[2]) / iterations
    }

    # Variables for ADAM update
    m.par <- numeric(length(par))
    m.A <- matrix(0, nrow=nrow(A), ncol=ncol(A))
    v.par <- numeric(length(par))
    v.A <- matrix(0, nrow=nrow(A), ncol=ncol(A))
  }

  if (optimization.method == SMD_) {
    if (length(optimization.method.pars$par.initial.step.size) == length(par)) {
      smd.par.a_i <- optimization.method.pars$par.initial.step.size
    } else {
      smd.par.a_i <- rep(optimization.method.pars$par.initial.step.size, length.out=length(par))
    }
    smd.par.a_i[par.fixed] <- 0
    smd.par.v_i <- rep(0, length(par))
    if (optimize.A) {
      if (length(optimization.method.pars$par.initial.step.size) == 1) {
        smd.A.a_i <- matrix(optimization.method.pars$A.initial.step.size, ncol=ncol(A), nrow=nrow(A))
      } else {
        stop("optimization.method.pars$A.initial.step.size should be a scalar value.")
      }
      smd.A.v_i <- matrix(0, ncol=ncol(A), nrow=nrow(A))
    }

    smd.mu <- optimization.method.pars$meta.step.size.mu
    smd.lambda <- optimization.method.pars$learning.rate.lambda
    if (smd.lambda < 0 | smd.lambda > 1) stop("optimization.method.pars$learning.rate.lambda must be between 0 and 1 (inclusive).")


  }

  par.hist <- matrix(0, ncol=length(par)+1, nrow=iterations + 1)
  if (optimize.A) {
    A.hist <- array(0, dim=c(iterations+1, dim(A)))
  }


  # Initialize some variables which we will use to avoid including NA entries
  X.na.indices <- which(is.na(X))
  num.points <- prod(dim(X)) - length(X.na.indices)
  X.num.not.na.per.row <- apply(X, 1, function(row) sum(!is.na(row)))
  while (iteration < iterations) {
    iteration <- iteration + 1
    if (optimize.A) {
      Z <- bc.Z(K.bc, A)
    }


    # Select a subset of points to use in this step
    if (ivm) {
      temp.sample.size <- ivm.selection.size
    } else {
      temp.sample.size <- points.in.approximation
    }

    if (optimizing.structure) {
      num.non.na.points <- 0
      randomly.sorted.indices <- sample(dim(X)[1], dim(X)[1], replace=FALSE)
      current.index <- 0
      while (num.non.na.points < temp.sample.size) {
        current.index <- current.index + 1
        num.non.na.points <- num.non.na.points + X.num.not.na.per.row[randomly.sorted.indices[current.index]]
      }
      preselect.index <- randomly.sorted.indices[seq(current.index)]
      # Use slice.index to slice X by row - "slice.index(X, 1) %in% preselect.index"
      # will return an array with the same dimensions as X which is TRUE that entry's
      # row number is in preselect.index, and FALSE otherwise. This lets us subset
      # based on the index of the first dimension
      preselect.non.na.indices <- which((!is.na(X)) & (slice.index(X, 1) %in% preselect.index))
      sample.points <- sample(preselect.non.na.indices, temp.sample.size, replace=FALSE)
      W.pre <- t(sapply(sample.points, vec_to_arr_index, limits=dim(X)))
    } else {
      sample.points <- sample.int(n=num.points, size=temp.sample.size, replace=FALSE)
      if (length(X.na.indices) > 0) {
        for (na.index in X.na.indices) {
          sample.points[sample.points >= na.index] <- sample.points[sample.points >= na.index] + 1
        }
      }
      W.pre <- t(sapply(sample.points, vec_to_arr_index, limits=dim(X)))
    }

    W.pre <- W.pre[order(W.pre[,1]),]
    W <- calc.W(W.pre=W.pre, Z=Z)

    if (ivm) {
      # select using IVM
      kernel.function <- function(x1, x2) LSA_BCSGPLVM.kernel(W=rbind(x1, x2), l = par[-(1:2)], alpha = par[1], 0)[1,2]

      ivm.sample.points <- IVM::IVM.regression(predictors = W,
                                               activeSetSize = points.in.approximation,
                                               variance = 1 / par[2],
                                               kernel.function = kernel.function)

      W.pre <- W.pre[ivm.sample.points$activeSet, ]
      W <- W[ivm.sample.points$activeSet, ]
    }

    X.sample <- as.matrix(X[W.pre])

    # Calculate the gradient w.r.t. the parameters
    K <- LSA_BCSGPLVM.kernel(W, par[-(1:2)], par[1], 1 / par[2])
    if (verbose) print(paste("Squared sum of K off-diag:", sum((K^2)) - sum(diag(K)^2)))
    dL.dK <- dL.dK(X.sample, K)
    L <- LSA_BCSGPLVM.L(W=W, X=X.sample, l=par[-(1:2)], alpha=par[1], sigma=1 / par[2], Z.prior=Z.prior, K=K, Z.prior.params=Z.prior.params)

    if (optimize.A) {
      A.hist[iteration, , ] <- A
    }
    par.hist[iteration, ] <- c(par, L)

    if (iteration > 1 & plot.freq != 0 & iteration %% plot.freq == 0) {
      LSA_BCSGPLVM.plot_iteration(A.hist, par.hist, optimize.A, iteration, classes, K.bc)
    }

    if (optimize.A) {
      dL.dA <- LSA_BCSGPLVM.dL.dA(W.pre=W.pre, X=X.sample, Z=Z,
                                  l=par[-(1:2)], alpha=par[1], sigma=1 / par[2],
                                  K.bc=K.bc, Z.prior=Z.prior, K=K, dL.dK=dL.dK,
                                  Z.prior.params=Z.prior.params)
    }

    dL.dpar <- LSA_BCSGPLVM.dL.dpar(W, X.sample, par[-(1:2)], par[1], 1 / par[2], K=K, dL.dK=dL.dK)
    dL.dpar[2] <- -dL.dpar[2] / par[2]^2

    if (optimization.method == ADAM_) {

      if (optimize.A) {
        m.A <- momentum.rate * m.A + (1 - momentum.rate) * dL.dA
        v.A <- learning.rate * v.A + (1 - learning.rate) * dL.dA^2

        m.A.hat <- m.A / (1 - momentum.rate^iteration)
        v.A.hat <- v.A / (1 - learning.rate^iteration)
        delta.A <- step.size / (sqrt(v.A.hat) + adam.epsilon) * m.A.hat
      }

      # Update the Adam variables
      m.par <- momentum.rate * m.par + (1 - momentum.rate) * dL.dpar
      v.par <- learning.rate * v.par + (1 - learning.rate) * dL.dpar^2

      m.par.hat <- m.par / (1 - momentum.rate^iteration)
      v.par.hat <- v.par / (1 - learning.rate^iteration)

      # Calculate the change (including momentum)
      delta.par <- par.step.size / (sqrt(v.par.hat) + adam.epsilon) * m.par.hat
    }

    if (optimization.method == SMD_) {
      if (optimize.A) {

        if (iteration != 1) {
          smd.A.a_i <- smd.A.a_i * pmax(matrix(0.5, ncol=ncol(A), nrow=nrow(A)),
                                        1 - smd.mu * smd.A.v_i * dL.dA)
        }
        delta.A <- smd.A.a_i * dL.dA
        g <- function(x) {
          par <- x[seq_along(par)]
          A <- matrix(x[-seq_along(par)], nrow=nrow(A), ncol=ncol(A))
          Z <- bc.Z(K.bc, A)
          W <- calc.W(W.pre=W.pre, Z=Z)
          K <- LSA_BCSGPLVM.kernel(W=W, l=par[-(1:2)], alpha=par[1], sigma=1 / par[2])
          dL.dK <- dL.dK(X.sample, K)
          dL.dA <- LSA_BCSGPLVM.dL.dA(W.pre=W.pre, X=X.sample, Z=Z,
                                      l=par[-(1:2)], alpha=par[1], sigma=1 / par[2],
                                      K.bc=K.bc, Z.prior=Z.prior,
                                      Z.prior.params=Z.prior.params,
                                      K=K, dL.dK=dL.dK)
          dL.dpar <- LSA_BCSGPLVM.dL.dpar(W=W, X=X.sample, l=par[-(1:2)], alpha=par[1], sigma=1 / par[2],
                                          K=K, dL.dK=dL.dK)
          dL.dpar[2] <- -dL.dpar[2] / par[2]^2
          return(c(dL.dpar, as.numeric(dL.dA)))
        }
        smd.H_i.v_i <- approx.hessian.vector.product(g,
                                                     c(par, as.numeric(A)),
                                                     c(smd.par.v_i, as.numeric(smd.A.v_i)))
        smd.par.H_i.v_i <- smd.H_i.v_i[seq_along(par)]
        smd.A.H_i.v_i <- matrix(smd.H_i.v_i[-seq_along(par)], nrow=nrow(A), ncol=ncol(A))
        smd.A.v_i <- smd.lambda * smd.A.v_i + smd.A.a_i * (smd.A.H_i.v_i - dL.dA)
      } else {
        g <- function(x) {
          out <- LSA_BCSGPLVM.dL.dpar(W, X.sample, x[-(1:2)], x[1], 1 / x[2])
          out[2] <- -out[2] / x[2]^2
          return(out)
        }
        smd.par.H_i.v_i <- approx.hessian.vector.product(g, par, smd.par.v_i)
      }

      if (iteration != 1) {
        smd.par.a_i <- smd.par.a_i * pmax(0.5, 1 - smd.mu * smd.par.v_i * dL.dpar)
      }
      if(verbose) print("dL.dpar")
      if(verbose) print(dL.dpar)
      if(verbose) print("a_i")
      if(verbose) print(smd.par.a_i)
      if(verbose) print("v_i")
      if(verbose) print(smd.par.v_i)

      delta.par <- smd.par.a_i * dL.dpar

      if(verbose) print("H_i.v_i")
      if(verbose) print(smd.par.H_i.v_i)
      smd.par.v_i <- smd.lambda * smd.par.v_i + smd.par.a_i * (smd.lambda * smd.par.H_i.v_i - dL.dpar)
    }

    # Update par and step size
    par[!par.fixed] <- par[!par.fixed] + delta.par[!par.fixed]
    if (optimize.A) {
      A <- A + delta.A
    }
    if (verbose) print(par)
    if (optimization.method == ADAM_) {
      par.step.size <- par.step.size - par.step.size.change
      step.size <- step.size - step.size.change
    }
  }

  L <- LSA_BCSGPLVM.L(W=W, X=X.sample, l=par[-(1:2)], alpha=par[1], sigma=1 / par[2], Z.prior=Z.prior, Z.prior.params=Z.prior.params)

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
#
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
#' @param Z.init Either a matrix of initial latent values, or "PCA" for PCA
#'   start values, or ISOMAP for ISOMAP start values.
#' @param A.init
#' @param K.bc.l lengthscale to use for backconstraints. Either a numeric value
#'   or "auto", which automatically determines an appropriate lengthscale.
#' @param K.bc.l.selection.params parameters for algorithm which automatically
#'   selects backconstraint lengthscale. See \link{select.bc.l.centile} for
#'   details.
#' @param par.init Vector of parameters: alpha, sigma, l_Z, followed by the
#'   lengthscales for the structural dimensions
#' @param points.in.approximation
#' @param Z.prior
#' @param parameter.opt.iterations
#' @param par.fixed.A.opt
#' @param verbose
#' @param subsample.flat.X
#' @param Z.prior.params
#' @param save.X
#' @param optimization.method
#' @param optimization.method.pars
#' @param ivm select points in each step using IVM
#' @param ivm.selection.size number of points to consider for each IVM
#'   selection, NULL for all points
#' @param par.fixed.par.opt
#' @param optimize.structure.params.first
#' @param optimize.all.params
#' @param K.bc.l.plot.graphs
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
                             K.bc.l.selection.params=NULL,
                             K.bc.l.plot.graphs=T,
                             Z.prior=c("normal", "uniform", "discriminative"),
                             par.init=NULL,
                             points.in.approximation=1024,
                             optimization.method=c("SMD", "ADAM"),
                             optimization.method.pars=NULL,
                             parameter.opt.iterations=300,
                             par.fixed.par.opt=NULL,
                             par.fixed.A.opt=NULL,
                             verbose=FALSE,
                             subsample.flat.X=NULL,
                             Z.prior.params=list(),
                             save.X=FALSE,
                             optimize.structure.params.first=TRUE,
                             optimize.all.params=FALSE,
                             ivm=FALSE,
                             ivm.selection.size=2048
) {
  optimization.method <- match.arg(optimization.method)

  if (is.null(optimization.method.pars)) {
    # Set up default optimization parameters
    optimization.method.pars <- list()
    if (optimization.method == ADAM_) {
      optimization.method.pars$learning.rate <- 0.9
      optimization.method.pars$momentum.rate <- 0.9
      optimization.method.pars$adam.epsilon <- 10^-8
      optimization.method.pars$step.size.range <- c(10^-1, 10^-3)
      optimization.method.pars$par.step.size.range <- c(10^-1, 10^-2)
    }

    if (optimization.method == SMD_) {
      optimization.method.pars$par.initial.step.size <- 10^-2
      optimization.method.pars$A.initial.step.size <- 10^-2
      optimization.method.pars$meta.step.size.mu <- 10^-2
      optimization.method.pars$learning.rate.lambda <- 0.9
    }
  }

  Z.prior <- match.arg(Z.prior, ZPriorOptions_)
  out <- list()
  out$Z.prior <- Z.prior
  out$Z.prior.params <- Z.prior.params

  if (!is.element("IVM", installed.packages()[,1]) & ivm) {
    stop("ivm=TRUE but IVM package is not installed. Rerun with ivm=FALSE.")
  }

  if (save.X) {
    out$X <- X
  }

  if (!(is.null(Z.init) | is.null(A.init))) {
    stop("Can't specify both Z.init and A.init")
  }

  out$dim.X <- dim(X)
  if (is.null(subsample.flat.X)) {
    X.unstructured <- t(apply(X, 1, as.numeric))
  } else {
    if (!is.numeric(subsample.flat.X)) stop("subsample.flat.X should be NULL (for no subsampling), or an integer.")
    ncols <- prod(dim(X)[-1])
    sampled.cols <- sample(ncols, subsample.flat.X)
    ind <- t(sapply(sampled.cols, vec_to_arr_index, limits=dim(X)[-1]))
    out$K.bc.X.ind <- ind
    X.unstructured <- sample.cols.from.array(X, ind)
  }

  #Deal with missing data in X - replace NAs with gaussian weighted average
  X.unstructured.na.indices <- which(is.na(X.unstructured), arr.ind=T)
  if (length(X.unstructured.na.indices) > 0) {
    warning("X contains NAs. Inferring missing entries using weighted average of surrounding pixels for back constraints and PCA initialisation.")
    if(any(apply(X, 1, function(x) all(is.na(x))))) stop("X contains an entry where all values are missing.")
    X.unstructured <-  t(apply(infer.missing.values(X), 1, as.numeric))
  }

  if (is.null(K.bc.l)) {
    #Turn off constraint
    K.bc <- diag(nrow(X.unstructured))
  } else {
    if (K.bc.l == "auto") {
      K.bc.l <- select.bc.l.centile(X.unstructured, K.bc.l.selection.params)
    }
    out$K.bc.l <- K.bc.l
    if(K.bc.l.plot.graphs) {
      bc.l.selection.plots(X.unstructured, chosen.lengthscale=K.bc.l)
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
      M <- matrix(M.1, nrow=nrow(A.init), ncol=ncol(A.init)) %*% diag(m, nrow=length(m), ncol=length(m))
      A.init <- A.init - M

      Z.sd <- apply(t(Z) - colMeans(Z), 1, sd)
      A.init <- A.init %*% diag(1 / Z.sd, nrow=length(Z.sd), ncol=length(Z.sd))
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

  # Optimize structure parameters assuming all samples are independent:
  if (parameter.opt.iterations > 0 & optimize.structure.params.first) {
    str.par <- par
    # tau
    str.par[2] <- 1
    # l_Z
    str.par[3] <- 10^-8
    if (!is.null(par.fixed.par.opt)) {
      str.fixed.par <- par.fixed.par.opt
    } else {
      str.fixed.par <- logical(length(str.par))
    }
    str.fixed.par[3] <- TRUE
    opt <- LSA_BCSGPLVM.sgdopt(X=X, A.init=A.init, par.init=str.par, K.bc=K.bc,
                               points.in.approximation=points.in.approximation,
                               iterations=parameter.opt.iterations,
                               optimize.A=FALSE, classes=classes, plot.freq=plot.freq,
                               par.fixed=str.fixed.par,
                               optimization.method=optimization.method,
                               optimization.method.pars=optimization.method.pars,
                               verbose=verbose, Z.prior=Z.prior,
                               Z.prior.params=Z.prior.params,
                               ivm=ivm, ivm.selection.size=ivm.selection.size,
                               optimizing.structure = TRUE)
    str.par <- opt$par
    out$str.par.opt <- opt

    par[-3] <- str.par[-3]
  }

  #optimize with fixed A (not recommended for SMD or other optimizations with unbounded step sizes)
  if (parameter.opt.iterations > 0 & optimize.all.params) {
    opt <- LSA_BCSGPLVM.sgdopt(X=X, A.init=A.init, par.init=par, K.bc=K.bc,
                               points.in.approximation=points.in.approximation,
                               iterations=parameter.opt.iterations,
                               optimize.A=FALSE, classes=classes, plot.freq=plot.freq,
                               par.fixed=par.fixed.par.opt,
                               optimization.method=optimization.method,
                               optimization.method.pars=optimization.method.pars,
                               verbose=verbose, Z.prior=Z.prior,
                               Z.prior.params=Z.prior.params,
                               ivm=ivm, ivm.selection.size=ivm.selection.size)
    par <- opt$par
    out$par.opt <- opt
  }

  # Optimize pars and A
  opt <- LSA_BCSGPLVM.sgdopt(X=X, A.init=A.init, par.init=par, K.bc=K.bc,
                             points.in.approximation=points.in.approximation,
                             iterations=iterations,
                             optimize.A=TRUE, classes=classes, plot.freq=plot.freq,
                             par.fixed=par.fixed.A.opt,
                             optimization.method=optimization.method,
                             optimization.method.pars=optimization.method.pars,
                             verbose=verbose, Z.prior=Z.prior,
                             Z.prior.params=Z.prior.params,
                             ivm=ivm, ivm.selection.size=ivm.selection.size)


  out$A.opt <- opt

  out$final.A <- opt$A
  out$final.Z <- opt$Z

  class(out) <- "LSA_BCSGPLVM"

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
      LSA_BCSGPLVM.plot_iteration(A.hist, par.hist, TRUE, iteration, classes, K.bc)
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

#' Predict data for an LSA_BCSGPLVM model
#'
#' @param object
#' @param data
#' @param training.data
#'
#' @return
#' @export
#'
#' @examples
predict.LSA_BCSGPLVM <- function(object, data, training.data=NULL) {
  if (is.null(object$K.bc.l)) {
    stop("Provided LSA_BCSGPLVM is unconstrained. Prediction for unconstrained models is not yet implemented.")
  }
  if (!identical(dim(data)[-1], object$dim.X[-1])) {
    stop(paste("data has the wrong dimensions. Expecting,", paste(object$dim.X, collapse=", ")))
  }
  if (is.null(training.data) & is.null(object$X)) {
    stop("The training data must either be saved in the LSA_BCSGPLVM object or provided using training.data")
  } else {
    if (!is.null(training.data) & !is.null(object$X)) {
      warning("Training data provided but also saved in LSA_BCSGPLVM object. Using data from the object.")
    }
    if (is.null(object$X)) {
      if (!identical(dim(training.data), object$dim.X)) {
        stop("Provided training data dimensions do not match expected dimensions for training data")
      }
      object$X <- training.data
    }
  }
  if (is.null(object$K.bc.X.ind)) {
    object$X <- t(apply(object$X, 1, as.numeric))
    data <- t(apply(data, 1, as.numeric))
  } else {
    object$X <- sample.cols.from.array(object$X, object$K.bc.X.ind)
    data <- sample.cols.from.array(data, object$K.bc.X.ind)
  }

  out <- list()
  out$K.star <- matrix(0, nrow=nrow(data), ncol=nrow(object$X))

  for (i in 1:nrow(data)) {
    for (j in 1:nrow(object$X)) {
      out$K.star[i, j] <- sqrt(sum((data[i,] - object$X[j, ])^2))
    }
  }
  out$K.star <- gplvm.SE.dist(out$K.star, object$K.bc.l, 1, 0)

  out$predictions <- out$K.star %*% object$final.A

  return(out)
}
