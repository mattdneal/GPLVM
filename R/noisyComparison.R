#' Compare BCGPLVM to BCSGPLVM
#'
#' @param base.data
#' @param true.latent.vars
#' @param classes
#' @param noise.variances
#'
#' @return
#' @export
#'
#' @importFrom caret train trainControl
#' @importFrom vegan isomap
#' @import gaussianProcess
#'
#' @examples
noisyComparison <- function(base.data,
                            base.data.nrows,
                            true.latent.vars,
                            classes,
                            noise.variances,
                            PCA.init=FALSE,
                            maxit=1000) {
  q = ncol(true.latent.vars)

  fits <- list()
  fits[["PCA"]] <- list()
  fits[["ISOMAP"]] <- list()
  fits[["BCGPLVM"]] <- list()
  fits[["BCSGPLVM"]] <- list()


  regressions <- list()
  regressions[["PCA"]] <- list()
  regressions[["ISOMAP"]] <- list()
  regressions[["BCGPLVM"]] <- list()
  regressions[["BCSGPLVM"]] <- list()

  results_colnames <- c("noise",
                        "pca.accuracy", "pca.accuracy.sd",
                        "isomap.accuracy", "isomap.accuracy.sd",
                        "bc.accuracy", "bc.accuracy.sd",
                        "bcs.accuracy", "bcs.accuracy.sd")

  #Regression colnames

  results_colnames <- c(results_colnames,
                        "pca.likelihood",
                        "isomap.likelihood",
                        "bc.likelihood",
                        "bcs.likelihood")

  results <- matrix(0, nrow=length(noise.variances), ncol=length(results_colnames))
  colnames(results) <- results_colnames

  # Set up the GP kernel for regression
  mt <- create.model.tree.builtin()
  mt <- insert.kernel.instance(mt, 1, "squaredExponential", NULL, hyper.params=c(l=NULL))
  k <- create.kernel.object.from.model.tree(mt)

  rownum <- 0
  for (sigma in noise.variances) {
    noisy.data <- base.data + rnorm(length(base.data), sd=sigma)
    noisy.data <- scale(noisy.data, center=TRUE, scale=FALSE)
    noisy.data.sd <- sd(as.numeric(noisy.data))
    noisy.data <- noisy.data / noisy.data.sd
    Z.pca <- prcomp(noisy.data)$x[, 1:q]

    isomap.model <- isomap(dist(noisy.data), ndim=q, k=8)$points

    K.bc.l <- select.bc.l(noisy.data, 0.1, 0.9)
    Z.init <- NULL
    if (PCA.init) Z.init <- "PCA"

    bcgplvm.model <- fit.bcgplvm(noisy.data,
                                 q,
                                 iterations=maxit,
                                 plot.freq=min(1000, maxit),
                                 classes=classes,
                                 Z.init=Z.init,
                                 A.init=NULL,
                                 num.init.params=100,
                                 K.bc.l=K.bc.l,
                                 Z.normal.prior=TRUE)

    Z.init <- NULL
    if (PCA.init) Z.init <- "PCA"

    bcsgplvm.model <- fit.bcsgplvm(noisy.data,
                                   nrows.X=base.data.nrows,
                                   q=q,
                                   iterations=maxit,
                                   plot.freq=min(1000, maxit),
                                   classes=classes,
                                   Z.init=Z.init,
                                   A.init=NULL,
                                   K.bc.l=K.bc.l,
                                   Z.normal.prior=TRUE,
                                   structured.C.init=c(2, 2),
                                   structured.alpha.init=1,
                                   se.par.init=NULL,
                                   fixed.C=FALSE,
                                   structured.C.bounds=c(1, 1, 5, 5),
                                   noise.lower.bound=10^-4,
                                   probabilistic.trace.estimate=FALSE)

    pca.res <- train(Z.pca,
                     factor(make.names(classes)),
                     method="nnet",
                     trControl=trainControl(method="repeatedcv",
                                            number=2, repeats=5,
                                            classProbs=TRUE))
    pca.acc.ind <- which.max(pca.res$results$Accuracy)

    pca.regr <- true.latent.vars
    pca.regr[] <- NA
    pca.likelihood <- 0

    for (i in 1:ncol(true.latent.vars)) {
      pca.gp <- create.gaussian.process(Z.pca, true.latent.vars[, i], k)
      pca.gp <- fit.hyperparams(pca.gp)
      if (pca.gp$optimx.obj$convcode) warning(paste("PCA GP failed to converge, sigma:", sigma))
      pca.regr[, i] <- as.numeric(predict(pca.gp$gp, Z.pca)$mean)
      pca.likelihood <- pca.likelihood +
        gaussianProcess:::get.marginal.likelihood(pca.gp$gp,
                                                  pca.gp$gp$optimized.sigma.n,
                                                  pca.gp$gp$optimized.hyperparams)
    }

    isomap.res <- train(isomap.model,
                     factor(make.names(classes)),
                     method="nnet",
                     trControl=trainControl(method="repeatedcv",
                                            number=2, repeats=5,
                                            classProbs=TRUE))
    isomap.acc.ind <- which.max(isomap.res$results$Accuracy)

    isomap.regr <- true.latent.vars
    isomap.regr[] <- NA
    isomap.likelihood <- 0

    for (i in 1:ncol(true.latent.vars)) {
      isomap.gp <- create.gaussian.process(isomap.model, true.latent.vars[, i], k)
      isomap.gp <- fit.hyperparams(isomap.gp)
      if (isomap.gp$optimx.obj$convcode) warning(paste("isomap GP failed to converge, sigma:", sigma))
      isomap.regr[, i] <- as.numeric(predict(isomap.gp$gp, isomap.model)$mean)
      isomap.likelihood <- isomap.likelihood +
        gaussianProcess:::get.marginal.likelihood(isomap.gp$gp,
                                                  isomap.gp$gp$optimized.sigma.n,
                                                  isomap.gp$gp$optimized.hyperparams)
    }

    bc.res <- train(bcgplvm.model$Z,
                    factor(make.names(classes)),
                    method="nnet",
                    trControl=trainControl(method="repeatedcv",
                                           number=2, repeats=5,
                                           classProbs=TRUE))
    bc.acc.ind <- which.max(bc.res$results$Accuracy)

    bc.regr <- true.latent.vars
    bc.regr[] <- NA
    bc.likelihood <- 0

    for (i in 1:ncol(true.latent.vars)) {
      bc.gp <- create.gaussian.process(bcgplvm.model$Z, true.latent.vars[, i], k)
      bc.gp <- fit.hyperparams(bc.gp)
      if (bc.gp$optimx.obj$convcode) warning(paste("bc GP failed to converge, sigma:", sigma))
      bc.regr[, i] <- as.numeric(predict(bc.gp$gp, bcgplvm.model$Z)$mean)
      bc.likelihood <- bc.likelihood +
        gaussianProcess:::get.marginal.likelihood(bc.gp$gp,
                                                  bc.gp$gp$optimized.sigma.n,
                                                  bc.gp$gp$optimized.hyperparams)
    }

    bcs.res <- train(bcsgplvm.model$Z,
                    factor(make.names(classes)),
                    method="nnet",
                    trControl=trainControl(method="repeatedcv",
                                           number=2, repeats=5,
                                           classProbs=TRUE))
    bcs.acc.ind <- which.max(bcs.res$results$Accuracy)

    bcs.regr <- true.latent.vars
    bcs.regr[] <- NA
    bcs.likelihood <- 0

    for (i in 1:ncol(true.latent.vars)) {
      bcs.gp <- create.gaussian.process(bcsgplvm.model$Z, true.latent.vars[, i], k)
      bcs.gp <- fit.hyperparams(bcs.gp)
      if (bcs.gp$optimx.obj$convcode) warning(paste("bcs GP failed to converge, sigma:", sigma))
      bcs.regr[, i] <- as.numeric(predict(bcs.gp$gp, bcsgplvm.model$Z)$mean)
      bcs.likelihood <- bcs.likelihood +
        gaussianProcess:::get.marginal.likelihood(bcs.gp$gp,
                                                  bcs.gp$gp$optimized.sigma.n,
                                                  bcs.gp$gp$optimized.hyperparams)
    }

    rownum <- rownum + 1
    temp.colnames <- c("Accuracy", "AccuracySD")
    results[rownum, ] <- as.numeric(c(sigma,
                           pca.res$results[pca.acc.ind, temp.colnames],
                           isomap.res$results[isomap.acc.ind, temp.colnames],
                           bc.res$results[bc.acc.ind, temp.colnames],
                           bcs.res$results[bcs.acc.ind, temp.colnames],
                           pca.likelihood,
                           isomap.likelihood,
                           bc.likelihood,
                           bcs.likelihood))

    plot(results[1:rownum, 1], results[1:rownum, 2], col=1, ylim=range(results[,c(2,4,6,8)]))
    points(results[1:rownum, 1], results[1:rownum, 4], col=2)
    points(results[1:rownum, 1], results[1:rownum, 6], col=3)
    points(results[1:rownum, 1], results[1:rownum, 8], col=4)


    fits[["PCA"]][[as.character(sigma)]] <- Z.pca
    fits[["ISOMAP"]][[as.character(sigma)]] <- isomap.model
    fits[["BCGPLVM"]][[as.character(sigma)]] <- bcgplvm.model
    fits[["BCSGPLVM"]][[as.character(sigma)]] <- bcsgplvm.model


    regressions[["PCA"]][[as.character(sigma)]] <- pca.regr
    regressions[["ISOMAP"]][[as.character(sigma)]] <- isomap.regr
    regressions[["BCGPLVM"]][[as.character(sigma)]] <- bc.regr
    regressions[["BCSGPLVM"]][[as.character(sigma)]] <- bcs.regr


  }
  return(list(results=results, fits=fits, regressions=regressions))
}
