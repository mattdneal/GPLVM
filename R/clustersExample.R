if (FALSE) {
  library(viridis)
  library(plot3D)
  library(GPLVM)
  library(ggplot2)
  library(cowplot)

  num_dim <- 2
  num_clusters <- 3
  num_samples <- 100

  centers <- matrix(rnorm(num_dim * num_clusters), nrow=num_clusters) * 10
  covariances <- array(dim=c(num_dim, num_dim, num_clusters))
  for (i in 1:num_clusters) {
    covariances[, , i] <- rWishart(1, num_dim + 1, diag(1, num_dim))
  }

  Z <- matrix(0, num_samples, num_dim)

  cluster <- numeric(num_samples)

  for (i in 1:num_samples) {
    cluster[i] <- sample(num_clusters, 1)
    Z[i, ] <- MASS::mvrnorm(1, centers[cluster[i], ], covariances[, , cluster[i]])
  }

  Z <- scale(Z)

  pairs(Z, col=cluster)

  S.mat <- matrix(0, 5, 5)
  for (i in 1:nrow(S.mat)) S.mat[i, ] <- i
  S <- cbind(as.numeric(S.mat), as.numeric(t(S.mat)))
  K_S <- GPLVM:::gplvm.SE(S, 3, 1, 1E-6)
  K_Z <- GPLVM:::gplvm.SE(Z, 0.4, 1, 1E-6)
  image(K_Z, useRaster=T)
  image(K_S)
  data <- sample.from.model(Z, 5, K_S, K_Z)

  plotData <- as.data.frame(Z)
  colnames(plotData) <- c("Latent var 1", "Latent var 2")
  plotData <- cbind.data.frame(plotData, "Cluster"=as.character(cluster))
  pointsToShow <- c()
  for (i in 1:3) {
    pointsToShow <- c(pointsToShow, sample(which(cluster==i), 2))
  }
  rightFigure <- function() {
    defMar <- par("mar")
    layout(matrix(1:6, 2))
    showTitle <- F
    topRow=T
    for (i in pointsToShow) {
      if (topRow) par(mar=c(3,2,10,2)) else par(mar=c(10,2,3,2))
      topRow <- !topRow
      if (showTitle) title <- paste("Cluster", cluster[i]) else title <- NULL
      #showTitle <- !showTitle
      image(matrix(data[i,], 5), zlim=range(data), col=viridis(256), main=title, xaxt="n", yaxt="n")
    }
    par(mar=defMar)
  }
  leftFigure <- ggplot(plotData, aes(y=`Latent var 2`, x=`Latent var 1`, color=Cluster)) +
    geom_point(alpha=0.8, shape=16) +
    scale_color_manual(values=viridis(3, end=0.8), guide=FALSE) +
    geom_point(data=plotData[pointsToShow,], size=1.2, shape=1, color="red") +
    theme(text = element_text(size=10))
  figureWidth <- 210 - 38*2
  png(filename="artificial_data_summary.png", width=figureWidth, height=floor(figureWidth*0.6), units="mm", res=600, pointsize=1)
  plot_grid(leftFigure, rightFigure)
  dev.off()
  data.n <- data + rnorm(prod(dim(data))) * 0.1

  data.pca <- prcomp(data)
  plot(data.pca)
  pairs(data.pca$x[,1:3], col=cluster)
  data.n.pca <- prcomp(data.n)
  plot(data.n.pca)
  pairs(data.n.pca$x[,1:3], col=cluster)
  data.array <- array(data, c(num_samples,5,5))
  data.n.array <- array(data.n, c(num_samples,5,5))
  image(data.n.array[2,,])
  dev.off()

  if (FALSE) {
    write.csv(as.numeric(data.n.array), file="synthetic_data_clusters.csv", row.names = F)
    write.csv(cluster, file="synth_data_clusters_class.csv", row.names=F)
  }

  data.partially_observed <- data.n
  hidden.indices <- sample(length(data.partially_observed), length(data.partially_observed)*0.1)
  data.partially_observed[hidden.indices] <- NA
  data.po.array <- array(data.partially_observed, c(num_samples,5,5))

  gplvm <- fit.lsa_bcsgplvm(data.po.array,
                            q=num_dim + 1,
                            iterations=10000,
                            plot.freq=50,
                            classes=cluster,
                            Z.init=NULL,
                            A.init=NULL,
                            K.bc.l="auto",
                            K.bc.target.median=0.00001,
                            Z.prior="normal",
                            par.init=NULL,
                            points.in.approximation=100,
                            optimization.method="ADAM",
                            optimization.method.pars=list(learning.rate=0.9, momentum.rate=0.9, adam.epsilon=10^-8, step.size.range=c(10^-2, 10^-2), par.step.size.range=c(10^-2, 10^-2)),
                            parameter.opt.iterations=1000,
                            par.fixed.par.opt=NULL,
                            par.fixed.A.opt=c(F,F,T,T,T),#c(F, F, F, F, F),
                            verbose=FALSE,
                            subsample.flat.X=NULL,
                            Z.prior.params=list(),
                            save.X=FALSE,
                            optimize.structure.params.first=TRUE,
                            optimize.all.params=T)


  gplvm.unstructured <- fit.lsa_bcsgplvm(data.n.array,
                                         q=num_dim + 1,
                                         iterations=1000,
                                         plot.freq=50,
                                         classes=cluster,
                                         Z.init=NULL,
                                         A.init=NULL,
                                         K.bc.l="auto",
                                         K.bc.target.median=0.00001,
                                         Z.prior="normal",
                                         par.init=c(1, 1, 1, 1e-8, 1e-8),
                                         points.in.approximation=100,
                                         optimization.method="ADAM",
                                         optimization.method.pars=list(learning.rate=0.9, momentum.rate=0.9, adam.epsilon=10^-8, step.size.range=c(10^-2, 10^-2), par.step.size.range=c(10^-2, 10^-2)),
                                         parameter.opt.iterations=1000,
                                         par.fixed.par.opt=c(F,F,T,T,T),
                                         par.fixed.A.opt=c(F,F,T,T,T),#c(F, F, F, F, F),
                                         verbose=FALSE,
                                         subsample.flat.X=NULL,
                                         Z.prior.params=list(),
                                         save.X=FALSE,
                                         optimize.structure.params.first=TRUE,
                                         optimize.all.params=F)

  pairs(cbind(prcomp(Z)$x, prcomp(gplvm$final.Z)$x, data.n.pca$x[,1:3]), col=cluster)
  pairs(cbind(prcomp(Z)$x, prcomp(gplvm$final.Z)$x, prcomp(gplvm.unstructured$final.Z)$x), col=cluster)

  pairs(cbind(gplvm$final.Z, gplvm.unstructured$final.Z), col=cluster)

  K_Z_estimate <- GPLVM:::gplvm.SE(gplvm$A.opt$Z,  gplvm$A.opt$par[3], gplvm$A.opt$par[1], 0)
  K_Z_estimate.unstructured <- GPLVM:::gplvm.SE(gplvm.unstructured$A.opt$Z,
                                                gplvm.unstructured$A.opt$par[3],
                                                gplvm.unstructured$A.opt$par[1],
                                                0)

  sum((K_Z - K_Z_estimate)^2)
  sum((K_Z - K_Z_estimate.unstructured)^2)

}

