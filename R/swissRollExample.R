if (FALSE) {
  library(viridis)
  library(plot3D)
  library(GPLVM)
  Z.temp <- matrix(rnorm(2000), ncol=2) + 4
  Z.temp[1:500, ] <- Z.temp[1:500, ] + 3
  Z.temp <- Z.temp[order(Z.temp[,1]), ]
  col_palette <- viridis(1000)
  col <- col_palette[rank(Z.temp[,1])]
  plot(Z.temp, col=col)
  r <- Z.temp[,1]
  Z <- cbind(r * cos(Z.temp[,1]), r * sin(Z.temp[,1]), Z.temp[,2])
  Z <- scale(Z)
  pairs(Z, col=col)

  S.mat <- matrix(0, 5, 5)
  for (i in 1:nrow(S.mat)) S.mat[i, ] <- i
  S <- cbind(as.numeric(S.mat), as.numeric(t(S.mat)))
  K_S <- GPLVM:::gplvm.SE(S, 3, 3, 1E-6)
  K_Z <- GPLVM:::gplvm.SE(Z.temp, 1, 1, 1E-6)
  image(K_Z, useRaster=T)
  data <- sample.from.model(Z.temp, 5, K_S, K_Z)

  data.n <- data + rnorm(prod(dim(data))) * 0.5
  data.pca <- prcomp(data)
  plot(data.pca)
  pairs(data.pca$x[,1:3], col=col)
  data.n.pca <- prcomp(data.n)
  plot(data.n.pca)
  pairs(data.n.pca$x[,1:3], col=col)
  data.array <- array(data, c(1000,5,5))
  data.n.array <- array(data.n, c(1000,5,5))
  image(data.n.array[2,,])
  dev.off()
  gplvm <- fit.lsa_bcsgplvm(data.n.array,
                            q=2,
                            iterations=10000,
                            plot.freq=100,
                            classes=col,
                            Z.init="PCA",
                            A.init=NULL,
                            K.bc.l=NULL,
                            K.bc.target.median=NULL,
                            Z.prior="normal",
                            par.init=NULL,
                            points.in.approximation=50,
                            optimization.method="ADAM",
                            optimization.method.pars=list(learning.rate=0.9, momentum.rate=0.9, adam.epsilon=10^-8, step.size.range=c(10^-1, 10^-1), par.step.size.range=c(10^-2, 10^-2)),
                            parameter.opt.iterations=1000,
                            par.fixed.par.opt=NULL,
                            par.fixed.A.opt=c(T,T,T,T,T),#c(F, F, F, F, F),
                            verbose=FALSE,
                            subsample.flat.X=NULL,
                            Z.prior.params=list(),
                            save.X=FALSE,
                            optimize.structure.params.first=TRUE,
                            optimize.all.params=T)

 g plvm.unstructured <- fit.lsa_bcsgplvm(data.n.array,
                            q=2,
                            iterations=10000,
                            plot.freq=10,
                            classes=col,
                            Z.init="PCA",
                            A.init=NULL,
                            K.bc.l=NULL,
                            K.bc.target.median=NULL,
                            Z.prior="normal",
                            par.init=c(1, 0.5, 1, 1e-8, 1e-8),
                            points.in.approximation=50,
                            optimization.method="ADAM",
                            optimization.method.pars=list(learning.rate=0.9, momentum.rate=0.9, adam.epsilon=10^-8, step.size.range=c(10^-1, 10^-1), par.step.size.range=c(10^-2, 10^-2)),
                            parameter.opt.iterations=1000,
                            par.fixed.par.opt=c(F, F, T, T, T),
                            par.fixed.A.opt=c(F, F, T, T, T),
                            verbose=FALSE,
                            subsample.flat.X=NULL,
                            Z.prior.params=list(),
                            save.X=FALSE,
                            optimize.structure.params.first=TRUE,
                            optimize.all.params=F)

  pairs(cbind(prcomp(Z.temp)$x, prcomp(gplvm$final.Z)$x, data.n.pca$x[,1:3]), col=col)
  pairs(cbind(gplvm$final.Z, gplvm.unstructured$final.Z), col=col)
  }

