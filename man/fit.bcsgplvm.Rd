% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BCSGPLVM.R
\name{fit.bcsgplvm}
\alias{fit.bcsgplvm}
\title{Fit a backconstrained structured GPLVM model}
\usage{
fit.bcsgplvm(X, nrows.X, q = 2, iterations = 1000, plot.freq = 100,
  classes = 1, Z.init = NULL, A.init = NULL, K.bc.l = 0.1,
  Z.normal.prior = TRUE, structured.C.init = c(2, 2),
  structured.alpha.init = 1, se.par.init = NULL, fixed.C = FALSE,
  structured.C.bounds = c(1, 1, 5, 5), noise.lower.bound = 10^-4,
  probabilistic.trace.estimate = TRUE)
}
\arguments{
\item{probabilistic.trace.estimate}{}
}
\description{
Fit a backconstrained structured GPLVM model
}
\examples{
Z.true <- cbind(rnorm(50, mean=-1), rnorm(50, mean=-1), rnorm(50, mean=-1))
Z.true <- rbind(Z.true, cbind(rnorm(50, mean=1), rnorm(50, mean=1), rnorm(50, mean=1)))
classes <- rep(1:2, each=50)
pairs(Z.true, col=classes, pty=".")
X <- generate.structured.dataset(100,20,20,num.latent.variables=NULL,5,0.1,nonlinear=TRUE, Z=Z.true)
X$data <- scale(X$data, center=T, scale=F)
X.sd <- sd(as.numeric(X$data))
X$data <- X$data / X.sd
X.noisy <- X$data + rnorm(length(X$data))
Z.init <- prcomp(X.noisy)$x[,1:3]
bcsgplvm.model <- fit.bcsgplvm(X=X.noisy, nrows.X=20,
                               q=3, iterations=1000, plot.freq=100,
                               classes=classes, Z.init="PCA", K.bc.l=50,
                               probabilistic.trace.estimate=FALSE)
}
