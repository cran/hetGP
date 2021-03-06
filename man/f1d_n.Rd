% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestFunctions.R
\name{f1d_n}
\alias{f1d_n}
\title{Noisy 1d test function (1)
Add Gaussian noise with variance r(x) = scale * (1.1 + sin(2 pi x))^2 to \code{\link[hetGP]{f1d}}}
\usage{
f1d_n(x, scale = 1)
}
\arguments{
\item{x}{scalar or matrix (size n x 1) in [0,1]}

\item{scale}{scalar in [0, Inf] to control the signal to noise ratio}
}
\description{
Noisy 1d test function (1)
Add Gaussian noise with variance r(x) = scale * (1.1 + sin(2 pi x))^2 to \code{\link[hetGP]{f1d}}
}
\examples{
X <- matrix(seq(0, 1, length.out = 101), ncol = 1)
Xr <- X[sort(sample(x = 1:101, size = 500, replace = TRUE)),, drop = FALSE]
plot(Xr, f1d_n(Xr))
lines(X, f1d(X), col = "red", lwd = 2)
}
