% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hetTP.R
\name{mleHomTP}
\alias{mleHomTP}
\title{Student-T process modeling with homoskedastic noise}
\usage{
mleHomTP(
  X,
  Z,
  lower = NULL,
  upper = NULL,
  known = list(beta0 = 0),
  noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), 10000), nu_bounds = c(2 +
    0.001, 30), sigma2_bounds = c(sqrt(.Machine$double.eps), 10000)),
  init = list(nu = 3),
  covtype = c("Gaussian", "Matern5_2", "Matern3_2"),
  maxit = 100,
  eps = sqrt(.Machine$double.eps),
  settings = list(return.Ki = TRUE, factr = 1e+09)
)
}
\arguments{
\item{X}{matrix of all designs, one per row, or list with elements:
\itemize{
  \item \code{X0} matrix of unique design locations, one point per row
  \item \code{Z0} vector of averaged observations, of length \code{nrow(X0)}
  \item \code{mult} number of replicates at designs in \code{X0}, of length \code{nrow(X0)}
}}

\item{Z}{vector of all observations. If using a list with \code{X}, \code{Z} has to be ordered with respect to \code{X0}, and of length \code{sum(mult)}}

\item{lower, upper}{bounds for the \code{theta} parameter (see \code{\link[hetGP]{cov_gen}} for the exact parameterization).
In the multivariate case, it is possible to give vectors for bounds (resp. scalars) for anisotropy (resp. isotropy)}

\item{known}{optional list of known parameters, e.g., \code{beta0} (default to \code{0}), \code{theta}, \code{g}, \code{sigma2} or \code{nu}}

\item{noiseControl}{list with element,
\itemize{
\item \code{g_bound}, vector providing minimal and maximal noise variance
\item \code{sigma2_bounds}, vector providing minimal and maximal signal variance
\item \code{nu_bounds}, vector providing minimal and maximal values for the degrees of freedom. 
The mininal value has to be stricly greater than 2. If the mle optimization gives a large value, e.g., 30,
considering a GP with \code{\link[hetGP]{mleHomGP}} may be better. 
}}

\item{init}{list specifying starting values for MLE optimization, with elements:
\itemize{
 \item \code{theta_init} initial value of the theta parameters to be optimized over (default to 10\% of the range determined with \code{lower} and \code{upper})
 \item \code{g_init} initial value of the nugget parameter to be optimized over (based on the variance at replicates if there are any, else 10\% of the variance)
 \item \code{sigma2} initial value of the variance paramter (default to \code{1})
 \item \code{nu} initial value of the degrees of freedom parameter (default to \code{3})
}}

\item{covtype}{covariance kernel type, either 'Gaussian', 'Matern5_2' or 'Matern3_2', see \code{\link[hetGP]{cov_gen}}}

\item{maxit}{maximum number of iteration for L-BFGS-B of \code{\link[stats]{optim}}}

\item{eps}{jitter used in the inversion of the covariance matrix for numerical stability}

\item{settings}{list with argument \code{return.Ki}, to include the inverse covariance matrix in the object for further use (e.g., prediction).
Arguments \code{factr} (default to 1e9) and \code{pgtol} are available to be passed to \code{control} for L-BFGS-B in \code{\link[stats]{optim}}.}
}
\value{
a list which is given the S3 class "\code{homGP}", with elements:
\itemize{
\item \code{theta}: maximum likelihood estimate of the lengthscale parameter(s),
\item \code{g}: maximum likelihood estimate of the nugget variance,
\item \code{trendtype}: either "\code{SK}" if \code{beta0} is given, else "\code{OK}" 
\item \code{beta0}: estimated trend unless given in input,
\item \code{sigma2}:  maximum likelihood estimate of the scale variance,
\item \code{nu2}:  maximum likelihood estimate of the degrees of freedom parameter,
\item \code{ll}: log-likelihood value,
\item \code{X0}, \code{Z0}, \code{Z}, \code{mult}, \code{eps}, \code{covtype}: values given in input,
\item \code{call}: user call of the function
\item \code{used_args}: list with arguments provided in the call
\item \code{nit_opt}, \code{msg}: \code{counts} and \code{msg} returned by \code{\link[stats]{optim}}
\item \code{Ki}, inverse covariance matrix (if \code{return.Ki} is \code{TRUE} in \code{settings})
\item \code{time}: time to train the model, in seconds.

}
}
\description{
Student-t process regression under homoskedastic noise based on maximum likelihood estimation of the 
hyperparameters. This function is enhanced to deal with replicated observations.
}
\details{
The global covariance matrix of the model is parameterized as \code{K = sigma2 * C + g * diag(1/mult)},
with \code{C} the correlation matrix between unique designs, depending on the family of kernel used (see \code{\link[hetGP]{cov_gen}} for available choices).

It is generally recommended to use \code{\link[hetGP]{find_reps}} to pre-process the data, to rescale the inputs to the unit cube and to normalize the outputs.
}
\examples{
##------------------------------------------------------------
## Example 1: Homoskedastic Student-t modeling on the motorcycle data
##------------------------------------------------------------
set.seed(32)

## motorcycle data
library(MASS)
X <- matrix(mcycle$times, ncol = 1)
Z <- mcycle$accel
plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")

noiseControl = list(g_bounds = c(1e-3, 1e4))
model <- mleHomTP(X = X, Z = Z, lower = 0.01, upper = 100, noiseControl = noiseControl)
summary(model)
  
## Display averaged observations
points(model$X0, model$Z0, pch = 20) 
xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
preds <- predict(x = xgrid, object =  model)

## Display mean prediction
lines(xgrid, preds$mean, col = 'red', lwd = 2)
## Display 95\% confidence intervals
lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.05, df = model$nu + nrow(X)), col = 2, lty = 2)
lines(xgrid, preds$mean + sqrt(preds$sd2) * qt(0.95, df = model$nu + nrow(X)), col = 2, lty = 2)
## Display 95\% prediction intervals
lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.05, df = model$nu + nrow(X)), 
  col = 3, lty = 2)
lines(xgrid, preds$mean + sqrt(preds$sd2 + preds$nugs) * qt(0.95, df = model$nu + nrow(X)), 
  col = 3, lty = 2)
}
\references{
M. Binois, Robert B. Gramacy, M. Ludkovski (2018), Practical heteroskedastic Gaussian process modeling for large simulation experiments,
Journal of Computational and Graphical Statistics, 27(4), 808--821.\cr 
Preprint available on arXiv:1611.05902.\cr \cr

A. Shah, A. Wilson, Z. Ghahramani (2014), Student-t processes as alternatives to Gaussian processes, Artificial Intelligence and Statistics, 877--885. \cr \cr

M. Chung, M. Binois, RB Gramacy, DJ Moquin, AP Smith, AM Smith (2019). 
Parameter and Uncertainty Estimation for Dynamical Systems Using Surrogate Stochastic Processes.
SIAM Journal on Scientific Computing, 41(4), 2212-2238.\cr
Preprint available on arXiv:1802.00852.
}
\seealso{
\code{\link[hetGP]{predict.homTP}} for predictions.
\code{summary} and \code{plot} functions are available as well.
}
