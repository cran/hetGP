% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim.R
\name{deriv_crit_EI}
\alias{deriv_crit_EI}
\title{Derivative of EI criterion for GP models}
\usage{
deriv_crit_EI(x, model, cst = NULL, preds = NULL)
}
\arguments{
\item{x}{matrix for the new design (size 1 x d)}

\item{model}{\code{homGP} or \code{hetGP} model}

\item{cst}{threshold for contour criteria}

\item{preds}{pre-computed preds for contour criteria}
}
\description{
Derivative of EI criterion for GP models
}
\references{
Ginsbourger, D. Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables Ecole Nationale Superieure des Mines de Saint-Etienne, Ecole Nationale Superieure des Mines de Saint-Etienne, 2009 \cr \cr
Roustant, O., Ginsbourger, D., DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization, Journal of Statistical Software, 2012
}
\seealso{
\code{\link[hetGP]{crit_EI}} for the criterion
}
