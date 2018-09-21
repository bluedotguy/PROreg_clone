\name{EffectsEst.BBNR}
\alias{EffectsEst.BBNR}

\title{
Specific Newton-Raphson based algorithm for the estimation of the fixed and random effects in a beta-binomial mixed effects model.
}

\description{
Specific Newton-Raphson based algorithm for the estimation of the fixed and random effects in a beta-binomial mixed effects model.
}

\usage{
EffectsEst.BBNR(y,m,beta,u,p,phi,D.,X,Z,maxiter)
}

\arguments{
\item{y}{dependent response variable in the model.}
\item{m}{maximum score number in each beta-binomial observation.}
\item{beta}{initial values of the fixed-effects.}
\item{u}{initial values of the random-effects.}
\item{p}{initial values of the probability parameter of the beta-binomial distribution.}
\item{phi}{estimated value of the dispersion parameter of the beta-binomial distribution.}
\item{D.}{estimated value of the variance-covariance matrix of the random effects.}
\item{X}{model matrix of the fixed effects.}
\item{Z}{model matrix of the random effects.}
\item{maxiter}{maximum number of iterations for the algorithm.}
}

\details{
\code{EffectsEst.BBNR} function performs a Newton-Raphson based mamximum likelihood estimation of the fixed and random effects in a beta-binomial mixed-effects models given some initial values.

See Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2018): A beta-binomial mixed-effects model approach for analysing longitudinal discrete and bounded outcomes, \emph{Biometrical Journal} for more information.
}

\value{
  \code{EffectsEst.BBNR} returns a list of the estimates and variances of the fixed and random effects.

  \item{fixed.est}{ estimated value of the fixed coefficients of the regression.}
  \item{random.est}{ estimated value of the random coefficients of the regression.}
  \item{vcov.fixed}{ variance-covariance matrix of the estiamtion of the fixed-effects.}
  \item{var.random}{ variance of the estimation of the random effects.}
  \item{iter.fixrand}{ numbero of iterations in the algorithm.}
  \item{conv.fixrand}{ convergence og the algorithm.}
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25.

Lee Y. & Nelder J. A. (1996): Hierarchical generalized linear models, \emph{Journal of the Royal Statistical
Society. Series B}, \bold{58}, 619-678.

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2018): A beta-binomial mixed-effects model approach for analysing longitudinal discrete and bounded outcomes, to appear in \emph{Biometrical Journal}.
}
