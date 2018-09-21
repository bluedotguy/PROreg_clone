\name{VarEst}
\alias{VarEst}

\title{
Estimation of the dispersion parameters in a beta-binomial mixed-effects model.
}

\description{
\code{VarEst} function performs a mamximum likelihood estimation of the dispersion parameters in a beta-binomial mixed-effects models given some initial values.
}

\usage{
VarEst(y,m,p,X,Z,u,nRand,nComp,nRandComp,OLDall.sigma,OLDphi,q)
}

\arguments{
\item{y}{ dependent response variable in the model.}
\item{m}{ maximum score number in each beta-binomial observation.}
\item{p}{ estimated values of the probability parameter.}
\item{X}{ model matrix of the fixed effects.}
\item{Z}{ model matrix of the random effects.}
\item{u}{ estimated values of the random-effects.}
\item{nRand}{ number of random effects in the model.}
\item{nComp}{ number of random effects components, i.e., number of different random effects.}
\item{nRandComp}{ the number of random effects in each random component of the model. It must be especified as a vector, where the 'i'th value must correspond with the number of random effects of the 'i'th random component}
\item{OLDall.sigma}{ initial values of the variance parameters of the random effects.}
\item{OLDphi}{ initial value of the dispersion parameter of the beta-binomial distribution.}
\item{q}{ number of independent variables or covariates in the fixed part of the model.}
}

\details{
\code{VarEst} function performs a mamximum likelihood estimation of the dispersion parameters in a beta-binomial mixed-effects models given some initial values.

It uses the \code{multiroot} and \code{uniroot} functions of the \code{rootSolve} R-package.
}

\value{
  \code{VarEst} returns a list of estimates and variance of the dispersion parameters.

(phi=phi,all.sigma=all.sigma,psi=psi,psi.var=psi.var,all.sigma.var=sigma.var)

  \item{phi}{ estimated value of the dispersion parameter in the beta-binomial distribution.}
  \item{all.sigma}{ estimated value of the variances of the random effects.}
  \item{psi}{ estimated value of the logarithm of the dispersion parameter in the beta-binomial distribution. }
  \item{psi.var}{ variance of the estimated value of the logarithm of the dispersion parameter in the beta-binomial distribution.}
  \item{all.sigma.var}{ variance of the estimated value of the variances of the random effects.}
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25

Lee Y. & Nelder J. A. (1996): Hierarchical generalized linear models, \emph{Journal of the Royal Statistical
Society. Series B}, \bold{58}, 619-678

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2018): A beta-binomial mixed-effects model approach for analysing longitudinal discrete and bounded outcomes, \emph{Biometrical Journal}.
}

\keyword{rootSolve}
\seealso{
The \code{\link{multiroot}} and \code{\link{uniroot}} functions of the R-package \code{rootSolve} for the general Newton-Raphson algorithm.
}
