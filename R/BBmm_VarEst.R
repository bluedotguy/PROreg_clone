VarEst <- function(y,m,p,X,Z,u,nRand,nComp,nRandComp,OLDall.sigma,OLDphi,q){

  # Number of observations
  nObs <- length(y)

  s <- c(p*(1-p))

  # Defining the variance of the random effects
  d <- d. <- NULL
  for (i in 1:nComp){
    d <- c(d,rep(OLDall.sigma[i]^2,nRandComp[i]))
    d. <- c(d.,rep(1/(OLDall.sigma[i])^2,nRandComp[i]))
  }
  D <- diag(d)
  D. <- diag(d.)

  ### ESTIMATION OF PSI
  function.psi <- function(psi){

    # Compute j and v
    j <- NULL
    v <- NULL
    hpsi <- NULL
    for (l in 1:nObs){
      j1 <- 0
      j2 <- 0
      v1 <- 0
      v2 <- 0
      hpsi1 <- 0
      hpsi2 <- 0
      hpsi3 <- 0

      if (y[l]==0){
        j1 <- 0
        v1 <- 0
        hpsi1 <- 0
      }else{
        for (k in 0:(y[l]-1)){
          j1 <- j1+k*exp(psi)/(p[l]+k*exp(psi))^3
          v1 <- v1+1/(p[l]+k*exp(psi))^2
          hpsi1 <- hpsi1+k*exp(psi)/(p[l]+k*exp(psi))
        }
      }
      if (y[l]==m[l]){
        j2 <- 0
        v2 <- 0
        hpsi2 <- 0
      }else{
        for (k in 0:(m[l]-y[l]-1)){
          j2 <- j2+k*exp(psi)/(1-p[l]+k*exp(psi))^3
          v2 <- v2+1/(1-p[l]+k*exp(psi))^2
          hpsi2 <- hpsi2+k*exp(psi)/(1-p[l]+k*exp(psi))
        }
      }
      for (k in 0:(m[l]-1)){
        hpsi3 <- hpsi3+k*exp(psi)/(1+k*exp(psi))
      }
      j <- c(j,j1+j2)
      v <- c(v,v1+v2)
      hpsi <- c(hpsi,hpsi1+hpsi2-hpsi3)
    }

    SVS <- diag(s*v*s)
    SJS <- diag(s*j*s)

    # Compute H (Expected Hessian Matrix)

    H1 <- cbind(t(X)%*%SVS%*%X,t(X)%*%SVS%*%Z)
    H1psi <- cbind(t(X)%*%SJS%*%X,t(X)%*%SJS%*%Z)
    H2 <- cbind(t(Z)%*%SVS%*%X,t(Z)%*%SVS%*%Z+D.)
    H2psi <- cbind(t(Z)%*%SJS%*%X,t(Z)%*%SJS%*%Z)
    H <- rbind(H1,H2)
    Hpsi <- rbind(H1psi,H2psi)

    trace <- sum(diag(solve(H)%*%Hpsi))

    out <- sum(hpsi)+trace

    return(out)

  }

  if ((function.psi(-10)*function.psi(4))>0){
    psi <- runif(1,-10,4)
    phi <- exp(psi)
    cat("-Phi estimation has not converged-\n")

  } else{
    mle.psi <- uniroot(function.psi,lower=-10,upper=4,tol=0.01)
    psi <- mle.psi$root
    phi <- exp(psi)

  }

  ##########

  # ESTIMATION OF SIGMA

  # Compute v because it is the same for all sigma score equation
  v <- NULL
  for (l in 1:nObs){
    v1 <- 0
    v2 <- 0
    if (y[l]==0){
      v1 <- 0
    }else{
      for (k in 0:(y[l]-1)){
        v1 <- v1+1/(p[l]+k*phi)^2
      }
    }
    if (y[l]==m[l]){
      v2 <- 0
    }else{
      for (k in 0:(m[l]-y[l]-1)){
        v2 <- v2+1/(1-p[l]+k*phi)^2
      }
    }
    v <- c(v,v1+v2)
  }

  # We get the score equation terms where are the same for all sigma
  SVS <- diag(s*v*s)
  K. <- t(Z)%*%SVS%*%Z-t(Z)%*%SVS%*%X%*%solve(t(X)%*%SVS%*%X)%*%t(X)%*%SVS%*%Z


  function.sigma <- function(sigma){

    nRand.previous <- 0
    d. <- NULL
    out <- NULL

    for (i in 1:nComp){
      d. <- c(d.,rep(1/(sigma[i])^2,nRandComp[i]))
    }

    # Compute the trace
    D. <- diag(c(d.))
    K <- solve(K.+D.)
    trace <- diag(K)

    for(i in 1:nComp){
      out <- c(out,-nRandComp[i]*(sigma[i])^2+t(u[seq(nRand.previous+1,nRand.previous+nRandComp[i])])%*%u[seq(nRand.previous+1,nRand.previous+nRandComp[i])]+sum(trace[seq(nRand.previous+1,nRand.previous+nRandComp[i])]))

      nRand.previous <- nRand.previous+nRandComp[i]
    }
    out
  }

  mle.sigma <- multiroot(function.sigma,start=OLDall.sigma,atol=0.01,ctol=0.01,rtol=0.01)
  all.sigma <- mle.sigma$root

  # The derivate in the point, the variance:
  psi.var <- -1/grad(function.psi,psi)
  sigma.var <- -1/grad(function.sigma,all.sigma)

  out <- list(phi=phi,all.sigma=all.sigma,psi=psi,psi.var=psi.var,all.sigma.var=sigma.var)
  return(out)
}
