# function computing varphi1(z,t) by central finite differences

varphi1 <- function(z,t,Nmax=12,eps=10^-6){
  
  fmid <- phiLogisticNormal(z, t, Nmax)
  fup <- phiLogisticNormal(z+eps, t, Nmax)
  fdown <- phiLogisticNormal(z-eps, t, Nmax)
  
  phid <- (fup-fdown)/(2*eps)
  phi1 <- t*phid + z*fmid
  return (phi1)
}

# function computing E[X^2] for X ~ logisticNorm(mu,sig)

momX2 <- function(mu,sig,Nmax=12,eps=10^-6){
  
  aux <- -1/sig^2*varphi1(-mu-sig^2,sig^2,Nmax,eps) 
  aux <- aux - (1 + mu/sig^2)*phiLogisticNormal(-mu-sig^2,sig^2,Nmax)
  aux <- exp(mu+0.5*sig^2)*aux
  
  return(aux)
}

# function computing the first two moments of X ~ logisticNorm(mu,sig)

momentsLogisticNorm <- function(mu,sig,Nmax=12){
  
  mean <- phiLogisticNormal(-mu,sig^2,Nmax)
  mom2 <- momX2(mu,sig,Nmax)
  var <- mom2 - mean^2
  
  meanvar <- c(mean,var)
  return (meanvar)
}

#momentsLogisticNorm(0.6,0.5) returns [1] 0.6381209 0.0120817


