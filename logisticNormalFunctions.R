#' Logistic normal integral
#'
#' @param z 
#' @param t 
#' @param Nmax 
#'
#' @return Returns a real value
#' @export
#'
#' @examples
phiLogisticNormal <- function(z,t,Nmax=5){
  
  phi <- 1/2*exp(-z/2+t/8)*gRec(z-t/2,t,Nmax)
  return(phi)
}


#-------------------------------------------------------------------------
#' Computes g(z,t) using the theta2 Poisson series (27)
#' 
#' @param Nmax, the truncation order of the series. Defaults to 5.
#' @export
#' @examples 
#' gPoisson(0.0,1.0) should return 0.901925

gPoisson <- function(x,t,Nmax=5){
  PI <- 4.0*atan(1.0) 
  
  q <- exp(-0.5*t)
  q1 <- exp(-2.0*PI*PI/t)
  
  x <- abs(x)
  
  # S1 = theta2(i*z/2, exp(-t/2)) 
  S1 <- exp(-1/2*x)*q^0.25 
  S2 <- 0.0
  S3 <- 0.0
  for (k in 1:Nmax){
    fact <- q^((k+0.5)^2)*exp(-(k+0.5)*x) + q^((-k+0.5)^2)*exp((k-0.5)*x)
    S1 <- S1 + fact
  }

  
  for (j in 0:2*Nmax-1){
    k <- j-Nmax+1
    fact <- q1^(k*k-0.25)
    den <- 1.0 - q1^(2*k-1)
    S2 <- S2 + (-1)^k*cos(2*PI/t*(k-0.5)*x)*fact/den
  }
  
  
  # S3 is called S4 in the paper
  for (k in -Nmax:Nmax){
    fact = q^(k*k+k)
    den = 1.0 + q^(2*k)
    S3 <- S3 + exp(k*x)*fact/den
  }
  
  y <- (4.0*PI/t*S2 + 2.0*S3)/(S1)
  
  return(y)
  
}

#------------------------------------------------------------------------------------
#' Computes the g(z,t) function by recursion
#' Does not work when x is a vector
#' 
#' The function g(z,t) is computed by recursion from its values in the primitive cell (-1/2 t,1/2 t)
#' 
#' @param Nmax is the truncation order in the evaluation of the Appell sums. Defaults to 5
#' @export
#' @examples 
#' gRec(0.0,1.0) 

gRec <- function(x,t,Nmax=5) {
  x <- abs(x)
  k <- floor(x/t)
  
  x0 <- x - k*t

  if (x0 > (0.5*t)) {
    x0 <- x0 - t
    k <- k+1
  }
  
  g0 <- gPoisson(x0,t,Nmax)
  z <- x0
  g <- g0
  
  if (k > 0 )
  {
    for (j in 1:k){
      g <- 2*exp(-1/2*z-3/8*t)- exp(-z-1/2*t)*g
      z <- z + t
    }
  }
  
  
  
  return (g)
}


#-------------------------------------------------------------------------------
#' exact evaluation of phi(k*t,t) with integer k, using (20) in JCAM

phiLNexact <- function(t, k){
  
  ks <- k
  k <- abs(k)
  #  print (k)
  
  if (k==0) sum = 1/2
  else if (k==1) sum = 1/2*exp(-1/2*t)
  else
  {
    sum <- 0
    for (j in 1:(k-1)){
      term <- (-1)^(j+1)*exp(j*(1/2*j-k)*t)
      sum <- sum + term
      #    print(j)
    }
    lterm <- 1/2*(-1)^(k-1)*exp(-1/2*k^2*t)
    sum <- sum + lterm
  }
  
  
  sumn <- 1 - sum
  
  if (ks>=0) return (sum)
  if (ks <0) return (sumn)
  
  
}
#-------------------------------------------------
