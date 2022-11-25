#' Computes g(z,t) using the theta2 Poisson series (27)
#' 
#' Uses the Poisson series (27)
#' @param Nmax, the truncation order of the series. Defaults to 5.
#' @export
#' @examples 
#' gPoisson(0.0,1.0) should return 0.901925

gPoisson <- function(x,t,Nmax=5){
  PI <- 4.0*atan(1.0) 
  
  q <- exp(-0.5*t)
  q1 <- exp(-2.0*PI*PI/t)
  
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