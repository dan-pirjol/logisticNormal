#' Computes the g(z,t) function by recursion
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
  if (x0 > 0.5*t) {
    x0 <- x0 - t
    k <- k+1
  }
  
  g0 <- gPoisson(x0,t,Nmax)
  z <- x0
  g <- g0
  
  if (k > 0 )
  {
    for (j in 1:k){
      g <- 2*exp(-0.5*z-0.375*t)- exp(-z-0.5*t)*g
      z <- z + t
    }
  }
  
  
  
  return (g)
}