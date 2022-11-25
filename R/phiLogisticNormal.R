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