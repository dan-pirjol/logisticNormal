#' exact evaluation of phi(k*t,t) with integer k, using (20) in JCAM
#'
#' @param t 
#' @param k 
phiLNexact <- function(t, k){
  
  ks <- k
  k <- abs(k)
  print (k)
  
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

phiLNexact(1.0,0)
phiLNexact(1.0,1)
1/2*exp(-1/2)

phiLNexact(1.0,-1)
1- 1/2*exp(-1/2)

phiLNexact(1.0,2)
exp(-1.5) - 1/2*exp(-2.0)

phiv <- c()

for (i in 1:20) phiv <- c(phiv,phiLNexact(1.0,i))

plot(1:20,log(phiv))
