# logisticNormal
R package with functions for the evaluation of the logistic-normal integral using the methods proposed in 

>D.Pirjol, The logistic-normal integral and its generalizations, Journal of Computational and Applied Mathematics 237, 460-469 (2012)

## **Sample usage**

# plot the logistic-normal integral
x <- seq(-5,5, 0.1)
n <- length(x)

d <- c()
for (i in 1:n) d <- c(d,phiLogisticNormal(x[i], 1.0, 10))

plot(d~x, type ="l",main="Logistic-normal integral (t=1)")





