# logisticNormal
R package with functions for the evaluation of the logistic-normal integral using the methods proposed in 

>D.Pirjol, The logistic-normal integral and its generalizations, Journal of Computational and Applied Mathematics 237, 460-469 (2012)

## **Sample usage**

```
# plot the logistic-normal integral 
x <- seq(-5,5, 0.1)
n <- length(x)

d <- c()
for (i in 1:n) d <- c(d,phiLogisticNormal(x[i], 1.0, 10))
plot(d~x, type ="l",main="Logistic-normal integral (t=1)")

y <- seq(-5,5,1)

ex <- c()
for (i in 1:11) ex <- c(ex, phiLNexact(1.0,y[i]))
                        
points(ex~y, type = "p", col="blue")
```

Numerical evaluation comparing with exact results (blue dots)

![test1](https://user-images.githubusercontent.com/60016102/204109670-daa6e5b0-8561-481a-9f80-549608d81698.png)

A more detailed comparison

```
# compare the Poisson summation formula with the exact evaluation (t=1.0 - 50.0)
t <- 1.0
v1 <- c()
v2 <- c()
for (k in 1:50) v1 <- c(v1,phiLNexact(t, k))
for (k in 1:50) v2 <- c(v2,phiLogisticNormal(k*t, t, 10))
v3 <- v1 - v2
vd <- data.frame("exact"=v1, "Poisson" = v2, "diff"=v3)
vd
```

Keeping N=10 terms in the Poisson sum, the error is reasonably small over the entire range of t

<img width="312" alt="test2" src="https://user-images.githubusercontent.com/60016102/204109785-b503943a-b339-4984-831d-ab8b58c1115b.png">




