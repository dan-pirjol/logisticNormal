# logisticNormal
R functions for the evaluation of the logistic-normal integral using the methods proposed in 

>D.Pirjol, The logistic-normal integral and its generalizations, Journal of Computational and Applied Mathematics 237, 460-469 (2012)

In this paper the logistic-normal integral is defined as $$\varphi(z,t)=\int_{-\infty}^\infty \frac{1}{1+e^x} e^{-\frac{1}{2t}(x-z)^2} \frac{dx}{\sqrt{2\pi t}}$$

The function $\varphi(z,t)$ has the properties:

* Its values for positive and negative argument $z$ are related as $\varphi(z,t) + \varphi(-z,t)=1$

* For any $z\in \mathbb{R}$, we have $\varphi(z+t,t)= e^{-z-\frac12 t} (1 - \varphi(z,t))$

The recursion relation implies that knowledge of $\varphi(z,t)$ in an interval in $z$ of length $t$ is sufficient to determine it anywhere else.
In particular, using $\varphi(0,t)=\frac12$ (which follows by taking $z=0$ in the first relation), gives exact evaluations at $z=kt$ with $k\in \mathbb{Z}$. The recursion relation gives also a Poisson summation formula for $\varphi(z,t)$. 

$\varphi(z,t)$ is related to the Mordell integral. See:

>D. Pirjol, Addendum: The logistic-normal integral and its generalizations, Journal of Computational and Applied Mathematics 237, 460-469 (2013)

The function **phiLogisticNormal(z,t,Nmax)** evaluates $\varphi(z,t)$ by Poisson summation truncating the sum to $Nmax$ terms. 

The function **phiLNexact(t,k)** with $k\in \mathbb{Z}$ gives the exact value of $\varphi(kt, t)$ 

The function **phiLinearInterp(x,t)** computes $\varphi(x,t)$ by linear interpolation from the exact values at the points $(kt,(k+1)t)$ bracketing $x$

**Relation to logistic-normal random variables**

The logistic-normal random variable $X \sim logitnorm(\mu,\sigma)$ is defined as $X=\frac{1}{1+e^{-Z}}$ with $Z\sim N(\mu,\sigma)$. The logistic-normal integral is related to the expectation of such a random variable as $$\mathbb{E}[X]=\varphi(-\mu,\sigma^2)$$

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

Keeping N=10 terms in the Poisson sum, the error is reasonably small for t sufficiently large (above 0.5), and all x.

<img width="312" alt="test2" src="https://user-images.githubusercontent.com/60016102/204109785-b503943a-b339-4984-831d-ab8b58c1115b.png">

For small t, the exact evaluation gives the values of $\varphi(z,t)$ on a fine grid with step $t$. Its values at other points can be obtained by interpolation. In the example below the blue dots show the exact values for $\varphi(k t, t)$ with $t=0.1$.

<img width="697" alt="varphit0p1" src="https://user-images.githubusercontent.com/60016102/204840327-ea06d9eb-72ba-4068-8590-438c3b5df258.png">



