library("tidyverse")
library("mgcov")
library("mvtnorm")

# Set Up Parameters -------------------------------------------------------
n <- 1000

S0 <- diag(3)
S1 <- matrix(c(1, .1, .1,
               .1, 1, .1, 
               .1, .1, 1), byrow = TRUE, nrow = 3)
S2 <- matrix(c(1, .2, .2,
               .2, 1, .2, 
               .2, .2, 1), byrow = TRUE, nrow = 3)

X0 <- rmvnorm(n, sigma = S0)
X1 <- rmvnorm(n, sigma = S1)
X2 <- rmvnorm(n, sigma = S2)

# Test of Different Data/Initial Matrix -----------------------------------
covchaud(S0, X0)
covchaud(S0, X1)
covchaud(S0, X2)

covchaud(S1, X0)
covchaud(S1, X1)
covchaud(S1, X2)

covchaud(S2, X0)
covchaud(S2, X1)
covchaud(S2, X2)

covchaud
