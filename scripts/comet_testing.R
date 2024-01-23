devtools::install_github("rakheon/mgcov", force = TRUE)

library("mgcov")
library("mvtnorm")

p <-10; n<- 100
Sigma = diag(p); diag(Sigma[-p,-1])=0.5; diag(Sigma[-1,-p])=0.5
dat <- rmvnorm(n,mean=rep(0,ncol(Sigma)),sigma=Sigma)
res_uni = COmet(dat, lambda = seq(0,1,0.01))
res_uni$cov_list[[which.min(res_uni$aic)]]; res_uni$cov_list[[which.min(res_uni$bic)]]
res_ada = COmet(dat, mul = 3)
res_ada$cov_list[[which.min(res_ada$aic)]]; res_ada$cov_list[[which.min(res_ada$bic)]]

res_uni$cov_list |> str()

S_comet <- COmet(dat, lambda = .2)$cov_list[[1]]
S_spcov <- spcov(cov(dat), cov(dat), step.size = 100, lambda = .2)$Sigma
