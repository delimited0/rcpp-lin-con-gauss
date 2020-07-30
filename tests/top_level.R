d <- 1
mu <- rep(0, d)
Sigma <- diag(d)
lb <- rep(0, d)
ub <- rep(Inf, d)

lincongauss::ptmvn(mu, Sigma, lb, ub, 
                   n_sub_samples = 16, n_sub_skip = 3, 
                   n_hdr_samples = 512, n_hdr_skip = 9)