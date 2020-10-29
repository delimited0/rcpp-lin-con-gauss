d <- 30
mu <- rep(0, d)
Sigma <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb <- rep(0, d)
ub <- rep(Inf, d)

samples <- lincongauss::rtmvn(1000, mu, Sigma, lb, ub, x_init = rep(10,  d))
plot(samples[, 1:2])

set.seed(0)
lincongauss::pmvn(mu, Sigma, lb, ub, 
                  n_sub_samples = 160, n_sub_skip = 10, 
                  n_hdr_samples = 128, n_hdr_skip = 2)


foo <- Rcpp::cppFunction(
  "
  arma::vec foo() {
    arma::vec x = {1,2,3,4,5,6};
    return x.head(x.n_elem - 1);
  }                         
  ", depends = "RcppArmadillo"
)
