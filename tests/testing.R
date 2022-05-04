A <- matrix(1:6, byrow = TRUE, ncol = 2)
b <- rep(1, 3)

X <- matrix(c(-4, 1, 2, 3, 2, 3, -1, -1, 1), byrow = TRUE, nrow = 3, ncol = 3)

# evaled <- evaluate(A, b, X)

shiftvals <- c(4, -5, -14)
lincongauss::test_update_fix_shift(A, b, .5, 5, shiftvals)


sample_ex <- Rcpp::cppFunction("
  arma::uvec sample_index(const int &size) {
    arma::uvec sequence = arma::linspace<arma::uvec>(0, size-1, size);
    arma::uvec out = Rcpp::RcppArmadillo::sample(sequence, size, false);
    return out; 
  }",
  depends = c("RcppArmadillo"),
  includes = "RcppArmadilloExtensions/sample.h")


sample_ex <- Rcpp::cppFunction("
  IntegerVector sample_dbl(IntegerVector x, int sz, bool rep = false, sugar::probs_t p = R_NilValue)
  {
    return sample(x, sz, rep, p);
  }",
  )

sample2 <- Rcpp::cppFunction("
  int sample_int(int n, int sz) {
    IntegerVector sequence = seq(0, n-1);
    return sample(sequence, sz, false, R_NilValue)[0];
  }")

regspace <- Rcpp::cppFunction("
  arma::uvec regspace(int n, int s) {
    return arma::regspace<arma::uvec>(0, s, n-1);
  }", depends  = "RcppArmadillo")


lincongauss::hdr_prob(A, b, TRUE, .5, n_sub_samples = 16, 1, 16, 1)



# compound symmetric covariance orthant -----------------------------------

d = 30
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
mu = rep(0, d)
lb = rep(0, d)
ub = rep(Inf, d)

lincongauss::pmvn(mu, Sigma, lb, ub, 
                  n_sub_samples = 16, n_hdr_samples = 420,
                  n_sub_skip = 10, n_hdr_skip = 2, n_est = 10)



