func <- Rcpp::cppFunction("arma::mat foo() {
  arma::vec theta(6, arma::fill::zeros);
  arma::vec a1 = {1, 0, 3};
  return a1 * arma::cos(theta);
}", depends = "RcppArmadillo")

reshape_func <- Rcpp::cppFunction("arma::mat baz() {
  arma::vec slices = {1, 2, 3, 4, 6, 7};
  return arma::reshape(slices, slices.size() / 2, 2);
}", depends = "RcppArmadillo")

foo <- lincongauss::rtmvn(1000, 0, 1, -1, 1, 0)


mu <- c(0, 0)
Sigma <- diag(2)
lb <- c(0, 0)
ub <- c(1, 1)
baz <- lincongauss::rtmvn(1000, mu, Sigma, lb, ub, c(.5, .5))