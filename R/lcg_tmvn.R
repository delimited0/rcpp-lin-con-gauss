#' @export
rtmvn <- function(n, mu, Sigma, lb, ub, x_init = NULL,
                  mode = "intersection", verbose = FALSE) {
  U <- chol(Sigma)
  d <- length(mu)
  A <- rbind(U, -U)
  b <- c(mu - lb, -mu + ub)
  
  if (mode == "intersection")
    mode_bool = TRUE
  else if (mode == "union")
    mode_bool = FALSE
  else
    stop("Invalid mode, must be intersection or union")
  
  if (is.null(x_init)) {
    tmvtnorm::rtmvnorm(d, mean = mu, sigma = Sigma, lower = lb)
  }
  
  std_samples <- ess(n, A, b, x_init, mode_bool, verbose)
  samples <- std_samples %*% t(U) + matrix(mu, nrow = d, byrow = TRUE)
  return(samples)
}