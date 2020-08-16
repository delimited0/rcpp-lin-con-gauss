get_polytope_constraints <- function(lb, ub, mu, U) {
  d <- length(lb)
  
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  
  A <- rbind(diag(d), -diag(d))[!inf_idx, ]
  A <- A %*% U
  
  b <- c(mu - lb, -mu + ub)[!inf_idx]
  
  return(list(A = A, b = b))
}

#' @export
rtmvn <- function(n, mu, Sigma, lb, ub, x_init = NULL,
                  mode = "intersection", verbose = FALSE) {
  # U <- chol(Sigma)
  L <- t(chol(Sigma))
  # d <- length(mu)
  # A <- rbind(U, -U)
  # b <- c(mu - lb, -mu + ub)
  
  pc <- get_polytope_constraints(lb, ub, mu, L)
  A <- pc$A
  b <- pc$b
  
  if (mode == "intersection")
    mode_bool = TRUE
  else if (mode == "union")
    mode_bool = FALSE
  else
    stop("Invalid mode, must be intersection or union")
  
  if (is.null(x_init)) {
    x_init <- tmvtnorm::rtmvnorm(1, mean = mu, sigma = Sigma, lower = lb)
  }
  
  std_samples <- t(ess(n, A, b, x_init, mode_bool, verbose))
  # samples <- std_samples %*% t(U) + matrix(rep(mu, n), nrow = n, byrow = TRUE)
  samples <- std_samples %*% t(L) + matrix(rep(mu, n), nrow = n, byrow = TRUE)
  return(samples)
}

#' @export
ptmvn <- function(mu, Sigma, lb, ub, mode = "intersection",
                  n_sub_samples = 10, domain_fraction = .5, n_sub_skip = 0,
                  n_hdr_samples = 10, n_hdr_skip = 0,
                  log = FALSE) {
  L <- t(chol(Sigma))
  # d <- length(mu)
  # A <- rbind(U, -U)
  # b <- c(mu - lb, -mu + ub)
  
  pc <- get_polytope_constraints(lb, ub, mu, L)
  A <- pc$A
  b <- pc$b
  
  if (mode == "intersection")
    mode_bool = TRUE
  else if (mode == "union")
    mode_bool = FALSE
  else
    stop("Invalid mode, must be intersection or union")
    
  logprob <- hdr_prob(A, b, mode_bool, 
                      n_sub_samples, domain_fraction, n_sub_skip,
                      n_hdr_samples, n_hdr_skip) 
  if (log)
    return(sum(logprob))
  else
    return(prod(exp(logprob)))
}