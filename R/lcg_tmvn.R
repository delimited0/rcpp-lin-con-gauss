get_polytope_constraints <- function(lb, ub, A, mu, L) {
  d <- length(lb)
  
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  
  Amat <- rbind(A, -A)[!inf_idx, ]
  Amat <- Amat %*% L
  
  b <- c(A %*% mu - lb,  -A %*% mu + ub)[!inf_idx]
  
  return(list(A = Amat, b = b))
}

#' @export
rtmvn <- function(n, mu, Sigma, lb, ub, A = NULL, x_init = NULL,
                  mode = "intersection", verbose = FALSE) {
  
  L <- t(chol(Sigma))
  d <- length(mu)
  # A <- rbind(U, -U)
  # b <- c(mu - lb, -mu + ub)
  
  if (is.null(A)) A = diag(d)
    
  pc <- get_polytope_constraints(lb, ub, A, mu, L)
  A <- pc$A
  b <- pc$b  
  
  if (mode == "intersection")
    mode_bool = TRUE
  else if (mode == "union")
    mode_bool = FALSE
  else
    stop("Invalid mode, must be intersection or union")
  
  # if (is.null(x_init)) {
  #   x_init <- ifelse(is.finite(lb), lb + 1e-12, ifelse(is.finite(ub), ub - 1e-12, 0))
  # }
  
  std_samples <- t(ess(n, A, b, x_init, mode_bool, verbose))
  # samples <- std_samples %*% t(U) + matrix(rep(mu, n), nrow = n, byrow = TRUE)
  samples <- std_samples %*% t(L) + matrix(rep(mu, n), nrow = n, byrow = TRUE)
  return(samples)
}

#' @export
pmvn <- function(mu, Sigma, lb, ub, A = NULL, mode = "intersection",
                  n_sub_samples = 10, domain_fraction = .5, n_sub_skip = 1,
                  n_hdr_samples = 10, n_hdr_skip = 1,
                  log = FALSE, n_est = 10) {
  
  L <- t(chol(Sigma))
  d <- length(mu)
  # A <- rbind(U, -U)
  # b <- c(mu - lb, -mu + ub)
  
  if (is.null(A)) A = diag(d)
  
  pc <- get_polytope_constraints(lb, ub, A, mu, L)
  A <- pc$A
  b <- pc$b
  
  if (mode == "intersection")
    mode_bool = TRUE
  else if (mode == "union")
    mode_bool = FALSE
  else
    stop("Invalid mode, must be intersection or union")
  
  ests = rep(NA, n_est)
  for (i in 1:n_est) {
    logprob <- hdr_prob(A, b, mode_bool, 
                        n_sub_samples, domain_fraction, n_sub_skip,
                        n_hdr_samples, n_hdr_skip) 
    ests[i] = sum(logprob)
  }
  
  log_est = matrixStats::logSumExp(ests) - log(n_est)
  pmu = exp(ests) - exp(log_est)
  
  if (log) {
    result = log_est
  }
  else {
    result = exp(log_est)
  }
  attr(result, "error") = sqrt( sum(pmu^2) / (n_est - 1) )
  
  return(result)
}

#' @param A m x d matrix of constraints, 
#' @export
pmvn2 <- function(mu, Sigma, A, lb, ub, mode = "intersection",
                  n_sub_samples = 10, domain_fraction = .5, n_sub_skip = 1,
                  n_hdr_samples = 10, n_hdr_skip = 1, log = FALSE) {
  
  L <- t(chol(Sigma))
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  A <- rbind(A, -A)[!inf_idx, , drop=FALSE] 
  A <- A %*% L
  b <- c(mu - lb, -mu + ub)[!inf_idx]
  
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