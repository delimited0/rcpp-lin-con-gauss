# was not aware of this implementation on CRAN
# now there's three of them
# my implementation seems to be slightly (?) faster
# not much of a difference

n = 1000
d = 64

mu = rep(0, d)
Sigma = .5*diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)
x_init = rep(1, d)

A = diag(d)
b = -lb
constr = lincongauss:::get_polytope_constraints(lb, ub, diag(d), mu, t(chol(Sigma)))

microbenchmark::microbenchmark(
  lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = x_init),
  linconGaussR::linconGauss(n, A, b, Sigma, mu, x_init = x_init,
                            nskp = 0),
  LinConGauss::LinESSFast(A, b, N = n, x0 = x_init, nskip = 0),
  times = 20  
)

my_samples = lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = x_init)
their_samples = linconGaussR::linconGauss(n, A, b, Sigma, mu, 
                                          x_init = x_init, nskp = 0)
other_samples = LinConGauss::LinESSFast(constr$A, constr$b, N = n, x0 = x_init, nskip = 0)

plot(my_samples[, c(10, 30)])
points(their_samples[, c(10, 30)], col = "blue")
points(other_samples[, c(10, 30)], col = "red")

# probability estimation
sub_ret <- LinConGauss::SubsetSimFast(A, b, n, .5, nskip = 0)
hdr_ret <- LinConGauss::HDR_algo(A, b, c(sub_ret$shift_seq[1:(length(sub_ret$shift_seq)-1)],0), n )

lincongauss::pmvn(mu, Sigma, lb, ub, )
