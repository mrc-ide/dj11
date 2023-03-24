
devtools::load_all()
par(mfrow = c(1, 1))

# Create data
true_theta <- c(0, 1)
data <- rnorm(100, true_theta[1], true_theta[2])

# Create DrJacoby style inputs
data_list <- list(x = data)
df_params <- drjacoby::define_params(name = "mu", min = -Inf, max = Inf, init = 0.1, block = 1,
                                     name = "sigma", min = 0, max = Inf, init = 1, block = 1)

# Likelihood and prior
ll <- function(params, data, block){
  # calculate log-probability of data
  sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
}
lp <- function(params, misc = NULL){
  return(0)
}

# Example run
o1r <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  burnin = 5000L,
  samples = 5000L,
  chains = 2
)
head(o1r$output)
class(o1r) <- 'drjacoby_output'
drjacoby::plot_par(o1r, "sigma")

o1r$rung_index
o1r$swap_acceptance / 10000

plot(o1r$output$mu, t= "l")
plot(o1r$output$sigma, t=  "l")

# Load c++ log likelihoods
cpp11::cpp_source("test_cpp11_ll.cpp")
Rcpp::sourceCpp("test_rcpp_ll.cpp")

o1c <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = cpp11_ll,
  logprior = lp,
  burnin = 5000L,
  samples = 5000L
)

o2r <- drjacoby::run_mcmc(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  chains = 1,
  burnin = 5000,
  samples = 5000,
  silent = TRUE
)

o2c <- drjacoby::run_mcmc(
  data = data_list,
  df_params = df_params,
  loglike = "loglike",
  logprior = lp,
  chains = 1,
  burnin = 5000,
  samples = 5000,
  silent = TRUE
)

par(mfrow = c(1, 3))
plot(density(o1r$output$mu), ylim = c(0, 5), main = "mu")
lines(density(o2r$output$mu), col = "red")
lines(density(o1c$ou$mu), col = "blue")
lines(density(o2c$output$mu), col = "orange")

plot(density(o1r$output$sigma), main = "sigma")
lines(density(o2r$output$sigma), col = "red")
lines(density(o1c$ou$sigma), col = "blue")
lines(density(o2c$output$sigma), col = "orange")

bm <- microbenchmark::microbenchmark(
  new_r = run_dj11(
    data = data_list,
    df_params = df_params,
    loglike = ll,
    logprior = lp,
    burnin = 5000L,
    samples = 5000L
  ),
  new_c = run_dj11(
    data = data_list,
    df_params = df_params,
    loglike = cpp11_ll,
    logprior = lp,
    burnin = 5000L,
    samples = 5000L
  ),
  old_r = drjacoby::run_mcmc(
    data = data_list,
    df_params = df_params,
    loglike = ll,
    logprior = lp,
    chains = 1,
    burnin = 5000,
    samples = 5000,
    silent = TRUE
  ),
  old_c = drjacoby::run_mcmc(
    data = data_list,
    df_params = df_params,
    loglike = "loglike",
    logprior = lp,
    chains = 1,
    burnin = 5000,
    samples = 5000,
    silent = TRUE
  ),
  times = 10
)
bm
plot(bm)

