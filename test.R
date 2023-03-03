
devtools::load_all()

true_theta <- c(0, 1)
data <- rnorm(100, true_theta[1], true_theta[2])
ll <- function(theta, data, block){
  sum(dnorm(data, theta[1], theta[2], log = TRUE))
}
lp <- function(theta){
  return(0)
}
theta <- c(0.5, 2)
blocks <- c(1L, 1L)
n_unique_blocks = length(unique(blocks))

system.time({
  o1 <- mcmc(theta, c(0L, 2L),  c(-Inf, 0),  c(Inf, Inf), blocks, n_unique_blocks, data, 5000L, 5000L, ll, lp)
})

o1$acceptance / 10000

o1$proposal_sd

head(o1$out)
head(o1$log_likelihood)
plot(o1$out[,1], t = "l")
plot(o1$out[,2], t = "l")

head(o1$log_likelihood)


# Test against dr J
data_list <- list(x = data)
df_params <- drjacoby::define_params(name = "mu", min = -Inf, max = Inf,
                           name = "sigma", min = 0, max = Inf)
r_loglike <- function(params, data, misc) {

  # extract parameter values
  mu <- params["mu"]
  sigma <- params["sigma"]

  # calculate log-probability of data
  ret <- sum(dnorm(data$x, mean = mu, sd = sigma, log = TRUE))

  # return
  return(ret)
}
r_logprior <- function(params, misc) {
  # return
  return(0)
}

drjacoby::run_mcmc(
  data = data_list,
  df_params = df_params,
  loglike = r_loglike,
  logprior = r_logprior,
  chains = 1,
  burnin = 5000,
  samples = 5000,
  silent = TRUE
  )

bm <- microbenchmark::microbenchmark(
  new = mcmc(theta, c(0L, 2L),  c(-Inf, 0),  c(Inf, Inf), blocks, n_unique_blocks, data, 5000L, 5000L, ll, lp),
  old = drjacoby::run_mcmc(
    data = data_list,
    df_params = df_params,
    loglike = r_loglike,
    logprior = r_logprior,
    chains = 1,
    burnin = 5000,
    samples = 5000,
    silent = TRUE
  ),
  times = 10
)
bm


