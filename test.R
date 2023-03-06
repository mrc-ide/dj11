
devtools::load_all()

true_theta <- c(0, 1)
data <- rnorm(100, true_theta[1], true_theta[2])

data_list <- list(x = data)
df_params <- drjacoby::define_params(name = "mu", min = -Inf, max = Inf, init = 1, block = 1,
                                     name = "sigma", min = 0, max = Inf, init = 1, block = 1)

ll <- function(params, data, block){
  # calculate log-probability of data
  sum(dnorm(data$x, mean = params["mu"], sd = params["sigma"], log = TRUE))
}

lp <- function(params, misc = NULL){
  return(0)
}

o1 <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  burnin = 5000L,
  samples = 5000L
)

o1$acceptance / 10000
plot(o1$out[,1], t = "l")
plot(o1$out[,2], t = "l")

o2 <- drjacoby::run_mcmc(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  chains = 1,
  burnin = 5000,
  samples = 5000,
  silent = TRUE
)

bm <- microbenchmark::microbenchmark(
  new = run_dj11(data = data_list, df_params = df_params, loglike = ll, logprior = lp, burnin = 10000L, samples = 10000L),
  old = drjacoby::run_mcmc(
    data = data_list,
    df_params = df_params,
    loglike = ll,
    logprior = lp,
    chains = 1,
    burnin = 10000,
    samples = 10000,
    silent = TRUE
  ),
  times = 10
)
bm


