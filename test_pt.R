
devtools::load_all()
par(mfrow = c(1, 1))

# Create data
true_theta <- c(0, 1)
data <- rnorm(100, true_theta[1], true_theta[2])

# Create DrJacoby style inputs
data_list <- list(x = rnorm(100, mean = 10))
# define parameters dataframe
df_params <- drjacoby::define_params(
  name = "alpha", min = -10, max = 10, init = 5, block = 1,
  name = "beta", min = 0, max = 10, init = 5, block = 1,
  name = "epsilon", min = -Inf, max = Inf, init = 0, block = 1)

# Likelihood and prior
ll <- function(params, data, block){
  mean = params["alpha"] * params["alpha"] * params["beta"] + params["epsilon"];
  # calculate log-probability of data
  sum(dnorm(data$x, mean = mean, sd = 1, log = TRUE))
}
lp <- function(params, misc = NULL){
  -log(20.0) - log(10.0) + dnorm(params["epsilon"], 0.0, 1.0, TRUE)
}

# Example run
pt1 <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  burnin = 10000L,
  samples = 10000L,
  n_rungs = 1
)

plot(pt1$out[15000:20000,1:2], xlab = "alpha", ylab = "beta")

pt2 <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = ll,
  logprior = lp,
  burnin = 1e3,
  samples = 1e4,
  n_rungs = 20
)
pt2$swap_acceptance / 10000
plot(pt2$out[,1:2], xlab = "alpha", ylab = "beta")

microbenchmark::microbenchmark(
  pt2 = run_dj11(
    data = data_list,
    df_params = df_params,
    loglike = ll,
    logprior = lp,
    burnin = 1e3,
    samples = 1e4,
    n_rungs = 20
  ),
  times = 10
)
