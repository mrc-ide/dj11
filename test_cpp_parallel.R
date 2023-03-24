devtools::load_all()

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

cpp11::cpp_source("test_cpp11_ll.cpp")
t2 <- run_dj11(
  data = data_list,
  df_params = df_params,
  loglike = "cpp11_ll",
  logprior = lp,
  burnin = 5000L,
  samples = 5000L,
  chains = 1
)

cl <- parallel::makeCluster(4)
add_cl <- parallel::clusterEvalQ(cl = cl, cpp11::cpp_source("test_cpp11_ll.cpp"))
t2 <- run_dj11(
    data = data_list,
    df_params = df_params,
    loglike = "cpp11_ll",
    logprior = lp,
    burnin = 5000L,
    samples = 5000L,
    chains = 4,
    cluster = cl
  )
parallel::stopCluster(cl)
