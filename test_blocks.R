devtools::load_all()

g <- 5
mu <- rnorm(g)

n <- 5
data_list <- list()
for (i in 1:g) {
  data_list[[i]] <- rnorm(n, mean = mu[i])
}
names(data_list) <- sprintf("group%s", 1:g)

L <- 5
df_params <- drjacoby::define_params(name = "mu_1", min = -L, max = L, init = 0.1, block = c(1, 6),
                           name = "mu_2", min = -L, max = L, init = 0.1, block = c(2, 6),
                           name = "mu_3", min = -L, max = L, init = 0.1, block = c(3, 6),
                           name = "mu_4", min = -L, max = L, init = 0.1, block = c(4, 6),
                           name = "mu_5", min = -L, max = L, init = 0.1, block = c(5, 6))


ll_block <- function(params, data, misc){
  block <- misc[["block"]]
  if(block == 6){
    out <- sum(dnorm(params[1:5], mean = 0, sd = 1, log = TRUE))
  } else {
    x <- data[[block]]
    out <- sum(dnorm(x, mean = params[block], sd = 1, log = TRUE))
  }
  return(out)
}

lp <- function(params, misc = NULL){
  return(0)
}

mcmc <- run_dj11(data = data_list,
                df_params = df_params,
                loglike = ll_block,
                logprior = lp,
                burnin = 1e3L,
                samples = 1e5L,
                chains = 1)

output_sub <- subset(mcmc$output, phase == "sampling")
for (i in 1:5) {
  # get posterior draws
  mu_draws <- output_sub[[sprintf("mu_%s", i)]]

  # get analytical solution for this group
  x <- seq(-L, L, l = 1001)
  m <- mean(data_list[[i]])
  fx <- dnorm(x, mean = m * n/(n + 1), sd = sqrt(1/(n + 1)))

  # get analytical solution if no multi-level model
  fx2 <- dnorm(x, mean = m, sd = sqrt(1/n))

  # overlay plots
  hist(mu_draws, breaks = seq(-L, L, l = 1001), probability = TRUE, col = "black",
       main = sprintf("mu_%s", i))
  lines(x, fx, col = 2, lwd = 4)
  lines(x, fx2, col = 3, lwd = 4)
}


# extract sampling draws
output_sub <- subset(mcmc$output, phase == "sampling")

for (i in 1:5) {
  # get posterior draws
  mu_draws <- output_sub[[sprintf("mu_%s", i)]]

  # get analytical solution for this group
  x <- seq(-L, L, l = 1001)
  m <- mean(data_list[[i]])
  fx <- dnorm(x, mean = m * n/(n + 1), sd = sqrt(1/(n + 1)))

  # get analytical solution if no multi-level model
  fx2 <- dnorm(x, mean = m, sd = sqrt(1/n))

  # overlay plots
  hist(mu_draws, breaks = seq(-L, L, l = 1001), probability = TRUE, col = "black",
       main = sprintf("mu_%s", i))
  lines(x, fx, col = 2, lwd = 4)
  lines(x, fx2, col = 3, lwd = 4)
}
