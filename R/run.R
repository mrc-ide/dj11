run_dj11 <- function(data, df_params, loglike, logprior, burnin, samples, target_acceptance = 0.44, misc = list(),
                     n_rungs = 1L, swap = TRUE, chains = 1){

  theta_init <- unlist(df_params$init)
  theta_names <- unlist(df_params$name)
  theta_min <-  unlist(df_params$min)
  theta_max <-  unlist(df_params$max)
  theta_transform_type <- get_transform_type(theta_min, theta_max)
  blocks_list <- lapply(df_params$block, as.integer)
  n_unique_blocks <- length(unique(unlist(blocks_list)))
  beta_init <- seq(1, 0, length.out = n_rungs)

  stopifnot(is.integer(burnin))
  stopifnot(is.integer(samples))
  stopifnot(is.integer(n_unique_blocks))
  stopifnot(is.integer(n_rungs))

  out <- mcmc(theta_init, theta_names, theta_transform_type,  theta_min,  theta_max,
       blocks_list, n_unique_blocks, data, burnin, samples, loglike, logprior,
       target_acceptance, misc, n_rungs, beta_init, swap)
  out$output <- cbind(
    data.frame(chain = 1, phase = rep(c("burnin", "sampling"), c(burnin, samples))),
    out$output
    )
  names(out$output) <- c("chain", "phase", "iteration", theta_names, "logprior", "loglikelihood")

  # Diagnostics
  # DIC
  deviance <- -2 * out$output[out$output$phase == "sampling", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * var(deviance)
  out$diagnostics$DIC_Gelman <- dic
  # MC acceptance
  out$diagnostics$mc_accept <- out$swap_acceptance / (burnin + samples)
  # ESS
  out$diagnostics$ess <- apply(out$output[out$output$phase == "sampling", theta_names], 2, coda::effectiveSize)
  # Rung index
  out$diagnostics$rung_index <- out$rung_index

  # Parameters
  out$parameters <- list()
  out$parameters$data <- data
  out$parameters$df_params <- df_params
  out$parameters$loglike <- loglike
  out$parameters$logprior <- logprior
  out$parameters$burnin <- burnin
  out$parameters$samples <- samples
  out$parameters$rungs <- n_rungs
  out$parameters$chains <- 1
  out$parameters$coupling_on <- swap

  out <- out[c("output", "diagnostics", "parameters")]

  return(out)
}


get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}
