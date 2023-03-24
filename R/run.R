run_dj11 <- function(data, df_params, loglike, logprior, burnin, samples, target_acceptance = 0.44, misc = list(),
                     n_rungs = 1L, swap = TRUE, chains = 1, cluster = NULL){

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

  if(is.null(cluster)){
    mcmc_runs <- lapply(1:chains, function(x){
      mcmc_out <- mcmc(theta_init, theta_names, theta_transform_type,  theta_min,  theta_max,
                       blocks_list, n_unique_blocks, data, burnin, samples, loglike, logprior,
                       target_acceptance, misc, n_rungs, beta_init, swap)
      mcmc_out$output <- cbind(
        data.frame(chain = x, phase = rep(c("burnin", "sampling"), c(burnin, samples))),
        mcmc_out$output
      )
      names(mcmc_out$output) <- c("chain", "phase", "iteration", theta_names, "logprior", "loglikelihood")
      return(mcmc_out)
    })
  }

  out <- list()
  out$output <- dplyr::bind_rows(sapply(mcmc_runs, '[', 'output'))
  # Diagnostics
  # DIC
  deviance <- -2 * out$output[out$output$phase == "sampling", "loglikelihood"]
  dic  <- mean(deviance) + 0.5 * var(deviance)
  out$diagnostics$DIC_Gelman <- dic
  if(n_rungs > 1){
    # MC acceptance
    out$diagnostics$mc_accept <- dplyr::bind_rows(sapply(mcmc_runs, '[', 'swap_acceptance'))
    # Rung index
    out$diagnostics$rung_index <- dplyr::bind_rows(sapply(mcmc_runs, '[', 'rung_index'))
  } else {
    out$diagnostics$mc_accept <- NULL
    out$diagnostics$rung_index <- NULL
  }
  # ESS
  out$diagnostics$ess <- apply(out$output[out$output$phase == "sampling", theta_names], 2, coda::effectiveSize)
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in seq_along(theta_names)) {
      rhat_est[p] <- out$output[out$output$phase == "sampling", c("chain", theta_names[p])] |>
        drjacoby:::gelman_rubin(chains = chains, samples = samples)
    }
    names(rhat_est) <- theta_names
    out$diagnostics$rhat <- rhat_est
  }


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
