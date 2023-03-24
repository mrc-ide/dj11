run_dj11 <- function(data, df_params, loglike, logprior, burnin, samples, target_acceptance = 0.44, misc = list(),
                     n_rungs = 1L, swap = TRUE, chains = 1, cluster = NULL){

  # Save inputs
  # TODO: these will need to be chain-specific with different starting vals
  input <- list()
  input$theta_init <- unlist(df_params$init)
  input$theta_names <- unlist(df_params$name)
  input$theta_min <-  unlist(df_params$min)
  input$theta_max <-  unlist(df_params$max)
  input$theta_transform_type <- get_transform_type(input$theta_min, input$theta_max)
  input$blocks_list <- lapply(df_params$block, as.integer)
  input$n_unique_blocks <- length(unique(unlist(input$blocks_list)))
  input$data <- data
  input$burnin <- burnin
  input$samples <- samples
  input$loglike <- loglike
  input$logprior <- logprior
  input$target_acceptance <- target_acceptance
  input$misc <- misc
  input$n_rungs <- n_rungs
  input$beta_init <- seq(1, 0, length.out = input$n_rungs)
  input$swap <- swap
  input$chains <- chains

  stopifnot(is.integer(input$burnin))
  stopifnot(is.integer(input$samples))
  stopifnot(is.integer(input$n_unique_blocks))
  stopifnot(is.integer(input$n_rungs))

  if(is.null(cluster)){
    mcmc_runs <- lapply(1:chains, run_internal, input = input)
  } else {
    parallel::clusterEvalQ(cluster, library(dj11))
    mcmc_runs <- parallel::clusterApplyLB(cl = cluster, 1:chains, run_internal, input = input)
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
  out$diagnostics$ess <- apply(out$output[out$output$phase == "sampling", input$theta_names], 2, coda::effectiveSize)
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in seq_along(input$theta_names)) {
      rhat_est[p] <- out$output[out$output$phase == "sampling", c("chain", input$theta_names[p])] |>
        drjacoby:::gelman_rubin(chains = chains, samples = samples)
    }
    names(rhat_est) <- input$theta_names
    out$diagnostics$rhat <- rhat_est
  }


  # Parameters
  out$parameters <- input[c("data",
                            "df_params",
                            "loglike",
                            "logprior",
                            "burnin",
                            "samples",
                            "n_rungs",
                            "chains",
                            "swap")]

  out <- out[c("output", "diagnostics", "parameters")]

  return(out)
}

run_internal <- function(x, input){
  if(is.character(input$loglike)){
    input$loglike <- get(input$loglike)
  }
  if(is.character(input$loglike)){
    input$loglike <- get(input$logprior)
  }
  mcmc_out <- mcmc(input$theta_init, input$theta_names, input$theta_transform_type,  input$theta_min,  input$theta_max,
                   input$blocks_list, input$n_unique_blocks, input$data, input$burnin, input$samples, input$loglike, input$logprior,
                   input$target_acceptance, input$misc, input$n_rungs, input$beta_init, input$swap)
  mcmc_out$output <- cbind(
    data.frame(chain = x, phase = rep(c("burnin", "sampling"), c(input$burnin, input$samples))),
    mcmc_out$output
  )
  names(mcmc_out$output) <- c("chain", "phase", "iteration", input$theta_names, "logprior", "loglikelihood")
  return(mcmc_out)
}

get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}
