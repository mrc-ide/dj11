run_dj11 <- function(data, df_params, loglike, logprior, burnin, samples, target_acceptance = 0.44, misc = list(),
                     n_rungs = 1L, swap = TRUE){

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
  out$out <- cbind(
    data.frame(chain = 1, phase = rep(c("burnin", "sampling"), c(burnin, samples))),
    out$out
    )
  names(out$out) <- c("chain", "phase", "iteration", theta_names, "logprior", "loglikelihood")
  return(out)

}


get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}
