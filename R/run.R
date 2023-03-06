run_dj11 <- function(data, df_params, loglike, logprior, burnin, samples){

  theta_init <- unlist(df_params$init)
  theta_names <- unlist(df_params$name)
  theta_min <-  unlist(df_params$min)
  theta_max <-  unlist(df_params$max)
  theta_transform_type <- get_transform_type(theta_min, theta_max)
  blocks <- as.integer(unlist(df_params$block))
  n_unique_blocks <- length(unique(blocks))

  mcmc(theta_init, theta_names, theta_transform_type,  theta_min,  theta_max, blocks, n_unique_blocks, data, burnin, samples, loglike, logprior)
}


get_transform_type <- function(theta_min, theta_max){
  as.integer(2 * is.finite(theta_min) + is.finite(theta_max))
}
