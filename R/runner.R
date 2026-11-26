runner_serial <- function(progress = NULL) {
  run <- function(pars, model, control, steps, rng) {
    n_chains <- length(rng)
    pb <- progress_bar(n_chains, steps$total, progress, show_overall = TRUE)
    with_progress_fail_on_error(
      pb,
      lapply(
        seq_along(rng),
        function(i) {
          run_chain(i, pars[, i], model, control, steps,
                    pb$update, rng[[i]])
        }))
  }
  
  list(run = run)
}


runner_parallel <- function(n_workers) {
  
  run <- function(pars, model, control, steps, rng) {
    n_chains <- length(rng)
    cl <- parallel::makeCluster(min(n_chains, n_workers))
    on.exit(parallel::stopCluster(cl))
    pars_list <- asplit(pars, MARGIN = 2)
    rng_state <- lapply(rng, function(r) monty::monty_rng_state(r))
    
    args <- list(model = model,
                 control = control,
                 steps = steps)
    
    ## To debug issues in the parallel sampler, it's most efficient to
    ## replace this call with `Map` and drop the `cl` argument, then
    ## you get nice stack traces.
    parallel::clusterMap(
      cl,
      run_chain_parallel,
      pars = pars_list,
      rng = rng_state,
      MoreArgs = args)
  }
  
  list(run = run)
}


run_chain_parallel <- function(chain_id, pars, model, control, steps,
                               rng) {
  rng <- monty::monty_rng_create(seed = rng)
  progress <- progress_bar_none()$update
  run_chain(chain_id, pars, model, control, steps, progress, rng)
}


run_chain <- function(chain_id, pars, model, control, steps,
                      progress, rng) {
  r_rng_state <- get_r_rng_state()
  
  chain_state <- initialise_state(pars, model, rng)
  
  run_chain2(chain_id, chain_state, model, control,
             steps, progress, rng, r_rng_state)
}


run_chain2 <- function(chain_id, chain_state, model,
                       control, steps, progress, rng, r_rng_state) {
  initial <- chain_state$pars
  n_pars <- length(model$parameters)
  
  burnin <- steps$burnin
  thinning_factor <- steps$thinning_factor
  n_steps <- steps$total
  n_steps_record <- ceiling((steps$total - burnin) / thinning_factor)
  
  history_pars <- matrix(NA_real_, n_pars, n_steps_record)
  history_density <- rep(NA_real_, n_steps_record)
  
  
  dim_est_dates <- dim(chain_state$augmented_data$estimated_dates)
  dim_err_inds <- dim(chain_state$augmented_data$error_indicators)
  history_augmented_data <- 
    list(estimated_dates = array(NA_real_, c(dim_est_dates, n_steps_record)),
         error_indicators = array(NA, c(dim_err_inds, n_steps_record)))
  
  history_full_chains <- list(pars = matrix(NA_real_, n_pars, n_steps),
                              density = rep(NA_real_, n_steps))
  
  j <- 1L
  for (i in seq_len(n_steps)) {
    chain_state <- model$augmented_data_update(chain_state, model, rng)
    chain_state <- update_pars(chain_state, control, model, rng)
    if (i > burnin && i %% thinning_factor == 0) {
      history_pars[, j] <- chain_state$pars
      history_density[[j]] <- chain_state$density
      history_augmented_data$estimated_dates[, , j] <-
        chain_state$augmented_data$estimated_dates
      history_augmented_data$error_indicators[, , j] <-
        chain_state$augmented_data$error_indicators
      j <- j + 1L
    }
    history_full_chains$pars[, i] <- chain_state$pars
    history_full_chains$density[i] <- chain_state$density
    progress(chain_id, i)
  }
  
  ## Pop the parameter names on last
  rownames(history_pars) <- model$parameters
  rownames(history_full_chains$pars) <- model$parameters
  
  ## This list will hold things that we'll use internally but not
  ## surface to the user in the final object (or summarise them in
  ## some particular way with no guarantees about the format).  We
  ## might hold things like start and stop times here in future.
  history <- list(
    pars = history_pars,
    density = history_density,
    augmented_data = history_augmented_data,
    full_chains = history_full_chains)
  
  list(
    history = history,
    initial = initial)
}


initialise_state <- function(pars, model, rng) {
  chain_state <- list(pars = pars)
  
  chain_state <- model$augmented_data_update(chain_state, model, rng)
  
  chain_state
}
