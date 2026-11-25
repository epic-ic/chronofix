##' Run MCMC
##'
##' @title Run MCMC
##'
##' @param data Observed data
##'
##' @param delay_map Delays information
##' 
##' @param hyperparameters List of hyperparameters
##'
##' @param initial Initial values
##' 
##' @param control List of control parameters
##'
##' @return Output
##'
##' @export
chronofix_mcmc <- function(data,
                           delay_map,
                           hyperparameters = chronofix_hyperparameters(),
                           initial = chronofix_mcmc_initial(),
                           control = chronofix_mcmc_control()) {
  
  data <- chronofix_prepare_data(data)
  
  model <- chronofix_model(data, delay_map, hyperparameters, control)
  sampler <- chronofix_sampler(control)
  
  pars_initial <- mcmc_pars_initialise(model, initial)
  pars <- matrix(pars_initial, length(model$parameters), 
                 control$n_chains)
  
  parameters <- model$parameters
  
  runner <- 
    if (control$parallel) monty::monty_runner_callr(control$n_workers) else
      monty::monty_runner_serial()
  
  rng <- initial_rng(control$n_chains)
  steps <- mcmc_steps(control$n_steps, control$burnin, 
                      control$thinning_factor,
                      save_full_chains = TRUE)
  
  res <- runner$run(pars, model, sampler, steps, rng)
  
  samples <- combine_chains(res, data, delay_map)
  
  samples <- flatten_chains(samples)
  
  ## Unpack augmented data
  samples$augmented_data <- 
    unpack_augmented_data(samples$augmented_data, model$data_packer)
  id <- attr(data, "id")
  group <- attr(data, "group")
  rownames(samples$augmented_data$estimated_dates) <- data[[id]]
  colnames(samples$augmented_data$estimated_dates) <- 
    setdiff(names(data), c(id, group))
  rownames(samples$augmented_data$error_indicators) <- data[[id]]
  colnames(samples$augmented_data$error_indicators) <- 
    setdiff(names(data), c(id, group))
  
  samples
  
}

##' Create control parameters
##'
##' @title Create control parameters
##'
##' @param n_steps The number of steps to run in each MCMC chain
##' 
##' @param burnin The number of steps at the beginning of each chain to discard
##'   as burnin
##'
##' @param thinning_factor A thinning factor applied to the chains. If given,
##'   every`thinning_factor`'th step is retained 
##'
##' @param n_chains The number of chains to run
##' 
##' @param parallel Logical, indicating whether or not to run chains in parallel
##'
##' @param n_workers Number of workers to use for parallelisation
##'
##' @param lower_quantile Lower quantile used for initialisation of true dates
##' 
##' @param upper_quantile Upper quantile used for initialisation of true dates
##' 
##' @param earliest_possible_date The earliest possible date (in YYYY-MM-DD 
##'   format). Dates will not be estimated as occurring earlier than this date.
##'   If `NULL` (the default), then it will be determined by the 
##'   earliest date in the data and `date_buffer`.
##'   
##' @param latest_possible_date The latest possible date (in YYYY-MM-DD 
##'   format). Dates will not be estimated as occurring later than this date.
##'   If `NULL` (the default), then it will be determined by the latest
##'   date in the data and `date_buffer`.
##'   
##' @param date_buffer The date buffer in terms of days to determine the
##'   earliest and/or latest possible dates from the data. If 
##'   `earliest_possible_date` is not specified then the earliest possible date
##'   will be taken as `date_buffer` days before the earliest date in the data.
##'   Similarly, if `latest_possible_date` is not specified then the latest
##'   possible date will be taken as `date_buffer` days after the latest date in
##'   the data.
##' 
##' @param prob_update_estimated_dates The probability of proposing an update
##'   to each estimated date at each iteration in the MCMC
##'
##' @param prob_update_error_indicators The probability of proposing an update
##'   to each error indicator (with the corresponding estimated date updated
##'   accordingly) at each iteration in the MCMC
##'   
##' @param prob_error_swap The probability of proposing to swap all errors to
##' non-errors and vice versa (excluding missing dates) for individuals with at
##' least one error and non-error at each iteration in the MCMC
##' 
##' @return List of control parameters
##'
##' @export
chronofix_mcmc_control <- function(n_steps = 1000,
                                   burnin = 0,
                                   thinning_factor = 1,
                                   n_chains = 1,
                                   parallel = FALSE,
                                   n_workers = 1,
                                   lower_quantile = 0.01,
                                   upper_quantile = 0.99,
                                   earliest_possible_date = NULL,
                                   latest_possible_date = NULL,
                                   date_buffer = 30,
                                   prob_update_estimated_dates = 1,
                                   prob_update_error_indicators = 1,
                                   prob_error_swap = 1) {
  
  list(n_steps = n_steps,
       burnin = burnin,
       thinning_factor = thinning_factor,
       n_chains = n_chains,
       parallel = parallel,
       n_workers = n_workers,
       lower_quantile = lower_quantile,
       upper_quantile = upper_quantile,
       earliest_possible_date = earliest_possible_date,
       latest_possible_date = latest_possible_date,
       date_buffer = date_buffer,
       prob_update_estimated_dates = prob_update_estimated_dates,
       prob_update_error_indicators = prob_update_error_indicators,
       prob_error_swap = prob_error_swap)
}

##' Create initial parameter values
##'
##' @title Create initial parameter values
##'
##' @param initial_delay_shape The initial value for the shape parameter of
##'   gamma-distributed delays
##' 
##' @param initial_delay_mean The initial value for the mean of
##'   gamma-distributed delays
##'
##' @param initial_delay_meanlog The initial value for the mean on the log-scale
##'   of log-normal-distributed delays
##'
##' @param initial_delay_precisionlog The initial value for the precision on the
##'  log-scale of log-normal-distributed delays
##'
##' @param initial_prob_error The initial value for the probability of error  
##' 
##' @return List of initial parameter values
##'
##' @export
chronofix_mcmc_initial <- function(initial_delay_shape = 1,
                                   initial_delay_mean = 5,
                                   initial_delay_meanlog = 1,
                                   initial_delay_precisionlog = 1,
                                   initial_prob_error = 0.1) {
  
  list(initial_delay_shape = initial_delay_shape,
       initial_delay_mean = initial_delay_mean,
       initial_delay_meanlog = initial_delay_meanlog,
       initial_delay_precisionlog = initial_delay_precisionlog,
       initial_prob_error = initial_prob_error)
}


mcmc_pars_initialise <- function(model, initial) {
  
  pars_initial <- numeric(length(model$parameters))
  pars_initial[model$parameters == "prob_error"] <- initial$initial_prob_error
  pars_initial[endsWith(model$parameters, "shape")] <- 
    initial$initial_delay_shape
  pars_initial[endsWith(model$parameters, "mean")] <- initial$initial_delay_mean
  pars_initial[endsWith(model$parameters, "meanlog")] <- 
    initial$initial_delay_meanlog
  pars_initial[endsWith(model$parameters, "precisionlog")] <- 
    initial$initial_delay_precisionlog
  
  pars_initial
}


mcmc_steps <- function(n_steps, burnin = NULL, thinning_factor = NULL,
                       save_full_chains = FALSE,
                       call = parent.frame()) {
  if (inherits(n_steps, "monty_sample_steps")) {
    return(n_steps)
  }
  assert_scalar_size(n_steps, call = call)
  if (is.null(burnin)) {
    burnin <- 0
  } else {
    assert_scalar_size(burnin, allow_zero = TRUE, call = call)
    if (burnin >= n_steps) {
      cli::cli_abort("'burnin' must be smaller than 'n_steps'",
                     arg = "burnin", call = call)
    }
  }
  if (is.null(thinning_factor)) {
    thinning_factor <- 1
  } else {
    assert_scalar_size(thinning_factor, allow_zero = FALSE, call = call)
    if (n_steps %% thinning_factor != 0) {
      cli::cli_abort(
        "'thinning_factor' must be a divisor of 'n_steps'",
        call = call)
    }
  }
  assert_logical(save_full_chains)
  ret <- list(total = n_steps,
              burnin = burnin,
              thinning_factor = thinning_factor,
              save_full_chains = save_full_chains)
  class(ret) <- "monty_sample_steps"
  ret
}


initial_rng <- function(n_chains, seed = NULL) {
  lapply(rng_distributed_state(n_nodes = n_chains, seed = seed),
         function(s) monty::monty_rng_create(seed = s))
}


rng_distributed_state <- function(seed = NULL, n_nodes = 1L) {
  rng <- monty::monty_rng_create(seed = seed)
  ret <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    ret[[i]] <- monty::monty_rng_state(rng)
    if (i < n_nodes) {
      monty::monty_rng_long_jump(rng)
    }
  }
  ret
}


##' @importFrom abind abind
combine_chains <- function(res, data, delay_map) {
  if (inherits(res, "chronofix_mcmc_samples")) {
    return(res)
  }
  
  history <- lapply(res, "[[", "history")
  pars <- abind::abind(lapply(history, "[[", "pars"), along = 3)
  density <- abind::abind(arrays = lapply(history, "[[", "density"), along = 2)
  
  initial <- abind::abind(arrays = lapply(res, "[[", "initial"), along = 2)
  augmented_data <- 
    abind::abind(lapply(history, "[[", "augmented_data"), along = 3)
  
  full_chains <- lapply(history, "[[", "full_chains")
  pars_full <- abind::abind(lapply(full_chains, "[[", "pars"), along = 3)
  density_full <- abind::abind(lapply(full_chains, "[[", "density"), along = 2)
  full_chains <- chronofix_samples(pars_full, density_full, initial, data,
                                   delay_map)
  
  
  chronofix_samples(pars, density, initial, data, delay_map, 
                    augmented_data, full_chains)
}


flatten_chains <- function(samples) {
  d <- dim(samples$pars)
  n_samples <- d[[2]]
  n_chains <- d[[3]]
  
  chain <- rep(seq_len(n_chains), each = n_samples)
  
  samples$pars <- array_flatten(samples$pars, c(2, 3))
  samples$density <- array_flatten(samples$density, c(1, 2))
  samples$augmented_data <- array_flatten(samples$augmented_data, c(2, 3))
  
  attr(samples, "chain") <- chain
  
  samples
}


chronofix_samples <- function(pars, density, initial, data, delay_map,
                              augmented_data = NULL, full_chains = NULL,
                              flatten_chains = FALSE) {
  rownames(initial) <- rownames(pars)
  samples <- list(pars = pars,
                  density = density,
                  initial = initial,
                  data = data,
                  delay_map = delay_map,
                  augmented_data = augmented_data,
                  full_chains = full_chains)
  class(samples) <- "chronofix_mcmc_samples"
  samples
}
