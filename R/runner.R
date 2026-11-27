runner_serial <- function(progress = NULL) {
  run <- function(pars, model, control, rng) {
    n_chains <- length(rng)
    pb <- progress_bar(n_chains, control$n_steps, progress, show_overall = TRUE)
    with_progress_fail_on_error(
      pb,
      lapply(
        seq_along(rng),
        function(i) {
          run_chain(i, pars[, i], model, control, pb$update, rng[[i]])
        }))
  }
  
  list(run = run)
}


runner_callr <- function(n_workers, progress = NULL) {
  assert_scalar_size(n_workers, allow_zero = FALSE)
  progress <- show_progress_bar(progress)
  env <- new.env(parent = emptyenv())
  
  launch <- function(session_id) {
    pending <- which(env$status == "pending")
    if (length(pending) == 0) {
      return()
    }
    chain_id <- pending[[1]]
    env$sessions[[session_id]] <- callr::r_bg(
      function(id, path) run_chain_manual(id, path, "simple"),
      list(chain_id, env$path),
      package = TRUE)
    env$target[[session_id]] <- chain_id
    env$status[[chain_id]] <- "running"
  }
  
  ## callr::poll will poll, with a timeout, all our processes.
  ## There's not much downside to a long poll because if they *are*
  ## ready then they will return instantly. However, the process
  ## will only be interruptable each time the timeout triggers, so
  ## use 1000 here (1s).
  step <- function(timeout = 1000) {
    res <- callr::poll(env$sessions, timeout)
    is_done <- vcapply(res, "[[", "process") == "ready"
    if (any(is_done)) {
      for (session_id in which(is_done)) {
        chain_id <- env$target[[session_id]]
        env$result_path[[chain_id]] <-
          callr_safe_result(env$sessions[[session_id]])
        env$status[[chain_id]] <- "done"
        env$n_steps_progress[[chain_id]] <- env$n_steps
        launch(session_id)
      }
    }
    
    if (progress != "none") {
      has_stderr <- vcapply(res, "[[", "error") == "ready" & !is_done
      if (any(has_stderr)) {
        for (session_id in which(has_stderr)) {
          stderr <- env$sessions[[session_id]]$read_error_lines()
          res <- parse_progress_bar_simple(stderr)
          if (!is.null(res)) {
            env$n_steps_progress[[res$chain_id]] <- res$step
          }
        }
      }
      env$progress$update(env$chain_ids, env$n_steps_progress)
    }
    
    all(env$status == "done")
  }
  
  loop <- function(path, n_workers, n_chains, n_steps, progress) {
    pb <- progress_bar(n_chains, n_steps, progress, show_overall = TRUE)
    n_workers <- min(n_chains, n_workers)
    env$path <- path
    env$sessions <- vector("list", n_workers)
    env$target <- rep(NA_integer_, n_workers)
    env$status <- rep("pending", n_chains)
    env$chain_ids <- seq_len(n_chains)
    env$result_path <- rep(NA_character_, n_chains)
    env$n_steps <- n_steps
    env$n_steps_progress <- rep(0, n_chains)
    env$progress <- pb
    for (session_id in seq_len(n_workers)) {
      launch(session_id)
    }
    with_progress_fail_on_error(
      pb,
      while (!step()) {
      })
    
    res <- lapply(env$result_path, readRDS)
    unlink(env$path, recursive = TRUE)
    res
  }
  
  run <- function(pars, model, control, rng) {
    seed <- unlist(lapply(rng, function(r) monty::monty_rng_state(r)))
    n_chains <- length(rng)
    path <- tempfile()
    prepare_chain_manual(
      model = model, control = control, path = path, pars = pars, seed = seed)
    loop(path, n_workers, control$n_chains, control$n_steps, progress)
  }
  
  list(run = run)
}


run_chain <- function(chain_id, pars, model, control, progress, rng) {
  r_rng_state <- get_r_rng_state()
  
  chain_state <- initialise_state(pars, model, rng)
  
  run_chain2(chain_id, chain_state, model, control, progress, rng, r_rng_state)
}


run_chain2 <- function(chain_id, chain_state, model, control,
                       progress, rng, r_rng_state) {
  initial <- chain_state$pars
  n_pars <- length(model$parameters)
  
  burnin <- control$burnin
  thinning_factor <- control$thinning_factor
  n_steps <- control$n_steps
  n_steps_record <- ceiling((n_steps - burnin) / thinning_factor)
  
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


prepare_chain_manual <- function(model, control, path, pars,
                                 seed = NULL, call = parent.frame()) {
  
  rng <- initial_rng(control$n_chains, seed = seed)
  rng_state <- lapply(rng, function(r) monty::monty_rng_state(r))
  
  dat <- list(pars = pars,
              model = model,
              control = control,
              rng_state = rng_state)
  
  chain_manual_path_create(path, dat, call = call)
}


run_chain_manual <- function(chain_id, path, progress = NULL) {
  assert_scalar_size(chain_id, allow_zero = FALSE)
  assert_scalar_character(path)
  path <- chain_manual_path(path, chain_id)
  
  inputs <- readRDS(path$inputs)
  
  control <- inputs$control
  if (chain_id > control$n_chains) {
    cli::cli_abort("'chain_id' must be an integer in 1..{inputs$n_chains}")
  }
  
  pb <- progress_bar(control$n_chains, control$n_steps, progress,
                     show_overall = FALSE, single_chain = TRUE)
  
  with_progress_fail_on_error(
    pb,
    {
      pars <- inputs$pars
      model <- inputs$model
      rng <- monty_rng_create(seed = inputs$rng_state[[chain_id]])
      res <- monty_run_chain(chain_id, pars[, chain_id], model, control,
                             pb$update, rng)
    }
  )
  
  saveRDS(res, path$results)
  invisible(path$results)
}


chain_manual_path_create <- function(path, dat, call = parent.frame()) {
  assert_scalar_character(path, call = call)
  exists <- file.exists(path) && !is_empty_directory(path)
  dir_create(path)
  path <- chain_manual_path(path, check = FALSE)
  
  ## NOTE: The 'suppressWarnings' avoids warnings about monty not
  ## being available on load.
  suppressWarnings(saveRDS(dat, path$inputs))
  invisible(path$root)
}


chain_manual_path <- function(path, chain_id = NULL, check = TRUE) {
  path_inputs <- file.path(path, "inputs.rds")
  
  list(root = path,
       inputs = path_inputs,
       results = file.path(path, sprintf("results_%d.rds", chain_id)))
}


callr_safe_result <- function(rs, grace = 2, dt = 0.1) {
  for (i in seq(0, grace, by = dt)) {
    if (!rs$is_alive()) {
      break
    }
    Sys.sleep(dt)
  }
  ## Failure here will throw and propagate nicely to the controlling
  ## process.
  rs$get_result()
}
