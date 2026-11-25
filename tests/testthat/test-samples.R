test_that("can print samples", {
  set.seed(1)
  toy <- toy_data()
  data <- toy$data$observed_data
  delay_map <- toy$delay_map
  
  hyperparameters <- chronofix_hyperparameters()
  initial <- chronofix_mcmc_initial()
  control <- chronofix_mcmc_control(n_steps = 50, n_chains = 3)
  
  samples <- chronofix_mcmc(data, delay_map, hyperparameters, initial, control)
  
  res <- evaluate_promise(withVisible(print(samples)))
  expect_mapequal(res$result, list(value = samples, visible = FALSE))
  expect_match(
    res$messages,
    "<chronofix_mcmc_samples: 13 parameters x 50 samples x 3 chains>",
    fixed = TRUE, all = FALSE)
  
})


test_that("can convert to posterior draws types", {
  skip_if_not_installed("posterior")
  set.seed(1)
  toy <- toy_data()
  data <- toy$data$observed_data
  delay_map <- toy$delay_map
  
  hyperparameters <- chronofix_hyperparameters()
  initial <- chronofix_mcmc_initial()
  control <- chronofix_mcmc_control(n_steps = 50, n_chains = 3)
  
  samples <- chronofix_mcmc(data, delay_map, hyperparameters, initial, control)
  
  samples_df <- posterior::as_draws_df(samples)
  expect_s3_class(samples_df, "draws_df")
  expect_setequal(names(samples_df), 
                  c(rownames(samples$pars), ".iteration", ".chain", ".draw"))
  
  samples_arr <- posterior::as_draws_array(samples)
  expect_s3_class(samples_arr, "draws_array")
  expect_equal(dimnames(samples_arr),
               list(chain = NULL, iteration = NULL, 
                    variable = rownames(samples$pars)))
})


test_that("can convert to coda type", {
  skip_if_not_installed("coda")
  set.seed(1)
  toy <- toy_data()
  data <- toy$data$observed_data
  delay_map <- toy$delay_map
  
  hyperparameters <- chronofix_hyperparameters()
  initial <- chronofix_mcmc_initial()
  control <- chronofix_mcmc_control(n_steps = 50, n_chains = 3)
  
  samples <- chronofix_mcmc(data, delay_map, hyperparameters, initial, control)
  
  samples_coda <- coda::as.mcmc.list(samples)
  expect_s3_class(samples_coda, "mcmc.list")
  expect_length(samples_coda, 3)
  expect_s3_class(samples_coda[[1]], "mcmc")
})
