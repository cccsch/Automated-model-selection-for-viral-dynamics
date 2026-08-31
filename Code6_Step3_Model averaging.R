
## ==========================================================================================================================
## Step III: Model averaging across the statistical models retained for one candidate biological model, to account for the uncertainty associated with statistical model selection.
##
## ===============================================================================================================================

library(dplyr)
library(tidyr)
library(data.table)
library(mvtnorm)  # rmvnorm
library(Matrix)   # nearPD

## Code to have the estimation of parameters of the averaged model 
set.seed(123456)

CANDIDATE_ID <- 1  # <-- edit for each candidate biological model
nrep <- 1000

# All parameters that appear anywhere in the model (structural + IIV
## standard deviations). Only those actually estimated in a given retained
## model (see estimated_params below, derived per model) are used.
param_names <- c(
  "h_pop", "init_pop", "c_pop", "k_pop", "delta_pop", "tau_pop", "a_pop",
  "omega_h", "omega_init", "omega_k", "omega_delta", "omega_tau", "omega_a",
  "beta_pop", "p_pop", "rho_pop", "phi_pop", "delta_adap_pop"
)

retained <- read.csv(sprintf("Retained_models_Candidate%d.csv", CANDIDATE_ID))
fixed_values <- read.csv(sprintf("Values_ACO_Candidate%d.csv", CANDIDATE_ID))

models <- retained$model_id

fixed_by_model <- as.data.table(merge(retained, fixed_values, by.x = "model_id", by.y = "model"))
setnames(fixed_by_model, "model_id", "model")

non_pop_cols <- setdiff(names(fixed_by_model), c("model", "BICc"))
setnames(fixed_by_model, non_pop_cols, paste0(non_pop_cols, "_pop"))

path_models <- file.path("Model_By_Hand", sprintf("Model_ACO_Candidate%d_", CANDIDATE_ID))

################################################################################
## Helpers
################################################################################

make_psd <- function(Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  as.matrix(nearPD(Sigma, corr = FALSE)$mat)
}

read_model_outputs <- function(model_id,path_models) {
  
  result <- read.table(
    paste0(path_models, model_id, "/populationParameters.txt"),
    header = TRUE, sep = ",", dec = "."
  )
  
  var_covar <- read.table(
    paste0(path_models, model_id, "/FisherInformation/covarianceEstimatesLin.txt"),
    header = FALSE, sep = ",", dec = "."
  )
  
  var_covar <- var_covar[,-1]
  
  list(
    result = result,
    var_covar = as.matrix(var_covar)
  )
}

build_mu_phi <- function(result, estimated_params) {
  
  result <- result %>%
    filter(parameter %in% estimated_params)
  
  vals <- result$value
  names(vals) <- result$parameter
  
  mu_phi <- log(vals[estimated_params])
  names(mu_phi) <- estimated_params
  
  mu_phi
}

build_sigma_phi <- function(var_covar, estimated_params) {
  
  n <- length(estimated_params)
  Sigma_phi <- var_covar[seq_len(n), seq_len(n)]
  
  rownames(Sigma_phi) <- estimated_params
  colnames(Sigma_phi) <- estimated_params
  
  Sigma_phi <- (Sigma_phi + t(Sigma_phi)) / 2
  
  eig <- eigen(Sigma_phi, symmetric = TRUE)$values
  
  if (min(eig) < -1e-8) {
    warning("Sigma_phi is not PSD; applying nearPD().")
    Sigma_phi <- make_psd(Sigma_phi)
  }
  
  Sigma_phi
}

################################################################################
## Build model objects
################################################################################

model_objects <- vector("list", length(models))
names(model_objects) <- paste0("M", models)

for (m in models) {
  
  out <- read_model_outputs(m, path_models)
  
  estimated_params <- intersect(param_names, out$result$parameter)
  
  mu_phi <- build_mu_phi(
    result = out$result,
    estimated_params = estimated_params
  )
  
  Sigma_phi <- build_sigma_phi(
    var_covar = out$var_covar,
    estimated_params = estimated_params
  )
  
  temp <- fixed_by_model[model == m]
  
  model_objects[[paste0("M", m)]] <- list(
    model = m,
    estimated_params = estimated_params,
    mu_phi = mu_phi,
    Sigma_phi = Sigma_phi,
    fixed = temp,
    IC = fixed_by_model[model == m, BICc]
  )
}

models_all <- 1:9
objects_all <- model_objects[paste0("M", models_all)]

################################################################################
## Model weights within each structural family
################################################################################

get_model_weights <- function(model_objects_subset) {
  
  IC_values <- sapply(model_objects_subset, function(x) x$IC)
  delta_IC <- IC_values - min(IC_values)
  
  weights <- exp(-0.5 * delta_IC) / sum(exp(-0.5 * delta_IC))
  
  data.frame(
    model = sapply(model_objects_subset, function(x) x$model),
    IC = IC_values,
    delta_IC = delta_IC,
    weight = weights
  )
}

weight_table_ALL <- get_model_weights(objects_all)

print(weight_table_ALL)

weights_ALL <- weight_table_ALL$weight

################################################################################
## Sampling one model
################################################################################

sample_one_from_model <- function(obj, param_names) {
  
  estimated_params <- obj$estimated_params
  
  phi_draw <- as.numeric(
    rmvnorm(
      n = 1,
      mean = obj$mu_phi,
      sigma = obj$Sigma_phi
    )
  )
  
  names(phi_draw) <- estimated_params
  theta_est <- exp(phi_draw)
  
  fixed_row <- as.data.frame(obj$fixed)
  
  theta <- setNames(rep(NA_real_, length(param_names)), param_names)
  theta[estimated_params] <- theta_est[estimated_params]
  
  for (p in param_names) {
    
    if (p %in% colnames(fixed_row) == T) {
      
      theta[p] <- fixed_row[[p]]
    }
  }
  
  theta
}

################################################################################
## Sampling from a model-averaged distribution
################################################################################

sample_model_averaged_parameters <- function(model_objects_subset,
                                             weights,
                                             nrep,
                                             param_names,
                                             method_name) {
  
  chosen_models <- sample(
    x = seq_along(model_objects_subset),
    size = nrep,
    replace = TRUE,
    prob = weights
  )
  
  draws <- vector("list", nrep)
  
  for (r in seq_len(nrep)) {
    
    obj <- model_objects_subset[[chosen_models[r]]]
    
    theta_r <- sample_one_from_model(
      obj = obj,
      param_names = param_names
    )
    
    draws[[r]] <- data.frame(
      rep = r,
      method = method_name,
      model = obj$model,
      t(theta_r),
      check.names = FALSE
    )
  }
  
  bind_rows(draws)
}


theta_MA <- sample_model_averaged_parameters(
  model_objects_subset = objects_all,
  weights = weights_ALL,
  nrep = nrep,
  param_names = param_names,
  method_name = "MA"
)

################################################################################
## Confidence intervals
################################################################################

summarize_parameters <- function(theta_df, param_names) {
  
  theta_df %>%
    pivot_longer(
      cols = all_of(param_names),
      names_to = "param",
      values_to = "value"
    ) %>%
    group_by(method, param) %>%
    summarize(
      median = median(value, na.rm = TRUE),
      q025 = quantile(value, 0.025, na.rm = TRUE),
      q975 = quantile(value, 0.975, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
}

CI_MA <- summarize_parameters(theta_MA, param_names)
write.csv(CI_MA, sprintf("CI_ModelAveraged_Candidate%d.csv", CANDIDATE_ID), row.names = FALSE)