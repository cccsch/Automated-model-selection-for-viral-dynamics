## ============================================================================
## Intermediate Step II: from ACO's final ranking for one candidate biological model, select statistical models (DeltaBICc < 2), and refit
## each of them individually (needed for Step III, model averaging).
##
## ============================================================================

initializeLixoftConnectors(software="monolix",force=TRUE)

source("Helpers.R")

CANDIDATE_ID <- 1  # <-- edit for each candidate biological model

initializeLixoftConnectors(software = "monolix", force = TRUE)

## All parameters that can be given inter-individual variability 
iiv_params <- c("a", "h", "init", "beta", "p", "c", "delta", "k", "rho", "phi", "delta_adap", "tau")

MULTIPLIERS <- c(NA, 1, 2, 4, 8, 0.5, 0.25, 0.125)

result_files <- list.files(pattern = sprintf("^res_C%d_[0-9]+\\.csv$", CANDIDATE_ID))
stopifnot(length(result_files) > 0)
iterations <- as.integer(sub(sprintf("^res_C%d_([0-9]+)\\.csv$", CANDIDATE_ID), "\\1", result_files))
last_result_file <- result_files[which.max(iterations)]

Result <- read.csv(last_result_file)
Result$Diff <- Result$BICc - min(Result$BICc)
Result <- Result[Result$Diff <= 2, ]
Result <- Result[!duplicated(Result$BICc), ]

# ------------------------------------------------------------------------------------------------------------------------
## Non-identifiable parameters and their DT-estimated (base) values for this
## candidate, read from GH3's output rather than typed in by hand.
## -----------------------------------------------------------------------------------------------------------------------
unidentifiable <- read.csv("Unidentifiable_parameters.csv")
unidentifiable <- unidentifiable[unidentifiable$model_id == CANDIDATE_ID, ]
non_identifiables <- unidentifiable$parameter
base_values <- setNames(unidentifiable$estimate, unidentifiable$parameter)

## ------------------------------------------------------------------------------------------------------------------------

## parse_fixed_values(): from an ACO model fil path recover the fixed value used for each non-identifiable parameter (NA if the parameter is estimated)

parse_fixed_values <- function(path, non_identifiables, base_values) {
  chunks <- strsplit(basename(path), "-")[[1]]
  
  values <- sapply(non_identifiables, function(param) {
    chunk <- chunks[grepl(paste0("^", param, "_op"), chunks)]
    if (length(chunk) == 0) return(NA_real_)
    op <- as.integer(sub(paste0("^", param, "_op"), "", chunk[1]))
    base_values[[param]] * MULTIPLIERS[op]
  })
  names(values) <- non_identifiables
  values
}

values <- do.call(rbind, lapply(seq_len(nrow(Result)), function(i) {
  v <- parse_fixed_values(Result$model[i], non_identifiables, base_values)
  as.data.frame(as.list(c(v, model = i)))
}))
write.csv(values, sprintf("Values_ACO_Candidate%d.csv", CANDIDATE_ID), row.names = FALSE)

write.csv(
  data.frame(model_id = seq_len(nrow(Result)), BICc = Result$BICc),
  sprintf("Retained_models_Candidate%d.csv", CANDIDATE_ID),
  row.names = FALSE
)

structural_dir <- paste0("Candidate", CANDIDATE_ID)
model_out_dir  <- "Model_By_Hand"
dir.create(model_out_dir, showWarnings = FALSE)

for (i in seq_len(nrow(Result))) {
  
  df_i <- Result[i, ]
  structural_model <- basename(df_i$model[1])
  
  initializeLixoftConnectors()
  
  newProject(
    modelFile = file.path(structural_dir, structural_model),
    data = list(dataFile = "Data/Data2.csv",
                headerTypes = c("ignore","ignore","id", "time", "obsid", "observation", "cens"))
  )
  
  project <- file.path(model_out_dir, sprintf("Model_ACO_Candidate%d_%d.mlxtran", CANDIDATE_ID, i))
  saveProject(project)
  loadProject(project)
  
  init_param()
  
  setErrorModel(yV = "constant")
  setErrorModel(yVI = "constant")
  
  iiv_flags <- setNames(as.logical(unlist(df_i[iiv_params])), iiv_params)
  iiv_flags[is.na(iiv_flags)] <- FALSE
  
  fixed_vals_i <- parse_fixed_values(df_i$model[1], non_identifiables, base_values)
  is_fixed_i <- names(fixed_vals_i)[!is.na(fixed_vals_i)]
  iiv_flags <- iiv_flags[setdiff(names(iiv_flags), is_fixed_i)]
  
  setIndividualParameterModel(list(variability = list(id = as.list(iiv_flags))))
  
  scenario <- getScenario()
  
  scenario$tasks=c(populationParameterEstimation=T,conditionalDistributionSampling=F,conditionalModeEstimation=F,standardErrorEstimation=F,logLikelihoodEstimation=T)
  setGeneralSettings(autoChains = FALSE, nbchains = 7)
  setScenario(scenario)
  setPopulationParameterEstimationSettings(variability="firstStage", nbExploratoryIterations=1000, exploratoryAutoStop=F,
                                           nbSmoothingIterations = 1000,smoothingAutoStop = F )
  
  runPopulationParameterEstimation()
  runStandardErrorEstimation(linearization = TRUE)
  runLogLikelihoodEstimation(linearization = FALSE)
  saveProject(project)
}
