## ===========================================================================================================================
## Intermediate Step I: aggregate DT results across the 6 feature orders,
## select candidate biological models (DeltaBICc < 2), and refit each of
## them individually to identify poorly estimated parameters (RSE > 50%),
## needed to build the ACO search space in Step II.
##

## =========================================================================================================================

source("GH1_Biological_Search_Space.R")
source("Helpers.R")

FEATURES <- c("ecl", "Ref", "TG", "in_p", "adap_h", "adap_a", "adap_delta", "adap_c")

initializeLixoftConnectors(software = "monolix", force = TRUE)  

## ----------------------------------------------------------------------------
## Step A: aggregate the 6 per-order DT logs and select candidate models
## ----------------------------------------------------------------------------

ORDER_LABELS <- paste0("Order", 1:6)


results <- do.call(rbind, lapply(seq_along(ORDER_LABELS), function(i) {
  df <- read.csv(sprintf("DT_%s_full_log.csv", ORDER_LABELS[i]))
  df$Order <- i
  df
}))

results$Difference <- results$BICc - min(results$BICc)
retained <- results[results$Difference <= 2, ]

parse_model_filename <- function(path, features = FEATURES) {
  chunks <- strsplit(basename(path), "-")[[1]]
  opts <- sapply(features, function(feat) {
    chunk <- chunks[grepl(paste0("^", feat, "_"), chunks)]
    if (length(chunk) == 0) NA_character_ else sub(paste0("^", feat, "_"), "", chunk[1])
  })
  as.data.frame(as.list(opts), stringsAsFactors = FALSE)
}

feature_sets <- do.call(rbind, lapply(retained$model, parse_model_filename))
feature_sets <- unique(feature_sets)

write.csv(feature_sets, "Candidate_biological_models.csv", row.names = FALSE)

## ----------------------------------------------------------------------------
## Step B: refit each candidate model individually, and record which
## parameters are unidentifiable (RSE > 50%), for use as the Step II (ACO)
## search space.
## ----------------------------------------------------------------------------

unidentifiable_results <- list()
 
for (i in seq_len(nrow(feature_sets))) {
 
  features <- as.list(feature_sets[i, ])
 
  model_path <- generate_model(features, out_dir = "Structural_By_Hand")
 
  newProject(
    modelFile = model_path,
    data = list(dataFile = "Data/Data.csv",
                headerTypes = c("ignore","id", "time", "obsid", "observation", "cens"))
  )
  project <- file.path("Model_By_Hand", paste0("Model", i, ".mlxtran"))
  dir.create("Model_By_Hand", showWarnings = FALSE)
  saveProject(project)
  loadProject(project)
 
  init_param(model = features, project = project)
  init_dist(model = features, project = project)
  setErrorModel(yV = "constant")
  setErrorModel(yVI = "constant")
 
  scenario <- getScenario()
  scenario$tasks <- c(
    populationParameterEstimation   = TRUE,
    conditionalDistributionSampling = FALSE,
    conditionalModeEstimation       = FALSE,
    standardErrorEstimation         = FALSE,
    logLikelihoodEstimation         = TRUE
  )
  setGeneralSettings(autoChains = FALSE, nbchains = 7)
  setScenario(scenario)
  setPopulationParameterEstimationSettings(
    variability = "firstStage",
    nbExploratoryIterations = 1000, exploratoryAutoStop = FALSE,
    nbSmoothingIterations   = 1000, smoothingAutoStop   = FALSE
  )
  runPopulationParameterEstimation()
  runStandardErrorEstimation(linearization = TRUE)
  runLogLikelihoodEstimation(linearization = FALSE)
 
  estimates <- getEstimatedPopulationParameters()
  se <- as.data.frame(getEstimatedStandardErrors())
 
  unidentifiable_names <- se$linearization.parameter[se$linearization.rse > 50]
  
  unidentifiable_names <- unidentifiable_names[grepl("_pop$", unidentifiable_names)]
 
  unidentifiable_results[[i]] <- data.frame(
    model_id            = i,
    features             = paste(names(features), unlist(features), sep = "_", collapse = "-"),
    parameter            = gsub("_pop$", "", unidentifiable_names),
    estimate              = estimates[unidentifiable_names]
  )
 
  message(sprintf("Model %d/%d: %d unidentifiable parameter(s)",
                   i, nrow(feature_sets), length(unidentifiable_names)))
}
 
unidentifiable_df <- do.call(rbind, unidentifiable_results)
write.csv(unidentifiable_df, "Unidentifiable_parameters.csv", row.names = FALSE)