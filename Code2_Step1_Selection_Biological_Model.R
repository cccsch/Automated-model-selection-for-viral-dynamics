
## =========================================================================================================================
## Step I: Decision Tree (DT) for biological model selection
##
## Sequentially tests the inclusion of each biological feature, keeping the
## change only if it improves the BICc (see Methods). 
##
## To run the 6 feature orders described in the Methods copy the "USAGE" section at the bottom into 6 small runner scripts, each reordering FEATURES accordingly and setting a distinct RUN_LABEL so their outputs do not overwrite each other
## Requirements: see GH1_Biological_Search_Space.R
## =======================================================================================================================

source("GH1_Biological_Search_Space.R")
source("Helpers.R")

## Features to test, in the order they are evaluated by the DT algorithm.
## Reorder this vector to test a different feature order (see note above).
FEATURES <- c("ecl", "Ref", "TG", "in_p", "adap_h", "adap_a", "adap_delta", "adap_c")

## Candidate options available for each feature, read directly from the
## SNIPPETS names defined in GH1 (e.g. "ecl_no", "ecl_expo" -> "no", "expo").
get_feature_options <- function(feat) {
  matches <- grep(paste0("^", feat, "_"), names(SNIPPETS), value = TRUE)
  gsub(paste0("^", feat, "_"), "", matches)
}

library("lixoftConnectors")

## Initial parameter values, use to initialize every candidate model before estimation.

## Run SAEM + BICc. Returns a very large value if the BICc could not be computed, so that a failed fit is never mistakenly selected as the best model
getMetric<- function(){
  
  scenario=getScenario()
  
  scenario$tasks=c(populationParameterEstimation=T,conditionalDistributionSampling=F,conditionalModeEstimation=F,standardErrorEstimation=F,logLikelihoodEstimation=T)
  setGeneralSettings(autoChains = FALSE, nbchains = 7)
  setScenario(scenario)
  setPopulationParameterEstimationSettings(variability="firstStage", nbExploratoryIterations=1000, exploratoryAutoStop=F,
                                           nbSmoothingIterations = 1000,smoothingAutoStop = F )
  runScenario()

  BICc=getEstimatedLogLikelihood()$importanceSampling[4]
  if(is.null(BICc)){
    BICc=10^10
    #OFV=10^10
  }

  return(BICc)

}

## fit_candidate_model(): builds, loads, initalizes, and fits one candidate model
fit_candidate_model <- function(included, project, model_dir = "generated_models") {
  
  model_path <- generate_model(included, model_dir)
  
  print(model_path)
  
  loadProject(project)
  setStructuralModel(model_path)
  init_param()
  init_dist()
  setErrorModel(yV = "constant")
  setErrorModel(yVI = "constant")
  
  bicc <- getMetric()
  
  list(
    model_path = model_path,
    BICc       = bicc,
    param      = getPopulationParameterInformation()
  )
}


## project: path to the base .mlxtran project
## start:  logical vector giving the model to start from
## search_space: named list giving, for each feature, the candidate values to consider
## best_metrics, best_BICc, best_model: previous BICc used 
## State: null for the first DT run, or the list returned by a previous call, to continue from its result without refitting the starting model 


decision_tree_biological_selection <- function(project,
                                               start = NULL,
                                               search_space = NULL,
                                               state = NULL,
                                               checkpoint_prefix = "DT_results",
                                               verbose = TRUE) {
  
  if (is.null(start)) start <- setNames(rep("no", length(FEATURES)), FEATURES)
  if (is.null(search_space)) {
    search_space <- setNames(lapply(FEATURES, get_feature_options), FEATURES)
  }
  

  
  if (is.null(state)) {
    fit <- fit_candidate_model(as.list(start), project)
    best <- list(included = start, BICc = fit$BICc, param = fit$param)
    saveProject(project)
    log_df <- data.frame(model = fit$model_path, BICc = fit$BICc)
    run_id <- 1
  } else {
    best   <- state$best
    log_df <- state$log
    run_id <- state$run_id
  }
  
  for (feat in FEATURES) {
    
    candidates <- setdiff(search_space[[feat]], best$included[[feat]])
    
    for (opt in candidates) {
      
      cat("---\n")
      cat("feat =", feat, "| opt =", opt, "\n")
      cat("AVANT : best$included[[feat]] =", best$included[[feat]], "\n")
      
    
      run_id <- run_id + 1
      candidate_included <- best$included
      candidate_included[[feat]] <- opt
      
      cat("APRES : candidate_included[[feat]] =", candidate_included[[feat]], "\n")
      print(candidate_included)   # avec les noms, pour tout voir clairement
      
      fit <- fit_candidate_model(candidate_included, project)
      
      log_df <- rbind(log_df, data.frame(model = fit$model_path, BICc = fit$BICc))
      
      if (!is.null(checkpoint_prefix)) {
        write.csv(log_df, sprintf("%s_%d.csv", checkpoint_prefix, run_id), row.names = FALSE)
      }
      if (verbose) {
        message(sprintf("[run_id=%d] feature=%s -> %s | BICc = %.3f",
                        run_id, feat, fit$model_path, fit$BICc))
      }
      
      if (fit$BICc < best$BICc) {
        best <- list(included = candidate_included, BICc = fit$BICc, param = fit$param)
        saveProject(project)
      }
    }
  }
  
  list(best = best, log = log_df, run_id = run_id)
}
  
## ===========================================================================================================================
## USAGE

## RUN_LABEL identifies this feature order's outputs, so that running this section 6 times produces 6 distinct sets of output files
RUN_LABEL <- "Order1"  # <-- edit for each of the 6 feature orders

initializeLixoftConnectors(software="monolix", force=TRUE)

## Build and register the core model (starting point)
core_features <- setNames(as.list(rep("no", length(FEATURES))), FEATURES)
core_model_path <- generate_model(core_features, out_dir = "generated_models")


newProject(modelFile =core_model_path ,data=list(dataFile="Data/Data.csv",
                                       headerTypes=c("ignore","id","time","obsid","observation","cens")))

dir.create("Models_DT", showWarnings = FALSE)

saveProject(projectFile = file.path("Models_DT", paste0("Model0_", RUN_LABEL, ".mlxtran")))
project <- file.path("Models_DT", paste0("Model0_", RUN_LABEL, ".mlxtran"))

res_1 = decision_tree_biological_selection(project, checkpoint_prefix = paste0("DT_", RUN_LABEL, "_run1"))

res_2 = decision_tree_biological_selection(project, 
                                           start = res_1$best$included, 
                                           state=run1, 
                                           checkpoint_prefix = paste0("DT_", RUN_LABEL, "_run2"))


write.csv(data.frame(t(unlist(res_2$best$included)), BICc = res_2$best$BICc),
          paste0("DT_", RUN_LABEL, "_selected_model.csv"), row.names = FALSE)
write.csv(res_2$log, paste0("DT_", RUN_LABEL, "_full_log.csv"), row.names = FALSE)