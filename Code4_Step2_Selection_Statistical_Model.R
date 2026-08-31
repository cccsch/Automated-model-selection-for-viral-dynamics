
library(parallel)
library("lixoftConnectors")

rm(list=ls())
set.seed(123456)


## Replace this template with the one of the biological model found by DT and make as input the parameters that are well estimated
TEMPLATE <- "
[LONGITUDINAL]
input = { a, delta, init, beta,  c, k, tau, delta_adap, phi, rho, p{input_params}}

EQUATION:

T_0 = (8*10^7)/30
R_0 = 0
V_0 = 0
I_0 = 1/30
E_0 = 0

t0 = 0
odeType=stiff

{equations}

if (t < init) 
  start = 0
else
 start = 1
end

if ( t < tau)
 delta2= delta 
 else
delta2 = delta_adap
end

VI = a*V^h

ddt_T = (-beta*VI*T - phi*I*T + rho*R)*start
ddt_R = (phi*I*T - rho*R)*start
ddt_E = (beta*VI*T - k*E)*start
ddt_I = (k*E - delta2*I)*start
ddt_V = (p*I - c*V)*start

V_log = max(log10(V), 0.001)
VI_log = max(log10(VI), 0.001)

OUTPUT:
output = {V_log, VI_log}
"
## Make snippets 

make_snippets <- function(non_identifiables, base_values) {
  
  multiplicateurs <- c(NA, 1, 2, 4, 8, 0.5, 0.25, 0.125)
  snippets <- list()
  
  for (j in seq_along(non_identifiables)) {
    
    param <- non_identifiables[j]
    feat_name <- paste0("F", j)
    
    snippets[[paste0(feat_name, "_op1")]] <- list(
      input_params = paste0(", ", param)
    )
    
    for (op in 2:8) {
      
      value <- base_values[[param]] * multiplicateurs[op]
      
      snippets[[paste0(feat_name, "_op", op)]] <- list(
        equations = paste0(param, " = ", format(value, scientific = FALSE))
      )
    }
  }
  
  snippets
}

generate_model <- function(
    model_id,
    non_identifiables,
    base_values,
    out_dir = "Candidate1"
) {
  
  dir.create(out_dir, showWarnings = FALSE)
  
  SNIPPETS <- make_snippets(
    non_identifiables = non_identifiables,
    base_values = base_values
  )
  
  features <- setNames(
    as.list(paste0("op", model_id)),
    paste0("F", seq_along(model_id))
  )
  
  slots <- list(
    input_params = "",
    equations = ""
  )
  
  for (feat in names(features)) {
    
    opt <- features[[feat]]
    key <- paste0(feat, "_", opt)
    snip <- SNIPPETS[[key]]
    
    if (!is.null(snip)) {
      for (block in names(snip)) {
        
        if (block == "input_params") {
          slots[[block]] <- paste0(slots[[block]], snip[[block]])
        }
        
        if (block == "equations") {
          slots[[block]] <- paste(slots[[block]], snip[[block]], sep = "\n")
        }
      }
    }
  }
  
  model_text <- TEMPLATE
  
  for (slot in names(slots)) {
    model_text <- gsub(
      paste0("\\{", slot, "\\}"),
      slots[[slot]],
      model_text
    )
  }
  
  name_parts <- paste0(
    non_identifiables,
    "_op",
    model_id,
    collapse = "-"
  )
  
  fname <- file.path(out_dir, paste0(name_parts, "-Model.txt"))
  
  writeLines(model_text, fname)
  
  return(fname)
}



get_model <- function(model_id, non_identifiables, base_values) {
  
  model <- setNames(
    as.list(paste0("op", model_id)),
    paste0("F", seq_along(model_id))
  )
  
  model <- generate_model(model_id = as.numeric(sub("op", "", model)),
                          non_identifiables = non_identifiables,
                          base_values = base_values)
  
  return(model)
}


init_param<- function(){
  param=getPopulationParameterInformation()
  
  if("beta_pop" %in% param$name){
    j=which(param$name=="beta_pop")
    param$initialValue[j]= 0.0001
  }
  
  if("p_pop" %in% param$name){
    j=which(param$name=="p_pop")
    param$initialValue[j]= 1000
  }
  
  if("c_pop" %in% param$name){
    j=which(param$name=="c_pop")
    param$initialValue[j]= 10
  }
  
  if("r_pop" %in% param$name){
    j=which(param$name=="r_pop")
    param$initialValue[j]= 0.1
  }
  
  if("init_pop" %in% param$name){
    j=which(param$name=="init_pop")
    param$initialValue[j]= 2
  }
  
  if("delta_pop" %in% param$name){
    j=which(param$name=="delta_pop")
    param$initialValue[j]= 2
  }
  
  if("k_pop" %in% param$name){
    j=which(param$name=="k_pop")
    param$initialValue[j]= 4
  }
  
  if("rho_pop" %in% param$name){
    j=which(param$name=="rho_pop")
    param$initialValue[j]= 0.01
  }
  
  if("phi_pop" %in% param$name){
    j=which(param$name=="phi_pop")
    param$initialValue[j]= 0.001
  }
  
  if("K_p_pop" %in% param$name){
    j=which(param$name=="K_p_pop")
    param$initialValue[j]= 0.2
  }
  
  
  if("a_pop" %in% param$name){
    j=which(param$name=="a_pop")
    param$initialValue[j]= 0.001
  }
  
  if("h_pop" %in% param$name){
    j=which(param$name=="h_pop")
    param$initialValue[j]= 0.9
  }
  
  if("tau_pop" %in% param$name){
    j=which(param$name=="tau_pop")
    param$initialValue[j]= 7.5
  }
  
  if("delta_adap_pop" %in% param$name){
    j=which(param$name=="delta_adap_pop")
    param$initialValue[j]= 5
  }
  
  if("c_max_pop" %in% param$name){
    j=which(param$name=="c_max_pop")
    param$initialValue[j]= 2
  }
  
  if("C_c_pop" %in% param$name){
    j=which(param$name=="C_c_pop")
    param$initialValue[j]= 0.2
  }
  
  
  if("C_a_pop" %in% param$name){
    j=which(param$name=="C_a_pop")
    param$initialValue[j]= 0.2
  }
  
  if("C_h_pop" %in% param$name){
    j=which(param$name=="C_h_pop")
    param$initialValue[j]= 0.2
  }
  
  if("V0_pop" %in% param$name){
    j=which(param$name=="V0_pop")
    param$initialValue[j]= 200
  }
  
  setPopulationParameterInformation(param)
}

ACO<-function(project,n_ant,rho=0.4, n_feat, IIV_params, non_identifiables, base_values, candidate_id = 1){

  initializeLixoftConnectors(software="monolix", force=TRUE)

  ## All files produced by this ACO run for this candidate biological model
  ## are kept under their own directories, so that running ACO for several
  ## candidate models (as recommended, see Methods) never overwrites another
  ## candidate's files.
  
  model_dir <- paste0("Models_Candidate", candidate_id)
  structural_dir <- paste0("Structural_Candidate", candidate_id)
  dir.create(model_dir, showWarnings = FALSE)
  
  feat_names <- paste0("feat", seq_len(n_feat))
  iiv_names <- paste0("iiv_", iiv_params)
  
  search_space_feat <- as.data.frame(
    replicate(n_feat, 1:8)
  )
  names(search_space_feat) <- feat_names
  
  search_space_iiv <- as.data.frame(
    replicate(length(iiv_params), c(1, 2,rep(0, 6)))
  )
  names(search_space_iiv) <- iiv_names
  
  search_space <- cbind(
    search_space_feat,
    search_space_iiv
  )
  
  ## Weight 
  
  weight_feat <- as.data.frame(
    replicate(n_feat, rep(1,8))
  )
  names(weight_feat) <- feat_names
  
  weight_iiv <- as.data.frame(
    replicate(length(iiv_params), c(1, 1,rep(0, 6)))
  )
  names(weight_iiv) <- iiv_names
  
  weight <- cbind(
    weight_feat,
    weight_iiv
  )
  
  
  
  data=data.frame(model=c(),IIV_on=c(),OFV=c(),BICc=c(),metric=c(),time=c())
  
  rank <- data.frame(
    id = numeric(0)
  )
  
  for (i in seq_len(n_feat)) {
    rank[[paste0("feat", i)]] <- numeric(0)
  }
  
  rank$metric <- numeric(0)
  
  RSE_model <- c()
  
  parameters <- data.frame(id = numeric(0),
                           matrix(nrow=0,
                                  ncol=length(iiv_params)))
  
  names(parameters) <- c("id", iiv_params)

  parameters_all <- as.data.frame(
    as.list(rep(2, length(iiv_params)))
  )
  
  names(parameters_all) <- iiv_params

  
  l=0
  g=rank
  par_d=parameters
  
  converged <- FALSE
  
  while (!converged) {
    l=l+1
    # print(paste("iteration while",l))
    prob=weight
    print(weight)
    for (i in 1:ncol(weight)) {
      
      s <- sum(weight[, i], na.rm = TRUE)
      
      if (s == 0 || is.na(s)) {
        cat("Problème colonne", i, ":", names(weight)[i], "\n")
        print(weight[, i])
        stop("Somme des poids nulle ou NA")
      }
      
      prob[, i] <- weight[, i] / s
      
      if (any(prob[, i] > 0.9, na.rm = TRUE)) {
        
        j <- which(prob[, i] > 0.9)[1]
        
        prob[, i] <- 0
        prob[j, i] <- 1
      }
    }

    
    
    sampled <- as.data.frame(
      lapply(
        seq_len(ncol(search_space)),
        function(j)
          sample(
            search_space[[j]],
            size = n_ant,
            prob = prob[[j]],
            replace = TRUE
          )
      )
    )
    names(sampled) <- names(search_space)
    
    
    # check later if we need to increase l 
    if (l<=1){
      
      sampled[which(colnames(sampled) %in% iiv_names)] <- as.data.frame(
        lapply(
          seq_len(ncol(search_space_iiv)),
          function(j)
            sample(
              c(1),
              size = n_ant,
              prob = c(1),
              replace = TRUE
            )
        )
      )
    
    }
    
    model_to_run <- c()
    model_info <- data.frame()
    for (i in 1:n_ant){
      loadProject(project)
      saveProject(projectFile = paste("Models_Candidate1/Model",i,".mlxtran",sep=""))
      

      
      IIV_on="IIV"
      
      model_id <- as.numeric(
        sampled[i, feat_names]
      )
      
      model=get_model(model_id, non_identifiables = non_identifiables, base_values = base_values)
      
      setStructuralModel(model)
      
      param_pop=getPopulationParameterInformation()
      param_with_pop <- param_pop$name[grep("_pop", param_pop$name)]
      
      # Remove "_pop" suffix
      param_without_pop <- gsub("_pop$", "", param_with_pop)
      args <- list()
      
      for (p in param_without_pop) {
        
        iiv_col <- paste0("iiv_", p)
        
        if (iiv_col %in% names(sampled)) {
          p2 <- sampled[[iiv_col]][i]
        } else {
          p2 <- 2
        }
        
        args[[p]] <- p2 == 1
      }
      IIV_param=names(args)[args==TRUE]
      for (k in IIV_param){
        IIV_on=paste(IIV_on,k,sep="_")
      }
      
      #saveProject(projectFile = paste("Models_ACO_Parallel/Model",i,".mlxtran",sep=""))
      
      setErrorModel(yV="constant")
      setErrorModel(yVI="constant")
      do.call(setIndividualParameterVariability, args)
      init_param()
      
      
      scenario=getScenario()
      scenario$tasks=c(populationParameterEstimation=T,conditionalDistributionSampling=F,conditionalModeEstimation=F,standardErrorEstimation=T,logLikelihoodEstimation=T)
      setScenario(scenario)
      setGeneralSettings(autoChains = FALSE, nbchains = 7)
      #5setPopulationParameterEstimationSettings(variability="firstStage")
      setPopulationParameterEstimationSettings(nbExploratoryIterations=1000, exploratoryAutoStop=F,
                                               nbSmoothingIterations = 500,smoothingAutoStop = F )
      model1=gsub(":","_",model)
      
      saveProject(projectFile = paste("Models_Candidate1/Model",i,".mlxtran",sep=""))
      
      
      p7=which(model==data$model & IIV_on==data$IIV_on)
      
      if (length(p7) != 0) {
        
        rank_i <- data.frame(id = i + n_ant)
        
        for (feat in feat_names) {
          rank_i[[feat]] <- sampled[[feat]][i]
        }
        
        rank_i$metric <- data$metric[p7[1]]
        
        parameters_i <- data.frame(id = i + n_ant)
        
        for (p in iiv_params) {
          if (p %in% names(args)) {
            parameters_i[[p]] <- args[[p]]
          } else {
            parameters_i[[p]] <- FALSE
          }
        }
        
        rank <- rbind(rank, rank_i)
        parameters <- rbind(parameters, parameters_i)
        
      } else {
        
        model_file <- paste("Models_Candidate1/Model", i, ".mlxtran", sep = "")
        
        model_to_run <- c(model_to_run, model_file)
        
        model_info_i <- data.frame(
          model_file = model_file,
          ant_id = i
        )
        
        for (feat in feat_names) {
          model_info_i[[feat]] <- sampled[[feat]][i]
        }
        
        for (p in iiv_params) {
          if (p %in% names(args)) {
            model_info_i[[p]] <- args[[p]]
          } else {
            model_info_i[[p]] <- FALSE
          }
        }
        
        model_info <- rbind(model_info, model_info_i)
      }
      
    }
    
     cl <- makeCluster(detectCores() - 1)
     clusterEvalQ(cl, library("lixoftConnectors"))
     clusterExport(cl, c("model_to_run","l"), envir = environment())
    
    
    
    
    res <- parLapply(cl, 1:length(model_to_run), function(i) {
      logfile <- paste0("debug_worker_test_", i, ".txt")
      
      write(paste("Starting worker", i,l), logfile)
      
      initializeLixoftConnectors(software="monolix", force=TRUE)
      loadProject(model_to_run[i])
      runPopulationParameterEstimation()
      runStandardErrorEstimation(linearization = TRUE)
      runLogLikelihoodEstimation(linearization = FALSE)
      saveProject(model_to_run[i])
      
      write("Project finish", logfile, append=TRUE)
    })
    stopCluster(cl)
    
    for (i in 1:length(model_to_run)){
      start.time=Sys.time()
      loadProject(model_to_run[i])
      model=getStructuralModel()
      corr_monolix=getCorrelationOfEstimates()
      se_rse=getEstimatedStandardErrors()
      rse=se_rse$linearization$rse
      corr_matrix=corr_monolix$stochasticApproximation
      k=length ( corr_matrix[1,])
      param=getPopulationParameterInformation()
      
      
      info_i <- model_info[model_info$model_file == model_to_run[i], ]
      
      for (feat in feat_names) {
        assign(feat, info_i[[feat]])
      }
      
      #######################################
      
      
      n_nan=sum(is.nan(rse))
      n_sup=sum(rse>50 & !is.nan(rse))
      
      OFV=getEstimatedLogLikelihood()$importanceSampling[1]
      BICc=getEstimatedLogLikelihood()$importanceSampling[4]
      if(is.null(BICc)){
        BICc=10^10
        OFV=10^10
      }

      metrics=BICc
      
      
      RSE=getEstimatedStandardErrors()$linearization$rse
      RSE_model_i=list(RSE)

      end.time=Sys.time()
      time.taken=round(end.time-start.time,2)
      param_pop=getPopulationParameterInformation()
      param_with_pop <- param_pop$name[grep("_pop", param_pop$name)]
      
      # Remove "_pop" suffix
      param_without_pop <- gsub("_pop$", "", param_with_pop)
      args <- list()
      
      for (p in iiv_params) {
        args[[p]] <- info_i[[p]]
      }
      IIV_on="IIV"
      IIV_param=names(args)[args==TRUE]
      for (k in IIV_param){
        IIV_on=paste(IIV_on,k,sep="_")
      }

      data_i=data.frame(model=c(model),IIV_on=c(IIV_on),OFV=c(OFV),BICc=c(BICc),metric=c(metrics),time=c(time.taken))
 
      for (p in iiv_params) {
        data_i[[p]] <- args[[p]]
      }
      
      rank_i <- data.frame(id = i)
      
      for (feat in feat_names) {
        rank_i[[feat]] <- info_i[[feat]]
      }
      
      rank_i$metric <- metrics
      
      parameters_i=cbind(id=c(i),parameters_all)
      
      for (p in param_without_pop){
        p1=p
        parameters_i[[p1]]=args[[p]]
      }
      rank=rbind(rank,rank_i)
      parameters=rbind(parameters,parameters_i)
      data=rbind(data,data_i)
      RSE_model=append(RSE_model,RSE_model_i)
    }
    
    rank=rank[order(rank$metric),]
    parameters=parameters[match(rank$id,parameters$id),]
    
    print(rank)
    
    #update weight
    weight=(1-rho)*weight
    
    n = length(feat_names)
    
    for (i in 1:n_ant){
      # 1 to number of features 
      for (j in seq_along(feat_names)){
        weight[rank[i,j+1],j]=weight[rank[i,j+1],j]+1/i
      }
      if(l>1){
        # number of feature + 1 to number of features + number of parameters
        for (j in (length(feat_names) + 1):ncol(weight)){
          # number of features - 1
          if (parameters[i,j-(n-1)]==TRUE){
            weight[1,j]=weight[1,j]+1/i
          }
          if (parameters[i,j-(n-1)]==FALSE){
            weight[2,j]=weight[2,j]+1/i
          }
        }
      }
    }
    if (l<=1){
      # number of feature + 1 to number of features + number of parameters
      for (j in (n+1):ncol(weight)){
        weight[1,j]=weight[1,j]/(1-rho)
        weight[2,j]=weight[2,j]/(1-rho)
        
      }
    }
    print(weight)
    # 2 to number of features + 1
    g = rank[, feat_names, drop = FALSE]
    rank=rank[1,]
    rank$id=n_ant+1
    # 2 to number of parameters + 1
    par_d <- parameters[, iiv_params, drop = FALSE]
    print(par_d)
    parameters=parameters[1,]
    parameters$id=n_ant+1
    write.csv(data[order(data$metric),], paste0("res_C", candidate_id, "_", l, ".csv"))
    
    conv_g <- nrow(g) > 0 &&
      max(table(apply(g, 1, paste, collapse = ","))) >= 3 * n_ant / 4
    
    conv_par <- nrow(par_d) > 0 &&
      max(table(apply(par_d, 1, paste, collapse = ","))) >= 3 * n_ant / 4
    
    converged <- conv_g & conv_par
    
    cat("conv_g =", conv_g, "\n")
    cat("conv_par =", conv_par, "\n")
    cat("converged =", converged, "\n")
    
  }
  
  return(list(data[order(data$metric),]))
}

## ============================================================================
## Usage: run ACO for ONE candidate biological model.
##
## ACO must be run once per candidate biological model retained by DT
## re-run this whole "Usage" section once per candidate, changing
## CANDIDATE_ID each time. 
##
## TEMPLATE must be replaced with the resolved Mlxtran text of
## the specific biological model being processed here (as produced by Code3 in
## Structural_By_Hand/ for this candidate), with the non-identifiable
## parameters listed below left as placeholders in the {equations}/
## {input_params} slots 
## ============================================================================

CANDIDATE_ID <- 1  # <-- edit for each candidate biological model

initializeLixoftConnectors(software = "monolix", force = TRUE)

start.time_for=Sys.time()


## Non-identifiable parameters (and their DT-estimated values) for this candidate

unidentifiable <- read.csv("Unidentifiable_parameters.csv")
unidentifiable <- unidentifiable[unidentifiable$model_id == CANDIDATE_ID, ]

structural_dir <- paste0("Structural_Candidate", CANDIDATE_ID)
model_dir <- paste0("Models_Candidate", CANDIDATE_ID)
dir.create(model_dir, showWarnings = FALSE)

non_identifiables <- unidentifiable$parameter
base_values <- setNames(unidentifiable$estimate, unidentifiable$parameter)
n_feat <- length(non_identifiables)

model0 <- generate_model(
  model_id = rep(1, n_feat),
  non_identifiables = non_identifiables,
  base_values = base_values,
  out_dir = structural_dir
)

newProject(modelFile =model0,
           data=list(dataFile="Data/Data.csv",
                     headerTypes=c("ignore","id","time","obsid","observation","cens")))

saveProject(projectFile = file.path(model_dir, "Model0.mlxtran"))
project <- file.path(model_dir, "Model0.mlxtran")

## Define all the parameters that can have IIV
iiv_params = c("a", "h", "init", "beta", "p", "c", "delta", "k", "rho", "phi", "delta_adap", "tau")

res=ACO(project,n_ant=35, rho = 0.4, n_feat = n_feat,
        IIV_params = iiv_params, non_identifiables = non_identifiables, 
        base_values = base_values, 
        candidate_id = CANDIDATE_ID)

write.csv(res[[1]],paste0("Models_Candidate", CANDIDATE_ID, ".csv"), row.names = FALSE)
