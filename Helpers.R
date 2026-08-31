
## =========================================================================================================================
## Shared parameter initialization helpers
##

## Initial parameter values, use to initialize every candidate model before estimation.

init_param<- function(){
  param=getPopulationParameterInformation()
  
  if("beta_pop" %in% param$name){
    j=which(param$name=="beta_pop")
    param$initialValue[j]= 0.0034
  }
  
  if("p_pop" %in% param$name){
    j=which(param$name=="p_pop")
    param$initialValue[j]= 4695
  }
  
  if("c_pop" %in% param$name){
    j=which(param$name=="c_pop")
    param$initialValue[j]= 3.87
  }
  
  if("r_pop" %in% param$name){
    j=which(param$name=="r_pop")
    param$initialValue[j]= 0.1
  }
  
  if("init_pop" %in% param$name){
    j=which(param$name=="init_pop")
    param$initialValue[j]= 0.33
  }
  
  if("delta_pop" %in% param$name){
    j=which(param$name=="delta_pop")
    param$initialValue[j]= 1.41
  }
  
  if("k_pop" %in% param$name){
    j=which(param$name=="k_pop")
    param$initialValue[j]= 2.08
  }
  
  if("rho_pop" %in% param$name){
    j=which(param$name=="rho_pop")
    param$initialValue[j]= 0.035
  }
  
  if("phi_pop" %in% param$name){
    j=which(param$name=="phi_pop")
    param$initialValue[j]= 2.72*10^(-5)
  }
  
  if("K_p_pop" %in% param$name){
    j=which(param$name=="K_p_pop")
    param$initialValue[j]= 0.2
  }
  
  
  if("a_pop" %in% param$name){
    j=which(param$name=="a_pop")
    param$initialValue[j]= 1.968*10^(-5)
  }
  
  if("h_pop" %in% param$name){
    j=which(param$name=="h_pop")
    param$initialValue[j]= 0.93
  }
  
  if("tau_pop" %in% param$name){
    j=which(param$name=="tau_pop")
    param$initialValue[j]= 7.5
  }
  
  if("delta_adap_pop" %in% param$name){
    j=which(param$name=="delta_adap_pop")
    param$initialValue[j]= 38.25
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
  
  setPopulationParameterInformation(param)
}

## Parameters constrained to (0,1) by construction are given a logit-normal rather than the default log-normal distribution, when present in the candidate model.
init_dist<- function(){
 param= getPopulationParameterInformation()
  
  if("C_c_pop" %in% param$name){
    setIndividualParameterDistribution(C_c = "logitNormal")
  }
  
  
  if("C_a_pop" %in% param$name){
    setIndividualParameterDistribution(C_a = "logitNormal")
  }
  
  if("C_h_pop" %in% param$name){
    setIndividualParameterDistribution(C_h = "logitNormal")
  }
 
 if("K_p_pop" %in% param$name){
   setIndividualParameterDistribution(K_p = "logitNormal")
 }
  
  setPopulationParameterInformation(param)
}
