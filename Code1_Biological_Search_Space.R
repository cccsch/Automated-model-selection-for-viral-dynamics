## ======================================================================================================================
## Biological search space construction
##
## Builds Monolix (.txt) model files from a modular template, by combining
## "features" (biological mechanisms) that can each be turned on/off (or set
## to a specific variant), following the biological model bank described in
## the paper (Fig 2).
##

## =====================================================================================================================


set.seed(123456)

# Feature name -> biological mechanism (paper notation, see Fig 2)
## ecl          eclipse phase                          (TIVWE) "no"/"expo"
## Ref          refractory cells                        (TIVWR) "no"/"yes"
## TG           target cell regeneration                (TrIVW) "no"/"yes"
## in_p         interferon-mediated inhibition of p      (TIpVW) "no"/"Hill"
## adap_h       antibody effect on infectivity exponent h (TIVWh) "no"/"Hill"
## adap_a       antibody effect on infectivity scaling a  (TIVWa) "no"/"Hill"
## adap_delta   cell-mediated effect on death rate        (TIdeltaVW) "no"/"Hill"
## adap_c       antibody effect on viral clearance c      (TIVcW) "no"/"Hill"
##



## Template for the core model (TIVW), with placeholders filled in by the snippets selected for each feature

TEMPLATE <- "
[LONGITUDINAL]
input = {a, h, init, beta, p, c{input_params}}

EQUATION:

T_0 = (8*10^7)/30
V_0 = 0
VI_0 = 0
R_0 = 0
I_0 = 1/30
Z_0 = 0

{init_states}

d= 2

t0 = 0
odeType=stiff

beta_mod = beta
p_mod=p
c2=c

stop1 = 0
stop2 = 0

phi2 = 0
rho2 = 0

a2 = a
h2 = h


{equations1}

IFN = I/(T_0 + T_new)

if (t < init) 
  start = 0
else
 start = 1
end

{equations2}

{equations3}

VI=a2*V^h2

{derivatives_T}
ddt_R = (phi2*I*T*stop1 - rho2*R*stop1)*start
{derivatives}
ddt_V = (p_mod*I - c2*V)*start


V_log = max(log10(V),0.001)
VI_log = max(log10(VI),0.001)

OUTPUT:
output = {V_log,VI_log}
"

## Snippets defining, for each feature and each variant, the pieces of Mlxtran code to insert into the template

SNIPPETS <- list(
  
  ecl_no = list(input_params = ", delta",
                equations2 = "J = (D+I)/(T_0+ T_new)",
                equations3 = "delta2 = delta\n N = T + I",
                derivatives = "ddt_I = (beta_mod*VI*T  - delta2*I)*start\nddt_D = delta2*I*start"),
  
  ecl_expo = list(input_params = ", delta, k",
                  init_states = "E_0 = 0",
                  equations2 = "J = (D+I+E)/(T_0 + T_new)",
                  equations3 = "delta2=delta\nN = T + I + E",
                  derivatives = "ddt_E = (beta_mod*VI*T - k*E)*start\nddt_I = (k*E - delta2*I)*start\nddt_D = delta2*I*start"),
  
  Ref_yes = list(input_params = ", rho, phi",
                 equations1 = "stop1=1\nphi2=phi\nrho2=rho"),
  
  Ref_no= list(),
  
  TG_yes = list(input_params = ",r",
                equations1 = "T_new = Z",
                derivatives_T = "ddt_T = (-beta_mod*VI*T - phi2*I*T*stop1 + rho2*R*stop1 + r*T*(1-N/T_0))*start\nddt_Z =r*T*(1-N/T_0)*start"), 
  
  TG_no = list(equations1 = "T_new = 0",
               derivatives_T = "ddt_T = (-beta_mod*VI*T - phi2*I*T*stop1 + rho2*R*stop1 )*start"), 
  
  in_p_no = list(), 
  
  in_p_Hill = list(input_params = ", K_p",
                   equations3 = "p_mod = p*(1 - IFN/(IFN+K_p))"), 
  
  adap_h_no = list(), 
  
  adap_h_Hill = list(input_params = ", C_h",
                     equations3 = "h2 = h*(1 - J/(J + C_h))"),
  
  adap_a_no = list(),
  
  adap_a_Hill = list(input_params = ",  C_a",
                     equations3 = "a2 = a*(1 - J/(J + C_a))"),
  
  adap_delta_no = list(),
  
  adap_delta_Hill = list(input_params = ", delta_adap,tau",
                         equations3 = "if ( t < tau)\n delta2= delta \n else\ndelta2 = delta_adap\nend"),
  
  adap_c_no = list(),
  
  adap_c_Hill= list(input_params = ", c_max, C_c",
                    equations3 = "c2 = c*(1 + c_max*J/(J + C_c))")
)

## generate_model(): builds the Mlxtran text for the given feature combination and writes it to out_dir. 

generate_model <- function(features, out_dir) {
  
  dir.create(out_dir, showWarnings = FALSE)
  
  # Initialize slots
  slots <- list(
    input_params = "",
    init_states = "",
    equations1 = "",
    equations2 = "",
    equations3 = "",
    derivatives_T = "",
    derivatives = "")
  
  # Collect snippets
  for (feat in names(features)) {
    opt <- features[[feat]]
    key <- paste0(feat, "_", opt)
    snip <- SNIPPETS[[key]]
    
    if (is.null(snip) && !key %in% names(SNIPPETS)) {
      valid_opts <- gsub(paste0("^", feat, "_"), "",
                         grep(paste0("^", feat, "_"), names(SNIPPETS), value = TRUE))
      stop(sprintf("Unknown option '%s' for feature '%s'. Valid options: %s",
                   opt, feat, paste(valid_opts, collapse = ", ")))
    }
    
    if (!is.null(snip)) {
      for (block in names(snip)) {
        if (block == "input_params") {
          slots[[block]] <- paste0(slots[[block]], snip[[block]])
        } else {
          slots[[block]] <- paste(slots[[block]], snip[[block]], sep = "\n")
        }
      }
    }
  }
  
  # Fill template
  model_text <- TEMPLATE
  for (slot in names(slots)) {
    model_text <- gsub(paste0("\\{", slot, "\\}"), slots[[slot]], model_text)
  }
  
  # Build filename
  name_parts <- paste0(names(features), "_", unlist(features), collapse = "-")
  fname <- file.path(out_dir, paste0(name_parts, "-Model.txt"))
  
  # Write file
  writeLines(model_text, fname)
  
  return(fname)
}

## Examples
features = list(ecl = "no", 
                Ref = "no",
                TG = "no",  
                in_p = "no", 
                adap_h = "no", 
                adap_a = "no", 
                adap_delta = "no", 
                adap_c = "no" )

out_dir <- "generated_models" #" Change to your desired output directory

## Generation of a model
model = generate_model(features, out_dir)

