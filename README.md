Automated biological and statistical model selection for within-host viral dynamics

This repository implements the automated workflow for selection of viral dynamic model, combining a Decision Tree (DT) algorithm for biological model selection, an Ant Colony Optimization (ACO) algorithm for statistical model selection, and model averaging (MA) to account for model selection uncertainty.

Requirements
R (developed under R 4.x)
MonolixSuite, with the lixoftConnectors R package (bundled with MonolixSuite; initialize once via initializeLixoftConnectors() following Lixoft's documentation)
R packages: Rsmlx, parallel, dplyr, tidyr, data.table, mvtnorm, Matrix

The scripts are meant to be run in sequence; each stage reads files written by the previous one. Paths and identifiers below assume the working directory is the same throughout.

Code1 ──(sourced by)──> Code2 (x6, one per feature order) ──> Code3 ──> Code4 (x1 per candidate biological model) ──> Code5 ──> Code6

Step I — Biological model selection (Code2)

Because DT evaluates features sequentially, the feature order can influence the result (see Methods). 
Run Code2 6 times, once per permutation of the 3 feature groups (eclipse phase, refractory cells, and the 6 parameter-modifying features).
Edit FEATURES (top of Code2) to indicate the feature order for this run.

Intermediate Step I (Code3)

This file retains all models within ΔBICc < 2 of the best model among the models tested for the 6 orders, refits each distinct candidate individually, and writes:
Candidate_biological_models.csv: the distinct candidate biological models retained.
Unidentifiable_parameters.csv: one row per (candidate model, poorly estimated parameter), with its DT-estimated value. 

Step II — Statistical model selection (Code4)

ACO must be run once per candidate biological model retained by DT. For each candidate:

Set CANDIDATE_ID to that candidate's model_id (as it appears in Unidentifiable_parameters.csv).
Manual step: replace TEMPLATE (near the top of the file) with the resolved Mlxtran text of this specific candidate.

ACO can take a long time on large datasets (multiple candidates can be run in parallel, e.g. as separate HPC jobs, since their output files are isolated by CANDIDATE_ID). 

Intermediate Step II (GH5)

For a given CANDIDATE_ID, automatically locates ACO's last iteration, retains statistical models within ΔBICc < 2, and refits each individually.

Output files:
Values_ACO_Candidate<CANDIDATE_ID>.csv: the fixed value (or NA if left estimated) of each poorly estimated parameter, per retained statistical model.
Retained_models_Candidate<CANDIDATE_ID>.csv: the BICc of each retained statistical model.

Step III — Model averaging (GH6)

Computes BICc-based model weights for all the model selected by the workflow, and draws nrep parameter samples from the resulting mixture distribution. 
Compute median and 95% CI per parameter.

Data format

As an example to apply the workflow, we have put the SARS-CoV-2 data that can also be found in Killingley et al., PNAS,2022
Datasets are expected as Data/Data.csv. 

Configuration checklist

Before running the pipeline on a new machine or dataset, check every occurrence of:

MONOLIX_PATH / MonolixSuite install path: set in Code 2 to 5.
CANDIDATE_ID: set in Codes 4,5,6 consistent with the model_id values in Unidentifiable_parameters.csv.
RUN_LABEL: set in Code 2, one value per feature order.
The TEMPLATE variable in Code 4  must be replaced by hand for each candidate biological model (see Step II above).
