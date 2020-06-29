# Convenience

R package for convergence assessment for phylogenetic inference.

## Installation

You can install the package using devtools:
  
  > `install.packages("devtools")` <br />
  > `library(devtools)` <br />
  > `install_github("lfabreti/convenience")` <br />
  > `library(convenience)` <br />
  
  
---------------------------------------------------------

## Criteria for convergence

The output parameters from a phylogenetic inference can be divided in 2 categories:

 1. Phylogenetic tree parameters (discrete parameters);
 
 2. Evolutionary model parameters (continuous parameters);

The tree parameter we use in the convergence assessment is the frequency of splits (or bipartitions). The continuous parameters are the parameters from the evolutionary substitution model and the tree length.<br />
The criteria for each type of parameter is described below:

### Discrete parameter : splits

 1. Effective Sample Size;
 
 2. Expected difference in split frequencies within windows of the same run;
 
 3. Expected difference in split frequencies between different runs.
 
### Continuous parameters : substitution model parameters and tree length

 1. Effective Sample Size;
 
 2. Kolmogorov-Smirnov test within windows of the same run;
 
 3. Kolmogorov-Smirnov test between different runs.

## Example

 `library(convenience)`<br />
 `setwd(example/4_runs/)`<br />
 `checkConvergence(path = ".", list_files = NULL, control= makeControl())`<br />
 




 
