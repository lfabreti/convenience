# Convenience

R package for convergence assessment for phylogenetic inference.

## Installation

You can install the package using devtools:
  
  > install.packages("devtools") <br />
  > library(devtools) <br />
  > install_github("lfabreti/convenience") <br />
  > library(convenience) <br />
  
  
---------------------------------------------------------

## How to use

> 1. The first step is to load your output files from RevBayes. You can do it by giving the path to the output, all the files for all your runs will be loaded

> 2. The functions are divided to check for convergence in _**Continuous parameters**_ and _**Trees**_

> 3. You can see the output for each function separately or run one wrap up function that automatically checks for convergence 

## Functions

### For _**Trees**_

#### Load output

`loadFiles()`

### For _**Continuous Parameters**_

#### Difference between split frquencies

`splitFreq()`

#### ESS 

`essContParam()`

#### PSRF

`psrfContParam()`

#### Difference in measures of dispersion

`statsContParam()`

#### Standard error of the mean

`stderrContParam()`

### Wrap up function

`checkConvergence()`

## Example

To test the package, download the example output. There are two output examples provided, one with 2 runs and the other with 4 runs.
