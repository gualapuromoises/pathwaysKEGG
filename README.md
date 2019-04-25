# pathwaysKEGG
Extract pathways from KEGG database using EC entries. 

## Input files
Input files contain interaction information in three columns: 
1. R_HMR_#### : reaction name
2. v#### : edge name
3. #.#.#.#: EC number, enzyme nomenclature number
### Source of data
Information source are from: 

Tissue: Breast, Lung, Kidney or Renal, and Urothelial

Condition: Cancer, Healthy

Set: MDS, nonMDS, where MDS means Minimum Dominanting Set

Example of file name: 

*EC_BreastCancer_MDS*, 

*EC_UrothelialHealthy_nonDMS*

## Functions in the scripts
### First_KEGGPathways.R

In this script you have three functions: 

#### function Pathways
This function is for extracting the list of pathways set by set and save two results: 
- List of raw pathways
- List of pathways with the number of counts of that pathway. 
* Two options are available, run for all datasets or run files one by one 

#### function statistical_test
This function is to calculate the p-value for proportions differences using prop.test() or fisher.test(). It is recomended to use Fisher test. The results of this function are: 
- Files with the p-value for compared sets (120 in total)
* Two options are available run all combinations or run for a single pair of sets. 

#### function join-function
This function joins the results of proportional test or Fisher exact test, and returns three results.
- A raw file with all the combination name, list of pathways and the p-value, 
- A file with a list pathways whose p-values were less than 0.05, 
- A file with the count of the number of pathways significantly different on any of the set combinations. 
Join functios allows to entre "pt" or "ft" to run it. 

#### Running Time
1. MDS files: aproximately 10 minutes (08'-15') by file
2. nonMDS files: aproximately 42 minutes (38'-55') by file

### Second_fisherTest.R

In this script you have three functions: 

#### function Pathways
This function is for extracting the list of pathways set by set and save two results: 
- List of raw pathways
- List of pathways with the number of counts of that pathway. 
* Two options are available, run for all datasets or run files one by one 

#### function statistical_test
This function is to calculate the p-value for proportions differences using prop.test() or fisher.test(). It is recomended to use Fisher test. The results of this function are: 
- Files with the p-value for compared sets (120 in total)
* Two options are available run all combinations or run for a single pair of sets. 

#### function join-function
This function joins the results of proportional test or Fisher exact test, and returns three results.
- A raw file with all the combination name, list of pathways and the p-value, 
- A file with a list pathways whose p-values were less than 0.05, 
- A file with the count of the number of pathways significantly different on any of the set combinations. 
Join functios allows to entre "pt" or "ft" to run it. 

#### Running Time
1. MDS files: aproximately 10 minutes (08'-15') by file
2. nonMDS files: aproximately 42 minutes (38'-55') by file

