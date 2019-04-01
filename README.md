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

## Script Function


#### Running Time
1. MDS files: aproximately 10 minutes (08'-15') by file
2. nonMDS files: aproximately 42 minutes (38'-55') by file
