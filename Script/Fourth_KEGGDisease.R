# Packages
library(KEGGREST)
library(stringr)
library(tidyverse)
library(dplyr)
library(purrr)
library(magrittr)

###############################################################################################
# Temp data, to extract codes for pathways
input_PathClasst <- read.csv("dataPathways.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
listPathwayst <- c(unique(input_PathClass$pathway)) # select only enzymes and write in EC number format
write.csv(listPathwayst, file = "PathwayCodes_disease.csv")# Pathways with nnumber of ocurrences

# Input
input_PathClass <- read.csv("PathwayCodes_disease.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
listcodes <- c(unique(input_PathClass$Code)) # select pathway codes
datapath<- vector("list", length(listcodes))

# extract diseases information from KEGG
kegg_list <- vector("list",length(listcodes))
for(i in 1:length(listcodes)){
  kegg_list[i]=map(keggGet(listcodes[i]), extract, c("ENTRY", "NAME", "CLASS", "DISEASE"))
}
# Unlist and select disease records
library(dplyr)
kegg_df=bind_rows(lapply(kegg_list, function (x) data.frame(t(unlist(x)),stringsAsFactors = F)))
library(tidyr)
kegg_disease=na.omit(gather(kegg_df,key = "", value = DISEASE, -NAME))[,c(1,3)]



# Select only disease information, previously checked the rows to select
diseases <- kegg_disease[251:744,1:2]
write.csv(diseases, file = "ListDiseases_pathways.csv")# List of pathways and diseases
diseases_count <- diseases %>% group_by(DISEASE)%>% summarise(count=length(DISEASE))




another <- data.frame(entry = unlist(kegg_list[[1]]),  name= unlist(kegg_list[[2]]), disease =unlist(kegg_list[[3]]))



# extract trial KeggGet() function from KEGGREST package
listcodes_trial <- c("map00010", "map00020", "map095", "map01100","map01130", "map01110","map01120") # trial codes
kegg_listt=map(keggGet(listcodes_trial), extract, c("ENTRY", "NAME", "CLASS", "DISEASE"))
