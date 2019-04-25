# Packages
library(KEGGREST)
library(stringr)
library(tidyverse)
library(dplyr)
library(purrr)
library(magrittr)

###############################################################################################
# Temp data, to extract codes for pathways
setwd("Data/") # set working directory where input files are located
EC_files <- list.files() # list all input files 
EC_CombinedFiles <- do.call(rbind, lapply(EC_files, function(x) read.csv(x, header = FALSE, 
                                                                         sep = "\t",stringsAsFactors = FALSE)))
setwd("../") # return wd to parent directory

EC_records <- EC_CombinedFiles %>% select(V3) %>%
  mutate(EC_number= paste0("EC: ", V3)) # select enzyme records and transform in EC format
EC_unique <- unique (EC_records$EC_number) # select unique records only

EC_KeggPathways <- vector("list",length(EC_unique))
for(i in 1:length(EC_unique)){
  EC_KeggPathways[i]=tryCatch(map(keggGet(EC_unique[i]), extract, c("ENTRY", "NAME", "PATHWAY")),
                              error=function(e) NULL)
}
library(dplyr)
KeggPathways_df=bind_rows(lapply(EC_KeggPathways, function (x) data.frame(t(unlist(x)),stringsAsFactors = F)))
library(tidyr)
Pathways_Kegg_final= gather(KeggPathways_df,key = "", value = PATHWAY, -NAME)[,c(1,3)]


# list of enzymes to verify wich linked to KEGG database 
l_ecverify <- lapply (EC_KeggPathways, function(m) m["ENTRY"])
df_ecverify <- data.frame(ec_name= unlist(l_ecverify))
write.csv(df_ecverify, file = "ECnumbers_withLink_KEGG.csv")# enzymes KEGG

####################### TRIAL
EC_KeggPath <- vector("list",length(EC_unique))
for(i in 1:length(EC_unique)){
  EC_KeggPath[i]=tryCatch(map(keggGet(EC_unique[i]), extract, c("ENTRY", "NAME", "PATHWAY")),
                          error=function(e) list(ENTRY = EC_unique[i], NAME= "NULL", PATHWAY= "NULL"))
}
library(dplyr)
KeggPathways_df=bind_rows(lapply(EC_KeggPath, function (x) data.frame(t(unlist(x)),stringsAsFactors = F)))
library(tidyr)
Pathways_Kegg_final= gather(KeggPathways_df,key = "", value = PATHWAY, -ENTRY)[,c(1,3)]

trial_list <- list(ENTRY = NULL, NAME= NULL, PATHWAY= NULL)
trial_list
